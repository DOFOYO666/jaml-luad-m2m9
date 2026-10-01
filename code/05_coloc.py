# -*- coding: utf-8 -*-
"""
05_coloc.py - 贝叶斯共定位 (coloc.abf)
对 scMR 显著基因，用 eQTLGen (bulk blood) cis-eQTL 区域数据与 LUAD GWAS 做共定位，
评估共享因果变异证据（PP.H4 > 0.75 为强证据）。

coloc.abf 依据：Giambartolomei et al. 2014, PLoS Genet
eQTLGen beta 近似：beta=Zscore, varbeta=1（与 TwoSampleMR 一致）
"""
import os
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
RES = os.path.join(BASE, "results")
EQTLGEN = os.path.join(DATA, "eqtlgen", "cis-eQTL_significant_20181017.txt.gz")
GWAS = os.path.join(DATA, "luad_gwas", "GCST004744.h.tsv.gz")
SIG = os.path.join(RES, "scmr_luad_risk_sig.csv")
WINDOW = 500_000
N_EQTLGEN = 31684  # eQTLGen significant 子集
N_GWAS = 66756


def coloc_abf(beta1, varbeta1, beta2, varbeta2, p1=1e-4, p2=1e-4, p12=1e-5):
    beta1 = np.asarray(beta1, float); varbeta1 = np.asarray(varbeta1, float)
    beta2 = np.asarray(beta2, float); varbeta2 = np.asarray(varbeta2, float)
    if len(beta1) < 3:
        return None
    W = 0.15 ** 2
    lABF1 = 0.5 * (np.log(varbeta1 / (varbeta1 + W)) + (beta1 ** 2 / (varbeta1 + W)) * (W / varbeta1))
    lABF2 = 0.5 * (np.log(varbeta2 / (varbeta2 + W)) + (beta2 ** 2 / (varbeta2 + W)) * (W / varbeta2))
    l1 = lABF1 + np.log(p1) - np.log(1 - p1)
    l2 = lABF2 + np.log(p2) - np.log(1 - p2)
    lsum = l1 + l2 + np.log(p12) - np.log(p1 * p2)
    logH0 = 0.0
    logH1 = np.logaddexp.reduce(l1)
    logH2 = np.logaddexp.reduce(l2)
    logH3 = np.logaddexp.reduce(l1 + l2)
    logH4 = np.logaddexp.reduce(lsum)
    logall = np.logaddexp.reduce([logH0, logH1, logH2, logH3, logH4])
    pp = np.exp([logH0 - logall, logH1 - logall, logH2 - logall, logH3 - logall, logH4 - logall])
    return dict(zip(["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"], pp))


def main():
    if not os.path.exists(SIG):
        print("[跳过] 无 scMR 显著结果"); return
    sig = pd.read_csv(SIG)
    genes = sorted(set(sig["gene_symbol"]))
    print(f"[输入] coloc 候选基因 {len(genes)} 个")

    # GWAS 载入（chr/pos/rsid/beta/se）
    gcols = ["hm_chrom", "hm_pos", "hm_rsid", "hm_beta", "standard_error"]
    gchunks = []
    for ch in pd.read_csv(GWAS, sep="\t", compression="gzip", chunksize=2_000_000,
                          low_memory=False, usecols=gcols):
        gchunks.append(ch)
    g = pd.concat(gchunks, ignore_index=True)
    g = g.rename(columns={"hm_beta": "beta", "standard_error": "se"})
    g["chrom"] = pd.to_numeric(g["hm_chrom"], errors="coerce")
    g["pos"] = pd.to_numeric(g["hm_pos"], errors="coerce")
    g["beta"] = pd.to_numeric(g["beta"], errors="coerce")
    g["se"] = pd.to_numeric(g["se"], errors="coerce")
    g = g.dropna(subset=["chrom", "pos", "beta", "se"])
    g["varbeta"] = g["se"] ** 2
    g = g.drop_duplicates(subset=["hm_rsid"], keep="first")
    print(f"[GWAS] {g.shape}")

    # 流式提取显著基因的 eQTLGen 数据
    eqtl = []
    for ch in pd.read_csv(EQTLGEN, sep="\t", compression="gzip", chunksize=1_000_000,
                          low_memory=False,
                          usecols=["SNP", "SNPChr", "SNPPos", "Zscore", "AssessedAllele",
                                   "OtherAllele", "GeneSymbol", "GeneChr", "GenePos"]):
        ch = ch[ch["GeneSymbol"].isin(genes)]
        if not ch.empty:
            eqtl.append(ch)
    eq = pd.concat(eqtl, ignore_index=True)
    for c in ["Zscore", "SNPChr", "SNPPos", "GeneChr", "GenePos"]:
        eq[c] = pd.to_numeric(eq[c], errors="coerce")
    eq = eq.dropna(subset=["Zscore", "SNPChr", "SNPPos"])
    eq["ea"] = eq["AssessedAllele"].astype(str).str.upper()
    eq["oa"] = eq["OtherAllele"].astype(str).str.upper()
    print(f"[eQTLGen] 提取 {len(eq)} 行, {eq['GeneSymbol'].nunique()} 基因")

    results = []
    for gene in genes:
        ge = eq[eq["GeneSymbol"] == gene]
        if ge.empty:
            continue
        gpos = ge["GenePos"].median()
        gchr = ge["GeneChr"].median()
        # eQTL 区域（该基因 cis 窗口内的 SNP）
        ge = ge.drop_duplicates(subset=["SNP"], keep="first")
        # GWAS 区域
        region = g[(g["chrom"] == gchr) &
                   (g["pos"] >= gpos - WINDOW) & (g["pos"] <= gpos + WINDOW)]
        if len(region) < 3:
            continue
        # 按 chr:pos 合并
        region = region.rename(columns={"chrom": "chr_g", "pos": "pos_g"})
        ge = ge.rename(columns={"SNPChr": "chr_e", "SNPPos": "pos_e"})
        merged = region.merge(ge, left_on=["chr_g", "pos_g"], right_on=["chr_e", "pos_e"],
                              how="inner")
        merged = merged.drop_duplicates(subset=["SNP"], keep="first")
        if len(merged) < 3:
            continue
        # 方向对齐（简单版：eQTLGen Zscore 基于 AssessedAllele；GWAS hm_beta 基于 hm_effect_allele，
        # 此处未载入 hm_effect_allele，用绝对值近似——coloc 对方向敏感度有限，且主要看信号共定位）
        beta_e = merged["Zscore"].values
        varbeta_e = np.ones(len(merged))
        beta_g = merged["beta"].values
        varbeta_g = merged["varbeta"].values
        pp = coloc_abf(beta_e, varbeta_e, beta_g, varbeta_g)
        if pp:
            pp["gene"] = gene
            pp["n_snp"] = len(merged)
            results.append(pp)

    res = pd.DataFrame(results)
    if res.empty:
        print("无共定位结果"); return
    res = res.sort_values("PP.H4", ascending=False)
    res.to_csv(os.path.join(RES, "coloc_results.csv"), index=False)
    print(f"\n[完成] coloc 结果 {len(res)} 基因")
    print(f"PP.H4>0.75: {(res['PP.H4'] > 0.75).sum()} | PP.H4>0.5: {(res['PP.H4'] > 0.5).sum()}")
    print("\n结果（按 PP.H4 排序）:")
    print(res[["gene", "n_snp", "PP.H0", "PP.H3", "PP.H4"]].to_string(index=False))


if __name__ == "__main__":
    main()
