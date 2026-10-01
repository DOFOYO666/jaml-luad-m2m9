# -*- coding: utf-8 -*-
"""
04_bulk_validation.py - eQTLGen (bulk whole blood) 跨组织验证
对 scMR 显著（基因 × 细胞类型），用 eQTLGen 全血 eQTL 做 MR，
判断信号是"免疫细胞特异"还是"全血共享"。

eQTLGen full 格式：Pvalue SNP SNPChr SNPPos Zscore AssessedAllele OtherAllele Gene GeneSymbol ...
Zscore 基于 AssessedAllele；MR 近似 beta_exp=Zscore, se_exp=1（TwoSampleMR 标准做法）
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
RES = os.path.join(BASE, "results")
EQTLGEN = os.path.join(DATA, "eqtlgen", "cis-eQTL_significant_20181017.txt.gz")
GWAS = os.path.join(DATA, "luad_gwas", "GCST004744.h.tsv.gz")
SIG = os.path.join(RES, "scmr_luad_risk_sig.csv")


def main():
    if not os.path.exists(SIG):
        print("[跳过] 无 scMR 显著结果"); return
    sig = pd.read_csv(SIG)
    genes = sorted(set(sig["gene_symbol"]))
    print(f"[输入] scMR 显著基因 {len(genes)} 个")

    # GWAS 载入（仅 hm_rsid + beta/se/p + allele）
    gcols = ["hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "p_value"]
    gchunks = []
    for ch in pd.read_csv(GWAS, sep="\t", compression="gzip", chunksize=2_000_000,
                          low_memory=False, usecols=gcols):
        gchunks.append(ch)
    g = pd.concat(gchunks, ignore_index=True)
    g = g.rename(columns={"hm_beta": "beta", "standard_error": "se", "p_value": "p"})
    g["ea"] = g["hm_effect_allele"].astype(str).str.upper()
    g["oa"] = g["hm_other_allele"].astype(str).str.upper()
    for c in ["beta", "se", "p"]:
        g[c] = pd.to_numeric(g[c], errors="coerce")
    g = g[g["beta"].notna() & g["se"].notna() & (g["beta"] != 0)]
    g = g.drop_duplicates(subset=["hm_rsid"], keep="first")
    g_map = g.set_index("hm_rsid")
    print(f"[GWAS] {g.shape}")

    # 流式读取 eQTLGen，只保留显著基因
    eqtl = []
    for ch in pd.read_csv(EQTLGEN, sep="\t", compression="gzip", chunksize=1_000_000,
                          low_memory=False,
                          usecols=["SNP", "Zscore", "AssessedAllele", "OtherAllele",
                                   "Pvalue", "FDR", "GeneSymbol"]):
        ch = ch[ch["GeneSymbol"].isin(genes)]
        if not ch.empty:
            eqtl.append(ch)
    eq = pd.concat(eqtl, ignore_index=True)
    eq["Zscore"] = pd.to_numeric(eq["Zscore"], errors="coerce")
    eq["Pvalue"] = pd.to_numeric(eq["Pvalue"], errors="coerce")
    eq = eq.dropna(subset=["Zscore", "Pvalue"])
    eq["ea"] = eq["AssessedAllele"].astype(str).str.upper()
    eq["oa"] = eq["OtherAllele"].astype(str).str.upper()
    print(f"[eQTLGen] 显著基因命中 {len(eq)} 行, {eq['GeneSymbol'].nunique()} 基因")

    # 匹配（rsid） + 方向对齐 + Wald ratio
    g_sub = g_map.reset_index().rename(columns={"ea": "ea_g", "oa": "oa_g"})
    m = eq.merge(g_sub, left_on="SNP", right_on="hm_rsid", how="inner")
    same = (m["ea"] == m["ea_g"]) & (m["oa"] == m["oa_g"])
    flip = (m["ea"] == m["oa_g"]) & (m["oa"] == m["ea_g"])
    m["beta_g"] = np.where(same, m["beta"], np.where(flip, -m["beta"], np.nan))
    m = m.dropna(subset=["beta_g", "Zscore"])
    m = m[m["Zscore"] != 0]
    m["b_mr"] = m["beta_g"] / m["Zscore"]
    m["se_mr"] = np.abs(m["se"] / m["Zscore"])
    m["p_mr"] = 2 * (1 - stats.norm.cdf(np.abs(m["b_mr"] / m["se_mr"])))
    m["or"] = np.exp(m["b_mr"])
    # 每基因取最强信号
    m = m.sort_values("Pvalue").drop_duplicates(subset=["GeneSymbol"], keep="first")
    m["fdr"] = stats.false_discovery_control(m["p_mr"], method="bh")
    m = m.sort_values("p_mr")
    m.to_csv(os.path.join(RES, "bulk_eqtlgen_validation.csv"), index=False)

    print(f"\n[完成] eQTLGen 验证 {len(m)} 基因")
    print("eQTLGen FDR<0.05:", (m["fdr"] < 0.05).sum(), "| P<0.05:", (m["p_mr"] < 0.05).sum())
    print("\n结果（按 p 排序）:")
    print(m[["GeneSymbol", "SNP", "Zscore", "or", "p_mr", "fdr"]].to_string(index=False))


if __name__ == "__main__":
    main()
