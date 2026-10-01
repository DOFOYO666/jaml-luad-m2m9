# -*- coding: utf-8 -*-
"""05b_coloc_gtex.py - GTEx Lung allpairs 版贝叶斯共定位（正式版 coloc）
流程:
1) 从 gencode v38 提取 21 个显著基因坐标(b38)，cis 窗口 ±250kb
2) GTEx Lung 全量文件流式提取 21 基因区域 eQTL(beta/se/p)
3) LUAD GWAS(b37) 提取同区域 SNP(beta/se/p)
4) chr:ref:alt 合并 → Wakefield coloc.abf → PP.H0..H4
5) 输出 results/coloc_gtex_lung.csv
用法: python 05b_coloc_gtex.py <gtex_lung_path>
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")
GTF = os.path.join(DATA, "onek1k", "gencode.v38.annotation.gtf.gz")
GWAS = os.path.join(DATA, "luad_gwas", "GCST004744.h.tsv.gz")
SIG_CSV = os.path.join(RES, "scmr_luad_risk_sig.csv")
WINDOW = 250_000  # cis 窗口 ±250kb


def logsumexp(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


def coloc_abf(beta1, se1, beta2, se2, p1=1e-4, p2=1e-4, p12=1e-5, w=0.15):
    """Wakefield 近似 coloc.abf，返回 (PP.H0..H4, n_snp)。"""
    beta1 = np.asarray(beta1, float)
    se1 = np.asarray(se1, float)
    beta2 = np.asarray(beta2, float)
    se2 = np.asarray(se2, float)
    m = min(len(beta1), len(beta2), len(se1), len(se2))
    beta1, se1, beta2, se2 = beta1[:m], se1[:m], beta2[:m], se2[:m]
    if m == 0:
        return (np.nan,) * 5, 0
    v1, v2 = se1 ** 2, se2 ** 2
    z1, z2 = beta1 / se1, beta2 / se2
    lABF1 = 0.5 * (np.log(v1 / (v1 + w)) + z1 ** 2 * w / (v1 + w))
    lABF2 = 0.5 * (np.log(v2 / (v2 + w)) + z2 ** 2 * w / (v2 + w))
    lH0 = 0.0
    lH1 = logsumexp(lABF1)
    lH2 = logsumexp(lABF2)
    lH3 = logsumexp(lABF1) + logsumexp(lABF2)
    lH4 = logsumexp(lABF1 + lABF2) + np.log(p12) - np.log(p1) - np.log(p2)
    denom = logsumexp(np.array([lH0, lH1, lH2, lH3, lH4]))
    pp = np.exp(np.array([lH0, lH1, lH2, lH3, lH4]) - denom)
    return tuple(pp), m


def load_genes():
    sig = pd.read_csv(SIG_CSV)
    genes = sorted(sig["gene_symbol"].dropna().unique())
    ensg2sym = dict(zip(sig["gene_symbol"], sig["molecular_trait_id"]))
    # gencode v38 提取基因坐标
    coords = {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            attrs = {}
            for kv in p[8].rstrip(";").split(";"):
                kv = kv.strip()
                if " " in kv:
                    k, v = kv.split(" ", 1)
                    attrs[k] = v.strip('"')
            gid = attrs.get("gene_id", "").split(".")[0]
            sym = attrs.get("gene_name", "")
            if sym in genes:
                coords[sym] = (p[0].replace("chr", ""), int(p[3]), int(p[4]), gid)
    return genes, coords, ensg2sym


def extract_gtex(gtex_path, coords, ensg2sym):
    """流式扫描 GTEx Lung 全量文件，提取 21 基因 cis 区域 eQTL。"""
    # 列索引: variant(0) pvalue(2) gene_id(6) beta(8) se(9) chromosome(12) position(13) ref(14) alt(15)
    cols = {"variant": 0, "pvalue": 2, "gene_id": 6, "beta": 8, "se": 9,
            "chrom": 12, "pos": 13, "ref": 14, "alt": 15}
    out = {}
    n_total = 0
    n_hit = 0
    with gzip.open(gtex_path, "rt", errors="replace") as f:
        header = f.readline().split("\t")
        for line in f:
            n_total += 1
            if n_total % 5_000_000 == 0:
                print(f"  [GTEx 扫描] {n_total/1e6:.0f}M 行, 命中 {n_hit}", flush=True)
            c = line.split("\t")
            if len(c) < 17:
                continue
            try:
                chrom = c[cols["chrom"]]
                pos = int(c[cols["pos"]])
                gene = c[cols["gene_id"]]
                beta = float(c[cols["beta"]])
                se = float(c[cols["se"]])
                pval = float(c[cols["pvalue"]])
                ref = c[cols["ref"]]
                alt = c[cols["alt"]]
            except ValueError:
                continue
            # 目标基因的 cis 区域
            for sym, (gchr, gstart, gend, gid) in coords.items():
                if gene == gid and gchr == chrom and gstart - WINDOW <= pos <= gend + WINDOW:
                    out.setdefault(sym, []).append((chrom, pos, ref, alt, beta, se, pval))
                    n_hit += 1
                    break
    print(f"  [GTEx] 总行 {n_total/1e6:.0f}M, 命中 {n_hit}")
    return out


def load_gwas_region(coords):
    """提取 GWAS(b37) 区域 SNP。用基因 b38 坐标近似（局部平移 <10kb 可接受）。"""
    gwas = {}
    with gzip.open(GWAS, "rt", errors="replace") as f:
        header = f.readline().split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in f:
            c = line.split("\t")
            try:
                chrom = c[idx.get("chromosome", 0)]
                pos = int(c[idx.get("base_pair_location", 1)])
            except (ValueError, KeyError):
                continue
            for sym, (gchr, gstart, gend) in coords.items():
                if gchr == chrom and gstart - WINDOW <= pos <= gend + WINDOW:
                    gwas.setdefault(sym, []).append(c)
                    break
    return gwas


def main():
    gtex_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA, "gtex", "Lung.tsv.gz")
    genes, coords, ensg2sym = load_genes()
    print(f"[输入] 显著基因 {len(genes)} 个")
    for g in genes:
        print(f"  {g}: {coords.get(g)}")
    # 缺少坐标的基因
    missing = [g for g in genes if g not in coords]
    if missing:
        print("缺少 gencode 坐标:", missing)
        return
    print(f"\n[GTEx] 提取区域 eQTL ...")
    gtex = extract_gtex(gtex_path, coords, ensg2sym)
    print(f"[GWAS] 提取区域 SNP ...")
    # GWAS 列名探测
    with gzip.open(GWAS, "rt") as f:
        print("  GWAS 表头:", f.readline().strip()[:200])
    rows = []
    with gzip.open(GWAS, "rt", errors="replace") as f:
        header = f.readline().split("\t")
        idx = {h: i for i, h in enumerate(header)}
        # harmonized 列映射（hm_* 为 b37 标准方向）
        bcol = idx.get("hm_beta")
        scol = idx.get("standard_error")
        pcol = idx.get("p_value")
        ccol = idx.get("hm_chrom")
        pcol_n = idx.get("hm_pos")
        rcol = idx.get("hm_other_allele")
        acol = idx.get("hm_effect_allele")
        print(f"  列映射: beta={bcol} se={scol} p={pcol} chr={ccol} pos={pcol_n} ref={rcol} alt={acol}")
        for line in f:
            c = line.split("\t")
            try:
                chrom = c[ccol]
                pos = int(c[pcol_n])
                beta = float(c[bcol])
                se = float(c[scol])
                pval = float(c[pcol])
                ref = c[rcol].strip()
                alt = c[acol].strip()
            except (ValueError, IndexError, TypeError):
                continue
            for sym, (gchr, gstart, gend) in coords.items():
                if gchr == chrom and gstart - WINDOW <= pos <= gend + WINDOW:
                    rows.append({"gene": sym, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                                 "beta": beta, "se": se, "pvalue": pval})
                    break
    gwas_df = pd.DataFrame(rows)
    print(f"[GWAS] 区域 SNP 行数: {len(gwas_df)}")

    # 合并 coloc
    results = []
    for sym in genes:
        gt = gtex.get(sym, [])
        if not gt:
            print(f"  {sym}: GTEx 无区域 eQTL")
            continue
        gdf = pd.DataFrame(gt, columns=["chrom", "pos", "ref", "alt", "beta", "se", "pvalue"])
        gdf["key"] = gdf["chrom"] + ":" + gdf["ref"] + ":" + gdf["alt"]
        wdf = gwas_df[gwas_df["gene"] == sym].copy()
        wdf["key"] = wdf["chrom"] + ":" + wdf["ref"] + ":" + wdf["alt"]
        merged = gdf.merge(wdf[["key", "beta", "se"]], on="key", suffixes=("_eqtl", "_gwas"))
        if len(merged) < 3:
            print(f"  {sym}: 仅 {len(merged)} 个共享 SNP，跳过")
            results.append({"gene": sym, "n_snp": len(merged), "PP.H0": np.nan, "PP.H1": np.nan,
                            "PP.H2": np.nan, "PP.H3": np.nan, "PP.H4": np.nan})
            continue
        pp, n = coloc_abf(merged["beta_eqtl"], merged["se_eqtl"], merged["beta_gwas"], merged["se_gwas"])
        print(f"  {sym}: {n} SNP 共享, PP.H4={pp[4]:.4f}")
        results.append({"gene": sym, "n_snp": n, "PP.H0": pp[0], "PP.H1": pp[1],
                        "PP.H2": pp[2], "PP.H3": pp[3], "PP.H4": pp[4]})
    out = pd.DataFrame(results).sort_values("PP.H4", ascending=False)
    out.to_csv(os.path.join(RES, "coloc_gtex_lung.csv"), index=False)
    print("\n[完成] 结果已保存 results/coloc_gtex_lung.csv")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
