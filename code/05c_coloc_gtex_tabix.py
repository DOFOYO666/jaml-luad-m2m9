# -*- coding: utf-8 -*-
"""05c_coloc_gtex_tabix.py - GTEx Lung(eQTL Catalogue tabix QTD000271) × LUAD GWAS 贝叶斯共定位
正式版 coloc：用 eQTL Catalogue 标准 tabix 全量区域数据（含非显著 SNP），
替代此前 eQTLGen significant 子集（PP.H4 偏保守的局限）。
用法: python 05c_coloc_gtex_tabix.py
"""
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tools_tabix import TabixQuery

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")
GTF = os.path.join(DATA, "onek1k", "gencode.v38.annotation.gtf.gz")
GWAS = os.path.join(DATA, "luad_gwas", "GCST004744.h.tsv.gz")
SIG_CSV = os.path.join(RES, "scmr_luad_risk_sig.csv")
GTEX_URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000271/QTD000271.all.tsv.gz"
GTEX_TBI = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/gtex/GTEx_Lung_QTD000271.tbi"
WINDOW = 250_000
GWAS_WINDOW = 600_000  # b37/b38 在 chr15 等区域存在 ~400kb 偏移，GWAS 提取放宽


def logsumexp(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


def coloc_abf(beta1, se1, beta2, se2, p1=1e-4, p2=1e-4, p12=1e-5, w=0.15):
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
            sym = attrs.get("gene_name", "")
            if sym in genes:
                gchr = p[0].replace("chr", "")
                coords[sym] = (gchr, int(p[3]), int(p[4]), attrs.get("gene_id", "").split(".")[0])
    return genes, coords, ensg2sym


def main():
    genes, coords, ensg2sym = load_genes()
    missing = [g for g in genes if g not in coords]
    if missing:
        print("缺少 gencode 坐标:", missing)
        return
    print(f"[输入] 显著基因 {len(genes)} 个")
    tq = TabixQuery(GTEX_URL, GTEX_TBI)

    # 1) GTEx Lung 区域 eQTL（缓存避免重复提取）
    CACHE = os.path.join(RES, "gtex_lung_regions.pkl")
    if os.path.exists(CACHE):
        gtex_df = pd.read_pickle(CACHE)
        print(f"[GTEx] 从缓存加载 {len(gtex_df)} 行")
    else:
        gtex_rows = []
        for sym, (chrom, gstart, gend, gid) in coords.items():
            beg = max(1, gstart - WINDOW)
            end = gend + WINDOW
            rows = tq.fetch(chrom, beg, end, chrom_col=2, pos_col=3, header_rows=1)
            hit = 0
            for c in rows:
                if len(c) < 19:
                    continue
                try:
                    # QTD sumstats 列: 0=molecular_trait_id 1=chromosome 2=position 3=ref 4=alt
                    #   9=beta 10=se 16=gene_id 18=rsid
                    gene = c[16]
                    beta = float(c[9])
                    se = float(c[10])
                    ref = c[3]
                    alt = c[4]
                    pos = int(c[2])
                    rsid = c[18].strip()
                except (ValueError, IndexError):
                    continue
                if gene == gid:
                    gtex_rows.append({"gene": sym, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                                      "beta": beta, "se": se, "rsid": rsid})
                    hit += 1
            print(f"  {sym}: 区域 {chrom}:{beg}-{end} 原始行 {len(rows)}, 目标基因行 {hit}", flush=True)
        gtex_df = pd.DataFrame(gtex_rows)
        gtex_df.to_pickle(CACHE)
        print(f"[GTEx] 总行: {len(gtex_df)}")

    # 2) GWAS 区域（b37, harmonized 列；区域放宽 ±600kb 覆盖 b37/b38 偏移）
    gwas_rows = []
    with gzip.open(GWAS, "rt", errors="replace") as f:
        header = f.readline().split("\t")
        idx = {h: i for i, h in enumerate(header)}
        bcol, scol, pcol = idx["hm_beta"], idx["standard_error"], idx["p_value"]
        ccol, pcol_n, rcol, acol = idx["hm_chrom"], idx["hm_pos"], idx["hm_other_allele"], idx["hm_effect_allele"]
        rscol = idx["hm_rsid"]
        for line in f:
            c = line.split("\t")
            try:
                chrom = c[ccol]
                pos = int(c[pcol_n])
                beta = float(c[bcol])
                se = float(c[scol])
                pval = float(c[pcol])
                ref, alt = c[rcol].strip(), c[acol].strip()
                rsid = c[rscol].strip()
            except (ValueError, IndexError):
                continue
            for sym, (gchr, gstart, gend, gid) in coords.items():
                if gchr == chrom and gstart - GWAS_WINDOW <= pos <= gend + GWAS_WINDOW:
                    gwas_rows.append({"gene": sym, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                                      "beta": beta, "se": se, "pvalue": pval, "rsid": rsid})
    gwas_df = pd.DataFrame(gwas_rows)
    print(f"[GWAS] 区域行: {len(gwas_df)}")

    # 3) 合并（rsid 优先，fallback chr:pos:ref:alt）+ coloc
    results = []
    for sym in genes:
        gd = gtex_df[gtex_df["gene"] == sym].copy()
        wd = gwas_df[gwas_df["gene"] == sym].copy()
        if gd.empty or wd.empty:
            results.append({"gene": sym, "n_snp": 0, "PP.H0": np.nan, "PP.H1": np.nan,
                            "PP.H2": np.nan, "PP.H3": np.nan, "PP.H4": np.nan})
            continue
        # key1: rsid; key2: chrom:pos:ref:alt
        gd["k1"] = np.where(gd["rsid"].isin(["NA", "", "."]), None, gd["rsid"])
        wd["k1"] = np.where(wd["rsid"].isin(["NA", "", "."]), None, wd["rsid"])
        gd["k2"] = gd["chrom"] + ":" + gd["pos"].astype(str) + ":" + gd["ref"] + ":" + gd["alt"]
        wd["k2"] = wd["chrom"] + ":" + wd["pos"].astype(str) + ":" + wd["ref"] + ":" + wd["alt"]
        wd["gwas_effect"] = wd["alt"]  # hm_effect_allele
        m1 = gd.dropna(subset=["k1"]).merge(wd.dropna(subset=["k1"])[["k1", "beta", "se", "gwas_effect"]],
                                            left_on="k1", right_on="k1", suffixes=("_eqtl", "_gwas"))
        m2 = gd.merge(wd[["k2", "beta", "se", "gwas_effect"]], on="k2", suffixes=("_eqtl", "_gwas"))
        merged = pd.concat([m1, m2], ignore_index=True).drop_duplicates(subset=["k2"])
        if len(merged) >= 3:
            # 方向校正: GTEx beta 对应 alt allele；GWAS beta 对应 effect allele(hm_effect_allele)
            merged["beta_gwas"] = np.where(merged["alt"] == merged["gwas_effect"],
                                          merged["beta_gwas"], -merged["beta_gwas"])
        if len(merged) < 3:
            print(f"  {sym}: 仅 {len(merged)} 共享 SNP")
            results.append({"gene": sym, "n_snp": len(merged), "PP.H0": np.nan, "PP.H1": np.nan,
                            "PP.H2": np.nan, "PP.H3": np.nan, "PP.H4": np.nan})
            continue
        pp, n = coloc_abf(merged["beta_eqtl"], merged["se_eqtl"], merged["beta_gwas"], merged["se_gwas"])
        print(f"  {sym}: {n} SNP 共享, PP.H4={pp[4]:.4f}")
        results.append({"gene": sym, "n_snp": n, "PP.H0": pp[0], "PP.H1": pp[1],
                        "PP.H2": pp[2], "PP.H3": pp[3], "PP.H4": pp[4]})
    out = pd.DataFrame(results).sort_values("PP.H4", ascending=False)
    out.to_csv(os.path.join(RES, "coloc_gtex_lung.csv"), index=False)
    print("\n[完成] results/coloc_gtex_lung.csv")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
