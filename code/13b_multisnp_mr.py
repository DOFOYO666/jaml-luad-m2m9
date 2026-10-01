# -*- coding: utf-8 -*-
"""13b_multisnp_mr.py - M1 多工具变量 MR（IVW/MR-Egger/加权中位数），按 QTD 分组流式读取
OneK1K QTD .all（官方 beta/se/rsid, GRCh38）→ cis 区域显著 eQTL SNP → rsid 匹配 LUAD GWAS(b37)
→ 位置 clumping（500kb）→ 多工具 MR。LD 剪枝正式版待 1000G 参考到位。
"""
import gzip
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
ALL_DIR = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/onek1k/all"
GWAS = os.path.join(BASE, "data", "luad_gwas", "GCST004744.h.tsv.gz")
SIG = os.path.join(RES, "scmr_luad_risk_sig.csv")
WINDOW = 250_000


def load_gwas_rsid():
    keep = []
    for ch in pd.read_csv(GWAS, sep="\t", compression="gzip",
                          usecols=["hm_rsid", "hm_beta", "standard_error", "p_value"],
                          dtype={"hm_rsid": str}, chunksize=3_000_000):
        ch = ch.rename(columns={"hm_rsid": "rsid", "hm_beta": "beta",
                                "standard_error": "se", "p_value": "p"})
        for c in ["beta", "se", "p"]:
            ch[c] = pd.to_numeric(ch[c], errors="coerce")
        ch = ch.dropna(subset=["beta", "se"]).drop_duplicates("rsid")
        keep.append(ch[["rsid", "beta", "se", "p"]])
    g2 = pd.concat(keep, ignore_index=True).drop_duplicates("rsid", keep="first")
    return dict(zip(g2["rsid"], zip(g2["beta"], g2["se"], g2["p"])))


def ivw(df):
    w = 1 / df["se_gw"] ** 2
    b = (df["beta_gw"] * w).sum() / w.sum()
    se = np.sqrt(1 / w.sum())
    return b, se


def egger(df):
    x = df["beta"].values
    y = df["beta_gw"].values
    w = 1 / df["se_gw"].values ** 2
    X = np.column_stack([np.ones(len(x)), x])
    W = np.diag(w)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
    resid = y - X @ beta
    sigma2 = (resid ** 2 * w).sum() / (len(x) - 2)
    cov = sigma2 * np.linalg.inv(X.T @ W @ X)
    return beta[1], np.sqrt(cov[1, 1]), beta[0], np.sqrt(cov[0, 0])


def weighted_median(df):
    b = df["beta_gw"].values
    w = 1 / df["se_gw"].values ** 2
    w /= w.sum()
    order = np.argsort(b)
    cw = np.cumsum(w[order])
    return b[order][np.searchsorted(cw, 0.5)]


def clump_by_distance(df, kb=500):
    df = df.sort_values("pvalue").copy()
    kept = []
    for _, r in df.iterrows():
        if all(abs(r["pos"] - k["pos"]) > kb * 1000 for k in kept):
            kept.append(r)
    return pd.DataFrame(kept)


def main():
    sig = pd.read_csv(SIG)
    top = sig.sort_values("F", ascending=False).drop_duplicates("gene_symbol")
    gq = dict(zip(top["gene_symbol"], top["qtd"]))
    print(f"[输入] {len(gq)} 基因", flush=True)
    gwas = load_gwas_rsid()
    print(f"[GWAS] rsid 条目 {len(gwas)}", flush=True)

    coords = {}
    with gzip.open(os.path.join(BASE, "data", "onek1k", "gencode.v38.annotation.gtf.gz"), "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            m = re.search(r'gene_name "([^"]+)"', line)
            gid = re.search(r'gene_id "([^"]+)"', line)
            if not (m and gid and m.group(1) in gq):
                continue
            p = line.strip().split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            coords[m.group(1)] = (p[0].replace("chr", ""), int(p[3]), int(p[4]), gid.group(1).split(".")[0])
    print(f"[坐标] {len(coords)}", flush=True)

    by_qtd = defaultdict(list)
    for g, q in gq.items():
        if g in coords:
            by_qtd[q].append(g)

    results = []
    for qtd, genes in by_qtd.items():
        fp = os.path.join(ALL_DIR, f"{qtd}.all.tsv.gz")
        if not os.path.exists(fp):
            print(f"  QTD {qtd} 缺失", flush=True)
            continue
        targets = {coords[g][3]: g for g in genes}
        region = defaultdict(list)
        print(f"[读取] {qtd} ({len(genes)} 基因)...", flush=True)
        with gzip.open(fp, "rt", errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(hdr)}
            gi, ci, pi = idx["molecular_trait_id"], idx["chromosome"], idx["position"]
            ri, ai, vi = idx["ref"], idx["alt"], idx["variant"]
            bi, si, pvi, rsi = idx["beta"], idx["se"], idx["pvalue"], idx["rsid"]
            for line in f:
                c = line.split("\t")
                ensg = c[gi]
                if ensg not in targets:
                    continue
                try:
                    ch, pos = c[ci], int(c[pi])
                except ValueError:
                    continue
                g = targets[ensg]
                chrom, gs, ge, _ = coords[g]
                if ch == chrom and gs - WINDOW <= pos <= ge + WINDOW:
                    region[g].append((c[vi].strip(), pos, c[ri].strip(), c[ai].strip(),
                                      float(c[bi]), float(c[si]), float(c[pvi]), c[rsi].strip()))
        for gene in genes:
            if gene not in region:
                print(f"  {gene}: 无区域 SNP", flush=True)
                continue
            df = pd.DataFrame(region[gene], columns=["variant", "pos", "ref", "alt", "beta", "se", "pvalue", "rsid"])
            df = df[df["pvalue"] < 0.05].copy()
            df["rsid2"] = df["rsid"].where(df["rsid"] != ".", np.nan)
            m = df.dropna(subset=["rsid2"]).copy()
            m["gwas"] = m["rsid2"].map(gwas)
            m = m.dropna(subset=["gwas"]).copy()
            if len(m) < 2:
                print(f"  {gene}: rsid 匹配 <2", flush=True)
                continue
            m["beta_gw"] = [x[0] for x in m["gwas"]]
            m["se_gw"] = [x[1] for x in m["gwas"]]
            m = m.dropna(subset=["beta_gw", "se_gw"])
            m["b_mr"] = m["beta_gw"] / m["beta"]
            m["se_mr"] = m["se_gw"] / m["beta"].abs()
            cl = clump_by_distance(m)
            if len(cl) < 2:
                cl = m
            b_ivw, se_ivw = ivw(cl)
            p_ivw = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(b_ivw) / se_ivw))
            try:
                slope, se_slope, intercept, se_intercept = egger(cl)
                p_egger = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(slope) / se_slope))
                p_int = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(intercept) / se_intercept))
            except Exception:
                slope = se_slope = intercept = se_intercept = p_egger = p_int = np.nan
            wmed = weighted_median(cl)
            results.append({"gene": gene, "qtd": qtd, "n_snp_region": len(df), "n_clumped": len(cl),
                            "beta_ivw": b_ivw, "p_ivw": p_ivw, "or_ivw": np.exp(b_ivw),
                            "beta_egger": slope, "p_egger": p_egger,
                            "egger_intercept": intercept, "p_intercept": p_int,
                            "wmedian_beta": wmed})
            print(f"  {gene}: {len(cl)} SNP IVW OR={np.exp(b_ivw):.3f} P={p_ivw:.2g} EggerP={p_egger:.2g} 截距P={p_int:.2g}", flush=True)
    out = pd.DataFrame(results)
    out.to_csv(os.path.join(RES, "multisnp_mr_results.csv"), index=False)
    print(f"\n[完成] {len(out)} 基因; results/multisnp_mr_results.csv", flush=True)


if __name__ == "__main__":
    main()
