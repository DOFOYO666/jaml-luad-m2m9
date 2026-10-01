# -*- coding: utf-8 -*-
"""14_m2_cross_onek1k.py - M2 修正版：统一 OneK1K sc-eQTL 工具变量做跨病种 MR
之前 M2 失败原因：OneK1K permuted 无 rsid 且 GRCh38 vs GWAS GRCh37 坐标不匹配（cra 匹配误配邻近 SNP）。
现在 QTD .all 有 rsid 列 → 用 rsid 跨版本精确匹配 + 等位基因方向校正。
暴露：35 个显著基因-细胞对的 top cis-eQTL（与主分析完全相同的工具变量，官方 beta/se）
结局：LUAD（复核）+ IBD / BRCA / T2D / Asthma
"""
import gzip
import os
import pickle
import re

import numpy as np
import pandas as pd
from scipy import stats

BASE = r"C:/Users/86159/WorkBuddy/单细胞数据孟德尔随机化"
CACHE = os.path.join(BASE, ".workbuddy", "tmp", "qtd_cache")
ALL_DIR = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/onek1k/all"
SIG = os.path.join(BASE, "results", "scmr_luad_risk_sig.csv")
OUT = os.path.join(BASE, "results", "cross_disease_mr_onek1k.csv")
CROSS = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/cross_disease"

GWAS = [
    ("LUAD", os.path.join(BASE, "data", "luad_gwas", "GCST004744.h.tsv.gz"), "肺腺癌"),
    ("IBD", os.path.join(CROSS, "GCST004132.h.tsv.gz"), "炎症性肠病"),
    ("BRCA", os.path.join(CROSS, "GCST004988.h.tsv.gz"), "乳腺癌"),
    ("T2D", os.path.join(CROSS, "GCST006867.h.tsv.gz"), "2型糖尿病"),
    ("Asthma", os.path.join(CROSS, "GCST005038.h.tsv.gz"), "哮喘"),
]

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def parse_variant(v):
    m = re.match(r"chr(\w+)_(\d+)_(\w+)_(\w+)", str(v))
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3), m.group(4)


def allele_sign(ea, oa, ref, alt):
    if ea == alt: return 1
    if ea == ref: return -1
    if ea in COMP and COMP[ea] == alt: return 1
    if ea in COMP and COMP[ea] == ref: return -1
    return None


def load_gwas(path):
    keep = []
    for ch in pd.read_csv(path, sep="\t", compression="gzip",
                          usecols=["hm_rsid", "hm_beta", "standard_error", "p_value",
                                   "hm_effect_allele", "hm_other_allele"],
                          dtype={"hm_rsid": str}, chunksize=5_000_000):
        ch = ch.rename(columns={"hm_rsid": "rsid", "hm_beta": "gw_beta", "standard_error": "gw_se",
                                "p_value": "gw_p", "hm_effect_allele": "ea", "hm_other_allele": "oa"})
        for c in ["gw_beta", "gw_se", "gw_p"]:
            ch[c] = pd.to_numeric(ch[c], errors="coerce")
        ch = ch.dropna(subset=["gw_beta", "gw_se"]).drop_duplicates("rsid")
        keep.append(ch[["rsid", "gw_beta", "gw_se", "gw_p", "ea", "oa"]])
    g = pd.concat(keep, ignore_index=True).drop_duplicates("rsid", keep="first")
    return g.set_index("rsid")


def main():
    sig = pd.read_csv(SIG)
    print(f"[输入] {len(sig)} 显著对", flush=True)
    # 1) 从 QTD 缓存提取每对的 rsid + 官方 beta/se；未命中则从 QTD .all 全扫补全
    instr_rows = []
    missing = []
    for _, r in sig.iterrows():
        qtd, gene, var = r["qtd"], r["gene_symbol"], r["variant"]
        pv = parse_variant(var)
        if pv is None:
            continue
        chrom, pos, ref, alt = pv
        fp = os.path.join(CACHE, f"{qtd}_{gene}.pkl")
        hit_row = None
        if os.path.exists(fp):
            with open(fp, "rb") as f:
                df = pickle.load(f)
            h = df[(df["pos"] == pos) & (df["ref"] == ref) & (df["alt"] == alt)]
            if len(h) > 0:
                hit_row = h.iloc[0]
        if hit_row is None:
            missing.append({"qtd": qtd, "gene": gene, "variant": var,
                            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                            "cell_type": r["cell_type"]})
            continue
        instr_rows.append({"gene": gene, "cell_type": r["cell_type"], "qtd": qtd,
                           "variant": var, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                           "rsid": hit_row["rsid"], "eqtl_beta": hit_row["beta"], "eqtl_se": hit_row["se"],
                           "F_official": (hit_row["beta"] / hit_row["se"]) ** 2})
    # 补充提取缺失的工具（按 QTD 分组全扫，variant 精确匹配）
    if missing:
        from collections import defaultdict
        by_qtd = defaultdict(list)
        for m in missing:
            by_qtd[m["qtd"]].append(m)
        print(f"[补充] 全扫 {len(set(by_qtd.keys()))} 个 QTD 提取 {len(missing)} 对工具...", flush=True)
        extra_dir = os.path.join(BASE, ".workbuddy", "tmp", "qtd_cache_extra")
        os.makedirs(extra_dir, exist_ok=True)
        for qtd, rows in by_qtd.items():
            fp = os.path.join(ALL_DIR, f"{qtd}.all.tsv.gz")
            if not os.path.exists(fp):
                print(f"  QTD {qtd} 缺失", flush=True)
                continue
            targets = {r["variant"]: r for r in rows}
            found = {}
            with gzip.open(fp, "rt", errors="replace") as f:
                hdr = f.readline().rstrip("\n").split("\t")
                idx = {h: i for i, h in enumerate(hdr)}
                vi, gi = idx["variant"], idx["molecular_trait_id"]
                bi, si, rsi = idx["beta"], idx["se"], idx["rsid"]
                for line in f:
                    c = line.split("\t")
                    v = c[vi].strip()
                    if v in targets and v not in found:
                        found[v] = (c[rsi].strip(), float(c[bi]), float(c[si]))
            for r in rows:
                v = r["variant"]
                if v not in found:
                    print(f"  ! {r['gene']} {v}: QTD 无此行", flush=True)
                    continue
                rsid, ebeta, ese = found[v]
                instr_rows.append({"gene": r["gene"], "cell_type": r["cell_type"], "qtd": qtd,
                                   "variant": v, "chrom": r["chrom"], "pos": r["pos"],
                                   "ref": r["ref"], "alt": r["alt"], "rsid": rsid,
                                   "eqtl_beta": ebeta, "eqtl_se": ese,
                                   "F_official": (ebeta / ese) ** 2})
    instr = pd.DataFrame(instr_rows)
    print(f"[工具] 共 {len(instr)}/{len(sig)} 对", flush=True)
    # 2) 各病种 Wald MR
    all_res = []
    for dis, path, label in GWAS:
        g = load_gwas(path)
        print(f"[{dis}] GWAS rsid 条目 {len(g)}", flush=True)
        m = instr.merge(g, left_on="rsid", right_index=True, how="inner")
        m["sign"] = [allele_sign(str(e), str(o), ref, alt)
                     for e, o, ref, alt in zip(m["ea"], m["oa"], m["ref"], m["alt"])]
        m = m[m["sign"].notna()].copy()
        if len(m) == 0:
            print(f"  {dis}: 0 匹配", flush=True)
            continue
        m["beta_gw"] = m["gw_beta"] * m["sign"]
        m["b_mr"] = m["beta_gw"] / m["eqtl_beta"]
        m["se_mr"] = m["gw_se"] / m["eqtl_beta"].abs()
        m["or"] = np.exp(m["b_mr"])
        m["p_mr"] = 2 * (1 - stats.norm.cdf(m["b_mr"].abs() / m["se_mr"]))
        # 每病种 FDR
        m["fdr"] = stats.false_discovery_control(m["p_mr"].values)
        res = m[["gene", "cell_type", "qtd", "variant", "rsid", "eqtl_beta", "F_official",
                 "b_mr", "se_mr", "or", "p_mr", "fdr"]].copy()
        res.insert(1, "disease", dis)
        res.insert(2, "disease_label", label)
        all_res.append(res)
        n_sig = (m["fdr"] < 0.05).sum()
        print(f"  {dis}: {len(m)} 对匹配, FDR<0.05: {n_sig}", flush=True)
    out = pd.concat(all_res, ignore_index=True)
    out.to_csv(OUT, index=False)
    print(f"\n[完成] {len(out)} 行 -> results/cross_disease_mr_onek1k.csv", flush=True)


if __name__ == "__main__":
    main()
