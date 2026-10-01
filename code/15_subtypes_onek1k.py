# -*- coding: utf-8 -*-
"""15_subtypes_onek1k.py - 肺癌亚型 MR 统一 OneK1K 工具（M2 延伸）
18 个强工具基因（官方 F>10）× LUAD/SqCC/SCLC/Overall 4 结局，rsID 跨版本匹配 + 方向校正
"""
import os
import pandas as pd
import numpy as np
from scipy import stats

BASE = r"C:/Users/86159/WorkBuddy/单细胞数据孟德尔随机化"
RES = os.path.join(BASE, "results")
SUB = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/lung_subtypes"
LUAD_GWAS = os.path.join(BASE, "data", "luad_gwas", "GCST004744.h.tsv.gz")
INSTR = os.path.join(RES, "cross_disease_mr_onek1k_f10.csv")  # 18 强工具基因（含 rsid/ref/alt/eqtl_beta）

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
OUTCOMES = [
    ("LUAD", LUAD_GWAS, "肺腺癌", 11273),
    ("SqCC", os.path.join(SUB, "SqCC.h.tsv.gz"), "肺鳞癌", 7426),
    ("SCLC", os.path.join(SUB, "SCLC.h.tsv.gz"), "小细胞肺癌", 2664),
    ("Overall", os.path.join(SUB, "Overall.h.tsv.gz"), "肺癌总体", 29266),
]


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
    full = pd.read_csv(os.path.join(RES, "cross_disease_mr_onek1k.csv"))
    luad = full[full["disease"] == "LUAD"].sort_values("F_official", ascending=False).drop_duplicates("gene")
    instr = luad[luad["F_official"] > 10].copy()
    # 从 variant (chr11_118234507_C_T) 解析 ref/alt
    parts = instr["variant"].str.replace("chr", "").str.split("_", expand=True)
    instr["ref"] = parts[2].values
    instr["alt"] = parts[3].values
    print(f"[工具] {len(instr)} 个强工具基因（官方 F>10）", flush=True)
    all_res = []
    for dis, path, label, n_case in OUTCOMES:
        g = load_gwas(path)
        print(f"[{dis}] GWAS rsid {len(g)}", flush=True)
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
        m["fdr"] = stats.false_discovery_control(m["p_mr"].values)
        res = m[["gene", "cell_type", "variant", "rsid", "or", "p_mr", "fdr"]].copy()
        res.insert(1, "subtype", dis)
        res.insert(2, "subtype_label", label)
        res.insert(3, "n_case", n_case)
        all_res.append(res)
        print(f"  {dis}: {len(m)} 基因, FDR<0.05: {(m['fdr']<0.05).sum()}", flush=True)
    out = pd.concat(all_res, ignore_index=True)
    out.to_csv(os.path.join(RES, "lung_subtypes_mr_onek1k.csv"), index=False)
    print(f"\n[完成] {len(out)} 行 -> results/lung_subtypes_mr_onek1k.csv", flush=True)


if __name__ == "__main__":
    main()
