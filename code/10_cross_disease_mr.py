# -*- coding: utf-8 -*-
"""10_cross_disease_mr.py - 跨病种 MR（多效性筛查）
对 21 个显著基因的最强 cis-eQTL 工具变量，在 6 个病种（LUAD + IBD + 乳腺癌 + RA + T2D + 哮喘）GWAS
中做 Wald ratio MR，输出 gene × disease 效应矩阵（OR 与 P），用于多效性图谱。
用法: python 10_cross_disease_mr.py
"""
import gzip
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
CROSS = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/cross_disease"
LUAD = os.path.join(BASE, "data", "luad_gwas", "GCST004744.h.tsv.gz")

DISEASES = [
    ("LUAD", LUAD, "肺腺癌", "GCST004744"),
    ("IBD", os.path.join(CROSS, "GCST004132.h.tsv.gz"), "炎症性肠病", "GCST004132"),
    ("BRCA", os.path.join(CROSS, "GCST004988.h.tsv.gz"), "乳腺癌", "GCST004988"),
    ("RA", os.path.join(CROSS, "GCST005543.h.tsv.gz"), "类风湿关节炎", "GCST005543"),
    ("T2D", os.path.join(CROSS, "GCST006867.h.tsv.gz"), "2型糖尿病", "GCST006867"),
    ("Asthma", os.path.join(CROSS, "GCST005038.h.tsv.gz"), "哮喘", "GCST005038"),
]


def load_instruments():
    """工具变量：21 基因最强 cis-eQTL SNP（eQTLGen 16 + GTEx Lung 5）"""
    v = pd.read_csv(os.path.join(RES, "bulk_eqtlgen_validation.csv"))
    # eQTLGen: 每基因一个 SNP（AssessedAllele 为效应等位基因）
    eq = v[["GeneSymbol", "SNP", "AssessedAllele", "OtherAllele", "beta", "se", "p"]].drop_duplicates("GeneSymbol")
    eq = eq.rename(columns={"GeneSymbol": "gene", "SNP": "rsid", "AssessedAllele": "ea", "OtherAllele": "oa",
                            "beta": "beta_e", "se": "se_e", "p": "p_e"})
    # GTEx 补充
    gt = pd.read_pickle(os.path.join(RES, "gtex_lung_regions.pkl"))
    gt = gt[gt["rsid"] != "NA"].copy()
    gt["z"] = (gt["beta"] / gt["se"]).abs()
    gtop = gt.sort_values("z", ascending=False).drop_duplicates("gene")
    # 只取 eQTLGen 缺失的基因
    have = set(eq["gene"])
    gtop = gtop[~gtop["gene"].isin(have)]
    # GTEx 的 ref/alt: GTEx beta 对应 alt allele
    gtop = gtop.rename(columns={"rsid": "rsid", "beta": "beta_e", "se": "se_e"})
    gtop["ea"] = gtop["alt"]
    gtop["oa"] = gtop["ref"]
    gtop["p_e"] = np.nan
    inst = pd.concat([eq[["gene", "rsid", "ea", "oa", "beta_e", "se_e", "p_e"]],
                      gtop[["gene", "rsid", "ea", "oa", "beta_e", "se_e", "p_e"]]], ignore_index=True)
    inst = inst.dropna(subset=["rsid", "beta_e", "se_e"])
    return inst


def extract_gwas(path, rsids):
    """分块读取 GWAS（harmonised），只保留工具变量 rsid 的行"""
    cols = ["hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "p_value"]
    keep = []
    try:
        for ch in pd.read_csv(path, sep="\t", compression="gzip", usecols=cols, dtype={"hm_rsid": str},
                              chunksize=2_000_000):
            ch = ch.rename(columns={"hm_rsid": "rsid", "hm_effect_allele": "ea", "hm_other_allele": "oa",
                                    "hm_beta": "beta", "standard_error": "se", "p_value": "p"})
            ch = ch[ch["rsid"].isin(rsids)]
            if len(ch):
                keep.append(ch)
    except EOFError:
        print("  [警告] 文件不完整(EOF)，请检查下载", flush=True)
        return pd.DataFrame(columns=["rsid", "ea", "oa", "beta", "se", "p"])
    if not keep:
        return pd.DataFrame(columns=["rsid", "ea", "oa", "beta", "se", "p"])
    g = pd.concat(keep, ignore_index=True)
    g["beta"] = pd.to_numeric(g["beta"], errors="coerce")
    g["se"] = pd.to_numeric(g["se"], errors="coerce")
    g["p"] = pd.to_numeric(g["p"], errors="coerce")
    g["ea"] = g["ea"].astype(str).str.strip()
    g["oa"] = g["oa"].astype(str).str.strip()
    return g.dropna(subset=["beta", "se"])


def main():
    inst = load_instruments()
    print(f"[工具变量] {len(inst)} 基因: {sorted(inst['gene'])}")
    rsids = set(inst["rsid"])
    rows = []
    for code, path, label, gwas_id in DISEASES:
        if not os.path.exists(path):
            print(f"[{code}] 文件不存在，跳过: {path}")
            continue
        print(f"[{code}] 提取 {label} ...", flush=True)
        g = extract_gwas(path, rsids)
        print(f"  命中 {len(g)} 行")
        for _, iv in inst.iterrows():
            hit = g[g["rsid"] == iv["rsid"]]
            if hit.empty:
                continue
            h = hit.iloc[0]
            # 方向校正：eQTL 效应等位基因(ea) vs GWAS effect allele(h.ea)
            b_out = h["beta"] if iv["ea"] == h["ea"] else -h["beta"]
            b_mr = b_out / iv["beta_e"]
            se_mr = h["se"] / abs(iv["beta_e"])
            p_mr = 2 * (1 - norm.cdf(abs(b_mr) / se_mr))
            or_ = np.exp(b_mr)
            rows.append({"gene": iv["gene"], "disease": code, "disease_label": label,
                         "rsid": iv["rsid"], "beta_mr": b_mr, "se_mr": se_mr,
                         "or": or_, "p_mr": p_mr, "n_case_gwas": np.nan})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "cross_disease_mr.csv"), index=False)
    # 矩阵视图
    pivot_or = out.pivot_table(index="gene", columns="disease", values="or")
    pivot_p = out.pivot_table(index="gene", columns="disease", values="p_mr")
    print("\n[OR 矩阵]")
    print(pivot_or.round(2).to_string())
    print("\n[P 矩阵]（<0.05 标记*）")
    print(pivot_p.applymap(lambda p: f"{p:.2g}" + ("*" if p < 0.05 else "")).to_string())
    pivot_or.to_csv(os.path.join(RES, "cross_disease_or_matrix.csv"))
    pivot_p.to_csv(os.path.join(RES, "cross_disease_p_matrix.csv"))
    # 显著计数
    n_sig = (out["p_mr"] < 0.05).sum()
    print(f"\n[完成] {len(out)} 条记录，P<0.05 显著 {n_sig} 条")


if __name__ == "__main__":
    main()
