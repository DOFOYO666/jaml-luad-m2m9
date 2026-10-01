# -*- coding: utf-8 -*-
"""10b_lung_subtypes_mr.py - 肺癌组织学亚型扩展 MR
用与 LUAD 相同的 21 基因工具变量，在 SqCC（肺鳞癌）与 SCLC（小细胞肺癌）GWAS 中做 Wald ratio MR，
输出 基因 × 亚型 的 OR/P/FDR 对比矩阵，检验 LUAD 结论能否扩展至其他组织学亚型。
用法: python 10b_lung_subtypes_mr.py
"""
import gzip
import os

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
SUB = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/lung_subtypes"
LUAD = os.path.join(BASE, "data", "luad_gwas", "GCST004744.h.tsv.gz")

DISEASES = [
    ("LUAD", LUAD, "肺腺癌", "GCST004744", 11273),
    ("SqCC", os.path.join(SUB, "SqCC.h.tsv.gz"), "肺鳞癌", "GCST004750", 7426),
    ("SCLC", os.path.join(SUB, "SCLC.h.tsv.gz"), "小细胞肺癌", "GCST004746", 2664),
    ("Overall", os.path.join(SUB, "Overall.h.tsv.gz"), "肺癌总体", "GCST004748", 29266),
]


def load_instruments():
    """工具变量：21 基因最强 cis-eQTL SNP（eQTLGen 16 + GTEx Lung 5）"""
    v = pd.read_csv(os.path.join(RES, "bulk_eqtlgen_validation.csv"))
    eq = v[["GeneSymbol", "SNP", "AssessedAllele", "OtherAllele", "beta", "se", "p"]].drop_duplicates("GeneSymbol")
    eq = eq.rename(columns={"GeneSymbol": "gene", "SNP": "rsid", "AssessedAllele": "ea", "OtherAllele": "oa",
                            "beta": "beta_e", "se": "se_e", "p": "p_e"})
    gt = pd.read_pickle(os.path.join(RES, "gtex_lung_regions.pkl"))
    gt = gt[gt["rsid"] != "NA"].copy()
    gt["z"] = (gt["beta"] / gt["se"]).abs()
    gtop = gt.sort_values("z", ascending=False).drop_duplicates("gene")
    have = set(eq["gene"])
    gtop = gtop[~gtop["gene"].isin(have)]
    gtop = gtop.rename(columns={"rsid": "rsid", "beta": "beta_e", "se": "se_e"})
    gtop["ea"] = gtop["alt"]
    gtop["oa"] = gtop["ref"]
    gtop["p_e"] = np.nan
    inst = pd.concat([eq[["gene", "rsid", "ea", "oa", "beta_e", "se_e", "p_e"]],
                      gtop[["gene", "rsid", "ea", "oa", "beta_e", "se_e", "p_e"]]], ignore_index=True)
    return inst.dropna(subset=["rsid", "beta_e", "se_e"])


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
        print("  [警告] 文件不完整(EOF)", flush=True)
        return pd.DataFrame(columns=["rsid", "ea", "oa", "beta", "se", "p"])
    if not keep:
        return pd.DataFrame(columns=["rsid", "ea", "oa", "beta", "se", "p"])
    g = pd.concat(keep, ignore_index=True)
    for c in ["beta", "se", "p"]:
        g[c] = pd.to_numeric(g[c], errors="coerce")
    g["ea"] = g["ea"].astype(str).str.strip()
    g["oa"] = g["oa"].astype(str).str.strip()
    return g.dropna(subset=["beta", "se"])


def main():
    inst = load_instruments()
    print(f"[工具变量] {len(inst)} 基因")
    rsids = set(inst["rsid"])
    rows = []
    for code, path, label, gwas_id, n_case in DISEASES:
        if not os.path.exists(path):
            print(f"[{code}] 文件不存在，跳过: {path}")
            continue
        print(f"[{code}] 提取 {label} (n_case={n_case}) ...", flush=True)
        g = extract_gwas(path, rsids)
        print(f"  命中 {len(g)} 行")
        for _, iv in inst.iterrows():
            hit = g[g["rsid"] == iv["rsid"]]
            if hit.empty:
                continue
            h = hit.iloc[0]
            b_out = h["beta"] if iv["ea"] == h["ea"] else -h["beta"]
            b_mr = b_out / iv["beta_e"]
            se_mr = h["se"] / abs(iv["beta_e"])
            p_mr = 2 * (1 - norm.cdf(abs(b_mr) / se_mr))
            rows.append({"gene": iv["gene"], "subtype": code, "subtype_label": label,
                         "gwas_id": gwas_id, "n_case": n_case, "rsid": iv["rsid"],
                         "beta_mr": b_mr, "se_mr": se_mr, "or": np.exp(b_mr), "p_mr": p_mr})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "lung_subtypes_mr.csv"), index=False)

    # 矩阵视图（仅含≥1次显著或全部有值的基因）
    pivot_or = out.pivot_table(index="gene", columns="subtype", values="or")
    pivot_p = out.pivot_table(index="gene", columns="subtype", values="p_mr")
    # 3 亚型 × 21 基因 = 63 检验，BH-FDR
    fdr_flat = multipletests(out["p_mr"], method="fdr_bh")[1]
    out["fdr"] = fdr_flat

    print("\n=== OR 矩阵（对比） ===")
    print(pivot_or.round(2).to_string())
    print("\n=== P 矩阵（*: P<0.05, †: FDR<0.05） ===")
    sig = out.set_index(["gene", "subtype"])
    disp = pd.DataFrame(index=pivot_p.index, columns=pivot_p.columns, dtype=object)
    for gene in pivot_p.index:
        for sub in pivot_p.columns:
            p = pivot_p.loc[gene, sub]
            if pd.isna(p):
                disp.loc[gene, sub] = "NA"
                continue
            fdr = sig.loc[(gene, sub), "fdr"]
            star = "†" if fdr < 0.05 else ("*" if p < 0.05 else "")
            disp.loc[gene, sub] = f"{p:.2g}{star}"
    print(disp.to_string())

    # 汇总：每个基因在哪些亚型显著
    print("\n=== 各基因亚型显著汇总（P<0.05） ===")
    sig_genes = out[out["p_mr"] < 0.05]
    for gene in sorted(out["gene"].unique()):
        sub = sig_genes[sig_genes["gene"] == gene]
        if len(sub):
            detail = "; ".join(f"{r['subtype']} OR={r['or']:.2f} P={r['p_mr']:.2g}"
                               for _, r in sub.iterrows())
            print(f"  {gene}: {detail}")
        else:
            print(f"  {gene}: (三个亚型均不显著)")

    pivot_or.to_csv(os.path.join(RES, "lung_subtypes_or_matrix.csv"))
    pivot_p.to_csv(os.path.join(RES, "lung_subtypes_p_matrix.csv"))
    out.to_csv(os.path.join(RES, "lung_subtypes_mr.csv"), index=False)
    n_sig = (out["p_mr"] < 0.05).sum()
    n_fdr = (out["fdr"] < 0.05).sum()
    print(f"\n[完成] {len(out)} 条记录；P<0.05: {n_sig}；FDR<0.05: {n_fdr}")


if __name__ == "__main__":
    main()
