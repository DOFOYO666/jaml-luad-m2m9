# -*- coding: utf-8 -*-
"""08_sensitivity.py - 依据研究局限性补充的敏感性分析
数据基础：主分析结果 results/scmr_luad_risk_sig.csv（已含 MR 估计、eQTL beta/se、F 值）
1) Steiger 方向性检验（z_exposure vs z_outcome，尺度无关）
2) F 阈值敏感性（F>20 后显著基因）
3) 排除 MHC(chr6:28-34Mb) 与 15q25(chr15:78-79.5Mb) 后显著基因
4) GTEx Lung 多 SNP IVW MR（缓解"单工具变量"局限；未 LD 剪枝，仅方向一致性参考）
"""
import gzip
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
GWAS = os.path.join(BASE, "data", "luad_gwas", "GCST004744.h.tsv.gz")
MHC = ("6", 28_000_000, 34_000_000)
P15Q25 = ("15", 78_000_000, 79_500_000)


def parse_variant(v):
    s = str(v).replace("chr", "")
    parts = s.split("_")
    if len(parts) < 4:
        return None
    return parts[0], parts[1], parts[2], parts[3]


def main():
    sig = pd.read_csv(os.path.join(RES, "scmr_luad_risk_sig.csv"))
    sig["chrom"] = sig["variant"].apply(lambda v: (lambda p: p[0] if p else None)(parse_variant(v)))
    sig["pos"] = sig["variant"].apply(lambda v: (lambda p: int(p[1]) if p else None)(parse_variant(v)))
    sig = sig.dropna(subset=["chrom", "pos", "beta", "se", "b_mr", "se_mr"])
    print(f"[输入] 主分析显著对 {len(sig)}，去重基因 {sig['gene_symbol'].nunique()}")

    # ---- 1) Steiger 方向性检验 ----
    print("\n[1] Steiger 方向性检验（eQTL→GWAS 因果方向）")
    z_e = (sig["beta"] / sig["se"]).abs()
    z_g = (sig["b_mr"] / sig["se_mr"]).abs()
    support = (z_e > z_g).mean()
    print(f"  |z_expo|>|z_out| 比例: {support:.1%} ({(z_e > z_g).sum()}/{len(sig)})")
    print(f"  解释: 比例越高 → 工具变量主要解释暴露(eQTL)而非结局 → 支持 eQTL→LUAD 方向")
    # 反向（GWAS→eQTL）的 SNP 数
    rev = sig[z_g > z_e]
    print(f"  方向可疑（|z_out|>|z_expo|）: {len(rev)} 对: {sorted(rev['gene_symbol'].unique())[:15]}")

    # ---- 2) F 阈值敏感性 ----
    print("\n[2] F 阈值敏感性（F>20）")
    sig20 = sig[sig["F"] > 20]
    print(f"  F>10 显著基因: {sig['gene_symbol'].nunique()} | F>20 显著基因: {sig20['gene_symbol'].nunique()}")
    print(f"  F>20 丢失: {sorted(set(sig['gene_symbol']) - set(sig20['gene_symbol']))}")
    print(f"  F>20 保留: {sorted(sig20['gene_symbol'].unique())}")

    # ---- 3) 排除 MHC / 15q25 ----
    print("\n[3] 排除 MHC 与 15q25（易受复杂连锁/已知位点影响）")
    mask_mhc = (sig["chrom"] == MHC[0]) & (sig["pos"] >= MHC[1]) & (sig["pos"] <= MHC[2])
    mask_15q = (sig["chrom"] == P15Q25[0]) & (sig["pos"] >= P15Q25[1]) & (sig["pos"] <= P15Q25[2])
    print(f"  MHC 区 SNP: {mask_mhc.sum()} (基因: {sorted(sig.loc[mask_mhc,'gene_symbol'].unique())})")
    print(f"  15q25 SNP: {mask_15q.sum()} (基因: {sorted(sig.loc[mask_15q,'gene_symbol'].unique())})")
    keep = sig[~mask_mhc & ~mask_15q]
    print(f"  排除后显著基因: {keep['gene_symbol'].nunique()}: {sorted(keep['gene_symbol'].unique())}")

    # ---- 4) GTEx Lung 多 SNP IVW（方向一致性参考） ----
    print("\n[4] GTEx Lung 多 SNP IVW（未 LD 剪枝，方向一致性参考）")
    gt_path = os.path.join(RES, "gtex_lung_regions.pkl")
    if not os.path.exists(gt_path):
        print("  无 GTEx 区域缓存，跳过")
    else:
        gt = pd.read_pickle(gt_path)
        gt["rsid"] = gt["rsid"].astype(str)
        # 分批读 GWAS 取 rsid 对应 beta/se
        need_rs = set(gt["rsid"][gt["rsid"] != "NA"])
        g = pd.read_csv(GWAS, sep="\t", compression="gzip",
                        usecols=["hm_rsid", "hm_beta", "standard_error"],
                        dtype={"hm_rsid": str})
        g = g.rename(columns={"hm_rsid": "rsid", "hm_beta": "bg", "standard_error": "sg"})
        g["bg"] = pd.to_numeric(g["bg"], errors="coerce")
        g["sg"] = pd.to_numeric(g["sg"], errors="coerce")
        g = g[g["rsid"].isin(need_rs)].dropna(subset=["bg", "sg"])
        out = []
        for gene, sub in gt.groupby("gene"):
            wd = sub.merge(g, on="rsid", suffixes=("_e", ""))
            wd = wd.dropna(subset=["beta", "se", "bg", "sg"])
            if len(wd) < 3:
                continue
            b = wd["bg"] / wd["beta"]
            se = wd["sg"] / wd["beta"].abs()
            w = 1 / se ** 2
            ivw_b = np.sum(w * b) / np.sum(w)
            ivw_se = np.sqrt(1 / np.sum(w))
            p = 2 * (1 - norm.cdf(abs(ivw_b) / ivw_se))
            out.append({"gene": gene, "n_snp": len(wd), "ivw_beta": ivw_b, "ivw_se": ivw_se, "ivw_p": p})
        ivw = pd.DataFrame(out).sort_values("ivw_p")
        print(f"  IVW 完成（{len(ivw)} 基因），top15（按 p）：")
        print(ivw.head(15).to_string(index=False))
        ivw.to_csv(os.path.join(RES, "gtex_ivw_sensitivity.csv"), index=False)

    summary = {
        "steiger_direction_support": float(support),
        "sig_pairs": int(len(sig)),
        "sig_genes": sorted(sig["gene_symbol"].unique()),
        "sig_genes_F20": sorted(sig20["gene_symbol"].unique()),
        "sig_genes_noMHC_no15q25": sorted(keep["gene_symbol"].unique()),
    }
    pd.Series(summary).to_json(os.path.join(RES, "sensitivity_summary.json"), force_ascii=False, indent=1)
    print("\n[完成] results/sensitivity_summary.json + gtex_ivw_sensitivity.csv")


if __name__ == "__main__":
    main()
