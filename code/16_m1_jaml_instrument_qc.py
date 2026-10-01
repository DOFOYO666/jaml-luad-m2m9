"""M1 步骤 1：从 OneK1K 官方 significant cis-QTL 汇总中提取 JAML(AMICA1) 的工具变量 QC。

数据源：data/onek1k/extracted/OneK1K_TensorQTL_top_eQTL_summary/
        OneK1K_{celltype}.sig_cis_qtl_pairs.chr11.csv  （官方 TensorQTL 输出）
注意：OneK1K 原始发布（Yazar 2022）使用旧基因符号 **AMICA1** 指代 JAML。
字段：phenotype_id, variant_id, af(=alt 频率), ma_count, ma_samples,
      slope, slope_se, pval_nominal, pval_perm, pval_beta, qval, tss_distance
F 统计量 = (slope / slope_se)^2
"""
import glob
import os
import re

import numpy as np
import pandas as pd

SRC = "data/onek1k/extracted/OneK1K_TensorQTL_top_eQTL_summary"
OUT = "results/m1_jaml_instrument_qc.csv"

# 细胞类型归并（与项目 QTD_MAP 口径一致）
CT_CLASS = {
    "B_IN": "B cell", "B_MEM": "B cell",
    "CD4_ET": "CD4 T", "CD4_NC": "CD4 T", "CD4_SOX4": "CD4 T",
    "CD8_ET": "CD8 T", "CD8_NC": "CD8 T", "CD8_S100B": "CD8 T",
    "NK": "NK", "NK_R": "NK",
    "Mono_C": "Monocyte", "Mono_NC": "Monocyte",
    "DC": "DC", "Plasma": "Plasma",
}

rows = []
files = sorted(glob.glob(os.path.join(SRC, "*chr11.csv")))
print(f"扫描 {len(files)} 个细胞类型文件 ...")

for f in files:
    ct_file = re.sub(r"^OneK1K_", "", os.path.basename(f)).replace(".sig_cis_qtl_pairs.chr11.csv", "")
    d = pd.read_csv(f, sep="\t")
    pid = d["phenotype_id"].astype(str).str.upper()
    hit = d[pid.isin(["AMICA1", "JAML"])].copy()
    if hit.empty:
        continue
    hit.insert(0, "celltype_file", ct_file)
    hit.insert(1, "cell_class", CT_CLASS.get(ct_file, "?"))
    rows.append(hit)
    print(f"  {ct_file:10s} ({CT_CLASS.get(ct_file,'?'):9s}): {len(hit)} 行")

allhit = pd.concat(rows, ignore_index=True)
allhit["F"] = (allhit["slope"] / allhit["slope_se"]) ** 2
allhit["MAF"] = np.minimum(allhit["af"], 1 - allhit["af"])
allhit = allhit.sort_values("F", ascending=False).reset_index(drop=True)

cols = ["celltype_file", "cell_class", "phenotype_id", "variant_id", "tss_distance",
        "af", "MAF", "ma_count", "ma_samples", "slope", "slope_se", "F",
        "pval_nominal", "pval_perm", "pval_beta", "qval"]
cols = [c for c in cols if c in allhit.columns]
out = allhit[cols]
out.to_csv(OUT, index=False)

pd.set_option("display.width", 200)
print("\n" + "=" * 100)
print("JAML (AMICA1) 各细胞类型的显著 cis-QTL —— 官方 TensorQTL，按 F 降序")
print("=" * 100)
print(out.to_string(index=False))

print("\n--- 汇总 ---")
print(f"细胞类型数: {out['celltype_file'].nunique()}（全部 14 个细胞类型均有显著 cis-QTL）")
print(f"F 范围: {out['F'].min():.1f} – {out['F'].max():.1f}")
print(f"F > 10 的细胞类型数: {(out['F'] > 10).sum()}")
print(f"MAF 范围: {out['MAF'].min():.4f} – {out['MAF'].max():.4f}")
print(f"MAF < 0.05 的条数: {(out['MAF'] < 0.05).sum()}")
print(f"\n已写出: {OUT}")
