# -*- coding: utf-8 -*-
"""M3 前置：从 GSE131907 raw UMI 矩阵中流式抽取基因面板，计算各细胞类型的表达统计。

设计考虑（本机仅 8 GB RAM）：
  - **不**把 29,634 × 208,506 的矩阵读入内存；只保留面板基因所在行
  - 逐行流式解析；仅对命中的行做字段切分（避免 61.8 亿次转换）
  - 输出「细胞类型 × 基因」的均值 counts 与阳性细胞比例

注意：raw UMI 未做文库大小归一化，跨细胞类型比较存在文库深度混杂；
      % 阳性细胞比例对该混杂相对稳健。正式的归一化版本另用 log2TPM 矩阵复核。
"""
import gzip
import json

import numpy as np
import pandas as pd

BASE = "."
MAT = f"{BASE}/rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
ANN = f"{BASE}/rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"

PANEL = [
    # 目标轴
    "JAML", "CXADR",
    # T 细胞
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "CTLA4", "PDCD1", "LAG3", "HAVCR2",
    # 细胞毒 / NK
    "NKG7", "GNLY", "KLRD1", "NCR1", "GZMB", "PRF1", "GZMA",
    # B / 浆细胞
    "MS4A1", "CD19", "CD79A", "MZB1",
    # 髓系
    "CD14", "CD68", "LYZ", "ITGAM", "ITGAX", "FCGR3A", "APOE", "C1QA", "MARCO",
    # 上皮 / 肿瘤
    "EPCAM", "KRT18", "KRT19", "NKX2-1", "NAPSA", "SFTPC", "TP63",
    # 基质 / 内皮
    "COL1A1", "COL1A2", "DCN", "PECAM1", "VWF",
    # 参照
    "PTPRC", "ACTB", "GAPDH", "MKI67",
    # 项目其他显著基因（有表达者）
    "FUBP1", "RNASET2", "GALK2", "RPS6KA2", "MAP4K4", "IREB2", "EID1", "SECISBP2L",
    "RAB31", "RAB4B", "CMIP", "STMN3", "NUMBL", "PARVA", "DNAJA4", "HLA-C",
]
PANEL = sorted(set(PANEL))
print(f"面板基因数: {len(PANEL)}")

# ---------- 注释 ----------
ann = pd.read_csv(ANN, sep="\t")
print("注释列:", list(ann.columns))
key = "Index" if "Index" in ann.columns else ann.columns[0]
ann[key] = ann[key].astype(str)
b2ct = dict(zip(ann[key], ann["Cell_type"].astype(str)))
b2sub = dict(zip(ann[key], ann["Cell_subtype"].fillna("NA").astype(str)))
b2origin = dict(zip(ann[key], ann["Sample_Origin"].astype(str)))
ct_levels = sorted(set(b2ct.values()))
sub_levels = sorted(set(b2sub.values()))
print(f"注释细胞数 {len(ann)}；Cell_type {len(ct_levels)} 类；Cell_subtype {len(sub_levels)} 类")

# ---------- 流式抽取 ----------
acc = {g: {} for g in PANEL}          # gene -> celltype -> [sum, n_pos]
acc_sub = {g: {} for g in PANEL}      # gene -> subtype  -> [sum, n_pos]
counts_by_ct = {}
counts_by_sub = {}
celltypes, subtypes = [], []

with gzip.open(MAT, "rt", errors="replace") as f:
    header = f.readline().rstrip("\n").split("\t")
    cell_ids = header[1:]
    print(f"细胞数（矩阵列）: {len(cell_ids)}；样例 {cell_ids[:2]}")
    miss = 0
    for cid in cell_ids:
        ct = b2ct.get(cid)
        if ct is None:
            miss += 1
            ct = "UNMATCHED"
        sc = b2sub.get(cid, "NA")
        celltypes.append(ct)
        subtypes.append(sc)
        counts_by_ct[ct] = counts_by_ct.get(ct, 0) + 1
        counts_by_sub[sc] = counts_by_sub.get(sc, 0) + 1
    print(f"未匹配到注释的细胞: {miss}")
    celltypes = np.array(celltypes)
    subtypes = np.array(subtypes)
    # 预先算好分组掩码，避免在基因循环里反复做 20 万次比较
    ct_masks = {ct: (celltypes == ct) for ct in np.unique(celltypes)}
    sub_masks = {sc: (subtypes == sc) for sc in np.unique(subtypes)}

    hits = 0
    for line in f:
        p = line.find("\t")
        name = line[:p]
        if name not in PANEL:
            continue
        hits += 1
        vals = np.fromstring(line[p + 1:], dtype=np.float32, sep="\t")
        for ct, m in ct_masks.items():
            v = vals[m]
            acc[name][ct] = [float(v.sum()), int((v > 0).sum()), int(m.sum())]
        for sc, m in sub_masks.items():
            v = vals[m]
            acc_sub[name][sc] = [float(v.sum()), int((v > 0).sum()), int(m.sum())]
        if hits % 10 == 0:
            print(f"  已抽取 {hits} 个面板基因...", flush=True)

print(f"\n命中面板基因 {hits} / {len(PANEL)}；未命中: {sorted(set(PANEL) - set(acc) | {g for g in PANEL if not acc[g]})}")

# ---------- 整理输出 ----------
rows = []
for g in PANEL:
    for ct, (s, npos, n) in sorted(acc[g].items()):
        rows.append({"gene": g, "level": "Cell_type", "group": ct,
                     "n_cells": n, "mean_counts": s / n, "pct_positive": 100 * npos / n})
    for sc, (s, npos, n) in sorted(acc_sub[g].items()):
        rows.append({"gene": g, "level": "Cell_subtype", "group": sc,
                     "n_cells": n, "mean_counts": s / n, "pct_positive": 100 * npos / n})
res = pd.DataFrame(rows)
res.to_csv("results/m3_gse131907_panel_expression.csv", index=False)

json.dump({"cell_type_counts": counts_by_ct, "cell_subtype_counts": counts_by_sub,
           "panel": PANEL, "n_cells": len(cell_ids), "unmatched": miss},
          open("results/m3_gse131907_panel_meta.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)

print("\n" + "=" * 96)
print("JAML 与 CXADR 在各细胞类型的表达（raw UMI counts）")
print("=" * 96)
for g in ["JAML", "CXADR", "CD3D", "EPCAM", "LYZ", "PECAM1"]:
    sub = res[(res["gene"] == g) & (res["level"] == "Cell_type")].sort_values("mean_counts", ascending=False)
    print(f"\n[{g}]")
    print(sub[["group", "n_cells", "mean_counts", "pct_positive"]].to_string(index=False))

print("\n" + "=" * 96)
print("JAML / CXADR 在 CD4 相关 Cell_subtype 中的表达")
print("=" * 96)
cd4 = res[(res["gene"].isin(["JAML", "CXADR"])) & (res["level"] == "Cell_subtype") &
          (res["group"].str.contains("CD4|Treg", case=False, na=False))]
print(cd4.sort_values(["gene", "mean_counts"], ascending=[True, False]).to_string(index=False))

print("\n已写出: results/m3_gse131907_panel_expression.csv")
