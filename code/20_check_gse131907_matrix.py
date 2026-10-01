# -*- coding: utf-8 -*-
"""校验 GSE131907 raw UMI 矩阵：维度、基因符号、以及 JAML(AMICA1) 的表达情况。

矩阵为 基因 × 细胞 的制表符文本（gz）。逐行流式读取，只保留首个字段以控制内存。
"""
import gzip
import json
import sys

import numpy as np

PATH = "rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
TARGETS = {"AMICA1", "JAML", "JAM3", "CD3D", "CD4", "CXADR", "PTPRC"}

n_genes = 0
target_rows = {}
header = None
first_genes = []

with gzip.open(PATH, "rt", errors="replace") as f:
    header = f.readline().rstrip("\n").split("\t")
    n_cells = len(header) - 1   # 首列为行名（可能为空或 "gene"）
    print(f"表头字段数 {len(header)}；首个字段 = {header[0]!r}；推断细胞数 = {n_cells}")
    print(f"细胞 ID 样例: {header[1:4]}")

    for line in f:
        p = line.rstrip("\n").find("\t")
        name = line[:p] if p > 0 else line.strip()
        n_genes += 1
        if n_genes <= 5:
            first_genes.append(name)
        if name in TARGETS:
            vals = line[p + 1:].rstrip("\n").split("\t")
            arr = np.array(vals, dtype=np.float32)
            nz = int((arr > 0).sum())
            target_rows[name] = {
                "row_index": n_genes - 1,
                "n_cells": len(arr),
                "n_nonzero": nz,
                "pct_positive": round(100 * nz / len(arr), 2),
                "mean_counts": round(float(arr.mean()), 4),
                "total_counts": int(arr.sum()),
                "max_counts": int(arr.max()),
            }
            print(f"  命中 {name}: 阳性细胞 {nz}/{len(arr)} ({100*nz/len(arr):.2f}%) "
                  f"总 counts={int(arr.sum())} 均值={arr.mean():.4f}")

print(f"\n矩阵维度: {n_genes} 基因 × {n_cells} 细胞")
print(f"前 5 个基因: {first_genes}")
print(f"目标基因命中: {sorted(target_rows)}")

out = {
    "path": PATH,
    "n_genes": n_genes,
    "n_cells": n_cells,
    "header_first_field": header[0],
    "first_genes": first_genes,
    "targets": target_rows,
}
with open("results/m3_gse131907_matrix_check.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
print("\n已写出 results/m3_gse131907_matrix_check.json")
