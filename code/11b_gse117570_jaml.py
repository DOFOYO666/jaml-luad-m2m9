# -*- coding: utf-8 -*-
"""11b_gse117570_jaml.py - GSE117570 (NSCLC 免疫 scRNA, Song 2021) JAML 表达验证
每样本 10x UMI 矩阵（无官方注释）→ 免疫 marker 基因推断细胞类型 → 计算 JAML(AMICA1)
在肿瘤 vs 正常、各免疫亚群中的表达（均值 CPM + 阳性率）。
用法: python 11b_gse117570_jaml.py
"""
import gzip
import os
from collections import defaultdict

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/gse117570"
OUT = os.path.join(BASE, "results")

# 关注的基因（symbol, 类别）
MARKERS = {
    "CD3D": "T", "CD3E": "T", "CD4": "CD4T", "CD8A": "CD8T", "CD8B": "CD8T",
    "TRDC": "gdT", "TRGC1": "gdT", "TRGC2": "gdT", "FOXP3": "Treg",
    "NKG7": "NK", "GNLY": "NK", "KLRD1": "NK", "NCAM1": "NK",
    "CD14": "MonoMac", "LYZ": "MonoMac", "CSF1R": "MonoMac", "FCGR3A": "MonoMac",
    "MS4A1": "B", "CD79A": "B", "CD79B": "B",
    "CLEC9A": "DC", "CD1C": "DC", "LILRA4": "DC",
    "AMICA1": "JAML", "JAML": "JAML",
}
# 类别优先级（用于打分）
PRIORITY = ["gdT", "NK", "CD8T", "Treg", "CD4T", "T", "MonoMac", "DC", "B"]

SAMPLES = [
    ("GSM3304007_P1_Tumor_processed_data.txt.gz", "P1", "Tumor", 2781642),
    ("GSM3304008_P1_Normal_processed_data.txt.gz", "P1", "Normal", 3122858),
    ("GSM3304009_P2_Tumor_processed_data.txt.gz", "P2", "Tumor", 2910994),
    ("GSM3304010_P2_Normal_processed_data.txt.gz", "P2", "Normal", 3969768),
    ("GSM3304011_P3_Tumor_processed_data.txt.gz", "P3", "Tumor", 545417),
    ("GSM3304012_P3_Normal_processed_data.txt.gz", "P3", "Normal", 762105),
    ("GSM3304013_P4_Tumor_processed_data.txt.gz", "P4", "Tumor", 3043332),
    ("GSM3304014_P4_Normal_processed_data.txt.gz", "P4", "Normal", 2465794),
]


def read_matrix(path, genes_of_interest):
    """读取 UMI 矩阵，返回 {gene: np.array(per cell)} + 每细胞总 UMI + barcodes"""
    with gzip.open(path, "rt", errors="replace") as f:
        header = f.readline()
        cells = header.rstrip("\n").replace('"', "").split("\t")
        # 去掉纯引号空值（文件末尾可能有）
        cells = [c for c in cells if c]
        n_cells = len(cells)
        genes = {g: np.zeros(n_cells, dtype=float) for g in genes_of_interest}
        tot = np.zeros(n_cells, dtype=float)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            gname = parts[0].replace('"', "")
            vals = parts[1:]
            if gname in genes:
                for i, v in enumerate(vals[:n_cells]):
                    if v:
                        x = float(v)
                        genes[gname][i] = x
                        tot[i] += x
            # 全矩阵行数过多，只统计总 UMI 用抽样？这里逐行累加太慢。
            # 改用只累加感兴趣基因 + 每行总和（成本高）。优化：总 UMI 仅需 per-cell，用部分行估计？
        return cells, genes, tot


def read_matrix_fast(path, genes_of_interest):
    """高效版：只读取感兴趣基因行，总 UMI 用全矩阵行求和"""
    with gzip.open(path, "rt", errors="replace") as f:
        header = f.readline()
        cells = [c for c in header.rstrip("\n").replace('"', "").split("\t") if c]
        n_cells = len(cells)
        genes = {g: np.zeros(n_cells, dtype=float) for g in genes_of_interest}
        tot = np.zeros(n_cells, dtype=float)
        for line in f:
            parts = line.split("\t")
            gname = parts[0].replace('"', "")
            arr = np.array([float(x) if x else 0.0 for x in parts[1:n_cells + 1]])
            tot += arr
            if gname in genes:
                genes[gname] = arr
        return cells, genes, tot


def classify(cells, genes, tot, min_umi=500):
    """按 marker 打分分配细胞类型；低质量细胞剔除"""
    out_celltype = []
    out_umi = []
    for i in range(len(cells)):
        umi = tot[i]
        out_umi.append(umi)
        if umi < min_umi:
            out_celltype.append("LowQC")
            continue
        scores = defaultdict(float)
        for g, arr in genes.items():
            if g in ("AMICA1", "JAML"):
                continue
            if arr[i] > 0:
                scores[MARKERS[g]] += 1
        # 取得分最高且 >0 的类别；Tie 按优先级
        best = None
        best_s = 0
        for ct in PRIORITY:
            if scores[ct] > best_s:
                best_s = scores[ct]
                best = ct
        if best is None:
            best = "Other/Unknown"
        out_celltype.append(best)
    return np.array(out_celltype), np.array(out_umi)


def main():
    all_rows = []
    for fn, patient, tissue, target_size in SAMPLES:
        path = os.path.join(DATA, fn)
        if not os.path.exists(path):
            print(f"[跳过] {fn} 未下载", flush=True)
            continue
        if os.path.getsize(path) != target_size:
            print(f"[跳过] {fn} 文件不完整（{os.path.getsize(path)} != {target_size}）", flush=True)
            continue
        print(f"[读取] {fn} ...", flush=True)
        cells, genes, tot = read_matrix_fast(path, set(MARKERS) | {"AMICA1", "JAML"})
        ct, umi = classify(cells, genes, tot)
        jaml = (genes.get("AMICA1", np.zeros(len(cells))) + genes.get("JAML", np.zeros(len(cells))))
        # 归一化 CPM（每 100 万 UMI）
        cpm = jaml / np.maximum(tot, 1) * 1e6
        df = pd.DataFrame({"celltype": ct, "JAML_UMI": jaml, "total_UMI": umi, "JAML_CPM": cpm})
        df["patient"] = patient
        df["tissue"] = tissue
        all_rows.append(df)
        # 打印概览
        print(f"  细胞 {len(cells)}, JAML+ 细胞 {(jaml > 0).sum()}, 按类型: {dict(pd.Series(ct).value_counts())}")
    if not all_rows:
        print("无可用样本")
        return
    d = pd.concat(all_rows, ignore_index=True)
    # 汇总：celltype × tissue
    summary = d.groupby(["celltype", "tissue"]).agg(
        n=("JAML_CPM", "size"),
        mean_CPM=("JAML_CPM", "mean"),
        pct_pos=("JAML_UMI", lambda x: (x > 0).mean() * 100),
        mean_UMI=("JAML_UMI", "mean"),
    ).round(2)
    print("\n=== JAML 表达 × 细胞类型 × 肿瘤/正常 ===")
    print(summary.to_string())
    summary.to_csv(os.path.join(OUT, "jaml_gse117570_summary.csv"))
    # 免疫亚群聚焦（只免疫）
    print("\n=== 免疫细胞亚群聚焦（mean CPM） ===")
    imm = summary.loc[[i for i in summary.index if i[0] in PRIORITY]]
    print(imm.to_string())
    # 保存全部细胞级数据（供画图）
    d.to_csv(os.path.join(OUT, "jaml_gse117570_cells.csv"), index=False)
    print(f"\n[完成] 总细胞 {len(d)}; 输出 results/jaml_gse117570_summary.csv")


if __name__ == "__main__":
    main()
