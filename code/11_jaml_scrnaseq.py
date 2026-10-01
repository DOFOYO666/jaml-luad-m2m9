# -*- coding: utf-8 -*-
"""11_jaml_scrnaseq.py - JAML 在 LUAD 肿瘤 scRNA-seq (GSE131907, Kim 2020) 中的细胞类型表达验证
流式读取归一化 log2TPM 矩阵（3GB），提取 JAML(AMICA1) 行，按 Cell_type / Cell_subtype 计算
平均表达与阳性率，验证 JAML 在免疫细胞（尤其 T/NK/单核）中的表达特异性。
用法: python 11_jaml_scrnaseq.py
"""
import gzip
import os

import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/gse131907"
ANN = os.path.join(DATA, "cell_annotation.txt.gz")
MTX = os.path.join(DATA, "matrix_log2TPM.txt.gz")
OUT = os.path.join(BASE, "results")
GENE_ALIASES = ["AMICA1", "JAML"]


def load_annotation():
    ann = pd.read_csv(ANN, sep="\t")
    print(f"[注释] {len(ann)} 细胞")
    print(f"  Cell_type 分布: {ann['Cell_type'].value_counts().to_dict()}")
    return ann


def extract_gene_row(mtx_path, aliases):
    """流式读取矩阵，提取目标基因行（矩阵: 行=基因, 列=细胞; 行首=基因名）"""
    with gzip.open(mtx_path, "rt", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        print(f"[矩阵] {len(cell_ids)} 细胞, 表头前3: {cell_ids[:3]}")
        for line in f:
            c = line.rstrip("\n").split("\t")
            gname = c[0]
            if gname in aliases or gname.split(".")[0] in aliases:
                vals = [float(x) if x not in ("", "NA") else 0.0 for x in c[1:]]
                print(f"[矩阵] 找到基因行: {gname}, 非零细胞 {sum(1 for v in vals if v > 0)}/{len(vals)}")
                return gname, cell_ids, vals
        return None, cell_ids, None


def main():
    ann = load_annotation()
    gname, cell_ids, vals = extract_gene_row(MTX, GENE_ALIASES)
    if vals is None:
        print("[错误] 矩阵中未找到 JAML/AMICA1 行")
        return
    # 细胞顺序: 注释的 Index 列应与矩阵列一一对应（按出现顺序）
    idx_col = [str(x) for x in ann["Index"]]
    if len(idx_col) != len(cell_ids):
        print(f"[警告] 细胞数不匹配: 注释 {len(idx_col)} vs 矩阵 {len(cell_ids)}")
        # 尝试用顺序对齐
    expr = pd.Series(vals, index=cell_ids, name="JAML_log2TPM")
    # 对齐注释（若 Index 与细胞 ID 不同，则按顺序）
    if idx_col == cell_ids:
        ann2 = ann.copy()
    else:
        ann2 = ann.iloc[: len(vals)].copy()
        ann2["JAML_log2TPM"] = vals
    ann2["JAML_log2TPM"] = vals
    ann2["JAML_pos"] = ann2["JAML_log2TPM"] > 0

    # 汇总：Cell_type 水平
    ct = ann2.groupby("Cell_type").agg(
        n_cells=("JAML_pos", "size"),
        mean_expr=("JAML_log2TPM", "mean"),
        pct_pos=("JAML_pos", "mean"),
    )
    ct["mean_expr"] = ct["mean_expr"].round(3)
    ct["pct_pos"] = (ct["pct_pos"] * 100).round(1)
    print("\n=== JAML 表达 × Cell_type ===")
    print(ct.sort_values("mean_expr", ascending=False).to_string())

    # Cell_subtype 水平（免疫亚群聚焦）
    cs = ann2.groupby("Cell_subtype").agg(
        n_cells=("JAML_pos", "size"),
        mean_expr=("JAML_log2TPM", "mean"),
        pct_pos=("JAML_pos", "mean"),
    )
    cs = cs[cs["n_cells"] >= 50]  # 排除极小亚群
    cs["mean_expr"] = cs["mean_expr"].round(3)
    cs["pct_pos"] = (cs["pct_pos"] * 100).round(1)
    print("\n=== JAML 表达 × Cell_subtype (n>=50) ===")
    print(cs.sort_values("mean_expr", ascending=False).head(30).to_string())

    # 保存
    ct.to_csv(os.path.join(OUT, "jaml_scrnaseq_celltype.csv"))
    cs.to_csv(os.path.join(OUT, "jaml_scrnaseq_subtype.csv"))
    print(f"\n[完成] 输出: results/jaml_scrnaseq_celltype.csv, jaml_scrnaseq_subtype.csv")
    # 免疫相关亚群摘要
    print("\n=== 免疫细胞亚群聚焦（T/NK/髓系/B） ===")
    imm = cs[cs.index.isin(["Cytotoxic CD8+ T", "Exhausted CD8+ T", "Naive CD8+ T", "CD8 low T",
                            "Naive CD4+ T", "Treg", "CD4+ Th", "NK",
                            "mo-Mac", "Monocytes", "Alveolar Mac", "CD1c+ DCs", "CD141+ DCs", "pDCs",
                            "MALT B cells", "Follicular B cells"])]
    print(imm.sort_values("mean_expr", ascending=False).to_string())


if __name__ == "__main__":
    main()
