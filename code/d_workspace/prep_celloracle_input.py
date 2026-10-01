# -*- coding: utf-8 -*-
"""【大内存任务 · 步骤 1】从 GSE131907 raw UMI 文本矩阵构建 CellOracle 输入 .h5ad。

⚠️ 运行环境要求：**≥ 64 GB 内存**（本机 8 GB 无法容纳 29,634 × 208,506 稀疏矩阵）。
   在 8 GB 机器上请勿运行本脚本。

输入：
  rawdata/scRNA_GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
  rawdata/scRNA_GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz
输出：
  results/GSE131907_full.h5ad          （全部细胞，layers['raw_count']）
  results/GSE131907_TNK_epi.h5ad       （仅 T/NK + 上皮 + 髓系谱系，供 CellOracle）

关键实现说明：
  - 文本矩阵为「基因 × 细胞」，逐行流式解析后构建 **CSC** 稀疏矩阵（列=细胞）
  - 注释以 `Index` 列与矩阵列名精确对齐（本项目实测 0 例不匹配）
  - `.obs["grn_unit"]` 由 **`Cell_subtype`** 映射而来 —— 注意 GSE131907 的
    `Cell_type`（10 类）与 `Cell_type.refined`（7 类）**都不含 CD4 T**，
    CD4 亚型只存在于 `Cell_subtype`（50 类）
"""
import gzip
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import anndata as ad
except ImportError:
    sys.exit("需要 anndata：pip install anndata")

DATA = "rawdata/scRNA_GSE131907"
MAT = os.path.join(DATA, "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz")
ANN = os.path.join(DATA, "GSE131907_Lung_Cancer_cell_annotation.txt.gz")
OUT_ALL = "results/GSE131907_full.h5ad"
OUT_SUB = "results/GSE131907_TNK_epi.h5ad"
os.makedirs("results", exist_ok=True)

t0 = time.time()

# ---------- 1. 注释 ----------
ann = pd.read_csv(ANN, sep="\t")
ann["Index"] = ann["Index"].astype(str)
print(f"注释 {ann.shape}；Cell_type {ann['Cell_type'].nunique()} 类；"
      f"Cell_subtype {ann['Cell_subtype'].nunique()} 类")

# CD4 T 谱系映射（依据 GSE131907 实际取值）
CD4_SUB = {
    "Naive CD4+ T": "CD4 T",
    "CD4+ Th": "CD4 T",
    "CD8+/CD4+ Mixed Th": "CD4 T",
    "Treg": "CD4 T",
}
CD8_SUB = {"Naive CD8+ T": "CD8 T", "Cytotoxic CD8+ T": "CD8 T",
           "Exhausted CD8+ T": "CD8 T", "CD8 low T": "CD8 T"}
NK_SUB = {"NK": "NK"}
MYE_SUB = {"mo-Mac": "Myeloid", "Alveolar Mac": "Myeloid", "Monocytes": "Myeloid",
           "Pleural Mac": "Myeloid", "DC": "Myeloid"}


def to_grn_unit(sc):
    for m in (CD4_SUB, CD8_SUB, NK_SUB, MYE_SUB):
        if sc in m:
            return m[sc]
    return "Other"


ann["grn_unit"] = ann["Cell_subtype"].map(to_grn_unit)
print("grn_unit 分布:", ann["grn_unit"].value_counts().to_dict())


# ---------- 2. 流式读取矩阵 → CSC ----------
def read_matrix():
    rows_idx, cols_idx, vals = [], [], []
    n_genes = 0
    with gzip.open(MAT, "rt", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        cells = header[1:]
        print(f"矩阵列（细胞）数: {len(cells)}")
        for line in f:
            p = line.find("\t")
            v = np.fromstring(line[p + 1:], dtype=np.float32, sep="\t")
            nz = np.nonzero(v)[0]
            if nz.size:
                rows_idx.append(np.full(nz.size, n_genes, dtype=np.int32))
                cols_idx.append(nz.astype(np.int32))
                vals.append(v[nz])
            n_genes += 1
            if n_genes % 2000 == 0:
                print(f"  已读 {n_genes} 基因，nnz={sum(a.size for a in vals):,} "
                      f"({time.time()-t0:.0f}s)", flush=True)
    rows = np.concatenate(rows_idx); cols = np.concatenate(cols_idx)
    vals = np.concatenate(vals)
    print(f"矩阵完成: {n_genes} 基因 × {len(cells)} 细胞, nnz={vals.size:,} "
          f"({time.time()-t0:.0f}s)")
    # 基因名
    names = []
    with gzip.open(MAT, "rt", errors="replace") as f:
        f.readline()
        for line in f:
            names.append(line[:line.find("\t")])
    X = sp.csc_matrix((vals, (rows, cols)), shape=(n_genes, len(cells)), dtype=np.float32)
    return X, names, cells


X, gene_names, cell_ids = read_matrix()

# ---------- 3. AnnData ----------
obs = ann.set_index("Index").reindex(cell_ids)
miss = obs["Cell_type"].isna().sum()
print(f"未匹配注释的细胞: {miss}")
obs["grn_unit"] = obs["grn_unit"].fillna("Unmatched")

adata = ad.AnnData(X=sp.csr_matrix(X.T).astype(np.float32),
                   obs=obs, var=pd.DataFrame(index=pd.Index(gene_names, name="gene")))
adata.layers["raw_count"] = adata.X.copy()
print(adata)
adata.write_h5ad(OUT_ALL, compression="gzip")
print(f"已写出 {OUT_ALL}")

# ---------- 4. 子集（CellOracle 用） ----------
keep = adata.obs["grn_unit"].isin(["CD4 T", "CD8 T", "NK", "Myeloid"]).values
sub = adata[keep].copy()
print(f"子集: {sub.shape}；grn_unit 分布 {sub.obs['grn_unit'].value_counts().to_dict()}")
# 记录来源样本，便于按原发灶过滤
print("Sample_Origin 分布:", sub.obs["Sample_Origin"].value_counts().head(10).to_dict())
sub.write_h5ad(OUT_SUB, compression="gzip")
print(f"已写出 {OUT_SUB}")
print(f"总耗时 {time.time()-t0:.0f}s")
