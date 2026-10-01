# -*- coding: utf-8 -*-
"""M3（Python 路线）：用 liana-py 推断 GSE131907 的细胞通讯网络，重点关注 JAML–CXADR。

替代 CellChat 的理由：本机无 Rtools，CellChat（GitHub，含 C++）无法编译。

内存策略（本机 8 GB）：
  1. 先取 liana 资源（consensus / CellChatDB 等）中出现的**全部配体/受体基因**；
  2. 只保留这些基因行（其余行整行跳过，不做字段切分）；
  3. 每个分组**随机抽取至多 MAX_PER_GROUP 个细胞**（组均值与置换检验所需样本量绰绰有余）。

分组列使用 **`Cell_subtype`** —— GSE131907 的 `Cell_type` 与 `Cell_type.refined` 都不含 CD4 T。
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

MAT = "rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
ANN = "rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
OUTDIR = "results"
MAX_PER_GROUP = 800
RNG = np.random.default_rng(20260929)

try:
    import anndata as ad
    import liana as li
except ImportError as e:
    sys.exit(f"缺少依赖: {e}\n请先 pip install liana")

# ---------------- 1. 取 LR 资源基因 ----------------
try:
    res = li.resource.select_resource("consensus")
    print(f"[资源] consensus: {res.shape}")
except Exception as e:
    print(f"[资源] consensus 不可用({e})，退回 cellchatdb")
    res = li.resource.select_resource("cellchatdb")
lr_genes = sorted(set(res["ligand"]).union(set(res["receptor"])))
print(f"[资源] 配体+受体基因 {len(lr_genes)} 个；JAML 在内: {'JAML' in lr_genes}；CXADR 在内: {'CXADR' in lr_genes}")

# ---------------- 2. 注释与抽样 ----------------
ann = pd.read_csv(ANN, sep="\t")
ann["Index"] = ann["Index"].astype(str)
ann = ann.set_index("Index")

# 只取原发灶 / 正常肺（排除脑转移、淋巴结转移、胸水），避免区室混杂
KEEP_ORIGIN = ("tLung", "nLung", "tL/B", "NT")
mask_origin = ann["Sample_Origin"].astype(str).isin(KEEP_ORIGIN)
print(f"[过滤] Sample_Origin ∈ {KEEP_ORIGIN}: {mask_origin.sum()} / {len(ann)}")
ann = ann[mask_origin]

keep_idx = []
for sub, grp in ann.groupby("Cell_subtype"):
    if len(grp) <= MAX_PER_GROUP:
        keep_idx.extend(grp.index.tolist())
    else:
        keep_idx.extend(RNG.choice(grp.index.values, MAX_PER_GROUP, replace=False).tolist())
keep_set = set(keep_idx)
print(f"[抽样] 保留 {len(keep_set)} 细胞，覆盖 {ann.loc[list(keep_set), 'Cell_subtype'].nunique()} 个 Cell_subtype")

# ---------------- 3. 流式读取目标基因 ----------------
rows_i, cols_i, vals = [], [], []
gene_names = []
with gzip.open(MAT, "rt", errors="replace") as f:
    header = f.readline().rstrip("\n").split("\t")
    cells = header[1:]
    col_keep = np.array([c in keep_set for c in cells])
    n_keep = col_keep.sum()
    print(f"[矩阵] 细胞 {len(cells)}；抽样命中 {n_keep}")
    gi = -1
    for line in f:
        p = line.find("\t")
        name = line[:p]
        if name not in lr_genes:
            continue
        gi += 1
        gene_names.append(name)
        v = np.fromstring(line[p + 1:], dtype=np.float32, sep="\t")[col_keep]
        nz = np.nonzero(v)[0]
        if nz.size:
            rows_i.append(np.full(nz.size, gi, dtype=np.int32))
            cols_i.append(nz.astype(np.int32))
            vals.append(v[nz])
        if gi % 200 == 0 and gi > 0:
            print(f"  已读 {gi} 个目标基因...", flush=True)

if not rows_i:
    sys.exit("未取到任何目标基因的表达值")

X = sp.csr_matrix((np.concatenate(vals),
                   (np.concatenate(rows_i), np.concatenate(cols_i))),
                  shape=(len(gene_names), n_keep), dtype=np.float32).T
print(f"[矩阵] {X.shape[0]} 细胞 × {X.shape[1]} 基因, nnz={X.nnz:,}, 内存约 {X.data.nbytes/1e6:.0f} MB")

obs = ann.loc[[c for c, k in zip(cells, col_keep) if k]].copy()
adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(gene_names, name="gene")))
adata.layers["counts"] = adata.X.copy()

# CPM 归一化 + log1p（liana 的 CellPhoneDB 方法期望归一化数据）
sp_sum = np.asarray(adata.X.sum(axis=1)).ravel()
sp_sum[sp_sum == 0] = 1
adata.X = sp.csr_matrix(adata.X.multiply(1e6 / sp_sum[:, None]))
adata.X.data = np.log1p(adata.X.data)
print("[归一化] CPM + log1p 完成")

# ---------------- 4. liana 通讯推断 ----------------
MIN_CELLS = 20
vc = adata.obs["Cell_subtype"].value_counts()
dropped = vc[vc < MIN_CELLS].index.tolist()
if dropped:
    print(f"[过滤] 细胞数 < {MIN_CELLS} 的亚型被剔除: {dropped}")
    adata = adata[~adata.obs["Cell_subtype"].isin(dropped)].copy()
print(f"[输入] {adata.shape}；分组数 {adata.obs['Cell_subtype'].nunique()}")

li.mt.rank_aggregate(adata, groupby="Cell_subtype", resource_name="consensus",
                     expr_prop=0.1, verbose=True, n_perms=200, seed=20260929,
                     use_raw=False)
lr = adata.uns["liana_res"]
lr.to_csv(os.path.join(OUTDIR, "m3_liana_rank_aggregate.csv"), index=False)
print(f"[liana] rank_aggregate 完成，{lr.shape} 行")

# ---------------- 5. 抽取 JAML / CXADR 相关 ----------------
cols = [c for c in lr.columns]
print("[liana] 输出列:", cols)
sel = lr[(lr["ligand_complex"].astype(str).str.contains("JAML", na=False)) |
         (lr["receptor_complex"].astype(str).str.contains("CXADR", na=False)) |
         (lr["receptor_complex"].astype(str).str.contains("JAML", na=False)) |
         (lr["ligand_complex"].astype(str).str.contains("CXADR", na=False))]
sel.to_csv(os.path.join(OUTDIR, "m3_liana_jaml_cxadr.csv"), index=False)
print(f"\n[liana] JAML/CXADR 相关交互 {len(sel)} 条")
pd.set_option("display.width", 220)
show = [c for c in ["source", "target", "ligand_complex", "receptor_complex", "lr_means",
                    "cellphone_pvals", "magnitude_rank", "specificity_rank"] if c in sel.columns]
if len(sel):
    print(sel.sort_values("magnitude_rank").head(25)[show].to_string(index=False))

# 网络图中 JAML–CXADR 的具体方向
jc = lr[(lr["ligand_complex"].astype(str).str.upper() == "JAML") &
        (lr["receptor_complex"].astype(str).str.upper() == "CXADR")]
print(f"\n[liana] JAML → CXADR 方向共 {len(jc)} 条")
if len(jc):
    print(jc.sort_values("magnitude_rank")[show].to_string(index=False))
    jc.to_csv(os.path.join(OUTDIR, "m3_liana_JAML_CXADR_directed.csv"), index=False)

print("\n[完成] results/m3_liana_rank_aggregate.csv, m3_liana_jaml_cxadr.csv(, m3_liana_JAML_CXADR_directed.csv)")
