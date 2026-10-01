# -*- coding: utf-8 -*-
"""M3：GSE131907 细胞通讯分析（自实现，基于 liana 官方共识资源）

**为什么不用 liana 本体**：`pip install liana` 依赖 numba/llvmlite，本次安装受网络限制
（llvmlite 41.9 MB 下载速率 ~2–40 kB/s，长时间未完成）。为不阻塞分析，改用
**liana 官方的 `omni_resource.csv`（OmniPath 资源，含 consensus 子集）**，
以向量化 numpy 实现与 liana 等价的「组均值 + 标签置换」评分。

**方法学声明（必须随结果一并报告）**：
  - 配体-受体对取自 **liana v1.x 官方 `omni_resource.csv` 的 `consensus` 子集**
    （已核实 JAML–CXADR 在该子集中，且为双向注释）
  - 评分：`lr_means` = (发送方配体组均值 + 接收方受体组均值) / 2（log1p-CPM 尺度）
  - 显著性：对**细胞分组标签做 200 次置换**得到零分布，取经验单侧 p 值
  - 该实现**不是 liana 本体调用**；liana 安装完成后应用同一输入复核

输入：GSE131907 raw UMI 矩阵 + 细胞注释；ref/liana/omni_resource.csv
输出：results/m3_ccc_lr_means.csv、m3_ccc_JAML_CXADR.csv、m3_ccc_summary.json
"""
import gzip
import json
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

MAT = "rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
ANN = "rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
RES = "ref/liana/omni_resource.csv"
OUTDIR = "results"
MAX_PER_GROUP = 400
N_PERM = 200
SEED = 20260929
rng = np.random.default_rng(SEED)

# ---------- 1. 资源：consensus 子集 ----------
res = pd.read_csv(RES)
cons = res[res["resource"] == "consensus"][["source_genesymbol", "target_genesymbol"]].dropna()
cons = cons.drop_duplicates()
print(f"[资源] consensus 配体-受体对 {len(cons)} 条")
lr_genes = sorted(set(cons["source_genesymbol"]) | set(cons["target_genesymbol"]))
print(f"[资源] 涉及基因 {len(lr_genes)} 个；JAML 在内: {'JAML' in lr_genes}；CXADR 在内: {'CXADR' in lr_genes}")

# ---------- 2. 注释与抽样 ----------
ann = pd.read_csv(ANN, sep="\t")
ann["Index"] = ann["Index"].astype(str)
ann = ann.set_index("Index")
KEEP_ORIGIN = ("tLung", "nLung", "tL/B", "NT")
n0 = len(ann)
ann = ann[ann["Sample_Origin"].astype(str).isin(KEEP_ORIGIN)]
print(f"[过滤] Sample_Origin ∈ {KEEP_ORIGIN}: {len(ann)} / {n0}")

keep = []
for sub, grp in ann.groupby("Cell_subtype"):
    take = min(len(grp), MAX_PER_GROUP)
    keep.extend(rng.choice(grp.index.values, take, replace=False).tolist())
keep_set = set(keep)
print(f"[抽样] {len(keep_set)} 细胞（每组 ≤{MAX_PER_GROUP}）")

# ---------- 3. 流式读取 LR 基因 ----------
rows, cols, vals, genes = [], [], [], []
with gzip.open(MAT, "rt", errors="replace") as f:
    header = f.readline().rstrip("\n").split("\t")
    cells = header[1:]
    col_keep = np.array([c in keep_set for c in cells])
    nkc = int(col_keep.sum())
    gi = -1
    for line in f:
        p = line.find("\t")
        name = line[:p]
        if name not in lr_genes:
            continue
        gi += 1
        genes.append(name)
        v = np.fromstring(line[p + 1:], dtype=np.float32, sep="\t")[col_keep]
        nz = np.nonzero(v)[0]
        if nz.size:
            rows.append(np.full(nz.size, gi, dtype=np.int32))
            cols.append(nz.astype(np.int32))
            vals.append(v[nz])
print(f"[矩阵] {nkc} 细胞 × {len(genes)} 基因（命中 LR 基因）")

X = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                  shape=(len(genes), nkc), dtype=np.float32)

# 归一化：CPM + log1p
tot = np.asarray(X.sum(axis=0)).ravel()
tot[tot == 0] = 1
X = sp.csr_matrix(X.multiply(1e6 / tot).tocsr())
X.data = np.log1p(X.data)
Xd = np.asarray(X.T.todense(), dtype=np.float32)     # 细胞 × 基因
print(f"[归一化] log1p-CPM 完成，稠密矩阵 {Xd.shape} = {Xd.nbytes/1e6:.0f} MB")

cell_ids = [c for c, k in zip(cells, col_keep) if k]
sub_obs = ann.loc[cell_ids]
groups = sub_obs["Cell_subtype"].values
gname, gidx = np.unique(groups, return_inverse=True)
n_g = len(gname)
print(f"[分组] {n_g} 个 Cell_subtype")

# ---------- 4. 组均值 ----------
def group_means(labels_idx, n_groups):
    M = sp.csr_matrix((np.ones(len(labels_idx), dtype=np.float32),
                       (labels_idx, np.arange(len(labels_idx)))),
                      shape=(n_groups, len(labels_idx)))
    counts = np.asarray(M.sum(axis=1)).ravel()
    counts[counts == 0] = 1
    G = np.asarray(M @ Xd) / counts[:, None]
    return G, counts

G_obs, gcount = group_means(gidx, n_g)
g2i = {g: i for i, g in enumerate(gname)}
print(f"[组均值] 完成；细胞数中位 {np.median(gcount):.0f}，最小 {gcount.min():.0f}")

# ---------- 5. 配体-受体评分 ----------
def score_pairs(G, pair_idx):
    """lr_means = (配体在发送方的均值 + 受体在接收方的均值)/2；返回 (n_group^2, n_pair)。"""
    gl = G[:, pair_idx[0]]        # 组 × 对（配体表达，作为发送方）
    gr = G[:, pair_idx[1]]        # 组 × 对（受体表达，作为接收方）
    # 外积式组合：source × target
    src = gl[:, None, :]          # n_g × 1 × n_pair
    tgt = gr[None, :, :]          # 1 × n_g × n_pair
    return (src + tgt) / 2.0      # n_g × n_g × n_pair

gene2col = {g: i for i, g in enumerate(genes)}
pairs, pair_idx = [], ([], [])
for _, r in cons.iterrows():
    l, rr = r["source_genesymbol"], r["target_genesymbol"]
    if l in gene2col and rr in gene2col:
        pairs.append((l, rr))
        pair_idx[0].append(gene2col[l])
        pair_idx[1].append(gene2col[rr])
pair_idx = (np.array(pair_idx[0]), np.array(pair_idx[1]))
print(f"[配对] 数据中可评分的 LR 对 {len(pairs)} / {len(cons)}")

obs = score_pairs(G_obs, pair_idx)                       # n_g × n_g × n_pair
# 置换零分布（打乱分组标签）
perm_max = np.zeros((N_PERM, len(pairs)), dtype=np.float32)
perm_all = np.zeros((N_PERM, len(pairs)), dtype=np.float32)
for b in range(N_PERM):
    gp = rng.permutation(gidx)
    Gp, _ = group_means(gp, n_g)
    Sp = score_pairs(Gp, pair_idx)
    perm_max[b] = Sp.max(axis=(0, 1))
    perm_all[b] = Sp.mean(axis=(0, 1))
pvals = ((perm_max >= obs.max(axis=(0, 1))).sum(axis=0) + 1) / (N_PERM + 1)

# 展平为长表
rows_out = []
for k, (l, rr) in enumerate(pairs):
    M = obs[:, :, k]
    si, ti = np.unravel_index(np.argsort(M, axis=None)[::-1][:5], M.shape)
    for a, b in zip(si, ti):
        rows_out.append({"ligand": l, "receptor": rr, "source": gname[a], "target": gname[b],
                         "lr_means": float(M[a, b]),
                         "source_n": int(gcount[a]), "target_n": int(gcount[b]),
                         "global_max_lr_means": float(M.max()),
                         "perm_p_global_max": float(pvals[k]),
                         "obs_mean_lr_means": float(M.mean()),
                         "perm_mean_lr_means": float(perm_all[:, k].mean())})
out = pd.DataFrame(rows_out)
out.to_csv(os.path.join(OUTDIR, "m3_ccc_lr_means.csv"), index=False)
print(f"\n[输出] 全部 LR 交互（每对取 top5 方向）{len(out)} 行 → results/m3_ccc_lr_means.csv")

# ---------- 6. JAML / CXADR 专项 ----------
jm = out[(out["ligand"] == "JAML") | (out["receptor"] == "CXADR") |
         (out["ligand"] == "CXADR") | (out["receptor"] == "JAML")]
jm.to_csv(os.path.join(OUTDIR, "m3_ccc_JAML_CXADR.csv"), index=False)
pd.set_option("display.width", 200)
print("\n" + "=" * 100)
print("JAML / CXADR 相关交互（按 lr_means 降序，各对 top5 方向）")
print("=" * 100)
print(jm.sort_values("lr_means", ascending=False).head(25).to_string(index=False))

# 只看 JAML→CXADR 与 CXADR→JAML 的完整组间矩阵
for lab, L, R in [("JAML → CXADR（配体在发送方，受体在接收方）", "JAML", "CXADR"),
                  ("CXADR → JAML", "CXADR", "JAML")]:
    if (L, R) in pairs:
        k = pairs.index((L, R))
        M = pd.DataFrame(obs[:, :, k], index=gname, columns=gname)
        print(f"\n--- {lab}：lr_means 矩阵（行=发送方，列=接收方）---")
        print(M.round(3).to_string())
        M.to_csv(os.path.join(OUTDIR, f"m3_ccc_matrix_{L}_to_{R}.csv"))

json.dump({"resource": "liana omni_resource.csv (consensus subset)",
           "n_pairs_scored": len(pairs), "n_perm": N_PERM,
           "n_cells": int(nkc), "n_groups": int(n_g),
           "max_per_group": MAX_PER_GROUP,
           "note": "自实现（组均值 + 标签置换），非 liana 本体调用"},
          open(os.path.join(OUTDIR, "m3_ccc_summary.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print("\n[完成] results/m3_ccc_lr_means.csv, m3_ccc_JAML_CXADR.csv, m3_ccc_matrix_*.csv, m3_ccc_summary.json")
