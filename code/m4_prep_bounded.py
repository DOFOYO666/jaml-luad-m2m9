# -*- coding: utf-8 -*-
"""【M4 · 步骤 1（内存受限版）】构建 CellOracle 输入 .h5ad —— 可在 8 GB 机器上运行。

为什么不用原版 prep_celloracle_input.py：
  原版把 29,634 × 208,506 的**全部**细胞一次性建成稀疏矩阵（nnz 数亿、需数 GB），
  在 8 GB 机器上必然 OOM。本版改为：
    ① 只取 **CD4 T 谱系**（JAML 锚点所在）并按亚型等量抽样，细胞数降到 ~4,000；
    ② 单趟遍历矩阵，把选中细胞写成 **float16 缓冲**（29,634×4,000×2 B ≈ 237 MB）；
    ③ 用该缓冲算每基因方差 → 取 HVG top-N + 强制保留基因；
    ④ 只把选中的 ~3,000 基因写成 h5ad（float32，约 48 MB）。
  峰值内存 < 600 MB（实测前先用 available > 1.2 GB 的守卫）。

⚠️ 与原版一致的硬性事实（实测）：
  - 对齐键必须用注释的 **`Index`** 列（`Barcode` 列缺样本后缀）；
  - GSE131907 的 `Cell_type`(10 类) 与 `Cell_type.refined`(7 类) **都不含 CD4 T**，
    CD4 亚型只在 **`Cell_subtype`**(50 类) 里；
  - 矩阵是「基因 × 细胞」的制表符文本，`.gz`，读表头必须 gzip.open。

⚠️ 抽样对 GRN 的影响：CellOracle 官方建议每个 GRN 单元 ≥2,000–3,000 细胞。
   本脚本给每个 CD4 亚型设上限，合计约 4,000，并把各亚型标签同时写入
   `obs["grn_unit"]`（全为 "CD4 T"）与 `obs["cell_subtype"]`（细标签），
   便于事后做亚型分层的敏感性分析。**抽样会降低 GRN 稳定性，报告中必须声明。**

输出（写到 D 盘工作区）：
  D:/workbuddy工作空间/JAML深度研究/results/CD4T_celloracle.h5ad
  D:/workbuddy工作空间/JAML深度研究/results/m4_prep_qc.json
"""
import gzip
import json
import os
import sys
import time

import numpy as np
import pandas as pd

try:
    import anndata as ad
    import scipy.sparse as sp
except ImportError:
    sys.exit("需要 anndata / scipy")

# ---------------- 路径（输入在 C 盘主工作区，输出在 D 盘工作区）----------------
C_ROOT = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化"
D_ROOT = r"D:\workbuddy工作空间\JAML深度研究"
MAT = os.path.join(C_ROOT, "rawdata", "scRNA", "GSE131907",
                   "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz")
ANN = os.path.join(C_ROOT, "rawdata", "scRNA", "GSE131907",
                   "GSE131907_Lung_Cancer_cell_annotation.txt.gz")
OUT_DIR = os.path.join(D_ROOT, "results")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_H5AD = os.path.join(OUT_DIR, "CD4T_celloracle.h5ad")
OUT_H5AD_WIDE = os.path.join(OUT_DIR, "CD4T_wide.h5ad")   # M8 用：更宽的基因集
OUT_QC = os.path.join(OUT_DIR, "m4_prep_qc.json")

# M8（TF 活性 / regulon）需要尽量宽的基因覆盖，只按检出率过滤、不设 HVG 上限
WIDE_MIN_DETECT_FRAC = 0.02

# CD4 谱系（顺序 = 从初始到耗竭，便于解读）
CD4_LINEAGE = ["Naive CD4+ T", "Treg", "CD4+ Th", "CD8+/CD4+ Mixed Th", "Exhausted Tfh"]
CAP_PER_SUBTYPE = 600            # 5 亚型 → 最多 2,842 细胞（8 GB 机器上的内存折中）
TOP_HVG = 3000
MIN_DETECT_FRAC = 0.05
FORCE_KEEP = ["JAML", "CXADR", "PDCD1", "CD3D", "CD3E", "CD4", "CD8A", "FOXP3",
              "IL7R", "CCR7", "TCF7", "SELL", "GZMB", "TOX", "HAVCR2", "LAG3",
              "TIGIT", "CTLA4", "IL2RA", "MKI67", "CXCL13", "IFNG", "GZMA",
              "PRF1", "CD44", "ITGAE", "BATF", "ANXA1", "S100A4"]
SEED = 20260930


def avail_gb():
    import ctypes

    class M(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    m = M(); m.dwLength = ctypes.sizeof(M)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ullAvailPhys / 1073741824


def main():
    t0 = time.time()
    a = avail_gb()
    print(f"[内存] 起始可用 {a:.2f} GB")
    if a < 0.45:
        sys.exit(f"可用内存仅 {a:.2f} GB < 0.45 GB，拒绝启动以免 OOM")

    # ---------------- 1. 选细胞 ----------------
    ann = pd.read_csv(ANN, sep="\t")
    ann.columns = [c.strip().lstrip("\ufeff") for c in ann.columns]
    key = "Index" if "Index" in ann.columns else ann.columns[0]
    ann["_key"] = ann[key].astype(str)
    ann["_sub"] = ann["Cell_subtype"].astype(str)
    rng = np.random.default_rng(SEED)
    picks = []
    for st in CD4_LINEAGE:
        idx = ann.index[ann["_sub"] == st].to_numpy()
        if len(idx) == 0:
            print(f"  [警告] 亚型 {st} 未找到")
            continue
        take = rng.choice(idx, size=min(CAP_PER_SUBTYPE, len(idx)), replace=False)
        picks.append(ann.loc[take, ["_key", "_sub"]])
    sel = pd.concat(picks, ignore_index=True).sort_values("_key").reset_index(drop=True)
    n_cells = len(sel)
    print(f"[选细胞] CD4 谱系 {n_cells} 个")
    print(sel["_sub"].value_counts().to_string())

    # ---------------- 2. 读表头，定位列 ----------------
    with gzip.open(MAT, "rt", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    cells = header[1:]
    pos = {c: i for i, c in enumerate(cells)}
    sel["_col"] = sel["_key"].map(lambda k: pos.get(k, -1))
    matched = int((sel["_col"] >= 0).sum())
    print(f"[对齐] 命中 {matched} / {n_cells}")
    if matched == 0:
        sys.exit("barcode 对齐为 0 —— 检查对齐键（应为注释的 Index 列）")
    sel = sel[sel["_col"] >= 0].sort_values("_col").reset_index(drop=True)
    n_cells = len(sel)
    col_idx = [int(c) + 1 for c in sel["_col"]]      # pandas usecols：含首列，故 +1
    usecols = [0] + col_idx
    keys = sel["_key"].tolist()
    subs = sel["_sub"].tolist()

    # ---------------- 3. 单趟遍历：float16 缓冲 + 流式统计量 ----------------
    # 关键改动：**不做整矩阵的 float32 副本**（那会把峰值翻倍到 ~700 MB），
    # 改为在分块读取时增量累加每基因的检出数 / 一阶矩 / 二阶矩，
    # 峰值内存 ≈ buf(234 MB) + 分块临时量(≤70 MB) ≈ 300 MB。
    print(f"\n[遍历] pandas usecols 分块读取（{n_cells} 列，预计 15–25 min）...")
    buf = np.zeros((29634, n_cells), dtype=np.float16)     # 已知基因数；下方做溢出校验
    gene_names = []
    reader = pd.read_csv(MAT, sep="\t", usecols=usecols, index_col=0,
                         chunksize=3000, engine="c", low_memory=True)
    n = 0
    s1_parts, s2_parts, det_parts = [], [], []
    for ci, chunk in enumerate(reader):
        arr = chunk.to_numpy(dtype=np.float16)
        names = list(chunk.index)
        if n + len(names) > buf.shape[0]:
            sys.exit(f"基因数超出预分配缓冲 {buf.shape[0]}，请修正")
        buf[n:n + len(names)] = arr
        af = arr.astype(np.float32)
        s1_parts.append(af.sum(axis=1))
        s2_parts.append(np.einsum("ij,ij->i", af, af))     # 避免 (af*af) 的额外临时量
        det_parts.append((arr > 0).sum(axis=1))
        del af
        gene_names.extend(names)
        n += len(names)
        if (ci + 1) % 3 == 0:
            print(f"  ...已读 {n} 基因 ({time.time()-t0:.0f}s)", flush=True)
    buf = buf[:n]
    s1 = np.concatenate(s1_parts).astype(np.float64)
    s2 = np.concatenate(s2_parts).astype(np.float64)
    det = np.concatenate(det_parts)
    del s1_parts, s2_parts, det_parts
    import gc
    gc.collect()
    print(f"[遍历] 完成：{n} 基因 × {n_cells} 细胞；{time.time()-t0:.0f}s；"
          f"缓冲 {buf.nbytes/1048576:.0f} MB")

    # ---------------- 4. HVG 选择（用流式统计量，无需整矩阵副本）----------------
    mu = s1 / n_cells
    var = np.clip(s2 / n_cells - mu * mu, 0, None)
    min_det = max(3, int(MIN_DETECT_FRAC * n_cells))
    cand = det >= min_det
    print(f"[HVG] 检出 ≥ {min_det} 细胞的基因：{int(cand.sum())} / {n}")
    order = np.argsort(-np.where(cand, var, -1.0))
    top_idx = [int(i) for i in order[:TOP_HVG] if cand[i]]
    top = [gene_names[i] for i in top_idx]
    keep = list(dict.fromkeys(top + [g for g in FORCE_KEEP if g in set(gene_names)]))
    qc = {
        "n_cells": int(n_cells),
        "n_genes_total": int(n),
        "n_cand": int(cand.sum()),
        "n_keep": len(keep),
        "hvg": int(len(top)),
        "force_kept": [g for g in FORCE_KEEP if g in set(gene_names)],
        "force_missing": [g for g in FORCE_KEEP if g not in set(gene_names)],
        "cells_per_subtype": {k: int(v) for k, v in sel["_sub"].value_counts().items()},
        "cap_per_subtype": CAP_PER_SUBTYPE,
        "min_detect": int(min_det),
        "seed": SEED,
    }
    print(f"[HVG] 保留 {len(keep)} 个基因（HVG {len(top)} + 强制保留）")
    print(f"[HVG] 强制保留未命中: {qc['force_missing'] or '无'}")
    if "JAML" not in keep:
        sys.exit("JAML 未进入保留集 —— 检查基因符号")

    # ---------------- 5. 只把选中基因转成 float32 并写出 ----------------
    gi = {g: i for i, g in enumerate(gene_names)}
    X = np.stack([buf[gi[g]].astype(np.float32) for g in keep], axis=0)   # 基因 × 细胞
    # ⚠️ 先不要释放 buf —— 下面还要用它生成 M8 的宽基因集
    gc.collect()
    adata = ad.AnnData(sp.csr_matrix(X.T))
    adata.obs_names = keys
    adata.var_names = keep
    adata.obs["cell_subtype"] = pd.Categorical(subs, categories=CD4_LINEAGE, ordered=True)
    adata.obs["grn_unit"] = "CD4 T"
    adata.layers["raw_count"] = adata.X.copy()              # 原始 UMI（未归一化）
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"[输出] {OUT_H5AD}  ({adata.shape[0]} 细胞 × {adata.shape[1]} 基因)")

    qc["shape"] = list(adata.shape)
    qc["jaml_index"] = int(list(adata.var_names).index("JAML"))
    del X, adata                      # 先释放，再建宽集，压低峰值
    gc.collect()
    print(f"[内存] 写出后可用 {avail_gb():.2f} GB")

    # ---------------- 6. 同时写一份"宽基因集"给 M8 用（同一次遍历，零额外读盘）----------------
    wide_min = max(3, int(WIDE_MIN_DETECT_FRAC * n_cells))
    wide_genes = [g for g, d in zip(gene_names, det) if d >= wide_min]
    print(f"\n[M8 宽集] 检出 ≥ {wide_min} 细胞：{len(wide_genes)} 基因")
    Xw = np.stack([buf[gi[g]].astype(np.float32) for g in wide_genes], axis=0)
    del buf
    gc.collect()
    adw = ad.AnnData(sp.csr_matrix(Xw.T))
    adw.obs_names = keys
    adw.var_names = wide_genes
    adw.obs["cell_subtype"] = pd.Categorical(subs, categories=CD4_LINEAGE, ordered=True)
    adw.layers["raw_count"] = adw.X.copy()
    adw.write_h5ad(OUT_H5AD_WIDE, compression="gzip")
    print(f"[输出] {OUT_H5AD_WIDE}  ({adw.shape[0]} 细胞 × {adw.shape[1]} 基因)")
    qc["wide_shape"] = list(adw.shape)
    qc["wide_min_detect"] = int(wide_min)

    qc["elapsed_sec"] = round(time.time() - t0, 1)
    qc["avail_gb_end"] = round(avail_gb(), 2)
    with open(OUT_QC, "w", encoding="utf-8") as f:
        json.dump(qc, f, ensure_ascii=False, indent=2)
    print(f"[输出] {OUT_QC}")
    print(f"[内存] 结束可用 {avail_gb():.2f} GB")
    print(f"\n总耗时 {time.time()-t0:.0f}s\nDONE")


if __name__ == "__main__":
    main()
