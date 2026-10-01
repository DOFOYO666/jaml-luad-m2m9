# -*- coding: utf-8 -*-
"""
32_m7_pseudotime.py — M7 拟时序：JAML 沿 CD4⁺T 分化轨迹

问题：M3/M2 都已显示 JAML 在 CD4⁺T 内部并非均一（Naive CD4+ T 低、CD4+ Th 高、
      耗竭 Tfh 最高）。本脚本用扩散拟时序（diffusion pseudotime, DPT）把这条
      梯度放到一条连续轨迹上，检验"JAML 随 CD4⁺T 活化/耗竭而上调"。

⚠️ 性能教训（第一版跑挂了）：不要用 `np.fromstring(line, sep='\\t')` 逐行解析
   全矩阵——29,634 行 × 208,506 个值，实测 200 行就要 9.3 s，全矩阵约 23 min，
   而且原设计要跑三趟。改为 **pandas `read_csv(usecols=..., chunksize=...)`**
   只取需要的细胞列，单趟完成，约 10 min。
   另外**不要用预归一化的 log2TPM 再算文库大小**——官方 log2TPM 已是 log2(TPM+1)，
   单趟即可同时得到检出率、均值与方差，无需 libsize。

根节点（root）选取：
  在 Naive CD4+ T 亚群内，取 naive 评分最高（CCR7+IL7R+TCF7 − GZMB−PDCD1）
  的细胞作为轨迹起点。评分用 z 分数，公式在代码中显式写出，便于复核。

⚠️ 符号：GSE131907 用 JAML（非 AMICA1）；对齐键用注释的 Index 列。
输入：rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz
      rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz
输出：results/m7_pseudotime_cells.csv
      results/m7_pseudotime_jaml_binned.csv
      results/m7_pseudotime_summary.json
      results/m7_hvg_top50.csv
      results/figures/m7_pseudotime_jaml.{png,tif}
"""
import os
import json
import time
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "rawdata", "scRNA", "GSE131907")
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")
L2 = os.path.join(DATA, "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz")
ANN = os.path.join(DATA, "GSE131907_Lung_Cancer_cell_annotation.txt.gz")

CD4_LINEAGE = ["Naive CD4+ T", "Treg", "CD4+ Th", "CD8+/CD4+ Mixed Th", "Exhausted Tfh"]
CAP_PER_SUBTYPE = 500

FORCE_KEEP = ["JAML", "CXADR", "CCR7", "IL7R", "TCF7", "SELL", "GZMB", "PDCD1",
              "FOXP3", "CD4", "CD3D", "CD3E", "CD8A", "TRDC", "CD69", "HLA-DRA",
              "MKI67", "CXCL13", "IFNG", "GZMA", "PRF1", "ANXA1", "S100A4", "CD44",
              "ITGAE", "BATF", "TOX", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "IL2RA"]
NAIVE_POS = ["CCR7", "IL7R", "TCF7", "SELL"]
EFF_NEG = ["GZMB", "PDCD1"]

TOP_HVG = 2000
MIN_DETECT_FRAC = 0.05
SEED = 20260930


def get_selected_cells():
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
    sel = pd.concat(picks, ignore_index=True)

    # ⚠️ 文件是 .gz，读表头必须用 gzip.open（用普通 open 会读到二进制乱码 → 对齐为 0）
    import gzip
    with gzip.open(L2, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    cells = header[1:]
    pos = {c: i for i, c in enumerate(cells)}
    sel["_col"] = sel["_key"].map(lambda k: pos.get(k, -1))
    matched = int((sel["_col"] >= 0).sum())
    print(f"[选细胞] {len(sel)} 个；在矩阵中命中 {matched}")
    if matched == 0:
        raise RuntimeError("barcode 对齐为 0 —— 检查对齐键（应为注释的 Index 列）")
    sel = sel[sel["_col"] >= 0].reset_index(drop=True)
    print(sel["_sub"].value_counts().to_string())
    return sel


def main():
    os.makedirs(FIG, exist_ok=True)
    print("=" * 70)
    print("M7 拟时序：JAML 沿 CD4⁺T 分化轨迹")
    print("=" * 70)

    sel = get_selected_cells()
    # ⚠️ 顺序一致性：pandas 按**文件列顺序**返回，即 col_idx 升序。
    #    必须把 sel 也按 _col 升序排好，否则 keys 与 X 的列会错位（静默错误）。
    sel = sel.sort_values("_col").reset_index(drop=True)
    col_idx = sel["_col"].to_numpy()
    keys = sel["_key"].tolist()
    sub_of = dict(zip(sel["_key"], sel["_sub"]))
    assert list(col_idx) == sorted(col_idx), "col_idx 必须升序"
    assert len(keys) == len(set(keys)), "细胞键重复"
    n_cells = len(keys)
    min_det = max(3, int(MIN_DETECT_FRAC * n_cells))
    print(f"[阈值] 检出 ≥ {min_det} / {n_cells} 细胞的基因进入候选")

    # ---------------- 单趟抽取：检出率 + 均值 + 方差 + 保留向量 ----------------
    print("\n[抽取] pandas usecols 分块读取（约 10 min）...")
    t0 = time.time()
    # ⚠️ usecols 用的是**文件列位置**（含首列 Index），故需 +1；
    #    col_idx 是 cells 数组内的 0 基下标。
    usecols = [0] + [int(c) + 1 for c in col_idx]
    kept, stats = {}, []
    nline = 0
    reader = pd.read_csv(L2, sep="\t", usecols=usecols, index_col=0,
                         chunksize=3000, engine="c", low_memory=True)
    for ci, chunk in enumerate(reader):
        arr = chunk.to_numpy(dtype=np.float32)          # 行=基因，列=选中细胞
        names = list(chunk.index)
        det = (arr > 0).sum(axis=1)
        mu = arr.mean(axis=1)
        var = arr.var(axis=1)
        for i, g in enumerate(names):
            nline += 1
            if g in FORCE_KEEP or det[i] >= min_det:
                stats.append((g, mu[i], var[i], int(det[i])))
                kept[g] = arr[i].copy()
        if (ci + 1) % 2 == 0:
            print(f"  ...已读 {nline} 基因，已保留 {len(kept)}，用时 {time.time()-t0:.0f}s",
                  flush=True)
    print(f"[抽取] 完成：读取 {nline} 个基因，保留 {len(kept)} 个；"
          f"用时 {time.time()-t0:.0f}s")

    if "JAML" not in kept:
        raise RuntimeError("JAML 未进入保留集 —— 检查基因符号或在 FORCE_KEEP 中确认")

    st_df = pd.DataFrame(stats, columns=["gene", "mean", "var", "detect"])
    st_df = st_df.sort_values("var", ascending=False)
    st_df.head(50).to_csv(os.path.join(RES, "m7_hvg_top50.csv"), index=False)

    hv = [g for g in st_df["gene"].head(TOP_HVG) if g in kept]
    genes = list(dict.fromkeys(hv + [g for g in FORCE_KEEP if g in kept]))
    X = np.stack([kept[g] for g in genes], axis=0)       # 基因 × 细胞
    del kept
    print(f"[矩阵] {X.shape[0]} 基因 × {X.shape[1]} 细胞；"
          f"占用 {X.nbytes/1048576:.0f} MB")

    # ---------------- scanpy 流程 ----------------
    print("\n[scanpy] scale → PCA → neighbors → diffmap → DPT ...")
    import scanpy as sc
    adata = sc.AnnData(X.T.copy())
    adata.obs_names = keys
    adata.var_names = genes
    adata.obs["subtype"] = pd.Categorical([sub_of[k] for k in keys],
                                          categories=CD4_LINEAGE, ordered=True)

    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=30, svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.diffmap(adata, n_comps=15)

    def z(g):
        v = np.asarray(adata[:, g].X).ravel().astype(float)
        sd = v.std()
        return (v - v.mean()) / (sd if sd > 0 else 1.0)

    score = np.zeros(adata.n_obs)
    for g in NAIVE_POS:
        if g in adata.var_names:
            score += z(g)
    for g in EFF_NEG:
        if g in adata.var_names:
            score -= z(g)
    is_naive = (adata.obs["subtype"].astype(str) == "Naive CD4+ T").values
    cand_root = np.where(is_naive)[0]
    root = int(cand_root[np.argmax(score[cand_root])])
    print(f"[root] 第 {root} 个细胞；亚型 {adata.obs['subtype'].iloc[root]}；"
          f"naive 评分 {score[root]:.3f}")
    adata.uns["iroot"] = root
    sc.tl.dpt(adata)

    pt = np.asarray(adata.obs["dpt_pseudotime"], dtype=float)
    df = pd.DataFrame({
        "cell": adata.obs_names,
        "subtype": adata.obs["subtype"].astype(str),
        "dpt_pseudotime": pt,
        "JAML": np.asarray(adata[:, "JAML"].X).ravel(),
        "CXADR": np.asarray(adata[:, "CXADR"].X).ravel() if "CXADR" in adata.var_names else np.nan,
        "naive_score": score,
    })
    df.to_csv(os.path.join(RES, "m7_pseudotime_cells.csv"), index=False)

    order = df.groupby("subtype")["dpt_pseudotime"].median().sort_values()
    print("\n=== 各亚型在拟时序上的中位数位置 ===")
    print(order.to_string())

    # ---------------- JAML 沿拟时序分箱 ----------------
    v = df[np.isfinite(df["dpt_pseudotime"])].copy()
    v["bin"] = pd.qcut(v["dpt_pseudotime"], 20, labels=False, duplicates="drop")
    binned = (v.groupby("bin")
                .agg(pt_mid=("dpt_pseudotime", "median"),
                     JAML=("JAML", "mean"),
                     JAML_sem=("JAML", "sem"),
                     n=("JAML", "size"))
                .reset_index(drop=True))
    binned.to_csv(os.path.join(RES, "m7_pseudotime_jaml_binned.csv"), index=False)
    print("\n=== JAML 沿拟时序的分箱均值（20 箱，前 5 / 后 5）===")
    print(binned.round(4).head(5).to_string(index=False))
    print("...")
    print(binned.round(4).tail(5).to_string(index=False))

    from scipy import stats as st
    rho, pval = st.spearmanr(v["dpt_pseudotime"], v["JAML"])
    print(f"\n[相关] 拟时序 vs JAML  Spearman rho = {rho:.4f}, P = {pval:.3g}, n = {len(v)}")

    summ = {
        "n_cells": int(len(df)), "n_genes": int(len(genes)),
        "cap_per_subtype": CAP_PER_SUBTYPE,
        "subtypes": CD4_LINEAGE,
        "median_pseudotime_by_subtype": {k: float(x) for k, x in order.items()},
        "spearman_pseudotime_vs_JAML": {"rho": float(rho), "p": float(pval), "n": int(len(v))},
        "jaml_mean_by_subtype": {k: float(x) for k, x in df.groupby("subtype")["JAML"].mean().items()},
    }
    with open(os.path.join(RES, "m7_pseudotime_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=2)

    # ---------------- 图 ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.6))
    ax = axes[0]
    data = [df.loc[df["subtype"] == s, "dpt_pseudotime"].dropna().values for s in CD4_LINEAGE]
    data = [d[np.isfinite(d)] for d in data]
    bp = ax.boxplot(data, tick_labels=CD4_LINEAGE, patch_artist=True, widths=0.6,
                    medianprops=dict(color="#333", lw=1.1),
                    flierprops=dict(marker=".", markersize=2, alpha=0.35))
    cmap = plt.get_cmap("YlOrRd")
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor(cmap(0.25 + 0.6 * i / max(1, len(data) - 1)))
        b.set_edgecolor("#333")
    ax.set_ylabel("Diffusion pseudotime")
    ax.set_title("CD4\u207a T subsets along pseudotime", fontsize=9)
    ax.tick_params(axis="x", rotation=32, labelsize=7.5)

    ax = axes[1]
    ax.errorbar(binned["pt_mid"], binned["JAML"], yerr=binned["JAML_sem"],
                marker="o", ms=3.4, lw=1.3, color="#D73027", ecolor="#999", capsize=2)
    ax.set_xlabel("Diffusion pseudotime")
    ax.set_ylabel("JAML (scaled log2TPM)")
    ax.set_title(f"JAML along the trajectory\nSpearman \u03c1 = {rho:.3f}, P = {pval:.2g}", fontsize=9)

    ax = axes[2]
    rngp = np.random.default_rng(1)
    for s in CD4_LINEAGE:
        d = df.loc[df["subtype"] == s, "JAML"].dropna()
        if len(d):
            ax.scatter(rngp.normal(CD4_LINEAGE.index(s), 0.08, len(d)), d,
                       s=3, alpha=0.18, color="#555")
    mns = [df.loc[df["subtype"] == s, "JAML"].mean() for s in CD4_LINEAGE]
    ax.plot(range(len(CD4_LINEAGE)), mns, marker="s", ms=4.5, lw=1.4, color="#D73027")
    ax.set_xticks(range(len(CD4_LINEAGE)))
    ax.set_xticklabels(CD4_LINEAGE, rotation=32, ha="right", fontsize=7.5)
    ax.set_ylabel("JAML (scaled log2TPM)")
    ax.set_title("JAML by subset", fontsize=9)

    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "m7_pseudotime_jaml.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "m7_pseudotime_jaml.tif"), dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("\n[图] results/figures/m7_pseudotime_jaml.png / .tif")
    print("\nDONE")


if __name__ == "__main__":
    main()
