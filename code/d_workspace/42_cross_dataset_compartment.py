# -*- coding: utf-8 -*-
"""跨数据集汇总：JAML / CXADR 的区室分布（发现集 vs 独立验证集）

发现集   GSE131907（LUAD 组织，208,506 细胞；归一化 log2TPM 面板）
验证集   GSE127465（Zilionis 2019 NSCLC，54,773 细胞；normalized counts 面板）

两套数据的细胞类型命名不同，故按**区室**归并（加权按细胞数）：
  Myeloid/DC、T cells、B lineage、NK cells、Epithelial、Endothelial、Fibroblasts、Mast

⚠️ 阳性细胞的定义：两数据集均为"该基因在该细胞中计数 > 0"（GSE131907 用 log2TPM>0，
   等价于 UMI>0；GSE127465 用 normalized counts > 0）。归一化方式不同，故**只比较
   同一数据集内部与跨数据集的方向/量级**，不把百分比差异当作生物学差异。

输出：results/cross_dataset_compartment_pct.csv
      results/figures/cross_dataset_compartment.{png,tif}
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# GSE131907（发现集）cell_type -> 区室
MAP_131907 = {
    "Myeloid cells": "Myeloid/DC",
    "T lymphocytes": "T cells",
    "B lymphocytes": "B lineage",
    "NK cells": "NK cells",
    "Epithelial cells": "Epithelial",
    "Endothelial cells": "Endothelial",
    "Fibroblasts": "Fibroblasts",
    "MAST cells": "Mast",
    # Undetermined / Oligodendrocytes 排除（后者为脑转移样本成分）
}

# GSE127465（验证集）Major cell type -> 区室
MAP_127465 = {
    "tMoMacDC": "Myeloid/DC", "tpDC": "Myeloid/DC", "bpDC": "Myeloid/DC",
    "bMonocytes": "Myeloid/DC", "tNeutrophils": "Myeloid/DC",
    "bNeutrophils": "Myeloid/DC", "bMyeloid precursor-like": "Myeloid/DC",
    "tT cells": "T cells", "bT cells": "T cells",
    "tB cells": "B lineage", "bB cells": "B lineage",
    "tPlasma cells": "B lineage", "bPlasma cells": "B lineage",
    "tNK cells": "NK cells", "bNK cells": "NK cells",
    "Type I cells": "Epithelial", "Type II cells": "Epithelial",
    "Club cells": "Epithelial", "Ciliated cells": "Epithelial",
    "Endothelial cells": "Endothelial", "Fibroblasts": "Fibroblasts",
    "tMast cells": "Mast",
}

ORDER = ["Myeloid/DC", "Mast", "B lineage", "T cells", "NK cells",
         "Fibroblasts", "Endothelial", "Epithelial"]


def weighted(df, gcol, pcol, ncol, gmap, label):
    d = df.copy()
    d["compartment"] = d[gcol].map(gmap)
    d = d[d["compartment"].notna()]
    n_excl = int(df[~df[gcol].isin(gmap)].groupby(gcol)[ncol].sum().sum())
    out = []
    for c, sub in d.groupby("compartment"):
        w = sub[ncol].values.astype(float)
        out.append({"compartment": c, "n_cells": int(w.sum()),
                    "pct_positive": float(np.average(sub[pcol].values, weights=w)),
                    "mean_expr": float(np.average(sub["mean_expr"].values, weights=w))
                    if "mean_expr" in sub.columns else np.nan})
    r = pd.DataFrame(out)
    r["dataset"] = label
    print(f"[{label}] 归入区室 {int(r['n_cells'].sum()):,} 细胞；"
          f"未归入 {n_excl:,}")
    return r


def main():
    print("=" * 70)
    print("跨数据集区室比对：GSE131907（发现集）vs GSE127465（验证集）")
    print("=" * 70)

    a = pd.read_csv(os.path.join(RES, "m3_gse131907_l2tpm_panel_expression.csv"))
    b = pd.read_csv(os.path.join(RES, "m8_gse127465_panel_by_celltype.csv"))
    b = b[b["level"] == "Major"].copy()

    # ⚠️ 单位必须统一：GSE131907 的 `pos_frac` 是 0–1 的**分数**
    #（如 Myeloid 0.469973 = 47.0%），而 GSE127465 的 `pct_positive` 是 0–100。
    # 曾因混用导致发现集数值被整体缩小 100 倍（Myeloid 47% 被写成 0.47%）。
    a["pos_frac"] = a["pos_frac"] * 100.0
    print(f"[单位] GSE131907 pos_frac × 100 → 百分数（示例 Myeloid "
          f"{a[(a['gene']=='JAML') & (a['cell_type']=='Myeloid cells')]['pos_frac'].iloc[0]:.2f}%）")

    rows = []
    for gene in ["JAML", "CXADR"]:
        ra = weighted(a[a["gene"] == gene], "cell_type", "pos_frac", "n",
                      MAP_131907, "GSE131907")
        rb = weighted(b[b["gene"] == gene], "group", "pct_positive", "n_cells",
                      MAP_127465, "GSE127465")
        t = pd.concat([ra, rb], ignore_index=True)
        t["gene"] = gene
        rows.append(t)
    out = pd.concat(rows, ignore_index=True)
    out["compartment"] = pd.Categorical(out["compartment"], categories=ORDER, ordered=True)
    out = out.sort_values(["gene", "compartment"])
    out.to_csv(os.path.join(RES, "cross_dataset_compartment_pct.csv"), index=False)

    print("\n=== 区室阳性率（%）===")
    for gene in ["JAML", "CXADR"]:
        piv = (out[out["gene"] == gene]
               .pivot_table(index="compartment", columns="dataset", values="pct_positive",
                            observed=False))
        print(f"\n--- {gene} ---")
        print(piv.round(2).to_string())

    # 上皮/免疫比值（区室分离的量化）
    print("\n=== 区室分离量化：Epithelial / Myeloid 比值 ===")
    for gene in ["JAML", "CXADR"]:
        d = out[out["gene"] == gene]
        for ds in ["GSE131907", "GSE127465"]:
            s = d[d["dataset"] == ds].set_index("compartment")["pct_positive"]
            ep = float(s.get("Epithelial", np.nan))
            my = float(s.get("Myeloid/DC", np.nan))
            print(f"  {gene:5s} {ds:11s} 上皮 {ep:6.2f}%  髓系/DC {my:6.2f}%  "
                  f"比值 {ep/my if my else np.nan:5.2f}")

    # ---------------- 图 ----------------
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.8), sharey=True)
    handles = None
    for ax, gene, col in zip(axes, ["JAML", "CXADR"], ["#D73027", "#4575B4"]):
        d = out[out["gene"] == gene]
        piv = (d.pivot_table(index="compartment", columns="dataset",
                             values="pct_positive", observed=False)
               .reindex(ORDER))
        y = np.arange(len(piv))
        h = 0.38
        ax.barh(y + h / 2, piv["GSE131907"].values, height=h, color=col, alpha=0.95,
                edgecolor="#333", linewidth=0.5, label="GSE131907 (discovery, LUAD)")
        ax.barh(y - h / 2, piv["GSE127465"].values, height=h, color=col, alpha=0.45,
                edgecolor="#333", linewidth=0.5, hatch="//",
                label="GSE127465 (replication, NSCLC)")
        ax.set_yticks(y)
        ax.set_yticklabels(piv.index, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Positive cells (%)")
        ax.set_title(f"{gene}", fontsize=10.5)
        ax.spines[["top", "right"]].set_visible(False)
        if handles is None:
            handles, labels = ax.get_legend_handles_labels()
    fig.suptitle("JAML\u2013CXADR compartmental separation: discovery vs independent replication",
                 fontsize=10, y=0.98)
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=[0, 0.07, 1, 0.94])
    for ext, kw in [("png", {}), ("tif", dict(pil_kwargs={"compression": "tiff_lzw"}))]:
        fig.savefig(os.path.join(FIG, f"cross_dataset_compartment.{ext}"), dpi=300, **kw)
    plt.close(fig)
    print("\n[图] results/figures/cross_dataset_compartment.png / .tif")
    print("[输出] results/cross_dataset_compartment_pct.csv")
    print("DONE")


if __name__ == "__main__":
    main()
