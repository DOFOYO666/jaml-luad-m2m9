# -*- coding: utf-8 -*-
"""
33_m3_figures_normalized.py — M3 定稿图（改用官方 log2TPM 归一化均值）

为什么重做：scripts/23 用的是 raw UMI 的 mean_counts，会被文库大小混淆；
现用官方归一化矩阵（GSE131907 normalized log2TPM）的组均值，跨细胞类型可比。

两版互验：raw UMI 与 log2TPM 的阳性率在 260 个（基因 × 细胞类型）配对上
Pearson r = 1.0000、最大绝对差 0.0000 个百分点（scripts/31）——即解析无误。

产出：
  figures/m3_jaml_cxadr_normalized.{png,tif}      主图：JAML/CXADR 区室分离
  figures/m3_jaml_cd4_subtypes.{png,tif}          CD4 亚群的 JAML 梯度
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

ORDER = ["Myeloid cells", "T lymphocytes", "NK cells", "B lymphocytes",
         "MAST cells", "Endothelial cells", "Fibroblasts",
         "Epithelial cells", "Oligodendrocytes"]
LABEL = {"Myeloid cells": "Myeloid", "T lymphocytes": "T lymph.", "NK cells": "NK",
         "B lymphocytes": "B lymph.", "MAST cells": "MAST", "Endothelial cells": "Endothelial",
         "Fibroblasts": "Fibroblast", "Epithelial cells": "Epithelial",
         "Oligodendrocytes": "Oligodendro."}


def save(fig, name):
    fig.savefig(os.path.join(FIG, f"{name}.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, f"{name}.tif"), dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"[图] results/figures/{name}.png / .tif")


def main():
    c = pd.read_csv(os.path.join(RES, "m3_gse131907_l2tpm_panel_expression.csv"))
    s = pd.read_csv(os.path.join(RES, "m3_gse131907_l2tpm_panel_subtype.csv"))

    def series(gene):
        d = c[c["gene"] == gene].set_index("cell_type")["mean_log2tpm"]
        return d.reindex(ORDER)

    jaml, cxadr = series("JAML"), series("CXADR")

    # ---------------- 主图 ----------------
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7),
                             gridspec_kw={"width_ratios": [1.45, 1]})
    y = np.arange(len(ORDER))
    ax = axes[0]
    h = 0.38
    ax.barh(y + h / 2, cxadr.values, height=h, color="#4575B4",
            edgecolor="#333", linewidth=0.5, label="CXADR (receptor)")
    ax.barh(y - h / 2, jaml.values, height=h, color="#D73027",
            edgecolor="#333", linewidth=0.5, label="JAML (ligand)")
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[k] for k in ORDER], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("Mean normalized expression, log2(TPM+1)")
    ax.set_title("Compartmental separation of the JAML\u2013CXADR axis\nLUAD tissue, GSE131907 (n = 208,506 cells)",
                 fontsize=9)
    ax.set_xlim(0, 0.92)
    ax.legend(frameon=False, fontsize=8, loc="center right", bbox_to_anchor=(1.0, 0.62))
    ax.spines[["top", "right"]].set_visible(False)
    for yi, v in zip(y, jaml.values):
        ax.text(v + 0.015, yi - h / 2, f"{v:.2f}", va="center", fontsize=6.6, color="#7a1a12")
    for yi, v in zip(y, cxadr.values):
        ax.text(v + 0.015, yi + h / 2, f"{v:.2f}", va="center", fontsize=6.6, color="#1d3f6e")

    # 阳性率对照
    ax = axes[1]
    pj = c[c["gene"] == "JAML"].set_index("cell_type")["pos_frac"].reindex(ORDER) * 100
    pc = c[c["gene"] == "CXADR"].set_index("cell_type")["pos_frac"].reindex(ORDER) * 100
    ax.scatter(pj.values, pc.values, s=46, color="#333", zorder=3)
    _off = {"Epithelial cells": (6, -9), "Oligodendrocytes": (7, 4),
            "Myeloid cells": (-42, 6), "T lymphocytes": (6, 3)}
    for k, xv, yv in zip(ORDER, pj.values, pc.values):
        ax.annotate(LABEL[k], (xv, yv), textcoords="offset points",
                    xytext=_off.get(k, (5, 4)), fontsize=6.8)
    ax.axhline(5, ls=":", color="#aaa", lw=0.9)
    ax.axvline(5, ls=":", color="#aaa", lw=0.9)
    ax.set_xlabel("JAML-positive cells (%)")
    ax.set_ylabel("CXADR-positive cells (%)")
    ax.set_xlim(-3, 58)
    ax.set_ylim(-3, 48)
    ax.set_title("Positivity is mutually exclusive", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "m3_jaml_cxadr_normalized")

    # ---------------- CD4 亚群梯度 ----------------
    cd4_order = ["Naive CD4+ T", "Treg", "CD8+/CD4+ Mixed Th", "CD4+ Th", "Exhausted Tfh"]
    js = s[s["gene"] == "JAML"].set_index("cell_subtype")
    vals = np.array([js["mean_log2tpm"].get(k, np.nan) for k in cd4_order])
    pos = np.array([js["pos_frac"].get(k, np.nan) for k in cd4_order]) * 100
    nn = np.array([js["n"].get(k, 0) for k in cd4_order])

    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    cols = ["#FEE0B6", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"]
    bars = ax.bar(range(len(cd4_order)), vals, color=cols, edgecolor="#333", linewidth=0.5)
    for i, (b, p, n) in enumerate(zip(bars, pos, nn)):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.015,
                f"{vals[i]:.2f}\n{p:.1f}% pos\nn={int(n):,}", ha="center", va="bottom", fontsize=6.6)
    ax.set_xticks(range(len(cd4_order)))
    ax.set_xticklabels(cd4_order, rotation=22, ha="right", fontsize=8)
    ax.set_ylabel("JAML mean expression, log2(TPM+1)")
    ax.set_ylim(0, max(vals) * 1.42)
    ax.set_title("JAML increases along the CD4\u207a T activation axis\n(GSE131907)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "m3_jaml_cd4_subtypes")
    print("\nDONE")


if __name__ == "__main__":
    main()
