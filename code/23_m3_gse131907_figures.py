# -*- coding: utf-8 -*-
"""M3 图：GSE131907 (LUAD 组织) 中 JAML 与 CXADR 的细胞类型分布。

同时输出 JAML vs CXADR 的「区室分离」散点图，直观展示配体-受体分属免疫侧与上皮侧。
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RES = "results"
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

res = pd.read_csv(os.path.join(RES, "m3_gse131907_panel_expression.csv"))
ct = res[res["level"] == "Cell_type"]


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".tif"), dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  [图] {name}.png / .tif")


plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

order = ["Myeloid cells", "T lymphocytes", "NK cells", "B lymphocytes",
         "Epithelial cells", "Endothelial cells", "Fibroblasts", "MAST cells"]

fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6), sharex=True)
for ax, gene, color in zip(axes, ["JAML", "CXADR"], ["#D73027", "#4575B4"]):
    d = ct[ct["gene"] == gene].set_index("group").reindex(order)
    y = np.arange(len(order))
    ax.barh(y, d["pct_positive"], color=color, height=0.68)
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f"{gene}-positive cells (%)")
    ax.set_xlim(0, 70)
    for i, (v, n) in enumerate(zip(d["pct_positive"], d["n_cells"])):
        ax.text(v + 1.2, i, f"{v:.1f}%", va="center", fontsize=7.5)
    ax.set_title(f"{gene}  (n = 208,506 cells)", fontsize=9.5)
fig.suptitle("GSE131907 — ligand (JAML) is immune-restricted; receptor (CXADR) is epithelial",
             fontsize=10, y=1.02)
fig.tight_layout()
save(fig, "m3_gse131907_jaml_cxadr_by_celltype")

# ---- 区室分离散点 ----
j = ct[ct["gene"] == "JAML"].set_index("group")["pct_positive"]
c = ct[ct["gene"] == "CXADR"].set_index("group")["pct_positive"]
both = pd.DataFrame({"JAML": j, "CXADR": c}).reindex(order).dropna()

fig, ax = plt.subplots(figsize=(4.9, 4.1))
ax.scatter(both["JAML"], both["CXADR"], s=64, color="#666", zorder=3, alpha=0.85)
# 左下角 5 个点高度重叠 → 手工错开标签位置（offset points, 水平对齐）
OFFSET = {
    "Epithelial cells": (8, 2, "left"),
    "T lymphocytes": (8, -2, "left"),
    "Myeloid cells": (8, 2, "left"),
    "Endothelial cells": (-8, 14, "right"),
    "Fibroblasts": (8, 8, "left"),
    "MAST cells": (8, -10, "left"),
    "B lymphocytes": (-8, -14, "right"),
    "NK cells": (8, -2, "left"),
}
for g, r in both.iterrows():
    dx, dy, ha = OFFSET.get(g, (7, 4, "left"))
    ax.annotate(g, (r["JAML"], r["CXADR"]), textcoords="offset points",
                xytext=(dx, dy), ha=ha, fontsize=7.5)
ax.axhline(10, ls=":", color="#bbb", lw=0.9); ax.axvline(10, ls=":", color="#bbb", lw=0.9)
ax.text(11, 45, "JAML$^{+}$ immune\nCXADR$^{-}$", fontsize=7.5, color="#999", va="top")
ax.text(0.5, 10.8, "CXADR$^{+}$ epithelial / JAML$^{-}$", fontsize=7.5, color="#999")
ax.set_xlabel("JAML-positive cells (%)"); ax.set_ylabel("CXADR-positive cells (%)")
ax.set_xlim(-3, 58); ax.set_ylim(-3, 48)
ax.set_title("Compartmental separation of the JAML–CXADR axis\n(LUAD tissue, GSE131907)", fontsize=9)
save(fig, "m3_gse131907_jaml_cxadr_separation")

print("\nJAML / CXADR 阳性率（Cell_type 层）:")
print(both.round(2).to_string())
