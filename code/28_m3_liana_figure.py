# -*- coding: utf-8 -*-
"""M3 定稿图：用 **liana 本体**结果绘制 JAML(发送方) → CXADR(接收方) 通讯热图。"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RES, FIG = "results", "results/figures"
os.makedirs(FIG, exist_ok=True)


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".tif"), dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  [图] {name}.png / .tif")


d = pd.read_csv(os.path.join(RES, "m3_liana_JAML_CXADR_directed.csv"))
piv = d.pivot_table(index="source", columns="target", values="lr_means")
bg = float(pd.read_csv(os.path.join(RES, "m3_liana_rank_aggregate.csv"))["lr_means"].median())

senders = ["Exhausted CD8+ T", "Exhausted Tfh", "CD4+ Th", "CD8+/CD4+ Mixed Th",
           "Naive CD8+ T", "Naive CD4+ T", "Treg", "CD8 low T",
           "CD1c+ DCs", "Alveolar Mac", "mo-Mac"]
targets = ["AT1", "AT2", "Club", "Ciliated", "Malignant cells"]
senders = [s for s in senders if s in piv.index]
targets = [t for t in targets if t in piv.columns]
sub = piv.loc[senders, targets]

fig, ax = plt.subplots(figsize=(0.75 * len(targets) + 3.6, 0.40 * len(senders) + 2.0))
im = ax.imshow(sub.values, cmap="RdYlBu_r", aspect="auto", vmin=0, vmax=np.nanmax(sub.values))
ax.set_xticks(range(len(targets)))
ax.set_xticklabels(targets, rotation=35, ha="right", fontsize=8.5)
ax.set_yticks(range(len(senders)))
ax.set_yticklabels(senders, fontsize=8.5)
for i in range(sub.shape[0]):
    for j in range(sub.shape[1]):
        v = sub.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7.4,
                color="white" if v > 0.68 * np.nanmax(sub.values) else "#1a1a1a")
ax.axhline(7.5, color="#333", lw=1.1, ls="--")   # T 细胞 / 髓系 分隔
ax.set_xlabel("Receiver — CXADR-expressing cells", fontsize=9)
ax.set_ylabel("Sender — JAML-expressing cells", fontsize=9)
ax.set_title("JAML → CXADR ligand–receptor communication in LUAD tissue (GSE131907)\n"
             "liana rank_aggregate · consensus resource · log1p-CPM group means", fontsize=9.5)
cb = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02); cb.set_label("lr_means", fontsize=8)
ax.text(1.03, -0.13, f"global lr_means median = {bg:.2f}", transform=ax.transAxes,
        fontsize=7.5, color="#555")
save(fig, "m3_liana_JAML_to_CXADR_heatmap")
print(sub.round(3).to_string())
