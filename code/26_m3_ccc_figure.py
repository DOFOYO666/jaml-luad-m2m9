# -*- coding: utf-8 -*-
"""M3 图：JAML(发送方) → CXADR(接收方) 通讯强度热图（GSE131907，原发灶/正常肺）。"""
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


m = pd.read_csv(os.path.join(RES, "m3_ccc_matrix_JAML_to_CXADR.csv"), index_col=0)
senders = ["Exhausted CD8+ T", "Exhausted Tfh", "CD4+ Th", "CD8+/CD4+ Mixed Th",
           "Naive CD8+ T", "Naive CD4+ T", "Treg", "CD8 low T", "Cytotoxic CD8+ T",
           "Alveolar Mac", "CD1c+ DCs", "mo-Mac"]
targets = ["AT1", "AT2", "Club", "Ciliated", "Malignant cells", "Mesothelial cells",
           "Tumor ECs", "Alveolar Mac"]
senders = [s for s in senders if s in m.index]
targets = [t for t in targets if t in m.columns]
sub = m.loc[senders, targets]
bg = float(np.median(m.values))

fig, ax = plt.subplots(figsize=(0.62 * len(targets) + 3.4, 0.36 * len(senders) + 2.0))
im = ax.imshow(sub.values, cmap="RdYlBu_r", aspect="auto", vmin=0, vmax=np.nanmax(sub.values))
ax.set_xticks(range(len(targets))); ax.set_xticklabels(targets, rotation=40, ha="right", fontsize=8)
ax.set_yticks(range(len(senders))); ax.set_yticklabels(senders, fontsize=8)
for i in range(sub.shape[0]):
    for j in range(sub.shape[1]):
        v = sub.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.6,
                color="white" if v > 0.65 * np.nanmax(sub.values) else "#222")
ax.set_xlabel("Receiver (CXADR expression)", fontsize=9)
ax.set_ylabel("Sender (JAML expression)", fontsize=9)
ax.set_title("Ligand–receptor communication: JAML → CXADR\nGSE131907 (LUAD primary/normal lung, log1p-CPM group means)",
             fontsize=9.5)
cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02); cb.set_label("lr_means", fontsize=8)
ax.text(1.02, -0.16, f"whole-matrix median = {bg:.2f}", transform=ax.transAxes,
        fontsize=7.5, color="#555")
save(fig, "m3_ccc_JAML_to_CXADR_heatmap")

print(sub.round(3).to_string())
