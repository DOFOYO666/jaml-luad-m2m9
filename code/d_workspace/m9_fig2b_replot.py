# -*- coding: utf-8 -*-
"""Re-render the GSE127465 replication panel (JTM Figure 2b) at the JTM page width.

The original plotting script re-reads the 44-million-entry Matrix Market file, which is not
present on this machine; this script reproduces the *same* panel from the landed summary
table results/m8_gse127465_panel_by_celltype.csv, with the fonts scaled so every label
stays >= 6 pt when the figure is printed at 170 mm.

The seventeen-group selection, the key-compartment list, the colours, the ordering and the
bar annotations are identical to scripts/38_gse127465_replication.py, and the values drawn
are asserted against the numbers already reported in the manuscript.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")

SEL = ["tMoMacDC", "tpDC", "tNeutrophils", "bNeutrophils", "tT cells", "bT cells",
       "tNK cells", "bNK cells", "tB cells", "tPlasma cells", "tMast cells",
       "Fibroblasts", "Endothelial cells", "Type II cells", "Club cells",
       "Ciliated cells", "Type I cells"]
SPEC = [("JAML", [("tMoMacDC", "Myeloid/DC (t)"), ("tT cells", "T cells (t)"),
                  ("tNeutrophils", "Neutrophils (t)"), ("tB cells", "B cells (t)"),
                  ("tNK cells", "NK cells (t)"),
                  ("Type II cells", "Type II epi."), ("Club cells", "Club")]),
        ("CXADR", [("Type II cells", "Type II epi."), ("Club cells", "Club"),
                   ("Ciliated cells", "Ciliated"), ("Type I cells", "Type I epi."),
                   ("tMoMacDC", "Myeloid/DC (t)"), ("tT cells", "T cells (t)")])]
# the manuscript's own numbers, used as assertions
EXPECT = {("JAML", "tMoMacDC"): 23.1, ("JAML", "tT cells"): 7.4,
          ("JAML", "tNeutrophils"): 9.9, ("JAML", "tB cells"): 0.6,
          ("JAML", "tNK cells"): 4.9, ("JAML", "Type II cells"): 2.2,
          ("JAML", "Club cells"): 0.0, ("CXADR", "Type II cells"): 18.9,
          ("CXADR", "Club cells"): 25.0, ("CXADR", "Ciliated cells"): 14.3,
          ("CXADR", "Type I cells"): 21.7, ("CXADR", "tMoMacDC"): 1.0,
          ("CXADR", "tT cells"): 0.2}

data = pd.read_csv(os.path.join(RES, "m8_gse127465_panel_by_celltype.csv"))
maj = data[data["level"] == "Major"]
piv = maj.pivot_table(index="group", columns="gene", values="pct_positive")
n_by_group = maj.groupby("group")["n_cells"].first()
sel = [c for c in SEL if c in set(piv.index)]
assert len(sel) == 17, "expected 17 groups, found %d" % len(sel)
piv = piv.reindex(sel)

# statement 1: the replication header numbers
assert int(n_by_group.loc["tMoMacDC"]) == 9372
assert abs(piv.loc["tMoMacDC", "JAML"] - 23.14) < 0.01
assert abs(piv.loc["Type II cells", "CXADR"] - 18.86) < 0.01

fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.0), gridspec_kw={"width_ratios": [1.1, 1]})

ax = axes[0]
y = np.arange(len(piv))
h = 0.38
ax.barh(y + h / 2, piv["CXADR"].values, height=h, color="#4575B4",
        edgecolor="#333", linewidth=0.5, label="CXADR (receptor)")
ax.barh(y - h / 2, piv["JAML"].values, height=h, color="#D73027",
        edgecolor="#333", linewidth=0.5, label="JAML (ligand)")
ax.set_yticks(y)
ax.set_yticklabels([f"{c}  (n={int(n_by_group.loc[c]):,})" for c in piv.index], fontsize=11.6)
ax.invert_yaxis()
ax.set_xlabel("Positive cells (%)")
ax.set_title("Independent replication: JAML\u2013CXADR compartmental separation\n"
             "GSE127465 (Zilionis 2019, NSCLC, n = 54,773 cells)", fontsize=13.6)
ax.legend(frameon=False, fontsize=11.6, loc="center right")
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
labels, vals, cols = [], [], []
for gene, items in SPEC:
    for c, lab in items:
        d = maj[(maj["gene"] == gene) & (maj["group"] == c)]
        if len(d) == 0:
            continue
        v = float(d["pct_positive"].iloc[0])
        key = (gene, c)
        if key in EXPECT:
            assert abs(v - EXPECT[key]) <= 0.05, "%s %s: %.2f vs expected %.1f" % (gene, c, v, EXPECT[key])
        labels.append(lab)
        vals.append(v)
        cols.append("#D73027" if gene == "JAML" else "#4575B4")
yy = np.arange(len(labels))
ax.barh(yy, vals, color=cols, edgecolor="#333", linewidth=0.5)
for i, vv in enumerate(vals):
    ax.text(vv + 0.5, i, f"{vv:.1f}", va="center", fontsize=9.9)
ax.set_yticks(yy)
ax.set_yticklabels(labels, fontsize=11.6)
ax.invert_yaxis()
ax.set_xlabel("Positive cells (%)")
ax.set_title("Key compartments (red = JAML, blue = CXADR)", fontsize=13.6)
ax.spines[["top", "right"]].set_visible(False)
ax.text(0.98, 0.02, "few epithelial cells in this dataset",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=9.9, color="#777")

fig.tight_layout()
for ext, kw in (("png", {}), ("tif", {"pil_kwargs": {"compression": "tiff_lzw"}})):
    fig.savefig(os.path.join(FIG, "m8_gse127465_jaml_cxadr." + ext), dpi=300,
                bbox_inches="tight", **kw)
print("[re-rendered] results/figures/m8_gse127465_jaml_cxadr.png/.tif")
print("  17 groups, %d key-compartment bars, all values match the manuscript" % len(vals))
