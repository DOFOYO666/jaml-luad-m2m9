# -*- coding: utf-8 -*-
"""M4 结果图：TF 扰动对 JAML 表达的影响。

严格按用户约定：图片输出 300 dpi 的 TIFF（LZW 压缩）+ PNG 预览。
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D_ROOT = r"D:\workbuddy工作空间\JAML深度研究"
RES = os.path.join(D_ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)


def save(fig, name):
    for ext, kw in [("png", {}), ("tif", dict(pil_kwargs={"compression": "tiff_lzw"}))]:
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, **kw)
    plt.close(fig)
    print(f"[图] results/figures/{name}.png / .tif")


s = pd.read_csv(os.path.join(RES, "m4_tf_perturbation_summary.csv"))
rand = s[s["is_randomized_control"]]
real = s[~s["is_randomized_control"]].sort_values("JAML_mean_abs_delta", ascending=False)
base_max = float(rand["JAML_mean_abs_delta"].max())
base_mean = float(rand["JAML_mean_abs_delta"].mean())

# ---------------- 图 1：|ΔJAML| 排序条形图 ----------------
top = real.head(16).iloc[::-1]
fig, ax = plt.subplots(figsize=(6.4, 4.6))
colors = ["#C0392B" if t == "ID2" else "#7B8FA1" for t in top["TF"]]
ax.barh(top["TF"], top["JAML_mean_abs_delta"], color=colors, height=0.7)
ax.axvline(base_max, color="#D35400", ls="--", lw=1.0,
           label=f"randomized-GRN control (max = {base_max:.1e})")
ax.set_xlabel(r"mean $|\Delta$ JAML$|$ after TF knockout")
ax.set_title("Virtual knockout of candidate TFs: effect on JAML expression\n"
             "GSE131907 CD4$^+$T lineage, CellOracle (n = 2,842 cells)",
             fontsize=9)
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)
for y, v in enumerate(top["JAML_mean_abs_delta"]):
    ax.text(v + top["JAML_mean_abs_delta"].max() * 0.015, y, f"{v:.2e}",
            va="center", fontsize=6.4, color="#333")
fig.tight_layout()
save(fig, "m4_tf_ko_effect_on_jaml")

# ---------------- 图 2：全局效应 vs JAML 效应 ----------------
fig, ax = plt.subplots(figsize=(5.6, 4.6))
ax.scatter(real["mean_abs_delta_all_genes"], real["JAML_mean_abs_delta"],
           s=34, color="#7B8FA1", zorder=3, label="TF knockout")
ax.scatter(rand["mean_abs_delta_all_genes"], rand["JAML_mean_abs_delta"],
           s=44, marker="X", color="#D35400", zorder=4, label="randomized GRN (control)")
ax.axhline(base_max, color="#D35400", ls=":", lw=0.9)
# 标注
for _, r in real.iterrows():
    if r["JAML_mean_abs_delta"] > 0.0009 or r["mean_abs_delta_all_genes"] > 0.004:
        ax.annotate(r["TF"], (r["mean_abs_delta_all_genes"], r["JAML_mean_abs_delta"]),
                    textcoords="offset points", xytext=(5, 3), fontsize=7.2)
ax.set_xlabel(r"mean $|\Delta|$ across all genes (global effect)")
ax.set_ylabel(r"mean $|\Delta$ JAML$|$")
ax.set_title("TF knockout: global effect vs effect on JAML\n"
             "ID2 is the only TF whose JAML effect far exceeds the control range",
             fontsize=9)
ax.legend(frameon=False, fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "m4_tf_ko_global_vs_jaml")

# ---------------- 图 3：ID2 扰动下关键基因的 |Δ| ----------------
kt = pd.read_csv(os.path.join(RES, "m4_key_genes_shift.csv"))
col_ko = [c for c in kt.columns if c.endswith("_KO")][0]
kt2 = kt.sort_values(col_ko, ascending=False)
fig, ax = plt.subplots(figsize=(6.4, 3.8))
x = np.arange(len(kt2))
ax.bar(x, kt2[col_ko], color="#C0392B", width=0.66)
ax.set_xticks(x)
ax.set_xticklabels(kt2["gene"], rotation=45, ha="right", fontsize=7.5)
ax.set_ylabel(r"mean $|\Delta|$ after ID2 knockout")
ax.set_title("ID2 knockout: gene-level shifts (JAML ranks first)\n"
             "GSE131907 CD4$^+$T lineage", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
for xi, v in enumerate(kt2[col_ko]):
    ax.text(xi, v + kt2[col_ko].max() * 0.02, f"{v:.2e}", ha="center",
            fontsize=6.2, rotation=90, color="#333")
fig.tight_layout()
save(fig, "m4_id2_ko_gene_shifts")

print(f"\n噪声基线：mean={base_mean:.2e}, max={base_max:.2e}")
print(f"ID2 / 基线(max) = {real.iloc[0]['JAML_mean_abs_delta'] / base_max:.1f} 倍")
print(f"ID2 / 第二名   = {real.iloc[0]['JAML_mean_abs_delta'] / real.iloc[1]['JAML_mean_abs_delta']:.1f} 倍")
print("DONE")
