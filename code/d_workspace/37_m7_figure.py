# jtm-fontscale x1.60 (JTM 170 mm legibility)
# -*- coding: utf-8 -*-
"""
37_m7_figure.py — M7 定稿图（从已保存的 CSV 重绘，避免重跑 14 min 的矩阵抽取）

改动理由：原第 3 子图用散点展示 JAML，但 CD4⁺T 中 **>75% 细胞 JAML 表达为 0**
（零膨胀），散点图不可读、量纲也与前两图不一致。改为"均值 ± SEM 条形图"，
并把 y 轴统一标注为 scaled（z）单位。同时补标 JAML 阳性率，便于如实报告。

输入：results/m7_pseudotime_cells.csv
      results/m7_pseudotime_jaml_binned.csv
输出：results/figures/m7_pseudotime_jaml.{png,tif}
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

SUB = ["Naive CD4+ T", "Treg", "CD4+ Th", "CD8+/CD4+ Mixed Th", "Exhausted Tfh"]


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".tif"), dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"[图] results/figures/{name}.png / .tif")


def main():
    df = pd.read_csv(os.path.join(RES, "m7_pseudotime_cells.csv"))
    binned = pd.read_csv(os.path.join(RES, "m7_pseudotime_jaml_binned.csv"))
    with open(os.path.join(RES, "m7_pseudotime_summary.json"), encoding="utf-8") as f:
        summ = json.load(f)
    rho = summ["spearman_pseudotime_vs_JAML"]["rho"]
    pval = summ["spearman_pseudotime_vs_JAML"]["p"]
    # "阳性"定义为 scaled 值高于全部为零的基线（scaled 0 均值对应的常数）
    base = float(df["JAML"].quantile(0.25))
    df["pos"] = df["JAML"] > base + 1e-6

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.0))

    # --- (1) 亚型在拟时序上的位置 ---
    ax = axes[0]
    data = [df.loc[df["subtype"] == s, "dpt_pseudotime"].dropna().values for s in SUB]
    bp = ax.boxplot(data, tick_labels=SUB, patch_artist=True, widths=0.6,
                    medianprops=dict(color="#333", lw=1.1),
                    flierprops=dict(marker=".", markersize=2, alpha=0.35))
    cmap = plt.get_cmap("YlOrRd")
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor(cmap(0.25 + 0.6 * i / max(1, len(data) - 1)))
        b.set_edgecolor("#333")
    ax.set_ylabel("Diffusion pseudotime")
    ax.set_title("CD4\u207a T subsets along pseudotime", fontsize=14.4)
    ax.tick_params(axis="x", rotation=30, labelsize=7.5)

    # --- (2) JAML 沿拟时序的分箱均值 ---
    ax = axes[1]
    ax.errorbar(binned["pt_mid"], binned["JAML"], yerr=binned["JAML_sem"],
                marker="o", ms=3.4, lw=1.3, color="#D73027", ecolor="#999", capsize=2)
    ax.axhline(0, ls=":", color="#bbb", lw=0.8)
    ax.set_xlabel("Diffusion pseudotime")
    ax.set_ylabel("JAML, scaled (z units)")
    ax.set_title(f"JAML rises along the trajectory\nSpearman \u03c1 = {rho:.3f}, "
                 f"P = {pval:.1e}, n = {summ['spearman_pseudotime_vs_JAML']['n']:,}", fontsize=14.4)

    # --- (3) 亚型层面的均值 ± SEM + 阳性率（替代不可读的散点）---
    ax = axes[2]
    means, sems, poss = [], [], []
    for s in SUB:
        d = df.loc[df["subtype"] == s, "JAML"]
        means.append(d.mean())
        sems.append(d.sem())
        poss.append(100 * df.loc[df["subtype"] == s, "pos"].mean())
    cols = ["#FEE0B6", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"]
    ax.bar(range(len(SUB)), means, yerr=sems, color=cols, edgecolor="#333",
           linewidth=0.5, capsize=2.4, error_kw=dict(lw=0.7, ecolor="#555"))
    ax.axhline(0, ls=":", color="#bbb", lw=0.8)
    for i, (mn, po) in enumerate(zip(means, poss)):
        ax.text(i, mn + sems[i] + 0.012, f"{mn:+.2f}\n{po:.0f}% pos", ha="center",
                va="bottom", fontsize=10.6)
    ax.set_xticks(range(len(SUB)))
    ax.set_xticklabels(SUB, rotation=30, ha="right", fontsize=12.0)
    ax.set_ylabel("JAML, scaled (z units)")
    ax.set_ylim(min(-0.28, min(means) - 0.06), max(means) + 0.10)
    ax.set_title("JAML by subset (mean \u00b1 SEM)", fontsize=14.4)

    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "m7_pseudotime_jaml")

    # 补充：阳性率的说明
    print("\n=== 各亚型 JAML 阳性率（scaled 值高于零膨胀基线）===")
    for s, mn, po in zip(SUB, means, poss):
        n = (df["subtype"] == s).sum()
        print(f"  {s:22s} n={n:4d}  均值 {mn:+.3f}  阳性率 {po:5.1f}%")
    print(f"\n[注] 全表 JAML 中位数 {df['JAML'].median():.3f}，"
          f"75 百分位 {df['JAML'].quantile(0.75):.3f} → 零膨胀，"
          f"轨迹信号由少数阳性细胞驱动，报告时须说明。")
    print("\nDONE")


if __name__ == "__main__":
    main()
