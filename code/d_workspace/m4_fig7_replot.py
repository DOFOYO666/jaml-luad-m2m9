# -*- coding: utf-8 -*-
"""Figure 7a 重绘：改用**最终网络配置**（ext3074：3,074 基因，bagging = 20，71 个候选/70 个可扰动）。

为什么换配置（依据全部来自已落盘的结果文件，不重新推断）：
  1. ext3074 是三个配置里唯一非 test mode 的生产运行，基因覆盖面最大（3,074 vs 3,003/3,013），
     不可扰动因子最少（1 个 vs 15 个和 5 个）；bagging = 3 那次是 test_mode=True 的冒烟测试。
  2. **权威噪声校准是在 ext3074 的输入上跑的**（m4_noise_control_summary.json 的 input =
     CD4T_celloracle_ext2.h5ad）：全局噪声上界 9.5945e-3、ID2 配对信噪比 6.55x 均出自该配置。
     于是图内的横线与 Table 4、正文用的是同一个门槛，不再出现"图注写上界、图内画另一个数"。

本脚本**只读**已落盘的结果表重绘，不重跑 CellOracle / decoupleR——重跑会让稿件里的数字轻微漂移。
写图前先跑一致性断言，任何一项不符即中止且不写盘。
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

D_ROOT = r"D:\workbuddy工作空间\JAML深度研究"
RES = os.path.join(D_ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

SUM_EXT2 = os.path.join(RES, "m4_tf_perturbation_summary_ext2.csv")
META_EXT2 = os.path.join(RES, "m4_celloracle_summary_ext2.json")
CLASSIFY = os.path.join(RES, "m4_noise_control_classify.json")
PAIRED = os.path.join(RES, "m4_noise_control_summary.json")

COL = "JAML_mean_abs_delta"
EXPLICIT_XMIN = 1e-10       # 最小非零效应为 1.23e-10，故下限取 1e-10（8.4 个数量级）
ZERO_X = EXPLICIT_XMIN * 1.5   # Δ = 0 的因子画在此处（贴着轴内下限，避免标记被裁掉）
N_ZERO_EXPECTED = 11        # 效应恰为 0 的因子数（网络内无通路／检出率过低）

FAIL = []


def check(label, ok, detail=""):
    print(("  OK   " if ok else "  FAIL ") + label + ("  | " + detail if detail else ""))
    if not ok:
        FAIL.append(label)


def save(fig, name):
    for ext, kw in [("png", {}), ("tif", dict(pil_kwargs={"compression": "tiff_lzw"}))]:
        fig.savefig(os.path.join(FIG, "%s.%s" % (name, ext)), dpi=300, **kw)
    plt.close(fig)
    print("[图] results/figures/%s.png / .tif" % name)


def main():
    s = pd.read_csv(SUM_EXT2)
    rand = s[s["is_randomized_control"] == True]                     # noqa: E712
    real = (s[s["is_randomized_control"] != True]                    # noqa: E712
            .sort_values(COL, ascending=False).reset_index(drop=True))
    meta = json.load(open(META_EXT2, encoding="utf-8"))
    cls = json.load(open(CLASSIFY, encoding="utf-8"))
    paired = json.load(open(PAIRED, encoding="utf-8"))

    ceiling = float(cls["global_noise_ceiling_nondegenerate"])
    paired_rows = {r["TF"]: r for r in paired["results"]}
    id2_real = float(paired_rows["ID2"]["jaml_real"])
    calib = [r["TF"] for r in paired["results"]]

    print("[断言]")
    check("随机对照行 = 3", len(rand) == 3, "n=%d" % len(rand))
    check("可扰动 TF = 70", len(real) == 70, "n=%d" % len(real))
    check("ext3074 元数据一致（3,074 基因 / bagging 20 / ext3074）",
          meta["n_genes"] == 3074 and meta["bagging_number"] == 20
          and meta.get("input_tag") == "ext3074",
          "genes=%s bagging=%s tag=%s" % (meta["n_genes"], meta["bagging_number"],
                                          meta.get("input_tag")))
    check("仅 1 个候选不可扰动", meta["n_TF_skipped"] == 1, "skipped=%s" % meta["n_TF_skipped"])
    check("排序第一 = ID2", real.loc[0, "TF"] == "ID2", "top=%s" % real.loc[0, "TF"])
    check("图内 ID2 值 = 配对对照里的 ID2 真扰动值",
          abs(float(real.loc[0, COL]) - id2_real) < 1e-12,
          "%.15g vs %.15g" % (float(real.loc[0, COL]), id2_real))
    degen = set(cls.get("degenerate_controls", []))
    cand = [float(r["jaml_rand_max"]) for r in paired["results"] if r["TF"] not in degen]
    check("全局噪声上界 = 非退化因子中最大的随机读数",
          abs(ceiling - max(cand)) < 1e-15,
          "ceiling=%.6e max_rand=%.6e (n=%d, degenerate=%s)"
          % (ceiling, max(cand), len(cand), sorted(degen)))
    above = real[real[COL] >= ceiling]
    check("高于上界的因子恰好 1 个，且为 ID2",
          len(above) == 1 and above.iloc[0]["TF"] == "ID2",
          "above=%s" % list(above["TF"]))
    missing = [t for t in calib if t not in set(real["TF"])]
    check("9 个配对校准因子的真扰动值全部在表内", not missing, "missing=%s" % missing)
    check("ID2 = 第二名的 5 倍以上",
          float(real.loc[0, COL]) / float(real.loc[1, COL]) > 5,  # noqa: PLR2004
          "%.2fx" % (float(real.loc[0, COL]) / float(real.loc[1, COL])))
    n_zero = int((real[COL].astype(float) == 0).sum())
    check("效应恰为 0 的因子数 = %d" % N_ZERO_EXPECTED, n_zero == N_ZERO_EXPECTED, "n_zero=%d" % n_zero)
    if FAIL:
        raise SystemExit("断言未通过，未写盘：%s" % FAIL)

    # ---------------- 绘图：两栏棒棒糖图（全部 70 个因子都标名） ----------------
    left = real.iloc[:35]
    right = real.iloc[35:]
    calib_set = set(calib)

    fig, axes = plt.subplots(1, 2, figsize=(6.8, 5.0), sharex=True)
    fig.subplots_adjust(left=0.095, right=0.985, top=0.845, bottom=0.145, wspace=0.46)

    for ax, d, tag in ((axes[0], left, "Ranks 1-35"),
                       (axes[1], right, "Ranks 36-70")):
        dd = d.iloc[::-1].reset_index(drop=True)          # 第一名在顶部
        y = np.arange(len(dd))
        for yi, tf, v in zip(y, dd["TF"], dd[COL].astype(float)):
            if tf == "ID2":
                c, ms, lw = "#C0392B", 4.4, 1.5
            elif tf in calib_set:
                c, ms, lw = "#2C3E50", 3.3, 1.0
            else:
                c, ms, lw = "#9AA5B1", 2.9, 0.85
            if v <= 0:                       # log 轴无法表示 0：画空心标记在轴内下限处
                ax.plot(ZERO_X, yi, "o", ms=ms, mfc="none", mec=c, mew=0.85, zorder=3)
            else:
                ax.hlines(yi, EXPLICIT_XMIN, v, color=c, lw=lw, zorder=2)
                ax.plot(v, yi, "o", ms=ms, color=c, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(dd["TF"], fontsize=5.5)
        ax.set_ylim(-0.9, len(dd) - 0.15)
        ax.set_xscale("log")
        ax.set_xlim(EXPLICIT_XMIN, 2.8e-2)
        ax.axvline(ceiling, color="#D35400", ls="--", lw=1.0, zorder=1)
        ax.set_title(tag, fontsize=7.4, pad=3)
        ax.tick_params(axis="x", labelsize=5.8, length=2)
        ax.tick_params(axis="y", length=1.6)
        ax.grid(axis="x", lw=0.35, color="#DDDDDD", zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.spines["left"].set_linewidth(0.6)
        ax.spines["bottom"].set_linewidth(0.6)
        ax.set_xlabel(r"mean $|\Delta$ JAML$|$  (log scale)", fontsize=6.6)

    # 前三名标数值
    for rank in range(3):
        tf = real.loc[rank, "TF"]
        v = float(real.loc[rank, COL])
        ax = axes[0]
        yi = len(left) - 1 - rank
        ax.text(v * 1.22, yi, "%.2e" % v, va="center", fontsize=5.2,
                color="#C0392B" if tf == "ID2" else "#2C3E50")

    handles = [
        Line2D([], [], color="#C0392B", marker="o", ls="-", lw=1.5, ms=4.4, label="ID2"),
        Line2D([], [], color="#2C3E50", marker="o", ls="-", lw=1.0, ms=3.3,
               label="paired-calibrated factor"),
        Line2D([], [], color="#9AA5B1", marker="o", ls="-", lw=0.85, ms=2.9,
               label="other factor"),
        Line2D([], [], color="#9AA5B1", marker="o", ls="none", ms=3.0, mfc="none",
               label=r"$\Delta$ = 0 exactly (n = %d)" % n_zero),
        Line2D([], [], color="#D35400", ls="--", lw=1.0,
               label="global noise ceiling %.2e" % ceiling),
    ]
    axes[1].legend(handles=handles, frameon=False, fontsize=5.0, loc="lower right",
                   handlelength=1.7, borderpad=0.3, labelspacing=0.32)

    fig.suptitle("Virtual knockout of candidate transcription factors: effect on JAML expression\n"
                 "final network configuration (3,074 genes, bagging = 20) - GSE131907 CD4$^+$T lineage - "
                 "CellOracle 0.20.0 (n = 2,842 cells)",
                 fontsize=7.9, y=0.985)

    save(fig, "m4_tf_ko_effect_on_jaml")

    rec = {
        "figure": "Figure 7a (m4_tf_ko_effect_on_jaml)",
        "configuration": "ext3074 / bagging = 20 / alpha = 10 / n_propagation = 3",
        "source_summary": os.path.basename(SUM_EXT2),
        "source_meta": os.path.basename(META_EXT2),
        "n_genes": meta["n_genes"],
        "n_cells": meta["n_cells"],
        "n_tf_perturbed_done": len(real),
        "n_tf_skipped": meta["n_TF_skipped"],
        "top_tf": real.loc[0, "TF"],
        "top_value": float(real.loc[0, COL]),
        "second_tf": real.loc[1, "TF"],
        "second_value": float(real.loc[1, COL]),
        "min_value": float(real[COL].min()),
        "n_zero_effect": n_zero,
        "min_nonzero_value": float(real[COL].astype(float)[real[COL].astype(float) > 0].min()),
        "global_noise_ceiling": ceiling,
        "ceiling_source": "m4_noise_control_classify.json (paired randomized-network control "
                          "run on CD4T_celloracle_ext2.h5ad)",
        "n_above_ceiling": int(len(above)),
        "n_paired_calibrated": len(calib),
        "note_replaces": "previous version plotted the 3,013-gene configuration and showed only the "
                         "top 16 of 34 perturbed factors while the legend said all were shown",
    }
    with open(os.path.join(FIG, "m4_fig7_replot_record.json"), "w", encoding="utf-8") as fh:
        json.dump(rec, fh, ensure_ascii=False, indent=1)

    print()
    print("[事实表] 配置 %s / %d 细胞 / %d 基因" % (rec["configuration"], meta["n_cells"], meta["n_genes"]))
    print("   可扰动 %d 个，不可扰动 %d 个" % (len(real), meta["n_TF_skipped"]))
    print("   ID2 = %.6e（第二名 %s = %.6e，%.1f 倍）"
          % (rec["top_value"], rec["second_tf"], rec["second_value"],
             rec["top_value"] / rec["second_value"]))
    print("   效应恰为 0 者 %d 个；最小非零效应 = %.3e；全局噪声上界 = %.6e；高于上界者 %d 个"
          % (n_zero, rec["min_nonzero_value"], ceiling, len(above)))
    print("   9 个配对校准因子：%s" % ", ".join(calib))
    print("DONE")


if __name__ == "__main__":
    main()
