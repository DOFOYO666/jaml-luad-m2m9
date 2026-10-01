# -*- coding: utf-8 -*-
"""审稿意见 M3 的执行：把"与自身随机对照最大值之比"升级为"合并经验零分布 + 超越概率"。

审稿人指出：以 5 次打乱重复的最大值作为噪声上界，是一个方差很大的估计量；
且"比值"形式掩盖了绝对效应量。本脚本在不重跑 CellOracle 的前提下，
用已归档的每次打乱读数（m4_noise_control_classify.csv 的 jaml_rand_values 列）做：

  1. 合并 9 个 TF × 5 次 = 45 个打乱读数的经验零分布（均值/SD/中位数/p95/最大值）；
  2. 每个 TF 的经验超越概率：p_emp = (#{打乱读数 >= 观测效应} + 1) / (N + 1)；
  3. 最小可检出效应（MDE）= 该零分布的最大值——低于此值的效应无法与噪声区分；
  4. 三种上界口径下的判定一致性（自身最大值 / 合并 p95 / 合并最大值）。

输出：results/ 下 CSV + JSON + 图（PNG + 300 dpi LZW TIFF）。
用法：python 53_m4_null_pooled.py
"""
import io
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = r"D:\workbuddy工作空间\JAML深度研究\results"
OUT_CSV = os.path.join(D, "m4_null_pooled.csv")
OUT_JSON = os.path.join(D, "m4_null_pooled_summary.json")
FIGDIRS = [os.path.join(D, "figures"),
           r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化\results\figures"]

plt.rcParams["axes.unicode_minus"] = False


def main():
    src = os.path.join(D, "m4_noise_control_classify.csv")
    d = pd.read_csv(src)
    d = d[d["TF"].notna()].copy()

    # 展开每次打乱读数
    draws = []
    for _, r in d.iterrows():
        vals = [float(x) for x in str(r["jaml_rand_values"]).split(";") if x.strip() != ""]
        for v in vals:
            draws.append({"TF": r["TF"], "rand": v})
    dl = pd.DataFrame(draws)
    N = len(dl)
    null_max = float(dl["rand"].max())
    null_mean = float(dl["rand"].mean())
    null_sd = float(dl["rand"].std(ddof=1))
    null_p95 = float(np.percentile(dl["rand"], 95))
    null_med = float(dl["rand"].median())

    rows = []
    for _, r in d.iterrows():
        obs = float(r["jaml_real"])
        n_ge = int((dl["rand"] >= obs).sum())
        p_emp = (n_ge + 1) / (N + 1)
        rows.append({
            "TF": r["TF"],
            "observed_abs_delta_JAML": obs,
            "own_rand_max": float(r["jaml_rand_max"]),
            "own_ratio": obs / float(r["jaml_rand_max"]) if r["jaml_rand_max"] else np.nan,
            "n_randomized_draws_at_or_above": n_ge,
            "empirical_p_one_sided": p_emp,
            "degenerate_control": bool(r["degenerate_control"]),
            "above_pooled_p95": obs > null_p95,
            "above_pooled_max": obs > null_max,
        })
    out = pd.DataFrame(rows).sort_values("observed_abs_delta_JAML", ascending=False)
    out["mde_pooled_max"] = null_max
    out.to_csv(OUT_CSV, index=False)

    summary = {
        "n_randomized_draws_pooled": N,
        "n_TF": int(d["TF"].nunique()),
        "draws_per_TF": int(N / max(1, d["TF"].nunique())),
        "null_mean": null_mean, "null_sd": null_sd, "null_median": null_med,
        "null_p95": null_p95, "null_max": null_max,
        "minimum_detectable_effect": null_max,
        "MDE_note": ("低于该值的扰动效应无法与合并零分布区分；"
                     "该口径比各 TF 自身最大值更保守（后者随 TF 数增加而变小）"),
        "verdict_by_three_ceilings": {
            "above_own_max": out.loc[out["own_ratio"] > 1, "TF"].tolist(),
            "above_pooled_p95": out.loc[out["above_pooled_p95"], "TF"].tolist(),
            "above_pooled_max": out.loc[out["above_pooled_max"], "TF"].tolist(),
        },
    }
    with io.open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print("=== 合并经验零分布（全部打乱读数）===")
    print(f"  打乱读数总数 N = {N}（{summary['n_TF']} 个 TF × {summary['draws_per_TF']} 次）")
    print(f"  均值 {null_mean:.3e}  SD {null_sd:.3e}  中位 {null_med:.3e}")
    print(f"  p95 {null_p95:.3e}  最大 {null_max:.3e}  ← 最小可检出效应 (MDE)")
    print()
    print("=== 各 TF：三种上界口径下的判定 ===")
    print(out[["TF", "observed_abs_delta_JAML", "own_rand_max", "own_ratio",
               "n_randomized_draws_at_or_above", "empirical_p_one_sided",
               "above_pooled_p95", "above_pooled_max"]].round(6).to_string(index=False))
    print()
    print("  超过自身最大值的 TF：", summary["verdict_by_three_ceilings"]["above_own_max"])
    print("  超过合并 p95 的 TF：", summary["verdict_by_three_ceilings"]["above_pooled_p95"])
    print("  超过合并最大值的 TF：", summary["verdict_by_three_ceilings"]["above_pooled_max"])

    # ---------------- 图 ----------------
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6))
    ax = axes[0]
    ax.hist(dl["rand"], bins=22, color="#B0BEC5", edgecolor="#455A64", linewidth=0.5)
    ax.axvline(null_max, color="#D73027", ls="--", lw=1.1,
               label=f"pooled null max = {null_max:.2e}\n(MDE)")
    id2 = float(out.loc[out["TF"] == "ID2", "observed_abs_delta_JAML"].iloc[0])
    gata = float(out.loc[out["TF"] == "GATA3", "observed_abs_delta_JAML"].iloc[0])
    ax.axvline(id2, color="#1B5E20", lw=1.4, label=f"ID2 knockout = {id2:.2e}")
    ax.axvline(gata, color="#6A1B9A", lw=1.2, ls=":", label=f"GATA3 = {gata:.2e}")
    ax.set_xlabel(r"mean $|\Delta$JAML$|$ under randomized GRN")
    ax.set_ylabel(f"count (N = {N})")
    ax.set_title("Pooled empirical null across all randomized draws", fontsize=9)
    ax.legend(frameon=False, fontsize=6.6, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    s = out.copy().iloc[::-1]
    y = np.arange(len(s))
    ax.barh(y, s["observed_abs_delta_JAML"], height=0.6, color="#C0392B",
            edgecolor="#333", linewidth=0.4, label="real GRN knockout")
    ax.axvline(null_max, color="#D73027", ls="--", lw=1.1, label="pooled null max (MDE)")
    ax.axvline(null_p95, color="#EF9A9A", ls=":", lw=1.0, label="pooled null p95")
    ax.set_yticks(y)
    ax.set_yticklabels(s["TF"], fontsize=7.4)
    ax.set_xlabel(r"mean $|\Delta$JAML$|$")
    ax.set_title("Effect vs pooled null ceiling", fontsize=9)
    ax.legend(frameon=False, fontsize=6.6, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for d_ in FIGDIRS:
        os.makedirs(d_, exist_ok=True)
        fig.savefig(os.path.join(d_, "m4_null_pooled.png"), dpi=300)
        fig.savefig(os.path.join(d_, "m4_null_pooled.tiff"), dpi=300,
                    pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("\n[输出] m4_null_pooled.csv / _summary.json / figures/m4_null_pooled.{png,tiff}")
    print("DONE")


if __name__ == "__main__":
    main()
