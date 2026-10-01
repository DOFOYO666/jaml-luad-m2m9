# -*- coding: utf-8 -*-
"""任务 3：M8（decoupleR TF 活性）与 M4（CellOracle TF 扰动）的交叉比较。

两个模块回答的是**不同**问题，交叉比较时必须说清：

- **M8**：JAML 阳性 vs 阴性 CD4⁺T 细胞之间，各 TF 的**活性差异**（观察性关联）。
  指标：`diff` = mean(TF 活性 | JAML 阳性) − mean(TF 活性 | JAML 阴性)，另有 FDR。
- **M4**：把某个 TF 敲除后，JAML 表达的**变化量**（计算性干预）。
  指标：`JAML_mean_abs_delta`（另有有符号值，判断方向）。

两者一致 → 互相支持；不一致 → 需要解释（观察性关联可能来自共变或反向因果）。
**注意**：M8 用的是 7,757 基因的宽集，M4 用的是 3,013 基因集，
两个模块的可评估 TF 集合不同，交集才是可比的。

用法：python m4_m8_cross.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd

D_ROOT = r"D:\workbuddy工作空间\JAML深度研究"
RES = os.path.join(D_ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

# 后缀：用于区分不同轮的 M4 结果（ext3013 无后缀；ext3074 用 _ext2）
SUF = ""
if "--suf" in sys.argv:
    SUF = sys.argv[sys.argv.index("--suf") + 1]
print(f"[参数] 结果后缀 = '{SUF}'")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def save(fig, name):
    for ext, kw in [("png", {}), ("tif", dict(pil_kwargs={"compression": "tiff_lzw"}))]:
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, **kw)
    plt.close(fig)
    print(f"[图] results/figures/{name}.png / .tif")


# ---------------- 读入 ----------------
m4 = pd.read_csv(os.path.join(RES, f"m4_tf_perturbation_summary{SUF}.csv"))
m4 = m4[~m4["is_randomized_control"]][["TF", "JAML_mean_abs_delta",
                                       "JAML_mean_delta_signed",
                                       "mean_abs_delta_all_genes"]].copy()

m8 = pd.read_csv(os.path.join(RES, "m8_tf_group_diff.csv"))
m8j = m8[m8["gene"] == "JAML"].copy()
print("M8 中 JAML 的行数:", len(m8j), "| split:", m8j["split"].unique().tolist())

# 主指标取 positive_vs_negative（阳性 vs 阴性细胞，与 M8 的显著 TF 数一致）
m8p = m8j[m8j["split"] == "positive_vs_negative"][["tf", "diff", "fdr", "n_high", "n_low"]]
m8p = m8p.rename(columns={"tf": "TF", "diff": "M8_activity_diff",
                          "fdr": "M8_fdr"})

both = m4.merge(m8p, on="TF", how="inner")
print(f"\nM4 可评估 TF: {len(m4)}；M8 JAML 行: {len(m8p)}；**交集: {len(both)}**")
missing_in_m4 = sorted(set(m8p["TF"]) - set(m4["TF"]))
print(f"M8 有、M4 未能扰动（GRN/矩阵限制）的 TF 数: {len(missing_in_m4)}")

if len(both) >= 3:
    r_p = both["JAML_mean_abs_delta"].corr(both["M8_activity_diff"], method="pearson")
    r_s = both["JAML_mean_abs_delta"].corr(both["M8_activity_diff"], method="spearman")
    print(f"\n相关性（|ΔJAML| vs M8 活性差）：Pearson r = {r_p:.3f}，Spearman ρ = {r_s:.3f}")
    # 有符号方向的一致性：M4 的 signed Δ 与 M8 的活性差同号？
    sign_agree = int((np.sign(both["JAML_mean_delta_signed"]) ==
                      np.sign(both["M8_activity_diff"])).sum())
    print(f"方向同号数：{sign_agree} / {len(both)}"
          f"（预期随机约 {len(both)/2:.0f}）")

both = both.sort_values("JAML_mean_abs_delta", ascending=False)
both.to_csv(os.path.join(RES, f"m4_m8_cross_tf{SUF}.csv"), index=False)
print("\n=== 交集 TF 按 |ΔJAML|（M4）降序 ===")
print(both.round(5).to_string(index=False))

# ---------------- 图：散点 ----------------
fig, ax = plt.subplots(figsize=(6.0, 4.8))
x = both["M8_activity_diff"]
y = both["JAML_mean_abs_delta"]
ax.scatter(x, y, s=40, color="#7B8FA1", zorder=3)
for _, r in both.iterrows():
    ax.annotate(r["TF"], (r["M8_activity_diff"], r["JAML_mean_abs_delta"]),
                textcoords="offset points", xytext=(5, 3), fontsize=7)
ax.set_xlabel("M8: TF activity difference (JAML-positive − negative)")
ax.set_ylabel(r"M4: mean $|\Delta$ JAML$|$ after TF knockout")
ax.set_title("Cross-check: TF activity association (M8) vs\nvirtual-knockout effect on JAML (M4)",
             fontsize=9)
ax.axhline(0, color="#bbb", lw=0.8)
ax.axvline(0, color="#bbb", lw=0.8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, f"m4_m8_cross_tf{SUF}")

# ---------------- 汇总 ----------------
out = {
    "n_m4_evaluable": int(len(m4)),
    "n_m8_jaml_rows": int(len(m8p)),
    "n_intersection": int(len(both)),
    "n_in_m8_not_m4": int(len(missing_in_m4)),
    "tf_in_m8_not_m4": missing_in_m4[:40],
    "pearson_r_absdelta_vs_m8diff": None if len(both) < 3 else float(r_p),
    "spearman_rho": None if len(both) < 3 else float(r_s),
    "sign_agreement": None if len(both) < 3 else f"{sign_agree}/{len(both)}",
    "note": ("M8 为观察性关联（TF 活性差），M4 为计算性干预（敲除后 JAML 变化）；"
             "两者一致才构成互证，不一致需考虑共变/反向因果。M8 用 7757 基因宽集，"
             "M4 用 3003/3013/3074 基因集（视轮次），交集才是可比集合。"),
}
with open(os.path.join(RES, f"m4_m8_cross_summary{SUF}.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print(f"\n[输出] results/m4_m8_cross_tf{SUF}.csv / m4_m8_cross_summary{SUF}.json")
print("DONE")
