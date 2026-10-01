# jtm-fontscale x1.85 (JTM 170 mm legibility)
# -*- coding: utf-8 -*-
"""M4 三档比对 + 噪声标定。

三档：
  A = bagging=3,  基因集 3013（ext3013）—— 最初的测试档
  B = bagging=20, 基因集 3013（ext3013）—— 统计稳健档
  C = bagging=20, 基因集 3074（ext3074）—— 补齐 61 个 TF 后（闭合审计缺口）

两个不同的问题，必须分开谈：
  A↔B  衡量 **bagging 次数** 的影响（统计学稳健性）
  B↔C  衡量 **基因集扩大** 的影响（方法学副作用）

⚠️ 噪声基线不能用各轮自带的 random-GRN 对照：
   * `m4_celloracle_run.py` 取 `TF_CANDIDATES[:3]` 作打乱基准；在 C 轮里前三位恰是
     MYB/RORC/EOMES（检出率约 1%），其真扰动本身就无效应 → 基线被人为压到 0。
   * 且 CellOracle 的打乱系数矩阵默认缓存（seed=123），"3 次重复"并非独立抽样。
   因此本脚本改用 `m4_noise_control.csv` 的**配对**信噪比（同一 TF、5 次独立打乱）
   作为唯一权威的噪声标定，并以其中的最大随机读数作为全局噪声上界参考线。

用法：python m4_round_compare.py
输出：results/m4_round_compare.csv / m4_round_compare_summary.json
      results/figures/m4_round_compare.{png,tif}
"""
import json
import os

import numpy as np
import pandas as pd

D_ROOT = r"D:\workbuddy工作空间\JAML深度研究"
RES = os.path.join(D_ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FILES = {
    "A_bag3_ext3013": "m4_tf_perturbation_summary_bag3.csv",
    "B_bag20_ext3013": "m4_tf_perturbation_summary_bag20.csv",
    "C_bag20_ext3074": "m4_tf_perturbation_summary_ext2.csv",
}
REF = ["ID2", "GATA3", "STAT3", "MAF", "IRF1", "STAT1", "IRF4", "NFKB1", "BATF", "TCF7"]
NEW_TFS = ['MYB', 'RORC', 'EOMES', 'IRF6', 'IRF5', 'ZBTB7B', 'GFI1', 'BCL6', 'IRF8',
           'SPI1', 'CEBPA', 'ATF3', 'JDP2', 'MAFG', 'KLF4', 'KLF10', 'SP2', 'SP3', 'SP4',
           'EGR2', 'EGR3', 'RORB', 'NR1H3', 'STAT2', 'STAT5A', 'STAT6', 'SMAD2', 'SMAD3',
           'SMAD4', 'RUNX1', 'RUNX2', 'TCF3', 'TCF4', 'TCF12', 'TCF7L2', 'FOXO3', 'FOXO4',
           'IKZF4', 'ZEB1', 'SNAI1', 'SNAI2', 'TWIST1', 'ID1', 'ID4', 'HES1', 'HEY1',
           'NR2F6', 'NR2F1', 'NR2F2', 'ETS2', 'ELF4', 'FLI1', 'ERG', 'ETV6', 'ETV5',
           'GABPA', 'BHLHE41', 'ATF6', 'CREB1', 'ATF2', 'NFIL3']


def save(fig, name):
    for ext, kw in [("png", {}), ("tif", dict(pil_kwargs={"compression": "tiff_lzw"}))]:
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, **kw)
    plt.close(fig)
    print(f"[图] results/figures/{name}.png / .tif")


def load(fn):
    p = os.path.join(RES, fn)
    if not os.path.isfile(p):
        print(f"[跳过] 不存在: {fn}")
        return None
    d = pd.read_csv(p)
    rt = d[d["is_randomized_control"]]
    real = d[~d["is_randomized_control"]][
        ["TF", "JAML_mean_abs_delta", "JAML_mean_delta_signed"]].copy()
    return real.set_index("TF"), rt["JAML_mean_abs_delta"]


tabs = {}
for k, fn in FILES.items():
    r = load(fn)
    if r is not None:
        tabs[k] = r
keys = list(tabs.keys())
print(f"载入 {len(keys)} 档: {keys}")

print("\n=== 各轮自带的 random-GRN 对照（仅记录，不作为判据）===")
runctrl = {}
for k, (real, noise) in tabs.items():
    runctrl[k] = {"n": int(len(noise)), "mean": float(noise.mean()),
                  "max": float(noise.max())}
    print(f"  {k:18s} n={len(noise)}  均值 {noise.mean():.3e}  最大 {noise.max():.3e}")

# ---------------- 权威噪声：配对对照 ----------------
paired = pd.read_csv(os.path.join(RES, "m4_noise_control.csv"))
noise_ceiling = float(paired["jaml_rand_max"].max())
print(f"\n=== 配对噪声对照（权威）===  {len(paired)} 个 TF，每个 5 次独立打乱")
print(paired[["TF", "jaml_real", "jaml_rand_mean", "jaml_rand_max",
              "snr_vs_randmax"]].round(6).to_string(index=False))
print(f"全局噪声上界（配对打乱的最大读数）= {noise_ceiling:.6f}")

# ---------------- 表达量（解释效应为何为零）----------------
import anndata as ad
import scipy.sparse as sp
_a = ad.read_h5ad(os.path.join(RES, "CD4T_celloracle_ext2.h5ad"))
vn = list(map(str, _a.var_names))
RC = _a.layers["raw_count"]
expr = {}
for g in set(REF) | set(NEW_TFS):
    if g not in vn:
        continue
    col = RC[:, vn.index(g)]
    v = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
    expr[g] = 100.0 * float((v > 0).mean())
del _a, RC

# ---------------- 合并表 ----------------
union = sorted(set().union(*[set(t[0].index) for t in tabs.values()]))
cmp = pd.DataFrame(index=union)
for k in keys:
    cmp[f"absJAML_{k}"] = tabs[k][0]["JAML_mean_abs_delta"]
    cmp[f"signed_{k}"] = tabs[k][0]["JAML_mean_delta_signed"]
    cmp[f"rank_{k}"] = tabs[k][0]["JAML_mean_abs_delta"].rank(ascending=False)
cmp["pos_pct_ext2"] = pd.Series(expr)
ps = paired.set_index("TF")
cmp["paired_SNR"] = ps["snr_vs_randmax"]
cmp["paired_rand_max"] = ps["jaml_rand_max"]
cmp["is_new_tf"] = cmp.index.isin(NEW_TFS)
cmp = cmp.sort_values(f"absJAML_{keys[-1]}", ascending=False)
cmp.to_csv(os.path.join(RES, "m4_round_compare.csv"))

print("\n=== 参照 TF 的三档表现 ===")
cols = [c for c in cmp.columns if c.startswith("absJAML_") or c.startswith("rank_")] + \
       ["pos_pct_ext2", "paired_SNR"]
print(cmp.loc[[g for g in REF if g in cmp.index], cols].round(5).to_string())

# ---------------- 档间一致性 ----------------
stat = {}
print("\n=== 档间一致性 ===")
for i in range(len(keys) - 1):
    a, b = keys[i], keys[i + 1]
    common = [g for g in cmp.index
              if np.isfinite(cmp.loc[g, f"absJAML_{a}"])
              and np.isfinite(cmp.loc[g, f"absJAML_{b}"])]
    x = cmp.loc[common, f"absJAML_{a}"].astype(float)
    y = cmp.loc[common, f"absJAML_{b}"].astype(float)
    rp = float(np.corrcoef(x, y)[0, 1])
    rs = float(x.corr(y, method="spearman"))
    same = int((np.sign(x) == np.sign(y)).sum())
    stat[f"{a}_vs_{b}"] = {"n_common": len(common), "pearson_r": round(rp, 4),
                           "spearman_rho": round(rs, 4)}
    print(f"  {a} vs {b}: n={len(common)}  Pearson r = {rp:.3f}  Spearman ρ = {rs:.3f}")

# 排名稳定性（参照 TF）
print("\n=== 参照 TF 在三档中的排名 ===")
rk = cmp.loc[[g for g in REF if g in cmp.index],
             [f"rank_{k}" for k in keys]]
print(rk.apply(pd.to_numeric, errors="coerce").round(0).astype("Int64").to_string())

# ---------------- 新增 TF（仅在 C 中）----------------
new_in_c = [g for g in cmp.index if g not in tabs[keys[-2]][0].index]
sub = cmp.loc[new_in_c, [f"absJAML_{keys[-1]}", "pos_pct_ext2"]].sort_values(
    f"absJAML_{keys[-1]}", ascending=False)
print(f"\n=== 新增 TF（{len(new_in_c)} 个）：|ΔJAML| 与表达量 ===")
print(sub.round(6).to_string())
n_above = int((sub[f"absJAML_{keys[-1]}"] > noise_ceiling).sum())
print(f"  超出配对噪声上界（{noise_ceiling:.6f}）的: {n_above} / {len(new_in_c)}")

# ---------------- 图 ----------------
fig, axes = plt.subplots(1, 2, figsize=(10.8, 6.2),
                         gridspec_kw={"width_ratios": [1.05, 1]})

ax = axes[0]
a, b = keys[-2], keys[-1]
common = [g for g in cmp.index
          if np.isfinite(cmp.loc[g, f"absJAML_{a}"]) and np.isfinite(cmp.loc[g, f"absJAML_{b}"])]
x = cmp.loc[common, f"absJAML_{a}"].astype(float)
y = cmp.loc[common, f"absJAML_{b}"].astype(float)
isnew = cmp.loc[common, "is_new_tf"].astype(bool)
ax.scatter(x[~isnew], y[~isnew], s=34, color="#7B8FA1", zorder=3,
           label="TFs present in both runs")
ax.scatter(x[isnew], y[isnew], s=34, color="#E67E22", marker="^", zorder=3,
           label="TFs added in run C")
lim = max(x.max(), y.max()) * 1.15
ax.plot([0, lim], [0, lim], ls="--", color="#ccc", lw=0.8)
for g in [g for g in REF if g in x.index]:
    ax.annotate(g, (x[g], y[g]), textcoords="offset points", xytext=(5, 3), fontsize=13.0)
ax.set_xlabel(f"$|\\Delta$JAML$|$   {a}")
ax.set_ylabel(f"$|\\Delta$JAML$|$   {b}")
ax.set_title(f"(B) vs (C): gene-set enlargement\nSpearman ρ = "
             f"{stat[f'{a}_vs_{b}']['spearman_rho']:.3f}  (n = {len(common)})", fontsize=15.9)
ax.legend(frameon=False, fontsize=12.6, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
s2 = sub.iloc[::-1]
colors = ["#C0392B" if v > noise_ceiling else "#95A5A6"
          for v in s2[f"absJAML_{keys[-1]}"]]
ax.barh(range(len(s2)), s2[f"absJAML_{keys[-1]}"].values, color=colors, height=0.72)
ax.axvline(noise_ceiling, ls="--", color="#2C3E50", lw=0.9)
ax.text(noise_ceiling, len(s2) - 1, " paired randomization\n noise ceiling",
        color="#2C3E50", fontsize=11.8, va="top")
ax.set_yticks(range(len(s2)))
ax.set_yticklabels(s2.index, fontsize=10.0)
ax.tick_params(axis="y", pad=1)
ax.set_xlabel(r"$|\Delta$JAML$|$  (run C, bagging = 20)")
ax.set_title(f"Newly added TFs: {n_above}/{len(new_in_c)} exceed the noise ceiling",
             fontsize=15.9)
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
save(fig, "m4_round_compare")

out = {"run_own_random_controls": runctrl,
       "paired_noise_ceiling": noise_ceiling,
       "pairwise_agreement": stat,
       "n_new_tf": len(new_in_c),
       "n_new_tf_above_noise": n_above,
       "note": ("各轮自带的 random-GRN 基线不可作判据（基准 TF 表达量差异 + 打乱矩阵被缓存）；"
                "权威噪声来自 m4_noise_control.csv 的配对对照。")}
with open(os.path.join(RES, "m4_round_compare_summary.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n[输出] results/m4_round_compare.csv / m4_round_compare_summary.json")
print("DONE")
