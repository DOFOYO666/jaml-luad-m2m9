# -*- coding: utf-8 -*-
"""M2：TCGA-LUAD 临床关联（分期 + 生存 + 免疫 marker 相关性）

相对 03_tcga_survival.py 的改进：
  1. 生存分析改用**正规 Mantel-Cox log-rank 检验** + **一维 Cox 比例风险模型**
     （原脚本用"事件率之比 + 二项 z 检验"近似，方法学上站不住）
  2. 新增分期表达差异（Kruskal-Wallis + 趋势检验）
  3. 新增 JAML 与免疫 marker 基因的相关性（**明确标注为 marker 相关性，
     不是反卷积估计**）

输入：data/luad_survival/TCGA.LUAD.HiSeqV2.gz
      data/luad_survival/TCGA.LUAD.clinicalMatrix.txt
      results/scmr_luad_risk_sig.csv
输出：results/m2_tcga_survival_proper.csv
      results/m2_tcga_stage_jaml.csv
      results/m2_tcga_immune_marker_corr.csv
      results/figures/m2_tcga_*.png / .tif
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "luad_survival")
RES = os.path.join(BASE, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

IMMUNE_MARKERS = {
    "PTPRC": "CD45 (pan-leukocyte)", "CD3D": "CD3D (T cell)", "CD3E": "CD3E (T cell)",
    "CD4": "CD4 (T helper)", "CD8A": "CD8A (cytotoxic T)", "CD8B": "CD8B (cytotoxic T)",
    "FOXP3": "FOXP3 (Treg)", "IL2RA": "CD25 (Treg/activated)",
    "PDCD1": "PD-1 (exhaustion)", "CTLA4": "CTLA-4", "CD274": "PD-L1",
    "GZMB": "granzyme B", "PRF1": "perforin", "GZMA": "granzyme A",
    "NKG7": "NKG7 (NK/cytotoxic)", "KLRD1": "CD94 (NK)", "NCR1": "NKp46 (NK)",
    "CD14": "CD14 (monocyte)", "CD68": "CD68 (macrophage)", "ITGAM": "CD11b",
    "LYZ": "lysozyme (myeloid)", "FCGR3A": "CD16 (NK/mono)",
    "MS4A1": "CD20 (B cell)", "CD19": "CD19 (B cell)", "CD79A": "CD79A (B cell)",
    "MZB1": "MZB1 (plasma)", "CXCL9": "CXCL9 (IFN-γ-inducible)",
    "IFNG": "IFN-γ", "HLA-DRA": "HLA-DR", "ITGAX": "CD11c (DC/myeloid)",
}


# 基因符号别名（TCGA/Xena HiSeqV2 使用较老的注释）
GENE_ALIAS = {
    "JAML": ["JAML", "AMICA1"],
    "HYKK": ["HYKK", "AGPHD1"],
    "ZNRD1ASP": ["ZNRD1ASP", "ZNRD1-AS1", "ZNRD1AS1"],
    "CTC-490E21.14": ["CTC-490E21.14"],
    "RP11-514O12.4": ["RP11-514O12.4"],
    "RP1-167A14.3": ["RP1-167A14.3"],
}


def resolve_gene(sym, index):
    """返回数据矩阵中实际使用的符号；找不到返回 None。"""
    for cand in GENE_ALIAS.get(sym, [sym]):
        if cand in index:
            return cand
    if sym in index:
        return sym
    return None


# ---------- 统计工具 ----------
def logrank(time, event, group):
    """Mantel-Cox log-rank 检验。group 为 0/1。返回 (chi2, p)。

    正确做法：χ² = (Σ(O₁−E₁))² / ΣV，即先累积观测-期望差与方差，再取比值。
    （逐时点累加 (O−E)²/V 是错误写法，零假设下会产生严重反保守的 χ²。）
    """
    time = np.asarray(time, float)
    event = np.asarray(event, int)
    group = np.asarray(group, int)
    O1 = E1 = V = 0.0
    for t in np.unique(time[event == 1]):
        at = time >= t
        n = at.sum()
        n1 = (at & (group == 1)).sum()
        d = ((time == t) & (event == 1)).sum()
        d1 = ((time == t) & (event == 1) & (group == 1)).sum()
        if n <= 1:
            continue
        O1 += d1
        E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V <= 0:
        return np.nan, np.nan
    chi2 = (O1 - E1) ** 2 / V
    return float(chi2), float(stats.chi2.sf(chi2, 1))


def cox_hr(time, event, x):
    """一维 Cox 比例风险模型（Breslow 处理并列时间），返回 (HR, se, p)。"""
    time = np.asarray(time, float)
    event = np.asarray(event, int)
    x = np.asarray(x, float)
    x = (x - x.mean()) / x.std()

    def negll(b):
        b = np.atleast_1d(b)[0]
        ll = 0.0
        for t in np.unique(time[event == 1]):
            risk = time >= t
            d = (time == t) & (event == 1)
            eta = b * x[risk]
            m = eta.max()
            ll += b * x[d].sum() - (d.sum() * (m + np.log(np.exp(eta - m).sum())))
        return -ll

    r = optimize.minimize(negll, x0=[0.0], method="BFGS")
    b = r.x[0]
    # 数值 Hessian
    h = 1e-5
    f0 = negll([b])
    f1 = negll([b + h])
    f2 = negll([b - h])
    d2 = (f1 - 2 * f0 + f2) / h ** 2
    se = np.sqrt(1 / d2) if d2 > 0 else np.nan
    z = b / se if se > 0 else np.nan
    p = 2 * stats.norm.sf(abs(z)) if not np.isnan(z) else np.nan
    return float(np.exp(b)), float(se), float(p)


def save_fig(fig, name):
    png = os.path.join(FIG, name + ".png")
    tif = os.path.join(FIG, name + ".tif")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(tif, dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  [图] {name}.png / .tif")


# ---------- 载入 ----------
expr = pd.read_csv(os.path.join(DATA, "TCGA.LUAD.HiSeqV2.gz"), sep="\t",
                   compression="gzip", index_col=0)
clin = pd.read_csv(os.path.join(DATA, "TCGA.LUAD.clinicalMatrix.txt"), sep="\t", index_col=0)
print(f"[TCGA] 表达 {expr.shape}；临床 {clin.shape}")

vs = clin["vital_status"].astype(str).str.upper()
death = vs.str.contains("DECEASED", na=False)
dd = pd.to_numeric(clin["days_to_death"], errors="coerce")
dlf = pd.to_numeric(clin["days_to_last_followup"], errors="coerce")
clin = clin.assign(_time=np.where(death, dd, dlf), _event=death.astype(int))
clin["_time"] = pd.to_numeric(clin["_time"], errors="coerce")
surv = clin.dropna(subset=["_time", "_event"])
surv = surv[surv["_time"] > 0]
print(f"[TCGA] 可分析生存样本 {len(surv)}（事件 {int(surv['_event'].sum())}）")

# ---------- 1. 生存分析（正规 log-rank + Cox），全部 22 基因 ----------
sig = pd.read_csv(os.path.join(RES, "scmr_luad_risk_sig.csv"))
genes = sorted(sig["gene_symbol"].unique())
common = sorted(set(expr.columns) & set(surv.index))

rows = []
resolved, unavail = {}, []
for g in genes:
    s = resolve_gene(g, expr.index)
    if s is None:
        unavail.append(g)
    else:
        resolved[g] = s
print(f"[符号] 可分析 {len(resolved)} 基因；矩阵中不可得 {len(unavail)}: {unavail}")
aliased = {k: v for k, v in resolved.items() if k != v}
print(f"[符号] 经别名映射: {aliased}")
JAML_SYM = resolved.get("JAML")
print(f"[符号] JAML 在 TCGA 中的实际符号: {JAML_SYM}")

for g, sym in resolved.items():
    v = pd.to_numeric(expr.loc[sym, common], errors="coerce")
    d = surv.loc[common, ["_time", "_event"]].copy()
    d["expr"] = v.values
    d = d.dropna(subset=["expr"])
    if len(d) < 30 or d["_event"].sum() < 10:
        continue
    grp = (d["expr"] > d["expr"].median()).astype(int).values
    chi2, p_lr = logrank(d["_time"].values, d["_event"].values, grp)
    hr, se, p_cox = cox_hr(d["_time"].values, d["_event"].values, d["expr"].values)
    rows.append({"gene": g, "tcga_symbol": sym, "n": len(d),
                 "events": int(d["_event"].sum()),
                 "HR_per_SD_cox": hr, "se": se, "p_cox": p_cox,
                 "logrank_chi2": chi2, "p_logrank": p_lr})
res = pd.DataFrame(rows)
res["fdr_cox"] = stats.false_discovery_control(res["p_cox"], method="bh")
res["fdr_logrank"] = stats.false_discovery_control(res["p_logrank"], method="bh")
res = res.sort_values("p_cox")
res.to_csv(os.path.join(RES, "m2_tcga_survival_proper.csv"), index=False)
print(f"\n[M2-1] 生存分析（正规 log-rank + Cox）完成，{len(res)} 基因")
print(f"  最小 p_cox = {res['p_cox'].min():.4f}；FDR<0.05 的基因数 = {(res['fdr_cox'] < 0.05).sum()}")
print(res.head(8).to_string(index=False))

# ---------- 2. 分期 ----------
stage_map = {}
for s in clin["pathologic_stage"].dropna().astype(str).unique():
    u = s.upper()
    for r in ["IV", "III", "II", "I"]:
        if f"STAGE {r}" in u:
            stage_map[s] = r
            break
clin["_stage"] = clin["pathologic_stage"].map(stage_map)
print(f"\n[M2-2] 分期取值映射: {stage_map}")

jm = pd.DataFrame({
    "expr": pd.to_numeric(expr.loc[JAML_SYM, common], errors="coerce").values,
    "stage": clin.loc[common, "_stage"].values,
    "sample_type": clin.loc[common, "sample_type"].values,
}, index=common).dropna()
jm = jm[jm["sample_type"].astype(str).str.contains("Primary", case=False, na=False)]
order = ["I", "II", "III", "IV"]
groups = [jm.loc[jm["stage"] == s, "expr"].values for s in order]
groups = [g for g in groups if len(g) >= 5]
kw = stats.kruskal(*groups)
rho, p_rho = stats.spearmanr(jm["stage"].map({s: i + 1 for i, s in enumerate(order)}), jm["expr"])
stage_tbl = jm.groupby("stage")["expr"].agg(n="size", mean="mean", median="median", sd="std").reindex(order)
stage_tbl.to_csv(os.path.join(RES, "m2_tcga_stage_jaml.csv"))
print(f"  JAML 分期差异 Kruskal-Wallis H={kw.statistic:.3f}, P={kw.pvalue:.4f}")
print(f"  分期趋势 Spearman rho={rho:.3f}, P={p_rho:.4f}")
print(stage_tbl.to_string())

# ---------- 3. 免疫 marker 相关性（Spearman，肿瘤样本）----------
tum = [s for s in common if str(clin.loc[s, "sample_type"]).lower().find("primary") >= 0]
jm_expr = pd.to_numeric(expr.loc[JAML_SYM, tum], errors="coerce")
corr = []
for mg, lab in IMMUNE_MARKERS.items():
    if mg not in expr.index:
        continue
    mv = pd.to_numeric(expr.loc[mg, tum], errors="coerce")
    ok = jm_expr.notna() & mv.notna()
    if ok.sum() < 30:
        continue
    rho_, p_ = stats.spearmanr(jm_expr[ok], mv[ok])
    corr.append({"marker": mg, "label": lab, "rho": rho_, "p": p_, "n": int(ok.sum())})
ct = pd.DataFrame(corr).sort_values("rho", ascending=False)
ct["fdr"] = stats.false_discovery_control(ct["p"], method="bh")
ct.to_csv(os.path.join(RES, "m2_tcga_immune_marker_corr.csv"), index=False)
print(f"\n[M2-3] 免疫 marker 相关性（n_tumor={len(tum)}），Spearman")
print(ct.to_string(index=False))

# ---------- 4. 图 ----------
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300})

# 图 1：KM 曲线
d = surv.loc[common, ["_time", "_event"]].copy()
d["expr"] = pd.to_numeric(expr.loc[JAML_SYM, common], errors="coerce").values
d = d.dropna()
d["grp"] = (d["expr"] > d["expr"].median()).astype(int)
fig, ax = plt.subplots(figsize=(4.6, 3.6))
for g, col, lab in [(0, "#4575B4", "JAML low"), (1, "#D73027", "JAML high")]:
    sub = d[d["grp"] == g].sort_values("_time")
    t, e = sub["_time"].values, sub["_event"].values
    s = 1.0
    xs, ys = [0], [1.0]
    for i in range(len(t)):
        if e[i] == 1:
            s *= (1 - 1 / (len(t) - i))
        xs.append(t[i]); ys.append(s)
    ax.step(xs, ys, where="post", color=col, lw=1.4, label=f"{lab} (n={len(sub)})")
chi2, p_lr = logrank(d["_time"].values, d["_event"].values, d["grp"].values)
ax.set_xlabel("Overall survival (days)"); ax.set_ylabel("Survival probability")
ax.set_ylim(0, 1.02); ax.legend(frameon=False, fontsize=8)
ax.set_title(f"JAML expression and OS in TCGA-LUAD\nlog-rank P = {p_lr:.3f}", fontsize=9)
save_fig(fig, "m2_tcga_jaml_km")

# 图 2：分期箱线图
fig, ax = plt.subplots(figsize=(4.2, 3.6))
data = [jm.loc[jm["stage"] == s, "expr"].values for s in order]
bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.6,
                medianprops=dict(color="#333", lw=1.2))
for b in bp["boxes"]:
    b.set_facecolor("#4575B4"); b.set_alpha(0.75); b.set_edgecolor("#333")
ax.set_xlabel("Pathologic stage"); ax.set_ylabel("JAML expression (log2)")
ax.set_title(f"JAML by stage (Kruskal-Wallis P = {kw.pvalue:.3f})", fontsize=9)
for i, g in enumerate(data, start=1):
    ax.text(i, ax.get_ylim()[1], f"n={len(g)}", ha="center", va="bottom", fontsize=7)
save_fig(fig, "m2_tcga_jaml_stage")

# 图 3：免疫 marker 相关性
ct2 = ct.sort_values("rho")
fig, ax = plt.subplots(figsize=(5.2, 0.26 * len(ct2) + 1.2))
cols = ["#D73027" if r > 0 else "#4575B4" for r in ct2["rho"]]
ax.barh(range(len(ct2)), ct2["rho"], color=cols, height=0.72)
ax.set_yticks(range(len(ct2)))
ax.set_yticklabels([f"{m}  ({l})" for m, l in zip(ct2["marker"], ct2["label"])], fontsize=7)
ax.axvline(0, color="#333", lw=0.8)
ax.set_xlabel("Spearman ρ with JAML expression")
ax.set_title("JAML vs immune marker genes (TCGA-LUAD tumors)", fontsize=9)
save_fig(fig, "m2_tcga_jaml_immune_marker")

print("\n[完成] M2 输出：")
for f in ["m2_tcga_survival_proper.csv", "m2_tcga_stage_jaml.csv",
          "m2_tcga_immune_marker_corr.csv"]:
    print("  results/" + f)
