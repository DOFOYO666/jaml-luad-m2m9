# -*- coding: utf-8 -*-
"""审稿人要求补充的三项分析（作为"审稿意见执行"的计算部分）。

R1 —— ID2 扰动的**相对效应量**：把 |ΔJAML| 从 count 尺度换算为相对 JAML 基础表达的比例，
      并给出逐细胞分布，回答"0.0175 到底是多大"。
R2 —— TCGA 生存的**参数化敏感性**：同一批样本、同一批基因，比较
      ① Cox（连续，每 SD）② Cox（连续，每 log2 单位）③ Cox（中位数分组，分类）④ log-rank（中位数分组），
      用以解释"log-rank 显著而 Cox 不显著"是参数化差异还是假象。
R3 —— M5（GSE135222，n = 27）的**最小可检出效应量**（Schoenfeld 公式）与功效说明。

输出：results/reviewer_R1_ID2_relative_effect.csv
      results/reviewer_R2_TCGA_parameterisation.csv
      results/reviewer_R3_M5_power.json
      logs/reviewer_analyses.log（本脚本 stdout）
"""
import io
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
from scipy import optimize, stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
DATA = os.path.join(BASE, "data", "luad_survival")
D = r"D:\workbuddy工作空间\JAML深度研究"
DRES = os.path.join(D, "results")


def hdr(t):
    print("\n" + "=" * 22 + f" {t} " + "=" * 22)


# ==================================================================== R1
def r1_id2_relative_effect():
    hdr("R1  ID2 扰动的相对效应量")
    h5 = os.path.join(DRES, "CD4T_celloracle_ext2.h5ad")
    ko = os.path.join(DRES, "m4_deltaX_ID2_KO_ext2.npy")
    rnd = os.path.join(DRES, "m4_deltaX_randomGRN_ext2.npy")
    for p in (h5, ko, rnd):
        if not os.path.isfile(p):
            print(f"  [缺失] {p}")
            return None
    a = ad.read_h5ad(h5)
    names = list(map(str, a.var_names))
    ji = names.index("JAML")
    RC = a.layers["raw_count"]
    col = RC[:, ji]
    v = np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
    v = v.astype(float)
    d_ko = np.load(ko)
    d_rn = np.load(rnd)
    print(f"  delta_X 形状 {d_ko.shape}（细胞 × 基因）；JAML 列号 {ji}")
    if d_ko.shape[1] != len(names):
        print("  [警告] delta_X 基因数与 h5ad 不一致，按最末列处理前先核对")
    mean_cnt = float(v.mean())
    pos_frac = float((v > 0).mean())
    mean_cnt_pos = float(v[v > 0].mean()) if (v > 0).any() else np.nan
    dk = d_ko[:, ji].astype(float)
    dr = d_rn[:, ji].astype(float)
    abs_ko = float(np.abs(dk).mean())
    abs_rn = float(np.abs(dr).mean())
    # 相对效应：以阳性细胞的基础表达为分母（零膨胀下以全体均值为分母会低估）
    rel_all = abs_ko / mean_cnt if mean_cnt else np.nan
    rel_pos = abs_ko / mean_cnt_pos if mean_cnt_pos else np.nan
    rel_rn_all = abs_rn / mean_cnt if mean_cnt else np.nan
    print(f"  JAML 基础表达：均值 {mean_cnt:.4f} count；阳性率 {100*pos_frac:.1f}%；"
          f"阳性细胞均值 {mean_cnt_pos:.4f}")
    print(f"  真扰动 mean|ΔJAML| = {abs_ko:.5f}  →  相对全体均值 {100*rel_all:.2f}%  /  "
          f"相对阳性细胞均值 {100*rel_pos:.2f}%")
    print(f"  本次运行的打乱对照 mean|ΔJAML| = {abs_rn:.6f}"
          f"（→ 相对全体均值 {100*rel_rn_all:.4f}%）")
    # 该 npy 对应的单次打乱读数恰为 0（与"JAML 的打乱边被滤除"一致），比值不可计算；
    # 权威噪声来自配对对照的 5 次独立打乱（见 m4_noise_control_classify.csv）。
    pair_mean = pair_max = np.nan
    pclass = os.path.join(DRES, "m4_noise_control_classify.csv")
    if os.path.isfile(pclass):
        pc = pd.read_csv(pclass)
        row = pc[pc["TF"] == "ID2"]
        if len(row):
            pair_mean = float(row["jaml_rand_mean"].iloc[0])
            pair_max = float(row["jaml_rand_max"].iloc[0])
            print(f"  配对对照（5 次独立打乱，权威基线）：均值 {pair_mean:.6f}，最大 {pair_max:.6f}")
            print(f"  相对效应之比：真扰动 / 配对均值 = {abs_ko / pair_mean:.1f} 倍；"
                  f"真扰动 / 配对最大 = {abs_ko / pair_max:.1f} 倍")
    ratio = (abs_ko / pair_max) if (pair_max and pair_max > 0) else np.nan
    # 逐细胞分布
    q = np.percentile(np.abs(dk), [50, 75, 90, 95, 99, 100])
    print("  逐细胞 |ΔJAML| 分位（50/75/90/95/99/max）: " +
          " / ".join(f"{x:.4f}" for x in q))
    n_nonzero = int((np.abs(dk) > 1e-8).sum())
    print(f"  非零 |ΔJAML| 的细胞数：{n_nonzero} / {len(dk)}")
    # 方向
    signed = float(dk.mean())
    print(f"  签名均值 ΔJAML = {signed:+.5f}（负 = 敲除 ID2 使 JAML 下降）")
    # 受影响基因排名（相对各自基础表达）
    out = pd.DataFrame({
        "item": ["JAML_mean_count", "JAML_pos_frac", "JAML_mean_count_positive",
                 "mean_abs_delta_ID2_KO", "mean_abs_delta_randomized",
                 "relative_to_all_cells", "relative_to_positive_cells",
                 "relative_randomized", "ratio_vs_paired_max",
                 "signed_mean_delta", "cells_abs_gt_1e-8",
                 "pctl_50", "pctl_75", "pctl_90", "pctl_95", "pctl_99", "max"],
        "value": [mean_cnt, pos_frac, mean_cnt_pos, abs_ko, abs_rn, rel_all, rel_pos,
                  rel_rn_all, ratio, signed, n_nonzero, *q],
    })
    out.to_csv(os.path.join(RES, "reviewer_R1_ID2_relative_effect.csv"), index=False)

    # 全局：ID2 敲除影响最大的基因（按 |Δ| 相对该基因自身均值）
    RCd = RC.todense() if hasattr(RC, "todense") else np.asarray(RC)
    RCd = np.asarray(RCd, dtype=float)
    gmean = RCd.mean(axis=0)
    mag = np.abs(d_ko.astype(float)).mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(gmean > 0, mag / gmean, np.nan)
    order = np.argsort(-np.nan_to_num(rel))
    print("\n  ID2 敲除后相对效应最大的 12 个基因（|Δ| / 该基因均值）：")
    for i in order[:12]:
        print(f"    {names[i]:10s} 相对 {100*rel[i]:7.2f}%   mean={gmean[i]:8.3f}  |Δ|={mag[i]:.4f}")
    r = names.index("JAML")
    print(f"    … JAML 排名 {int((rel > rel[r]).sum())+1} / {len(names)}"
          f"（相对 {100*rel[r]:.2f}%）")
    return out


# ==================================================================== R2
def _cox(time, event, x, standardize=True):
    """一维 Cox（Breslow），返回 (coef, se, p, hr)。x 可为连续或 0/1。"""
    time = np.asarray(time, float); event = np.asarray(event, int)
    x = np.asarray(x, float)
    sd = x.std()
    xs = (x - x.mean()) / sd if standardize else x.copy()

    def negll(b):
        b = float(np.atleast_1d(b)[0])
        ll = 0.0
        for t in np.unique(time[event == 1]):
            risk = time >= t
            d = (time == t) & (event == 1)
            eta = b * xs[risk]
            m = eta.max()
            ll += b * xs[d].sum() - d.sum() * (m + np.log(np.exp(eta - m).sum()))
        return -ll

    res = optimize.minimize(negll, x0=[0.0], method="BFGS")
    b = float(res.x[0])
    h = 1e-5
    d2 = (negll([b + h]) - 2 * negll([b]) + negll([b - h])) / h ** 2
    se = float(np.sqrt(1 / d2)) if d2 > 0 else np.nan
    z = b / se if se > 0 else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if not np.isnan(z) else np.nan
    return b, se, p, float(np.exp(b))


def _logrank(time, event, group):
    time = np.asarray(time, float); event = np.asarray(event, int)
    group = np.asarray(group, int)
    O1 = E1 = V = 0.0
    for t in np.unique(time[event == 1]):
        at = time >= t
        n = at.sum(); n1 = (at & (group == 1)).sum()
        d = ((time == t) & (event == 1)).sum()
        d1 = ((time == t) & (event == 1) & (group == 1)).sum()
        if n <= 1:
            continue
        O1 += d1; E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    chi2 = (O1 - E1) ** 2 / V if V > 0 else np.nan
    return float(chi2), float(stats.chi2.sf(chi2, 1))


def r2_tcga_parameterisation():
    hdr("R2  TCGA 生存的参数化敏感性")
    expr = pd.read_csv(os.path.join(DATA, "TCGA.LUAD.HiSeqV2.gz"), sep="\t",
                       compression="gzip", index_col=0)
    clin = pd.read_csv(os.path.join(DATA, "TCGA.LUAD.clinicalMatrix.txt"), sep="\t", index_col=0)
    vs = clin["vital_status"].astype(str).str.upper()
    death = vs.str.contains("DECEASED", na=False)
    dd = pd.to_numeric(clin["days_to_death"], errors="coerce")
    dlf = pd.to_numeric(clin["days_to_last_followup"], errors="coerce")
    clin = clin.assign(_time=np.where(death, dd, dlf), _event=death.astype(int))
    surv = clin.dropna(subset=["_time", "_event"])
    surv = surv[surv["_time"] > 0]
    print(f"  表达 {expr.shape}；可分析生存样本 {len(surv)}（事件 {int(surv['_event'].sum())}）")

    sig = pd.read_csv(os.path.join(RES, "scmr_luad_risk_sig.csv"))
    genes = sorted(sig["gene_symbol"].unique())
    alias = {"JAML": ["JAML", "AMICA1"], "HYKK": ["HYKK", "AGPHD1"]}
    common = sorted(set(expr.columns) & set(surv.index))
    rows = []
    for g in genes:
        used = next((c for c in alias.get(g, [g]) if c in expr.index), None)
        if used is None or used not in expr.index:
            continue
        v = pd.to_numeric(expr.loc[used, common], errors="coerce")
        t = pd.to_numeric(surv.loc[common, "_time"], errors="coerce")
        e = pd.to_numeric(surv.loc[common, "_event"], errors="coerce")
        ok = v.notna() & t.notna() & e.notna()
        v, t, e = v[ok].values, t[ok].values, e[ok].values.astype(int)
        if len(v) < 50:
            continue
        b1, se1, p_c_sd, hr_c_sd = _cox(t, e, v, standardize=True)
        b2, se2, p_c_u, hr_c_u = _cox(t, e, v, standardize=False)   # 每 1 log2 单位
        grp = (v > np.median(v)).astype(int)
        b3, se3, p_c_cat, hr_c_cat = _cox(t, e, grp, standardize=False)
        chi2, p_lr = _logrank(t, e, grp)
        rows.append(dict(gene=g, symbol_used=used, n=len(v), events=int(e.sum()),
                         HR_cox_per_SD=hr_c_sd, p_cox_per_SD=p_c_sd,
                         HR_cox_per_log2unit=hr_c_u, p_cox_per_log2unit=p_c_u,
                         HR_cox_median_split=hr_c_cat, p_cox_median_split=p_c_cat,
                         chi2_logrank=chi2, p_logrank=p_lr))
    d = pd.DataFrame(rows)
    for col, name in [("p_cox_per_SD", "fdr_cox_per_SD"),
                      ("p_cox_per_log2unit", "fdr_cox_per_log2unit"),
                      ("p_cox_median_split", "fdr_cox_median_split"),
                      ("p_logrank", "fdr_logrank")]:
        d[name] = stats.false_discovery_control(d[col].values, method="bh")
    d = d.sort_values("p_logrank")
    d.to_csv(os.path.join(RES, "reviewer_R2_TCGA_parameterisation.csv"), index=False)
    show = d[["gene", "n", "HR_cox_per_SD", "p_cox_per_SD", "fdr_cox_per_SD",
              "HR_cox_median_split", "p_cox_median_split", "fdr_cox_median_split",
              "p_logrank", "fdr_logrank"]]
    print("\n  按 log-rank P 排序（前 8 行）：")
    print(show.head(8).round(4).to_string(index=False))
    j = d[d["gene"] == "JAML"]
    if len(j):
        r = j.iloc[0]
        print(f"\n  JAML 四种参数化：")
        print(f"    Cox 每 SD       HR = {r['HR_cox_per_SD']:.3f}  P = {r['p_cox_per_SD']:.4f}  "
              f"FDR = {r['fdr_cox_per_SD']:.4f}")
        print(f"    Cox 每 log2 单位 HR = {r['HR_cox_per_log2unit']:.3f}  "
              f"P = {r['p_cox_per_log2unit']:.4f}  FDR = {r['fdr_cox_per_log2unit']:.4f}")
        print(f"    Cox 中位数分组   HR = {r['HR_cox_median_split']:.3f}  "
              f"P = {r['p_cox_median_split']:.4f}  FDR = {r['fdr_cox_median_split']:.4f}")
        print(f"    log-rank         chi2 = {r['chi2_logrank']:.3f}  "
              f"P = {r['p_logrank']:.5f}  FDR = {r['fdr_logrank']:.4f}")
    n_sig = {c: int((d[c] < 0.05).sum()) for c in
             ["fdr_cox_per_SD", "fdr_cox_median_split", "fdr_logrank"]}
    print(f"\n  FDR < 0.05 的基因数：Cox 每 SD = {n_sig['fdr_cox_per_SD']}；"
          f"Cox 中位数分组 = {n_sig['fdr_cox_median_split']}；log-rank = {n_sig['fdr_logrank']}"
          f"（共 {len(d)} 个基因）")
    for c in ["fdr_cox_per_SD", "fdr_cox_median_split", "fdr_logrank"]:
        g = d.loc[d[c] < 0.05, "gene"].tolist()
        print(f"    {c} 通过的基因：{g if g else '无'}")
    both = d[(d["fdr_cox_per_SD"] < 0.05) & (d["fdr_cox_median_split"] < 0.05)]
    print(f"  两种 Cox 参数化下均通过 FDR 的基因：{both['gene'].tolist() if len(both) else '无'}")
    return d


# ==================================================================== R3
def r3_m5_power():
    hdr("R3  M5（GSE135222）最小可检出效应量")
    d = pd.read_csv(os.path.join(RES, "m5_gse135222_jaml_pfs.csv"))
    n = len(d); n_ev = int(d["pfs"].sum())
    grp = (d["expr"] > d["expr"].median()).astype(int).values
    p_hi = float(grp.mean())
    z_a = stats.norm.ppf(0.975); z_b = stats.norm.ppf(0.80)
    # Schoenfeld：所需事件数 d = (z_a + z_b)^2 / (p(1-p) * ln(HR)^2)  →  反解 HR
    ln_hr = (z_a + z_b) / np.sqrt(n_ev * p_hi * (1 - p_hi))
    hr_mde = float(np.exp(ln_hr))
    # 实际观测 HR（中位数分组，Cox）
    b, se, p, hr = _cox(d["pfs_time"].values, d["pfs"].values, grp, standardize=False)
    chi2, p_lr = _logrank(d["pfs_time"].values, d["pfs"].values, grp)
    print(f"  n = {n}，事件 {n_ev}；高表达组 {int(grp.sum())} 例、低表达组 {int((1-grp).sum())} 例")
    print(f"  观测：HR = {hr:.3f}（95% CI {np.exp(b-1.96*se):.3f}–{np.exp(b+1.96*se):.3f}），"
          f"Cox P = {p:.3f}；log-rank P = {p_lr:.3f}")
    print(f"  **在 80% 功效、α = 0.05（双侧）下，本样本可检出的最小 HR ≈ {hr_mde:.2f}**")
    print(f"  → 该队列只能排除「很强」的关联；对 HR 在 1.0–{hr_mde:.1f} 之间的效应功效不足，"
          f"阴性结果不能解释为「无关联」。")
    out = dict(n=n, events=n_ev, group_high=int(grp.sum()), group_low=int((1 - grp).sum()),
               observed_HR=hr, observed_HR_CI_low=float(np.exp(b - 1.96 * se)),
               observed_HR_CI_high=float(np.exp(b + 1.96 * se)),
               cox_p=float(p), logrank_p=float(p_lr),
               min_detectable_HR_80power=hr_mde,
               method="Schoenfeld 公式反解 HR；事件数 %d、分组比例 %.3f" % (n_ev, p_hi))
    with io.open(os.path.join(RES, "reviewer_R3_M5_power.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    return out


if __name__ == "__main__":
    r1_id2_relative_effect()
    r2_tcga_parameterisation()
    r3_m5_power()
    print("\nDONE")
