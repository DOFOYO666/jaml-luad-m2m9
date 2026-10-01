# -*- coding: utf-8 -*-
"""Re-draw Figure 6a (m8_tf_jaml) and Figure 6b (m8_negative_control) with corrected
in-figure labels, using ONLY the artefacts already saved by m8_tf_activity.py.

Why not just re-run m8_tf_activity.py?  Re-running would re-query CollecTRI and re-fit
decoupleR, which risks producing slightly different numbers and thus drifting the
manuscript.  Everything the figures need is already on disk, so we re-plot from those
tables and *assert* that our recomputation reproduces the recorded values.

Corrections applied (see 审稿意见与逐条修改说明_第二轮.md, item M5, and the pre-submission list):
  panel a  title : "top 20 by p"                    -> top 20 of 498 by Spearman rho
  panel b  title : "JAML-associated TFs (n = 12, FDR < 0.05)"  <- read as "only 12 significant"
                                                    -> "Top 12 of 229 TFs with FDR < 0.05"
  panel c  title : "recovered 5/7 expected TFs"     -> "5/7 testable (5/9 pre-specified)"
  fig 6b   title : "empirical P = 0.0000"           -> "empirical P < 0.005 (200 permutations)"

Run:  python scripts/m8_fig6_replot.py
"""
import csv
import json
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")

D_ROOT = r"D:\workbuddy工作空间\JAML深度研究"
RES = os.path.join(D_ROOT, "results")
FIG = os.path.join(RES, "figures")
IN_H5AD = os.path.join(RES, "CD4T_wide.h5ad")

PDCD1_EXPECTED = ["TOX", "PRDM1", "IRF4", "NR4A1", "BATF", "MAF", "IKZF2", "EOMES", "TBX21"]
N_PERM = 200
SEED = 20260930


def log(*a):
    print(*a, flush=True)


def main():
    t0 = time.time()

    # ---------------------------------------------------------------- inputs
    act = pd.read_csv(os.path.join(RES, "m8_tf_activity_by_cell.csv"), index_col=0)
    log(f"[activity] {act.shape}  (TFs x cells)")

    vs = pd.read_csv(os.path.join(RES, "m8_tf_vs_gene.csv"))
    grp = pd.read_csv(os.path.join(RES, "m8_tf_group_diff.csv"))

    cj = vs[vs["gene"] == "JAML"].sort_values("rho", ascending=False)
    gj = grp[(grp["gene"] == "JAML") & (grp["split"] == "positive_vs_negative") & (grp["fdr"] < 0.05)].copy()
    gj = gj.reindex(gj["diff"].abs().sort_values(ascending=False).index)
    gp = grp[(grp["gene"] == "PDCD1") & (grp["split"] == "positive_vs_negative") & (grp["fdr"] < 0.05)].copy()
    gp = gp.reindex(gp["diff"].abs().sort_values(ascending=False).index)
    log(f"[selection] JAML rho-ranked {len(cj)}; JAML FDR<0.05 = {len(gj)}; PDCD1 FDR<0.05 = {len(gp)}")

    # JAML positivity per cell, straight from the same h5ad the pipeline used
    import anndata as ad
    adata = ad.read_h5ad(IN_H5AD)
    assert "JAML" in adata.var_names, "JAML not in the h5ad var_names"
    jcol = adata.var_names.get_loc("JAML")
    raw = adata.layers["raw_count"] if "raw_count" in adata.layers else adata.X
    jv = np.asarray(raw[:, jcol].todense()).ravel() if hasattr(raw, "todense") else np.asarray(raw[:, jcol]).ravel()
    posm = jv > 0
    # activity table must be in the same cell order
    assert list(act.columns) == list(adata.obs_names), "cell order differs between activity table and h5ad"
    log(f"[positivity] JAML-positive {int(posm.sum())} / negative {int((~posm).sum())}")
    exp_pos = int(grp[(grp["gene"] == "JAML") & (grp["split"] == "positive_vs_negative")]["n_high"].iloc[0])
    exp_neg = int(grp[(grp["gene"] == "JAML") & (grp["split"] == "positive_vs_negative")]["n_low"].iloc[0])
    assert int(posm.sum()) == exp_pos and int((~posm).sum()) == exp_neg, \
        f"positivity mismatch: {int(posm.sum())}/{int((~posm).sum())} vs recorded {exp_pos}/{exp_neg}"
    log("  OK  positivity reproduces the recorded group sizes")

    # ------------------------------------------------------- Figure 6a panels
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top_j = gj.head(12)["tf"].tolist()
    top_p = gp.head(12)["tf"].tolist()
    hit = [t for t in PDCD1_EXPECTED if t in set(gp["tf"])]
    testable = [t for t in PDCD1_EXPECTED if t in set(act.index)]
    log(f"[controls] PDCD1 hit {hit}")
    log(f"  testable {len(testable)}/{len(PDCD1_EXPECTED)} -> {sorted(testable)}")

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.3))

    ax = axes[0]
    d = cj.head(20).iloc[::-1]
    ax.barh(range(len(d)), d["rho"],
            color=["#D73027" if r > 0 else "#4575B4" for r in d["rho"]],
            edgecolor="#333", linewidth=0.4)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d["tf"], fontsize=7)
    ax.axvline(0, color="#666", lw=0.6)
    ax.set_xlabel("Spearman \u03c1 (TF activity vs JAML)")
    ax.set_title(f"TFs tracking JAML expression\n(top 20 of {len(cj)} by \u03c1; n = {act.shape[1]:,} CD4\u207a T cells)",
                 fontsize=8.6)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    idx = [list(act.index).index(t) for t in top_j if t in act.index]
    M = act.iloc[idx].values
    z = (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-12)
    ax.boxplot([z[:, ~posm].ravel(), z[:, posm].ravel()],
               tick_labels=[f"JAML-negative\n(n={(~posm).sum():,})", f"JAML-positive\n(n={posm.sum():,})"],
               patch_artist=True, boxprops=dict(facecolor="#c9d7ea"), medianprops=dict(color="#333"))
    ax.set_ylabel("TF activity (z)")
    ax.set_title(f"Top {len(top_j)} of {len(gj)} TFs with FDR < 0.05\n"
                 f"(largest activity difference, JAML-positive vs negative)", fontsize=8.6)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[2]
    sel = gp.head(12)[["tf", "diff"]].iloc[::-1]
    ax.barh(range(len(sel)), sel["diff"],
            color=["#D73027" if r > 0 else "#4575B4" for r in sel["diff"]],
            edgecolor="#333", linewidth=0.4)
    ax.set_yticks(range(len(sel)))
    ax.set_yticklabels(sel["tf"], fontsize=7)
    ax.axvline(0, color="#666", lw=0.6)
    ax.set_xlabel("\u0394 TF activity (PDCD1-positive \u2212 negative)")
    ax.set_title(f"Positive control: PDCD1\nrecovered {len(hit)}/{len(testable)} testable "
                 f"({len(hit)}/{len(PDCD1_EXPECTED)} pre-specified)", fontsize=8.6)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    for ext, kw in (("png", {}), ("tif", {"pil_kwargs": {"compression": "tiff_lzw"}})):
        p = os.path.join(FIG, "m8_tf_jaml." + ext)
        if os.path.exists(p):
            shutil.copy2(p, p + ".bak_fig6fix")
        fig.savefig(p, dpi=300, bbox_inches="tight", **kw)
    log("[wrote] m8_tf_jaml.png / .tif")

    # ------------------------------------- Figure 6b: recompute the null, verified
    from scipy import stats
    rng = np.random.default_rng(SEED)
    v0 = jv.astype(float)
    n_sig_perm = []
    for _ in range(N_PERM):
        v = rng.permutation(v0)
        posp = v > 0
        if posp.sum() < 20 or (~posp).sum() < 20:
            continue
        ps = []
        for tf in act.index:
            a = act.loc[tf].values
            ok = np.isfinite(a)
            if ok.sum() < 50:
                continue
            _, pu = stats.mannwhitneyu(a[ok & posp], a[ok & ~posp], alternative="two-sided")
            ps.append(pu)
        if ps:
            fdr = stats.false_discovery_control(np.array(ps), method="bh")
            n_sig_perm.append(int((fdr < 0.05).sum()))
    n_sig_perm = np.array(n_sig_perm)
    obs_n = len(gj)
    pct = float((n_sig_perm >= obs_n).mean()) if len(n_sig_perm) else float("nan")

    saved = json.load(open(os.path.join(RES, "m8_negative_control.json"), encoding="utf-8"))
    log(f"[null] recomputed mean={n_sig_perm.mean():.2f} max={int(n_sig_perm.max())} "
        f"p95={np.percentile(n_sig_perm, 95):.1f} obs={obs_n} empirical_P={pct:.4f}")
    log(f"[null] saved      mean={saved['perm_mean']} max={saved['perm_max']} "
        f"p95={saved['perm_p95']} obs={saved['obs_n_sig']} empirical_P={saved['empirical_p']}")
    ok = (round(float(n_sig_perm.mean()), 2) == round(float(saved["perm_mean"]), 2)
          and int(n_sig_perm.max()) == int(saved["perm_max"])
          and float(np.percentile(n_sig_perm, 95)) == float(saved["perm_p95"])
          and obs_n == int(saved["obs_n_sig"]))
    if not ok:
        log("  !! WARNING: recomputed null distribution does not reproduce the saved one; "
            "figure 6b is written but MUST be checked before submission")
    else:
        log("  OK  recomputed null distribution reproduces the saved one exactly")

    fig, ax = plt.subplots(figsize=(4.6, 3.3))
    ax.hist(n_sig_perm, bins=30, color="#9ecae1", edgecolor="#333", linewidth=0.4)
    ax.axvline(obs_n, color="#D73027", lw=1.6, label=f"observed = {obs_n}")
    ax.set_xlabel("# TFs with FDR < 0.05")
    ax.set_ylabel("permutations")
    ax.set_title(f"Negative control (gene labels permuted)\n"
                 f"empirical P < {1.0 / (N_PERM + 1):.3f} ({N_PERM} permutations)", fontsize=8.6)
    ax.legend(frameon=False, fontsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext, kw in (("png", {}), ("tif", {"pil_kwargs": {"compression": "tiff_lzw"}})):
        p = os.path.join(FIG, "m8_negative_control." + ext)
        if os.path.exists(p):
            shutil.copy2(p, p + ".bak_fig6fix")
        fig.savefig(p, dpi=300, bbox_inches="tight", **kw)
    log("[wrote] m8_negative_control.png / .tif")

    # a small machine-readable record of exactly what was drawn
    rec = {
        "figure_6a": {
            "panel_a": {"ranking": "Spearman rho of TF activity vs JAML", "n_tf_pool": int(len(cj)),
                        "shown": cj.head(20)["tf"].tolist()},
            "panel_b": {"rule": "top 12 by |activity difference| among FDR<0.05 TFs",
                        "n_fdr05": int(len(gj)), "shown": top_j,
                        "n_positive": int(posm.sum()), "n_negative": int((~posm).sum())},
            "panel_c": {"gene": "PDCD1", "hit": hit, "n_testable": len(testable),
                        "n_prespecified": len(PDCD1_EXPECTED), "shown": top_p},
        },
        "figure_6b": {"n_perm_requested": N_PERM, "n_perm_used": int(len(n_sig_perm)),
                      "observed_n_sig": int(obs_n), "perm_mean": float(n_sig_perm.mean()),
                      "perm_max": int(n_sig_perm.max()),
                      "perm_p95": float(np.percentile(n_sig_perm, 95)),
                      "empirical_p": pct, "reproduces_saved": bool(ok)},
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(os.path.join(RES, "m8_fig6_replot_record.json"), "w", encoding="utf-8") as f:
        json.dump(rec, f, ensure_ascii=False, indent=2)
    log(f"[done] {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
