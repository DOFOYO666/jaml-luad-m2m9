# -*- coding: utf-8 -*-
"""【M8 替代路线】CD4⁺T 的转录因子活性分析（decoupleR + CollecTRI）

⚠️ 与规划里 M8（pySCENIC/SCENIC）的关系 —— **必须如实声明**
  pySCENIC 本体在本机 **不可行**：它需要 cisTarget motif 数据库
  （hg38 的 `hg38__refseq-r80__*_tss.mc9nr.feather` 合计约 1.6 GB），
  本机带宽实测 ~40 kB/s，下载需 10 小时以上。
  因此改用 **decoupleR + CollecTRI**（Badia-i-Mompel 2022, Bioinformatics Advances）：
    - CollecTRI 是**人工审编的 TF–靶基因**资源（约 1,100 个 TF），**不需要 motif 数据库**；
    - 用 `ulm`（单变量线性模型）从表达矩阵反推 **TF 活性**。
  **这不是 SCENIC，也不给出 SCENIC 的 regulon 特异度分数**；
  它回答的是"哪些 TF 的活性与 JAML 表达/细胞状态同步"，属于**关联性**证据。

⚠️ 同时声明：CellOracle 本体也跑不了（见 `docs/执行报告_M4M8_D盘尝试.md`）：
  celloracle 0.20.0 硬依赖 `velocyto`，而 velocyto 在 PyPI 上**全部版本都只有 sdist、
  0 个 wheel**，Windows 无编译器装不上；且该依赖位于转移概率核心（`velocyto.diffusion`），
  不可用替身绕过。故 M4 的 CellOracle 版本只能在大内存 + 有编译器的机器上做。

分析内容（含规划要求的正/负对照）：
  ① 每个 TF 的活性 与 **JAML 表达** 的 Spearman；
  ② **JAML 高/低分组**（中位数切分）的 TF 活性差异（Mann-Whitney + BH FDR）；
  ③ TF 活性沿 CD4 亚型（Naive→Treg→Th→Mixed Th→Exhausted Tfh）的变化趋势（Spearman）；
  ④ **阳性对照**：同样流程对 **PDCD1**（耗竭标志）做一遍——应富集耗竭相关 TF
     （TOX/PRDM1/IRF4/NR4A1/BATF 等）；若对照通不过，说明流程无效，结果不可信；
  ⑤ **阴性对照**：对**随机基因表达**做主分析同款检验（200 次置换）→ 应无可重复的 TF。

输入：D:/workbuddy工作空间/JAML深度研究/results/CD4T_wide.h5ad
输出：results/m8_tf_activity_by_cell.csv
      results/m8_tf_vs_gene.csv               （TF 活性 vs JAML / PDCD1）
      results/m8_tf_by_subtype.csv
      results/m8_tf_group_diff.csv            （高/低分组差异 + FDR）
      results/m8_negative_control.json
      results/figures/m8_tf_jaml.*
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

try:
    import anndata as ad
    import decoupler as dc
    import scanpy as sc
    from scipy import stats
except ImportError as e:
    sys.exit(f"缺少依赖: {e}")

D_ROOT = r"D:\workbuddy工作空间\JAML深度研究"
RES = os.path.join(D_ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)
IN_H5AD = os.path.join(RES, "CD4T_wide.h5ad")

SUBTYPE_ORDER = ["Naive CD4+ T", "Treg", "CD4+ Th", "CD8+/CD4+ Mixed Th", "Exhausted Tfh"]
# 阳性对照的"期望 TF"（用于验证流程能否恢复已知生物学）
PDCD1_EXPECTED = ["TOX", "PRDM1", "IRF4", "NR4A1", "BATF", "MAF", "IKZF2", "EOMES", "TBX21"]
N_PERM = 200
SEED = 20260930


def log(*a):
    print(*a, flush=True)


def main():
    t0 = time.time()
    log("=" * 70)
    log("M8 替代路线：CD4⁺T 的 TF 活性分析（decoupleR + CollecTRI）")
    log("=" * 70)

    if not os.path.exists(IN_H5AD):
        sys.exit(f"找不到输入 {IN_H5AD}（需先跑 scripts/m4_prep_bounded.py）")
    adata = sc.read_h5ad(IN_H5AD)
    log(f"[载入] {adata.shape}；亚型 {adata.obs['cell_subtype'].value_counts().to_dict()}")

    # 若已有归一化 X 就跳过；否则由 raw_count 归一化
    raw = adata.layers["raw_count"] if "raw_count" in adata.layers else adata.X
    adata.X = raw.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    log(f"[归一化] CPM(1e4) + log1p；基因数 {adata.n_vars}")

    # ---------------- 1. CollecTRI 网络 ----------------
    log("\n[网络] 获取 CollecTRI ...")
    net = None
    for name, fn in [("collectri", lambda: dc.op.collectri(organism="human")),
                     ("dorothea", lambda: dc.op.dorothea(organism="human"))]:
        try:
            net = fn()
            log(f"  {name}: {net.shape}（sources {net['source'].nunique()}）")
            net.to_csv(os.path.join(RES, f"m8_net_{name}.csv"), index=False)
            break
        except Exception as e:
            log(f"  {name} 获取失败: {type(e).__name__}: {str(e)[:120]}")
    if net is None:
        sys.exit("无法获取 TF 网络资源（网络受限）—— M8 无法进行")

    # 只保留在数据中存在的靶基因
    genes = set(adata.var_names)
    net = net[net["target"].isin(genes)].copy()
    log(f"[网络] 与数据交集后: {net.shape}；覆盖 TF {net['source'].nunique()} 个")
    if len(net) < 100:
        sys.exit("交集过少，无法做 TF 活性分析")

    # ---------------- 2. TF 活性 ----------------
    log("\n[活性] decoupleR ulm ...")
    acts = dc.mt.ulm(data=adata, net=net, tmin=5, verbose=False)
    log(f"  ulm 返回类型: {type(acts)}"
        + (f"，列={list(acts.columns)}，形状={acts.shape}"
           if isinstance(acts, pd.DataFrame) else ""))
    A = None
    if isinstance(acts, pd.DataFrame):
        if {"source", "score"} <= set(acts.columns):
            if "condition" in acts.columns:
                A = acts.pivot(index="source", columns="condition", values="score")
        else:
            A = acts                       # 已是 宽表 sources × samples
    if A is None:
        # 退回 obsm（decoupler 对 AnnData 输入会把结果写进 obsm，返回 None）
        for k in list(getattr(adata, "obsm", {}).keys()):
            if "ulm" in k.lower() or "score" in k.lower():
                obj = adata.obsm[k]
                if isinstance(obj, pd.DataFrame):
                    A = obj.T.copy()          # 保留 TF 名（index）
                    log(f"  obsm['{k}'] 是 DataFrame，已取 .T 保留 TF 名")
                else:
                    A = pd.DataFrame(np.asarray(obj).T)
                    log(f"  obsm['{k}'] 不是 DataFrame（无 TF 名），需注意")
                break
    if A is None:
        raise RuntimeError("无法解析 ulm 输出形态")
    # 列名对齐到细胞名
    if len(A.columns) == adata.n_obs:
        A = A.copy()
        A.columns = list(adata.obs_names)
    elif len(set(A.columns) & set(adata.obs_names)) > 0:
        A = A.reindex(columns=list(adata.obs_names))
    else:
        raise RuntimeError("ulm 输出的细胞名与 adata.obs_names 对不上")
    A = A.astype(np.float64)
    A = A.loc[A.notna().any(axis=1), A.notna().any(axis=0)]
    # ⚠️ 守卫：TF 名不能是数字索引，否则后续所有输出都不可解读
    if all(isinstance(x, (int, np.integer)) for x in A.index[:5]):
        raise RuntimeError(
            f"TF 名疑似位置索引而非符号（前 5 个: {list(A.index[:5])}）——"
            "说明 ulm 结果取错了对象，请检查 decoupler 版本与 obsm 键")
    log(f"  [活性] 矩阵 {A.shape}；列名与细胞名对齐: {set(A.columns) <= set(adata.obs_names)}")
    log(f"  [活性] TF 名示例: {list(A.index[:8])}")
    log(f"[活性] 矩阵 {A.shape}（TF × 细胞）")
    A.to_csv(os.path.join(RES, "m8_tf_activity_by_cell.csv"))

    X = pd.DataFrame(np.asarray(adata.X.todense() if hasattr(adata.X, "todense") else adata.X),
                     index=adata.obs_names, columns=adata.var_names)
    sub = pd.Series(adata.obs["cell_subtype"].astype(str).values, index=adata.obs_names)

    # ---------------- 3. 通用分析函数 ----------------
    def analyse(gene):
        """对某个基因做：TF 活性 vs 表达、以及分组差异。

        ⚠️ 关键：JAML（以及多数免疫基因）在 CD4⁺T 中**零膨胀**（>70% 细胞为 0）。
           若用"中位数切分"，因为中位数就是 0，会退化成"全部为高"，
           分组差异检验直接失效（这正是第一版跑出空表的原因）。
           因此本函数改用**表达阳性 vs 阴性**作为主分组，并在阳性细胞内部
           再做一次上/下半分组的次级检验。
        """
        if gene not in X.columns:
            return None, None
        v = X[gene].values.astype(float)
        pos = v > 0
        frac_pos = float(pos.mean())
        rows_c, rows_g = [], []
        for tf in A.index:
            a = A.loc[tf].values
            ok = np.isfinite(a) & np.isfinite(v)
            if ok.sum() < 50:
                continue
            rho, p = stats.spearmanr(a[ok], v[ok])
            rows_c.append({"tf": tf, "gene": gene, "rho": rho, "p": p, "n": int(ok.sum())})

            # 主分组：阳性 vs 阴性
            ga, gb = a[ok & pos], a[ok & ~pos]
            if len(ga) >= 20 and len(gb) >= 20:
                _, pu = stats.mannwhitneyu(ga, gb, alternative="two-sided")
                rows_g.append({"tf": tf, "gene": gene, "split": "positive_vs_negative",
                               "n_high": len(ga), "n_low": len(gb),
                               "mean_high": float(np.mean(ga)), "mean_low": float(np.mean(gb)),
                               "diff": float(np.mean(ga) - np.mean(gb)), "p": pu})

            # 次级分组：仅在阳性细胞内按中位数切
            sub_a = a[ok & pos]
            if len(sub_a) >= 40:
                med = np.median(sub_a)
                ga2, gb2 = sub_a[sub_a >= med], sub_a[sub_a < med]
                if len(ga2) >= 20 and len(gb2) >= 20:
                    _, pu2 = stats.mannwhitneyu(ga2, gb2, alternative="two-sided")
                    rows_g.append({"tf": tf, "gene": gene, "split": "within_positive_median",
                                   "n_high": len(ga2), "n_low": len(gb2),
                                   "mean_high": float(np.mean(ga2)),
                                   "mean_low": float(np.mean(gb2)),
                                   "diff": float(np.mean(ga2) - np.mean(gb2)), "p": pu2})
        c = pd.DataFrame(rows_c)
        g = pd.DataFrame(rows_g)
        if len(c):
            c["fdr"] = stats.false_discovery_control(c["p"], method="bh")
        if len(g):
            g["fdr"] = stats.false_discovery_control(g["p"], method="bh")
        log(f"  [分组] {gene}: 阳性细胞 {frac_pos:.1%}；生成 {len(g)} 条分组检验")
        return c, g

    corr_j, grp_j = analyse("JAML")
    corr_p, grp_p = analyse("PDCD1")
    if corr_j is None or corr_p is None:
        sys.exit("JAML 或 PDCD1 不在基因集里")

    both = pd.concat([corr_j, corr_p], ignore_index=True)
    both.to_csv(os.path.join(RES, "m8_tf_vs_gene.csv"), index=False)
    both_g = pd.concat([grp_j, grp_p], ignore_index=True)
    both_g.to_csv(os.path.join(RES, "m8_tf_group_diff.csv"), index=False)

    log("\n=== JAML：TF 活性与 JAML 表达相关（前 15 正 / 前 15 负）===")
    cj = corr_j.sort_values("rho", ascending=False)
    log(cj.head(15)[["tf", "rho", "p", "fdr"]].round(4).to_string(index=False))
    log("...")
    log(cj.tail(15)[["tf", "rho", "p", "fdr"]].round(4).to_string(index=False))

    log("\n=== JAML 阳性 vs 阴性 细胞的 TF 活性差异（FDR<0.05，按 |diff| 降序，前 20）===")
    gj = both_g[(both_g["gene"] == "JAML") & (both_g["split"] == "positive_vs_negative")
                & (both_g["fdr"] < 0.05)].copy()
    gj = gj.reindex(gj["diff"].abs().sort_values(ascending=False).index)
    log(gj.head(20)[["tf", "mean_low", "mean_high", "diff", "p", "fdr"]].round(4)
        .to_string(index=False))
    n_j_all = int((both_g["gene"] == "JAML").sum())
    log(f"（JAML：FDR<0.05 的 TF 共 {len(gj)}；参与检验的 TF-分组对 {n_j_all}）")

    # ---------------- 4. 阳性对照：PDCD1 ----------------
    log("\n" + "=" * 70)
    log("阳性对照：PDCD1（应富集耗竭相关 TF）")
    log("=" * 70)
    gp = both_g[(both_g["gene"] == "PDCD1") & (both_g["split"] == "positive_vs_negative")
                & (both_g["fdr"] < 0.05)].copy()
    gp = gp.reindex(gp["diff"].abs().sort_values(ascending=False).index)
    log(gp.head(20)[["tf", "mean_low", "mean_high", "diff", "p", "fdr"]].round(4)
        .to_string(index=False))
    hit = [t for t in PDCD1_EXPECTED if t in set(gp["tf"])]
    log(f"\n[阳性对照判读] 期望 TF {PDCD1_EXPECTED}")
    log(f"  在 PDCD1 高/低分组差异中 FDR<0.05 命中的: {hit}")
    log(f"  命中率 {len(hit)}/{len([t for t in PDCD1_EXPECTED if t in set(A.index)])}"
        f"（分母=期望 TF 中实际出现在活性矩阵里的个数）")
    ctrl_ok = len(hit) >= 3
    log(f"  → 阳性对照{'通过' if ctrl_ok else '未通过'}")

    # JAML 与 PDCD1 的 TF 谱重叠（想说明 JAML 是否偏向耗竭程序）
    set_j = set(gj["tf"]); set_p = set(gp["tf"])
    ov = sorted(set_j & set_p)
    log(f"\n[重叠] JAML 与 PDCD1 的显著 TF 交集 {len(ov)} 个: {ov[:25]}")
    log(f"  JAML {len(set_j)} 个；PDCD1 {len(set_p)} 个")

    # ---------------- 5. 阴性对照：随机基因 ----------------
    log("\n" + "=" * 70)
    log("阴性对照：随机置换基因标签（应无可重复的 TF 数）")
    log("=" * 70)
    rng = np.random.default_rng(SEED)
    n_sig_perm = []
    v0 = X["JAML"].values.astype(float)
    for i in range(N_PERM):
        v = rng.permutation(v0)
        posp = v > 0                      # 与主分析同款分组（阳性 vs 阴性）
        if posp.sum() < 20 or (~posp).sum() < 20:
            continue
        ps = []
        for tf in A.index:
            a = A.loc[tf].values
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
    log(f"  置换 {N_PERM} 次，每次 FDR<0.05 的 TF 数：均值 {n_sig_perm.mean():.1f}，"
        f"最大 {n_sig_perm.max()}，95 分位 {np.percentile(n_sig_perm, 95):.1f}")
    log(f"  实测 JAML 的显著 TF 数 = {obs_n}；经验 P = {pct:.4f}")
    nc = {"n_perm": int(N_PERM), "obs_n_sig": int(obs_n),
          "perm_mean": float(n_sig_perm.mean()), "perm_max": int(n_sig_perm.max()),
          "perm_p95": float(np.percentile(n_sig_perm, 95)), "empirical_p": pct}
    with open(os.path.join(RES, "m8_negative_control.json"), "w", encoding="utf-8") as f:
        json.dump(nc, f, ensure_ascii=False, indent=2)

    # ---------------- 6. 沿亚型的趋势 ----------------
    log("\n=== TF 活性沿 CD4 亚型的趋势（Spearman，仅列 JAML 显著 TF 的前 25）===")
    sub_code = sub.map({s: i for i, s in enumerate(SUBTYPE_ORDER)})
    keep_sub = sub_code.notna().values
    rows = []
    for tf in (list(gj["tf"])[:200] if len(gj) else list(A.index[:200])):
        a = A.loc[tf].values[keep_sub]
        c = sub_code.values[keep_sub]
        if len(a) < 50 or np.nanstd(a) == 0:
            continue
        rho, p = stats.spearmanr(c, a)
        rows.append({"tf": tf, "rho_vs_subtype_order": rho, "p": p, "n": int(len(a))})
    trend = pd.DataFrame(rows)
    if len(trend):
        trend["fdr"] = stats.false_discovery_control(trend["p"], method="bh")
        trend = trend.sort_values("rho_vs_subtype_order", ascending=False)
        trend.to_csv(os.path.join(RES, "m8_tf_by_subtype.csv"), index=False)
        log(trend.head(25).round(4).to_string(index=False))

    # ---------------- 7. 图 ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def save(fig, name):
        fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(FIG, name + ".tif"), dpi=300, bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"})
        plt.close(fig)
        log(f"[图] results/figures/{name}.png / .tif")

    top_j = gj.head(12)["tf"].tolist() if len(gj) else []
    top_p = gp.head(12)["tf"].tolist() if len(gp) else []
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.2))

    ax = axes[0]
    d = cj.head(20).iloc[::-1]
    ax.barh(range(len(d)), d["rho"],
            color=["#D73027" if r > 0 else "#4575B4" for r in d["rho"]],
            edgecolor="#333", linewidth=0.4)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d["tf"], fontsize=7)
    ax.axvline(0, color="#666", lw=0.6)
    ax.set_xlabel("Spearman \u03c1 (TF activity vs JAML)")
    ax.set_title(f"TFs tracking JAML expression\n(top 20 by \u03c1; n = {A.shape[1]:,} CD4\u207a T cells)",
                 fontsize=8.6)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    if top_j:
        idx = [list(A.index).index(t) for t in top_j if t in A.index]
        M = A.iloc[idx].values
        posm = X["JAML"].values > 0            # 零膨胀 → 用阳性/阴性分组
        z = (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-12)
        ax.boxplot([z[:, ~posm].ravel(), z[:, posm].ravel()],
                   tick_labels=[f"JAML-negative\n(n={(~posm).sum():,})",
                                f"JAML-positive\n(n={posm.sum():,})"],
                   patch_artist=True, boxprops=dict(facecolor="#c9d7ea"),
                   medianprops=dict(color="#333"))
        ax.set_ylabel("TF activity (z)")
        ax.set_title(f"Top {len(top_j)} of {len(gj)} TFs with FDR < 0.05\n"
                     f"(largest activity difference, JAML-positive vs negative)", fontsize=8.6)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[2]
    if top_p:
        sel = gp.head(12)[["tf", "diff"]].iloc[::-1]
        ax.barh(range(len(sel)), sel["diff"],
                color=["#D73027" if r > 0 else "#4575B4" for r in sel["diff"]],
                edgecolor="#333", linewidth=0.4)
        ax.set_yticks(range(len(sel)))
        ax.set_yticklabels(sel["tf"], fontsize=7)
        ax.axvline(0, color="#666", lw=0.6)
        ax.set_xlabel("\u0394 TF activity (PDCD1-positive \u2212 negative)")
        ax.set_title(f"Positive control: PDCD1\nrecovered {len(hit)}/"
                     f"{len([t for t in PDCD1_EXPECTED if t in set(A.index)])} testable "
                     f"({len(hit)}/{len(PDCD1_EXPECTED)} pre-specified)",
                     fontsize=8.6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "m8_tf_jaml")

    # 阴性对照分布
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    ax.hist(n_sig_perm, bins=30, color="#9ecae1", edgecolor="#333", linewidth=0.4)
    ax.axvline(obs_n, color="#D73027", lw=1.6, label=f"observed = {obs_n}")
    ax.set_xlabel("# TFs with FDR < 0.05")
    ax.set_ylabel("permutations")
    ax.set_title(f"Negative control (gene labels permuted)\n"
                 f"empirical P < {1.0 / (N_PERM + 1):.3f} ({N_PERM} permutations)",
                 fontsize=8.6)
    ax.legend(frameon=False, fontsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "m8_negative_control")

    # ---------------- 8. 汇总 ----------------
    summ = {
        "note": "decoupleR + CollecTRI 替代 pySCENIC；非 SCENIC 本体，不给出 regulon 特异度分数",
        "n_cells": int(A.shape[1]), "n_tf": int(A.shape[0]),
        "n_genes_input": int(adata.n_vars),
        "network": "CollecTRI (fallback: DoRothEA)",
        "jaml": {"n_sig_tf_fdr05": int(obs_n),
                 "top_pos": cj.head(10)[["tf", "rho"]].to_dict("records"),
                 "top_neg": cj.tail(10)[["tf", "rho"]].to_dict("records")},
        "pdcd1_positive_control": {"expected": PDCD1_EXPECTED, "hit": hit,
                                   "n_sig_tf_fdr05": int(len(gp)), "passed": bool(ctrl_ok)},
        "negative_control": nc,
        "overlap_jaml_pdcd1_tfs": ov,
        "elapsed_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(RES, "m8_tf_activity_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=2)
    log(f"\n[输出] results/m8_tf_activity_summary.json")
    log(f"\n总耗时 {time.time()-t0:.0f}s\nDONE")


if __name__ == "__main__":
    main()
