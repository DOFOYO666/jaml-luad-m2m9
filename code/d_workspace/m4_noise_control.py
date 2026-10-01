# -*- coding: utf-8 -*-
"""M4 配套：配对噪声对照（proper randomized-GRN control）

动机（自查发现的设计缺陷）
--------------------------
`m4_celloracle_run.py` 的阴性对照取 `TF_CANDIDATES[:3]` 作为"被打乱 GRN"的基准 TF。
在 ext3074 那一轮里，这个名单前三位恰好是 MYB / RORC / EOMES —— 它们在本 CD4⁺T
子集中检出率仅约 1%，**其真扰动本身就无效应**，因此"打乱 GRN"的读数自然也是 0，
使信噪比被人为抬高、且与 ext3013 那一轮（基准为 TCF7/LEF1/KLF2，检出 25%+）不可比。

正确做法：**同一 TF 作配对比较** —— 对同一个（高表达的）TF 分别跑
  ① 真 GRN 扰动  ② 打乱 GRN 扰动（重复多次取分布）
两者之差才是该 TF 的真实效应超出噪声的部分。

用法：python m4_noise_control.py
输出：results/m4_noise_control.csv  /  m4_noise_control_summary.json
      results/figures/m4_noise_control.{png,tif}
"""
import os
import sys
import gc
import json
import time
import ctypes

import numpy as np
import pandas as pd

t0 = time.time()

D_ROOT = r"D:\workbuddy工作空间\JAML深度研究"
RES = os.path.join(D_ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

INPUT_FILE = "CD4T_celloracle_ext2.h5ad"
INPUT_TAG = "ext3074"
CLUSTER_COL = "cell_subtype"
ALPHA = 10
BAGGING = 20
N_PROP = 3
K_IMPUTE = 30
N_PCA_DIMS = 50
FILTER_P = 0.05
FILTER_TOPN = 20000

# 默认覆盖"曾被列为候选"的全部 TF（含 STAT3 —— 第五轮报告曾把它列为次级候选但未做配对检验）
TFS = ["ID2", "GATA3", "STAT3", "MAF", "IRF1", "STAT1", "NFKB1", "BATF", "IRF4"]
# 可用 `--tfs A,B,C` 覆盖
if "--tfs" in sys.argv:
    TFS = [t.strip() for t in sys.argv[sys.argv.index("--tfs") + 1].split(",") if t.strip()]
N_RAND = 5                                          # 每个 TF 的打乱-GRN 重复次数

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def log(*a):
    print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)


def avail_gb():
    class M(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    m = M(); m.dwLength = ctypes.sizeof(M)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ullAvailPhys / 1073741824


def guard(need, where):
    a = avail_gb()
    log(f"[内存] {where} 可用 {a:.2f} GB")
    if a < need:
        sys.exit(f"可用内存 {a:.2f} GB < {need} GB，拒绝继续（{where}）")


DEGENERATE_RAND_MAX = 1e-4      # 打乱读数低于此值 ⇒ 随机化对照本身退化（分母近 0）


def classify(out):
    """双判据分类：配对信噪比 + 绝对效应量，并标出「退化分母」。

    为什么不能只看信噪比：若某 TF 本身几乎不表达（如 IRF4 检出 3.7%），
    打乱其 GRN 后传播不出任何变化 → 随机读数被压到 ~1e-5，
    比值（信噪比）会被抬到 80+ —— 这是**与主脚本同一类**的假象。
    """
    d = out.copy()
    d["degenerate_control"] = d["jaml_rand_max"] < DEGENERATE_RAND_MAX
    d["passes_paired"] = d["jaml_real"] > d["jaml_rand_max"]
    ok = d[~d["degenerate_control"]]
    ceiling = float(ok["jaml_rand_max"].max()) if len(ok) else np.nan
    d["global_noise_ceiling"] = ceiling
    d["passes_global_ceiling"] = d["jaml_real"] > ceiling
    d = d.sort_values("jaml_real", ascending=False)

    print("\n=== 双判据分类（★ 表示通过）===")
    hdr = (f"{'TF':8s} {'|ΔJAML|':>10s} {'配对打乱max':>12s} {'配对SNR':>9s} "
           f"{'退化?':>6s} {'稳健':>6s} {'超全局上界':>10s}")
    print(hdr)
    for _, r in d.iterrows():
        print(f"{r['TF']:8s} {r['jaml_real']:10.6f} {r['jaml_rand_max']:12.6f} "
              f"{r['snr_vs_randmax']:9.2f} "
              f"{('是' if r['degenerate_control'] else '-'):>6s} "
              f"{('★' if r['passes_paired'] else '-'):>6s} "
              f"{('★' if r['passes_global_ceiling'] else '-'):>10s}")
    print(f"\n全局噪声上界（仅用非退化 TF 的打乱最大值）= {ceiling:.6f}")
    print(f"配对检验通过（非退化）: "
          f"{sorted(d[(d['passes_paired']) & (~d['degenerate_control'])]['TF'])}")
    print(f"配对检验通过但对照退化（信噪比不可信）: "
          f"{sorted(d[(d['passes_paired']) & (d['degenerate_control'])]['TF'])}")
    print(f"超过全局上界者: {sorted(d[d['passes_global_ceiling']]['TF'])}")
    d.to_csv(os.path.join(RES, "m4_noise_control_classify.csv"), index=False)
    with open(os.path.join(RES, "m4_noise_control_classify.json"), "w",
              encoding="utf-8") as f:
        json.dump({"degenerate_rand_max_threshold": DEGENERATE_RAND_MAX,
                   "global_noise_ceiling_nondegenerate": ceiling,
                   "pass_paired": sorted(d[d["passes_paired"]]["TF"]),
                   "pass_paired_nondegenerate":
                       sorted(d[(d["passes_paired"]) & (~d["degenerate_control"])]["TF"]),
                   "degenerate_controls": sorted(d[d["degenerate_control"]]["TF"]),
                   "pass_global_ceiling":
                       sorted(d[d["passes_global_ceiling"]]["TF"]),
                   "note": ("信噪比在随机化对照退化（打乱读数 < 1e-4）时不可用——"
                            "低表达 TF 的 GRN 打乱后传播不出变化，会使比值虚高。"
                            "必须同时看绝对效应量与全局噪声上界。")},
                  f, ensure_ascii=False, indent=2)
    print("[输出] results/m4_noise_control_classify.csv / .json")
    return d


def make_figure(out):
    """由结果表绘制配对对照图（与计算解耦，便于 --fig-only 重绘）。"""
    d = out.copy()
    if "degenerate_control" not in d.columns:
        d["degenerate_control"] = d["jaml_rand_max"] < DEGENERATE_RAND_MAX
    d = d.iloc[::-1]
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    y = np.arange(len(d))
    for i, (rv, dg) in enumerate(zip(d["jaml_real"], d["degenerate_control"])):
        ax.barh(i, rv, height=0.55, color="#C0392B",
                hatch="xx" if dg else None, edgecolor="#333", linewidth=0.4,
                alpha=0.55 if dg else 1.0, zorder=2)
    ax.errorbar(d["jaml_rand_mean"], y, xerr=d["jaml_rand_sd"], fmt="o", ms=4,
                color="#2C3E50", capsize=2.5, lw=0.9, zorder=3,
                label=f"randomized GRN (n={N_RAND}, mean ± SD)")
    xmax = max(float(d["jaml_real"].max()),
               float((d["jaml_rand_mean"] + d["jaml_rand_sd"]).max()))
    for i, (snr, dg) in enumerate(zip(d["snr_vs_randmax"], d["degenerate_control"])):
        lab = f"SNR {snr:.1f}×" + ("  (control degenerate)" if dg else "")
        ax.text(xmax * 1.03, i, lab, va="center", ha="left", fontsize=7,
                color="#B03A2E" if dg else "#222")
    ax.set_xlim(0, xmax * 1.34)
    ax.barh([], [], color="#C0392B", hatch="xx", edgecolor="#333",
            label="knockout of a poorly expressed TF (control degenerate)")
    ax.set_yticks(y)
    ax.set_yticklabels(d["TF"], fontsize=8)
    ax.set_xlabel(r"mean $|\Delta$ JAML$|$ after knockout")
    ax.set_title("Paired randomized-GRN control (same TF, same cells)\n"
                 "CD4$^+$ T lineage, GSE131907 subset (n = 2,842)", fontsize=8.8)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=7.2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=[0, 0.09, 1, 1])
    for ext, kw in [("png", {}), ("tif", dict(pil_kwargs={"compression": "tiff_lzw"}))]:
        fig.savefig(os.path.join(FIG, f"m4_noise_control.{ext}"), dpi=300, **kw)
    plt.close(fig)
    print("[图] results/figures/m4_noise_control.png / .tif")


def main():
    import scanpy as sc
    import anndata as ad
    import scipy.sparse as sp
    import celloracle as co
    from celloracle.network_analysis.links_object import load_links

    log(f"配对噪声对照（celloracle {co.__version__}）")
    guard(0.40, "起始")

    H5AD = os.path.join(RES, INPUT_FILE)
    log(f"[1] 载入 {H5AD}")
    adata = ad.read_h5ad(H5AD)
    X = adata.X
    if sp.issparse(X):
        adata.X = X.astype(np.int32)
    log(f"    {adata.shape}")

    log("[2] UMAP")
    adn = adata.copy()
    sc.pp.normalize_total(adn, target_sum=1e4)
    sc.pp.log1p(adn)
    sc.tl.pca(adn, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adn, n_neighbors=15)
    sc.tl.umap(adn, min_dist=0.3)
    adata.obsm["X_umap"] = np.asarray(adn.obsm["X_umap"]).copy()
    del adn
    gc.collect()
    guard(0.40, "UMAP 之后")

    log("[3] base GRN → TFdict")
    GRN_PATH = os.path.join(os.path.expanduser("~"), "celloracle_data", "promoter_base_GRN",
                            "hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet")
    _gl = list(adata.var_names)
    _sub = pd.read_parquet(GRN_PATH, filters=[("gene_short_name", "in", _gl)])
    _sub = _sub.drop(columns=["peak_id"]).groupby("gene_short_name").sum()
    TF_COLS = [str(c) for c in _sub.columns]
    TFdict = {str(g): list(np.asarray(TF_COLS)[r > 0])
              for g, r in zip(_sub.index.to_numpy(), _sub.to_numpy())}
    del _sub
    for _ in range(2):
        gc.collect()
    log(f"    TFdict {len(TFdict)} 基因")
    guard(0.25, "base GRN 之后")

    log("[4] Oracle 流程")
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name=CLUSTER_COL,
                                       embedding_name="X_umap")
    oracle.import_TF_data(TFdict=TFdict)
    oracle.perform_PCA()
    oracle.knn_imputation(k=K_IMPUTE, n_pca_dims=N_PCA_DIMS, balanced=True,
                          b_sight=500, b_maxl=10)
    gc.collect()
    guard(0.30, "knn_imputation 之后")

    CACHE = os.path.join(RES, f"m4_links_{INPUT_TAG}_alpha{ALPHA}_bag{BAGGING}.celloracle.links")
    if not os.path.isfile(CACHE):
        sys.exit(f"缺少 links 缓存 {CACHE} —— 请先运行 m4_celloracle_run.py --ext2")
    log(f"    载入 links 缓存 {CACHE}")
    links = load_links(CACHE)
    links.filter_links(p=FILTER_P, weight="coef_abs", threshold_number=FILTER_TOPN)
    log(f"    已 filter_links(p={FILTER_P}, top{FILTER_TOPN})")
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    # ⚠️ 必须有这一步：coef_matrix_per_cluster 由 fit_GRN_for_simulation() 设置，
    # simulate_shift 依赖它（缺则 AttributeError）。加载缓存的 links 不会自动补上。
    oracle.fit_GRN_for_simulation(alpha=ALPHA, use_cluster_specific_TFdict=True)
    log("    fit_GRN_for_simulation 完成")
    guard(0.30, "links 载入之后")

    GIDX = {g: i for i, g in enumerate(oracle.adata.var_names)}
    JI = GIDX["JAML"]

    def run(tf, randomized, seed=None):
        if randomized:
            # ⚠️ CellOracle 会把打乱后的系数矩阵**缓存**在
            # `coef_matrix_per_cluster_randomized` 上（`calculate_randomized_coef_table`
            # 默认 random_seed=123），若不主动换种子，多次"重复"其实是同一份随机矩阵，
            # 得到的只是同一个数。故每次显式指定新种子。
            oracle.calculate_randomized_coef_table(random_seed=seed)
        oracle.simulate_shift(perturb_condition={tf: 0.0}, GRN_unit="cluster",
                              n_propagation=N_PROP, use_randomized_GRN=randomized,
                              clip_delta_X=True)
        d = oracle.adata.layers["delta_X"]
        d = np.asarray(d.todense()) if sp.issparse(d) else np.asarray(d)
        return d

    rows = []
    log(f"[5] 配对对照：{len(TFS)} 个 TF × (1 真 + {N_RAND} 次独立打乱)")
    for tf in TFS:
        if tf not in GIDX:
            log(f"    {tf} 不在矩阵中 —— 跳过")
            continue
        real = run(tf, False)
        r_jaml = float(np.abs(real[:, JI]).mean())
        r_all = float(np.abs(real).mean())
        rand_j = []
        for k in range(N_RAND):
            dr = run(tf, True, seed=1000 + 97 * k)
            rand_j.append(float(np.abs(dr[:, JI]).mean()))
        rand_j = np.array(rand_j)
        snr = r_jaml / rand_j.max() if rand_j.max() > 0 else np.inf
        rows.append({"TF": tf, "jaml_real": r_jaml, "jaml_real_signed_all": float(real[:, JI].mean()),
                     "allgenes_real": r_all,
                     "jaml_rand_mean": float(rand_j.mean()), "jaml_rand_max": float(rand_j.max()),
                     "jaml_rand_sd": float(rand_j.std(ddof=1)) if len(rand_j) > 1 else 0.0,
                     "snr_vs_randmax": float(snr), "n_rand": int(len(rand_j)),
                     "jaml_rand_values": ";".join(f"{x:.6f}" for x in rand_j)})
        log(f"    {tf:8s} 真扰动 |ΔJAML| = {r_jaml:.6f}（有符号均值 {real[:, JI].mean():+.6f}）"
            f" | 打乱 GRN 均值 {rand_j.mean():.6f}（最大 {rand_j.max():.6f}）"
            f" → 信噪比 {snr:.2f}")

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "m4_noise_control.csv"), index=False)
    log("\n=== 配对噪声对照汇总 ===")
    log(out[["TF", "jaml_real", "jaml_rand_mean", "jaml_rand_max",
             "snr_vs_randmax"]].round(6).to_string(index=False))

    # ---------------- 分类 + 图 ----------------
    classify(out)
    make_figure(out)

    summary = {
        "input": INPUT_FILE, "n_tf": len(out), "n_rand_per_tf": N_RAND,
        "note": ("配对对照：同一 TF 的真扰动 vs 打乱 GRN 重复。修复了主脚本中"
                 "『以低表达 TF 为基准做打乱对照』导致信噪比被人为抬高的问题。"
                 "另注：CellOracle 的打乱系数矩阵默认会被缓存（seed=123），"
                 "故每次重复都显式换种子以避免『同一份随机矩阵重复 5 次』。"),
        "results": out.drop(columns=["jaml_rand_values"]).to_dict(orient="records"),
    }
    with open(os.path.join(RES, "m4_noise_control_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    log("[输出] results/m4_noise_control.csv / m4_noise_control_summary.json")
    log(f"总耗时 {time.time()-t0:.0f}s\nDONE")


if __name__ == "__main__":
    if "--fig-only" in sys.argv:
        # 仅从已有 CSV 做分类与重绘，避免重跑约 14 分钟的模拟
        _out = pd.read_csv(os.path.join(RES, "m4_noise_control.csv"))
        classify(_out)
        make_figure(_out)
        print("DONE (post-hoc)")
    else:
        main()
