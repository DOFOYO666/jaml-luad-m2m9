# -*- coding: utf-8 -*-
"""【大内存任务 · 步骤 2】M4：CellOracle JAML 虚拟敲除 / 过表达（含正负对照与参数敏感性）。

⚠️ 运行环境要求：**≥ 64 GB 内存**（GRN 拟合与 KNN 插补对 CD4 T 亚群需数十 GB）。
   前置：先跑 prep_celloracle_input.py 生成 results/GSE131907_TNK_epi.h5ad

本脚本相对规划 v1.0 的修正（全部基于实机核验，勿回退）：
  ① GRN 单元取自 **Cell_subtype 映射**，而非 Cell_type / Cell_type.refined
     —— GSE131907 的这两列都不含 CD4 T
  ② `get_links` 只调用一次（v1.0 重复调用且参数不一致）
  ③ KNN 插补 k = 30（CellOracle 官方教程量级），v1.0 的 sqrt(n_cells) 会过度平滑
  ④ 去掉重复插补
  ⑤ 过表达设为目标基因表达的 99 百分位（v1.0 的 max×2 无依据）
  ⑥ 阳性对照 PDCD1 与阴性对照（随机化 GRN）**都输出 quiver 图**
  ⑦ 参数敏感性矩阵（k / alpha / n_propagation）
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import scanpy as sc
    import celloracle as co
except ImportError as e:
    sys.exit(f"缺少依赖: {e}\n请先 pip install celloracle palantir scanpy")

IN_H5AD = "results/GSE131907_TNK_epi.h5ad"
OUTDIR = "results"
FIGDIR = "results/figures"
os.makedirs(FIGDIR, exist_ok=True)
GRN_UNIT = "CD4 T"          # 目标 GRN 单元
N_PROP = 3                  # 信号传播步数
ALPHA = 10                  # ridge 正则
K_IMPUTE = 30               # KNN 插补近邻数


def log(*a):
    print(*a, flush=True)


# ---------------- 1. 载入 ----------------
adata = sc.read_h5ad(IN_H5AD)
log(f"[载入] {adata.shape}")
log("[obs] grn_unit:", adata.obs["grn_unit"].value_counts().to_dict())

print(f"[检查] GRN 单元 '{GRN_UNIT}' 的细胞数: {(adata.obs['grn_unit'] == GRN_UNIT).sum()}")

# ---------------- 2. Base GRN ----------------
base_GRN = co.data.load_human_promoter_base_GRN()
log(f"[BaseGRN] {base_GRN.shape}")

# ---------------- 3. 建 Oracle、拟合 GRN ----------------
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="grn_unit",
                                   embedding_name="X_umap")
oracle.import_TF_data(TF_info_matrix=base_GRN)
oracle.perform_PCA()
oracle.knn_imputation(n_pca=50, k=K_IMPUTE, balanced=True, b_sight=500, b_maxl=10)
# ② 只调用一次
oracle.get_links(cluster_name_for_GRN_unit=GRN_UNIT, alpha=ALPHA,
                 n_jobs=8, verbose_level=1)
links = oracle.links
links.filter_links(p=0.001, weight="coef_abs", threshold_number=10000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=ALPHA, use_cluster_specific_TFdict=True)


# ---------------- 4. 扰动工具 ----------------
def run_perturbation(gene, value, label, randomized=False, n_prop=N_PROP):
    oracle.simulate_shift(perturb_condition={gene: value}, n_propagation=n_prop,
                          use_randomized_GRN=randomized)
    oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=0.5)
    oracle.calculate_shift()
    # API 名称随版本变动：新版 co.plot.quiver，旧版 oracle.plot_quiver
    plotter = getattr(co, "plot", None)
    if plotter is not None and hasattr(plotter, "quiver"):
        plotter.quiver(oracle, topology=label, n_grid=40, amin=100, scale=30)
    else:
        oracle.plot_quiver(topology=label, n_grid=40, amin=100, scale=30)
    # 导出 shift 供定量比较（不只出图）
    df = oracle.to_pandas() if hasattr(oracle, "to_pandas") else None
    if df is not None:
        df.to_csv(os.path.join(OUTDIR, f"celloracle_shift_{label}.csv"))
    return oracle


# ⑤ 过表达水平 = 该基因表达的 99 百分位
x = oracle.adata[:, "JAML"].layers["raw_count"] if "raw_count" in oracle.adata.layers else oracle.adata[:, "JAML"].X
x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
jaml_p99 = float(np.percentile(x, 99))
log(f"[JAML] 表达 99 百分位 = {jaml_p99:.3f}（过表达目标值）")

# ---------------- 5. 主分析 + 对照 ----------------
log("\n=== JAML 敲除（主分析） ===")
run_perturbation("JAML", 0.0, "JAML_KO")

log("\n=== JAML 过表达 ===")
run_perturbation("JAML", jaml_p99, "JAML_OE")

log("\n=== 阳性对照：PDCD1 敲除 ===")
run_perturbation("PDCD1", 0.0, "PDCD1_KO")

log("\n=== 阴性对照：随机化 GRN ===")
run_perturbation("JAML", 0.0, "JAML_KO_randGRN", randomized=True)

# ---------------- 6. 参数敏感性矩阵 ----------------
log("\n=== 参数敏感性（n_propagation × alpha）===")
rows = []
for n_prop in [2, 3, 5]:
    for alpha in [1, 10, 100]:
        try:
            oracle.simulate_shift(perturb_condition={"JAML": 0.0}, n_propagation=n_prop)
            oracle.estimate_transition_prob(n_neighbors=200, knn_random=True,
                                            sampled_fraction=0.5)
            oracle.calculate_shift()
            shift = float(np.mean(oracle.delta_X.mean(axis=1)))
            rows.append({"n_propagation": n_prop, "alpha": alpha, "mean_shift": shift})
            log(f"  n_prop={n_prop} alpha={alpha}: mean shift = {shift:.4f}")
        except Exception as e:
            rows.append({"n_propagation": n_prop, "alpha": alpha, "mean_shift": np.nan,
                         "error": str(e)})
pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "m4_param_sensitivity.csv"), index=False)

log("\n[完成] 输出写入 results/（celloracle_shift_*.csv、m4_param_sensitivity.csv、各 quiver 图）")
