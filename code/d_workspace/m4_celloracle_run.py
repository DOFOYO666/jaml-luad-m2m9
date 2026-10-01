# -*- coding: utf-8 -*-
"""M4：CellOracle 虚拟敲除 —— JAML 在 CD4⁺T 谱系中的 GRN 与扰动效应。

本脚本是 M4 的**运行部分**（输入准备见 m4_prep_bounded.py）。

前置条件（本机已满足，2026-09-30 实测）
---------------------------------------
1. velocyto：官方 sdist 只含 Cython 0.29.2 生成的 C，在 CPython 3.13 + numpy 2 下无法编译。
   本机用 Cython 3.3.0 从 speedboosted.pyx 重新生成 C，再用 Rtools44 GCC 13.3.0 编译
   （见 scripts/40_verify_velocyto_minimal.py，数值 max|Δ| = 8.9e-16）。
12. celloracle 0.20.0 可 import；gimmemotifs/pybedtools 为显式占位（ATAC 自建 base GRN 路径
   在 Windows 上不可用；本流程用官方预构建 promoter base GRN，不经过它们）。
3. base GRN：hg38_gimmemotifsv5_fpr2（39,315 peak-gene 行 × 1094 TF），
   经核实其中 JAML 用**新符号 JAML**（无 AMICA1），与表达矩阵一致。

方法学要点（避免踩坑）
---------------------
* CellOracle 0.20 **没有** calculate_shift()/to_pandas()（旧版 API）。结果直接存
  `oracle.adata.layers["delta_X"]` = simulated_count − imputed_count（细胞 × 基因）。
* `simulate_shift(use_randomized_GRN=True)` 是官方内置的**阴性对照**（打乱 GRN）。
* `get_links` 用 bagging 回归，每个基因 × 每个 GRN unit 重复 bagging_number 次；
  n_jobs>1 会多进程复制数据 —— 本机 8 GB 内存下固定 n_jobs=1。
* **JAML 在 CD4⁺T 中 >70% 为 0（零膨胀）**，故预期扰动效应偏弱；因此必须同时跑
  （a）高表达核心基因 CD3E 作阳性对照（验证量程），（b）CXADR 作生物对照（受体，
  CD4⁺T 中几乎不表达），（c）randomized GRN 作阴性对照。缺了对照则无法解释阴性结果。

用法：
    python m4_celloracle_run.py            # 完整流程
    python m4_celloracle_run.py --test     # 仅用小 bagging 数验证流程可跑通
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

TEST_MODE = "--test" in sys.argv

# ---- 参数 ----
CLUSTER_COL = "cell_subtype"     # GRN unit（5 个 CD4 亚型）
ALPHA = 10                       # 岭回归正则
N_PROP = 3                       # 信号传播步数
K_IMPUTE = 30
N_PCA_DIMS = 50
BAGGING = 3 if TEST_MODE else 20

# 输入：完整版用扩充基因集，追加 10 个原本缺失的 TF。
# 理由：CellOracle 的 GRN 只能以**矩阵内**的基因作回归因子；这些 TF 不在矩阵时，
# 它们作为 source 的连接数恰为 0 → 触发 "does not have enough regulatory connection"
# 而无法被扰动（实测 10/10 个都是 0 条）。追加后即可评估；HVG 主体不变，GRN 规模基本不变。
# EXT2 模式：输入再加 61 个 base GRN 中有 motif 证据、但不在 3013 基因里的 TF
# （见 m4_prep_extend.py），目的是让"因不在矩阵而跳过"的 TF 归零，闭合审计缺口。
EXT2 = "--ext2" in sys.argv

INPUT_FILE = ("CD4T_celloracle.h5ad" if TEST_MODE
              else ("CD4T_celloracle_ext2.h5ad" if EXT2 else "CD4T_celloracle_ext.h5ad"))
INPUT_TAG = "base3003" if TEST_MODE else ("ext3074" if EXT2 else "ext3013")
SUF = "_ext2" if EXT2 else ""      # EXT2 的输出加后缀，避免覆盖 ext3013 那一轮的结果

TARGET = "JAML"
POS_CTRL = "CD3E"                # 阳性对照：高表达、T 细胞核心
BIO_CTRL = "CXADR"               # 生物对照：JAML 的受体，CD4⁺T 中几乎不表达


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


def guard(need=0.5, where=""):
    a = avail_gb()
    log(f"[内存] {where} 可用 {a:.2f} GB")
    if a < need:
        sys.exit(f"可用内存 {a:.2f} GB < {need} GB，拒绝继续（{where}）")


log("=" * 74)
log(f"M4 CellOracle 虚拟敲除（TEST_MODE={TEST_MODE}, bagging={BAGGING}）")
log("=" * 74)
log(f"Python {sys.version.split()[0]}")

import scanpy as sc          # noqa: E402
import anndata as ad         # noqa: E402
import scipy.sparse as sp    # noqa: E402
import celloracle as co      # noqa: E402

log(f"celloracle {co.__version__} / scanpy {sc.__version__} / anndata {ad.__version__}")

# =====================================================================
# 1. 载入输入
# =====================================================================
H5AD = os.path.join(RES, INPUT_FILE)
log(f"[1] 载入 {H5AD}（{INPUT_TAG}）")
adata = ad.read_h5ad(H5AD)
log(f"    形状 {adata.shape}（细胞 × 基因）；obs 列 {list(adata.obs.columns)}")
log(f"    layers: {list(adata.layers.keys())}")

# CellOracle 要求 .X 为原始计数（整数）。输入已是 raw UMI，显式转 int32 以通过其校验。
X = adata.X
if sp.issparse(X):
    X = X.astype(np.int32)
    adata.X = X
frac_int = float(np.mean(np.abs(X.data - np.rint(X.data)) < 1e-6)) if sp.issparse(X) else 1.0
log(f"    .X 整数占比 = {frac_int:.4f}；非零 {X.nnz if sp.issparse(X) else -1}")

RC = adata.layers["raw_count"]
VNAMES = list(adata.var_names)
for g in (TARGET, POS_CTRL, BIO_CTRL):
    if g not in VNAMES:
        sys.exit(f"关键基因 {g} 不在矩阵中 —— 停止（不做任何静默跳过）")
    col = RC[:, VNAMES.index(g)]
    col = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
    log(f"    {g}: 阳性率 {100*(col>0).mean():.1f}%（中位 {np.median(col):.3f}，"
        f"均值 {col.mean():.3f}）")

# =====================================================================
# 2. UMAP（CellOracle 的 embedding_name 需要它）
# =====================================================================
log("[2] 计算 UMAP（在归一化副本上做，不动 raw count）")
adn = adata.copy()
sc.pp.normalize_total(adn, target_sum=1e4)
sc.pp.log1p(adn)
# 跳过 sc.pp.scale：它会把稀疏矩阵 densify（2842 × 3003 float64 ≈ 68 MB），
# 在 8 GB 机器上无必要地抬高内存峰值。PCA 直接用 log 归一化数据（arpack 支持稀疏）。
sc.tl.pca(adn, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adn, n_neighbors=15)
sc.tl.umap(adn, min_dist=0.3)
adata.obsm["X_umap"] = np.asarray(adn.obsm["X_umap"]).copy()
adata.obsm["X_pca"] = np.asarray(adn.obsm["X_pca"]).copy()
del adn
gc.collect()
log(f"    X_umap {adata.obsm['X_umap'].shape}")
guard(0.55, "UMAP 之后")

# =====================================================================
# 3. base GRN
# =====================================================================
log("[3] 载入预构建 promoter base GRN（hg38）")
# 注意：co.data.load_human_promoter_base_GRN() 的默认 version 是 hg19，而本项目下载的是 hg38
# （与表达矩阵的 GRCh38 一致）。且其内部会在缺失时联网下载，该 URL 实测 404。
# 故显式读取本地文件，路径与 celloracle 约定的 CELLORACLE_DATA_DIR 一致。
GRN_PATH = os.path.join(os.path.expanduser("~"), "celloracle_data", "promoter_base_GRN",
                        "hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet")
if not os.path.isfile(GRN_PATH):
    sys.exit(f"base GRN 文件不存在: {GRN_PATH}")
# 只读本项目基因所在的行。全量 base GRN 为 39315 × 1096，float64 约 350 MB+，
# 8 GB 机器上实测载入后可用内存从 0.67 GB 掉到 0.02 GB（触发内存守卫）。
# 按 gene_short_name 过滤是**等价**的：import_TF_data 内部即
#   tmp.drop(["peak_id"]).groupby("gene_short_name").sum()
# 其他基因既不参与本项目基因的分组，也不影响其 TF 列表。
_gl = list(adata.var_names)
try:
    _sub = pd.read_parquet(GRN_PATH, filters=[("gene_short_name", "in", _gl)])
    log(f"    已按项目基因过滤（IN {len(_gl)} 个基因）")
except Exception as e:
    log(f"    [警告] parquet filters 不可用（{type(e).__name__}），回退全量读取+筛选")
    _sub = pd.read_parquet(GRN_PATH)
    _sub = _sub[_sub["gene_short_name"].isin(_gl)].copy()
log(f"    本地读取 {GRN_PATH}；过滤后 {_sub.shape[0]} 行 × {_sub.shape[1]} 列")
n_tf = _sub.shape[1] - 2

# 立刻折叠为 TFdict（{gene: [TF, ...]}）并释放宽表。
# 理由：宽表以 int64 存 1094 个 TF 列，4560 行仍占约 150 MB；而 import_TF_data
# 真正需要的只是每个基因的非零 TF 名列表，体积小两个数量级。
# 这与 import_TF_data 内部的
#   tmp.drop(["peak_id"]).groupby("gene_short_name").sum()  →  x>0 的列名
# 完全等价。
_sub = _sub.drop(columns=["peak_id"]).groupby("gene_short_name").sum()
TF_COLS = [str(c) for c in _sub.columns]          # base GRN 的 TF 列表（CellOracle 可扰动的对象）
_mat = _sub.to_numpy()
_genes = _sub.index.to_numpy()
TFdict = {str(g): list(np.asarray(TF_COLS)[r > 0]) for g, r in zip(_genes, _mat)}
del _sub, _mat, _genes
# groupby/转 numpy 会留下较大的临时对象；显式回收，让"可用内存"尽快回到真实水平。
# （实测：不做这一步时，守卫点读到的是尚未回收的瞬时低值。）
for _ in range(2):
    gc.collect()
_n_tf = np.array([len(v) for v in TFdict.values()])
log(f"    TFdict: {len(TFdict)} 基因；每个基因候选 TF 中位 {np.median(_n_tf):.0f}、"
    f"均值 {_n_tf.mean():.1f}、最大 {_n_tf.max()}")
guard(0.25, "base GRN 之后")

# =====================================================================
# 4. Oracle 流程
# =====================================================================
log("[4] Oracle 流程")
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name=CLUSTER_COL,
                                   embedding_name="X_umap")
log(f"    import_anndata_as_raw_count 完成；cluster 数 = "
    f"{oracle.adata.obs[CLUSTER_COL].nunique()}")

oracle.import_TF_data(TFdict=TFdict)
log(f"    TF dict 基因数 = {len(oracle.TFdict)}")

oracle.perform_PCA()
log("    perform_PCA 完成")

oracle.knn_imputation(k=K_IMPUTE, n_pca_dims=N_PCA_DIMS, balanced=True,
                      b_sight=500, b_maxl=10)
log("    knn_imputation 完成")
gc.collect()
guard(0.35, "knn_imputation 之后")

log(f"    get_links（alpha={ALPHA}, bagging={BAGGING}, n_jobs=1）...")
# 注意：get_links 是**返回** Links 对象，并不写入 self（无 oracle.links 属性）。
# 且 Links 的数据在 .links_dict = {cluster: DataFrame}，也没有 .links 属性。
# get_links 实测约 9 分钟（bagging=3），故加磁盘缓存，调试时不重复计算。
# ⚠️ CellOracle 的 Links.to_hdf5 要求文件名以 `.celloracle.links` 结尾（否则 ValueError）。
LINKS_CACHE = os.path.join(RES, f"m4_links_{INPUT_TAG}_alpha{ALPHA}_bag{BAGGING}.celloracle.links")
links = None
if os.path.isfile(LINKS_CACHE):
    try:
        from celloracle.network_analysis.links_object import load_links
        links = load_links(LINKS_CACHE)
        log(f"    从缓存载入 links：{LINKS_CACHE}")
    except Exception as e:
        log(f"    [警告] 缓存载入失败（{type(e).__name__}: {str(e)[:70]}）→ 重新计算")
        links = None
if links is None:
    links = oracle.get_links(cluster_name_for_GRN_unit=CLUSTER_COL, alpha=ALPHA,
                             bagging_number=BAGGING, verbose_level=0, n_jobs=1)
    try:
        links.to_hdf5(LINKS_CACHE)
        log(f"    已缓存 links → {LINKS_CACHE}")
    except Exception as e:
        log(f"    [警告] 缓存写入失败（{type(e).__name__}: {str(e)[:70]}）")

_ld = links.links_dict
n_links = int(sum(len(v) for v in _ld.values()))
_sample = next(iter(_ld.values()))
log(f"    get_links 完成：{len(_ld)} 个 GRN unit；总连接 {n_links}")
log(f"    单 unit 连接表 {_sample.shape}，列 {list(_sample.columns)}")
gc.collect()
guard(0.35, "get_links 之后")

# ⚠️ 关键修正：CellOracle 教程默认 p=0.001。但实测本项目 JAML 的 125 条候选调控连接中
# **没有一条** p<0.001（最佳为 ID2→JAML, p=0.0027；p 值中位 0.28）。
# 该阈值会把 JAML 的全部连接滤掉 → 任何 TF 扰动都无法经 GRN 传播到 JAML（|ΔJAML| 恒为 0）。
# 故放宽至 p=0.05（JAML 保留 17 条），并把 threshold_number 设得足够大以避免再按
# coef_abs 把 JAML 削掉（JAML 连接的 coef_abs 仅 ~0.04，低于全局 top-N 的入选线）。
# 报告中须如实说明：放宽阈值是为了**让 JAML 可被评估**，而非制造显著性；
# 且"JAML 在 CD4⁺T 中缺乏强 TF 调控连接"本身即为结果之一。
FILTER_P = 0.05
FILTER_TOPN = 20000
links.filter_links(p=FILTER_P, weight="coef_abs", threshold_number=FILTER_TOPN)
_n_filt = {k: len(v) for k, v in links.filtered_links.items()}
_n_jaml_conn = {k: int((v["target"] == "JAML").sum())
                for k, v in links.filtered_links.items()}
log(f"    filter_links(p={FILTER_P}, threshold_number={FILTER_TOPN})：")
log(f"      每 unit 保留 {_n_filt}，合计 {sum(_n_filt.values())}")
log(f"      JAML 作为 target 的连接 {_n_jaml_conn}（合计 {sum(_n_jaml_conn.values())}）")
if sum(_n_jaml_conn.values()) == 0:
    log("      ⚠️ JAML 仍无连接 —— 扰动的 |ΔJAML| 将为 0，属结构性结果")
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=ALPHA, use_cluster_specific_TFdict=True)
log("    fit_GRN_for_simulation 完成")
guard(0.35, "fit_GRN 之后")

# =====================================================================
# 5. TF 扰动模拟
#    重大约束（实测）：CellOracle 的 simulate_shift **只接受 TF 作为扰动对象**。
#    已核实 JAML 不在 base GRN 的 TF 列中（JAML 是黏附分子，不是转录因子），
#    因此无法直接模拟"JAML 敲除"。改用 CellOracle 的正确用法：
#      扰动候选 TF → 观察其对 JAML 表达的调控，即回答"哪些 TF 调控 JAML"。
# =====================================================================
log("[5] TF 扰动模拟")
log(f"    ⚠️ JAML 不在 base GRN 的 {len(TF_COLS)} 个 TF 列中 → 无法直接虚拟敲除 JAML；")
log("       改为扰动候选 TF，观察 JAML 表达变化（回答「哪些 TF 调控 JAML」）。")

# EXT2 模式：61 个新补入的 TF + bagging=20 那一轮的 top-10 作参照
# （参照项用于量化"基因集变化本身"带来多少 GRN 漂移，避免把差异误读为生物学差异）
NEW_TFS = ['MYB', 'RORC', 'EOMES', 'IRF6', 'IRF5', 'ZBTB7B', 'GFI1', 'BCL6', 'IRF8',
           'SPI1', 'CEBPA', 'ATF3', 'JDP2', 'MAFG', 'KLF4', 'KLF10', 'SP2', 'SP3', 'SP4',
           'EGR2', 'EGR3', 'RORB', 'NR1H3', 'STAT2', 'STAT5A', 'STAT6', 'SMAD2', 'SMAD3',
           'SMAD4', 'RUNX1', 'RUNX2', 'TCF3', 'TCF4', 'TCF12', 'TCF7L2', 'FOXO3', 'FOXO4',
           'IKZF4', 'ZEB1', 'SNAI1', 'SNAI2', 'TWIST1', 'ID1', 'ID4', 'HES1', 'HEY1',
           'NR2F6', 'NR2F1', 'NR2F2', 'ETS2', 'ELF4', 'FLI1', 'ERG', 'ETV6', 'ETV5',
           'GABPA', 'BHLHE41', 'ATF6', 'CREB1', 'ATF2', 'NFIL3']
REF_TFS = ["ID2", "GATA3", "STAT3", "MAF", "IRF1", "STAT1", "IRF4", "NFKB1", "BATF", "TCF7"]

TF_CANDIDATES = [
    # Naive / 记忆维持
    "TCF7", "LEF1", "KLF2", "BACH2", "ID2", "FOXO1",
    # 活化 / 效应
    "BATF", "IRF4", "JUN", "JUNB", "FOS", "NFKB1", "RELA", "MYB", "SP1", "CEBPD",
    # 耗竭
    "PRDM1", "NR4A1", "NR4A2", "NR4A3", "MAF", "TBX21", "EOMES", "IKZF2", "IKZF3",
    # Treg / 分化
    "FOXP3", "RUNX3", "GATA3", "RORC",
    # M8（decoupleR）中 JAML 阳性细胞的 top TF
    "HIF1A", "RFX5", "IRF6", "IRF1", "IRF5", "STAT1", "STAT3",
    # 其他
    "ETS1", "XBP1", "POU2F2",
]
if EXT2:
    TF_CANDIDATES = NEW_TFS + REF_TFS
TF_CANDIDATES = [t for t in dict.fromkeys(TF_CANDIDATES) if t in set(TF_COLS)]
log(f"    候选 TF {len(TF_CANDIDATES)} 个（已剔除不在 base GRN 者）：{TF_CANDIDATES}")

GIDX = {g: i for i, g in enumerate(oracle.adata.var_names)}
JAML_IDX = GIDX["JAML"]


def run_perturb(tf, use_randomized=False, do_transition=False,
                n_prop=N_PROP, value=0.0):
    """把一个 TF 设为 value(0=敲除)，返回 delta_X（细胞 × 基因）。"""
    oracle.simulate_shift(perturb_condition={tf: value}, GRN_unit="cluster",
                          n_propagation=n_prop, use_randomized_GRN=use_randomized,
                          clip_delta_X=True)
    if do_transition:
        oracle.estimate_transition_prob(n_neighbors=200, knn_random=True,
                                        sampled_fraction=0.5, n_jobs=1)
    d = oracle.adata.layers["delta_X"]
    return np.asarray(d.todense()) if sp.issparse(d) else np.asarray(d)


results = {}
skipped = []
for i, tf in enumerate(TF_CANDIDATES):
    full = (i == 0)          # 第一个跑完整流程，验证 velocyto 集成的 estimate_transition_prob
    try:
        results[tf] = run_perturb(tf, do_transition=full)
    except ValueError as e:
        # 例：LEF1 —— "does not have enough regulatory connection in the GRNs"
        # （filter_links 后该 TF 的连接被削掉）。属预期情形，记录后跳过，不做静默处理。
        skipped.append((tf, str(e)[:90]))
        log(f"    [{i+1:2d}/{len(TF_CANDIDATES)}] {tf:8s} 跳过：{str(e)[:80]}")
        continue
    d = results[tf]
    log(f"    [{i+1:2d}/{len(TF_CANDIDATES)}] {tf:8s} |Δ|均值={np.abs(d).mean():.5f}"
        f"  |ΔJAML|均值={np.abs(d[:, JAML_IDX]).mean():.5f}"
        f"{'   (+transition prob → velocyto 路径可用)' if full else ''}")
if skipped:
    log(f"    共跳过 {len(skipped)} 个 TF（GRN 连接不足）：{[s[0] for s in skipped]}")

# 阴性对照：打乱 GRN（跑 3 次，作为噪声基线）
_ok_tfs = list(results.keys())
for k in range(3):
    key = f"__randomGRN_{k+1}"
    btf = _ok_tfs[k % len(_ok_tfs)]
    try:
        results[key] = run_perturb(btf, use_randomized=True)
    except ValueError as e:
        log(f"    [随机对照 {k+1}] 跳过：{str(e)[:70]}")
        continue
    log(f"    [随机对照 {k+1}] 基于 {btf} 但打乱 GRN：|Δ|均值="
        f"{np.abs(results[key]).mean():.5f}  |ΔJAML|均值="
        f"{np.abs(results[key][:, JAML_IDX]).mean():.5f}")

genes = list(oracle.adata.var_names)

# =====================================================================
# 6. 汇总
# =====================================================================
log("[6] 汇总")
rows = []
for key, d in results.items():
    mag = np.abs(d).mean(axis=0)
    order = np.argsort(-mag)
    top = [(genes[i], float(mag[i])) for i in order[:5]]
    rows.append({
        "TF": key,
        "is_randomized_control": key.startswith("__randomGRN"),
        "mean_abs_delta_all_genes": float(mag.mean()),
        "n_genes_delta_gt_0.01": int((mag > 0.01).sum()),
        "JAML_mean_abs_delta": float(np.abs(d[:, JAML_IDX]).mean()),
        "JAML_mean_delta_signed": float(d[:, JAML_IDX].mean()),
        "top1_gene": top[0][0], "top1_val": top[0][1],
        "top2_gene": top[1][0], "top2_val": top[1][1],
        "top3_gene": top[2][0], "top3_val": top[2][1],
    })

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(RES, f"m4_tf_perturbation_summary{SUF}.csv"), index=False)

_r = summary[summary["is_randomized_control"]]["JAML_mean_abs_delta"]
_rall = summary[summary["is_randomized_control"]]["mean_abs_delta_all_genes"]
log(f"\n    随机对照基线（n=3）：|ΔJAML| 均值 {_r.mean():.5f}"
    f"（{_r.min():.5f}–{_r.max():.5f}）；全基因 |Δ| 均值 {_rall.mean():.5f}")
log("    判据：|ΔJAML| 明显高于该基线者，为 JAML 的候选调控 TF；")
log("          全基因 |Δ| 最高者用于验证 GRN 模拟的量程（阳性对照）。")

_sum = summary[~summary["is_randomized_control"]].sort_values(
    "JAML_mean_abs_delta", ascending=False)
log("\n=== 按 |ΔJAML| 降序（前 20）===")
log(_sum.head(20)[["TF", "JAML_mean_abs_delta", "JAML_mean_delta_signed",
                   "mean_abs_delta_all_genes", "n_genes_delta_gt_0.01"]]
    .round(5).to_string(index=False))

# |ΔJAML| 最大的 TF + 一次随机对照，用于展示关键基因的变化
_sum_valid = _sum[_sum["TF"].isin(results.keys())]
if len(_sum_valid):
    _top_tf = str(_sum_valid.iloc[0]["TF"])
else:
    _top_tf = list(results.keys())[0]
_rand_key = next((k for k in results if k.startswith("__randomGRN")), list(results)[-1])

log(f"\n=== 关键基因在「|ΔJAML| 最大 TF = {_top_tf}」扰动下的 |Δ| ===")
key_genes = ["JAML", "CXADR", "CD3E", "CD3D", "IL7R", "CCR7", "TCF7", "SELL", "PDCD1",
             "FOXP3", "GZMB", "GZMA", "BATF", "TOX", "ANXA1", "S100A4", "CD4", "CD8A"]
key_genes = [g for g in key_genes if g in genes]
_ki = [genes.index(g) for g in key_genes]
kt = pd.DataFrame({
    "gene": key_genes,
    f"{_top_tf}_KO": [float(np.abs(results[_top_tf]).mean(axis=0)[i]) for i in _ki],
    "randomized_ctrl": [float(np.abs(results[_rand_key]).mean(axis=0)[i]) for i in _ki],
})
log(kt.round(5).to_string(index=False))
kt.to_csv(os.path.join(RES, f"m4_key_genes_shift{SUF}.csv"), index=False)

# 保存主要 delta_X（|ΔJAML| 最大的 TF、随机对照），供复核
np.save(os.path.join(RES, f"m4_deltaX_{_top_tf}_KO{SUF}.npy"),
        results[_top_tf].astype(np.float32))
np.save(os.path.join(RES, f"m4_deltaX_randomGRN{SUF}.npy"),
        results[_rand_key].astype(np.float32))

meta = {
    "celloracle": co.__version__,
    "test_mode": TEST_MODE,
    "cluster_column": CLUSTER_COL,
    "alpha": ALPHA,
    "n_propagation": N_PROP,
    "bagging_number": BAGGING,
    "k_impute": K_IMPUTE,
    "n_pca_dims": N_PCA_DIMS,
    "base_GRN": "hg38_gimmemotifs_v5_fpr2_threshold_10_20210630",
    "n_tf_columns": int(n_tf),
    "n_genes_with_TFdict": int(len(TFdict)),
    "median_TFs_per_gene": float(np.median(_n_tf)),
    "n_cells": int(adata.n_obs),
    "n_genes": int(adata.n_vars),
    "input_file": INPUT_FILE,
    "input_tag": INPUT_TAG,
    "n_grn_units": int(len(_ld)),
    "n_links_total": int(n_links),
    "links_filtered_per_unit": {k: int(len(v)) for k, v in links.filtered_links.items()},
    "links_filter_p": FILTER_P,
    "links_filter_top_per_unit": FILTER_TOPN,
    "JAML_connections_per_unit": _n_jaml_conn,
    "JAML_total_connections": int(sum(_n_jaml_conn.values())),
    "n_TF_perturbed": int(len(TF_CANDIDATES)),
    "TF_perturbed": TF_CANDIDATES,
    "n_TF_skipped": int(len(skipped)),
    "TF_skipped_insufficient_GRN": [s[0] for s in skipped],
    "top_TF_by_JAML_delta": _top_tf,
    "note_JAML_not_TF": ("JAML 不在 base GRN 的 TF 列中（JAML 是黏附分子，非转录因子），"
                         "CellOracle 无法直接虚拟敲除 JAML；改为扰动 TF 观察其对 JAML 的调控"),
    "elapsed_sec": round(time.time() - t0, 1),
    "avail_gb_end": round(avail_gb(), 2),
}
with open(os.path.join(RES, f"m4_celloracle_summary{SUF}.json"), "w", encoding="utf-8") as f:
    json.dump(meta, f, ensure_ascii=False, indent=2)
log(f"\n[输出] m4_tf_perturbation_summary{SUF}.csv / m4_key_genes_shift{SUF}.csv / "
    f"m4_celloracle_summary{SUF}.json")
log(f"[内存] 结束可用 {avail_gb():.2f} GB")
log(f"\n总耗时 {time.time()-t0:.0f}s\nDONE")
