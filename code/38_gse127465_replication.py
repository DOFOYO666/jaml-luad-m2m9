# -*- coding: utf-8 -*-
"""
38_gse127465_replication.py — 独立数据集复核：JAML/CXADR 的区室分离

为什么做：规划 §7.3 明确"纯计算研究必须有独立数据集的重复验证"。GSE131907 既是
发现集也是验证集，不足以服人。本脚本用**另一项独立研究**（GSE127465 / Zilionis 2019，
NSCLC，54,773 细胞，含肿瘤与血液）复核：
  ① JAML 是否仍在免疫侧（尤其髓系/DC）；
  ② CXADR 是否仍在上皮侧；
  ③ CD4⁺T 内部的 JAML 阳性率梯度是否可重复。

数据格式（实测，2026-10-01）
---------------------------
GSE127465_human_counts_normalized_54773x41861.mtx.gz
  MatrixMarket coordinate，头 `54773 41861 44663765`
  数据行 `i j v`（1 基）。**文件按列(基因)聚集**：前 30 万行 j 只覆盖 1–766，
  而 i 已覆盖全部 1–54773 → **i = 细胞行号，j = 基因列号**。
GSE127465_gene_names_human_41861.tsv.gz  → 无表头，按行序给列号（LF 换行）
GSE127465_human_cell_metadata_54773x25.tsv.gz → 细胞行号 → 注释（含 Major/Minor）

⚠️ 上一版的 bug（已修）：`j2k = {v: k for k, v in enumerate(panel_genes)}` 里 v 是
   **基因名**而不是索引，导致用整数 j 查询永远 miss → 面板命中 0 行。
   正确写法：`{sym2idx[g]: k for k, g in enumerate(panel_genes)}`。

⚠️ 对齐校验（本版重做）：
   * 硬门槛 = **自算 MT 比例 vs 元数据 MT 比例**（逐细胞比例，实测 Spearman rho = 1.0000，
     n = 54,773）——只有行序完全一致才可能达到，是本数据集最可靠的对齐证据。
   * `Total counts` 列与自算总量仅 rho = 0.014（原因未明，疑为该列语义/来源不同），
     故**不作为判据**，仅记录在 JSON 中。
   * marker 方向检查为**辅助**证据。注意 LYZ 在中性粒/DC 也高表达，因此其比较范围
     限定为非髓系区室（见 MARKER_CHECK 的 compare_keys）。

输出：results/m8_gse127465_panel_by_celltype.csv
      results/m8_gse127465_alignment_check.json
      results/m8_gse127465_vs_gse131907.csv
      results/figures/m8_gse127465_jaml_cxadr.{png,tif}
      results/figures/m8_gse127465_cd4_jaml.{png,tif}
"""
import os
import gzip
import json
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "rawdata", "scRNA", "GSE127465")
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")
MTX = os.path.join(DATA, "GSE127465_human_counts_normalized_54773x41861.mtx.gz")
GENES = os.path.join(DATA, "GSE127465_gene_names_human_41861.tsv.gz")
META = os.path.join(DATA, "GSE127465_human_cell_metadata_54773x25.tsv.gz")

PANEL = ["JAML", "CXADR",
         "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "FOXP3", "IL7R", "CCR7", "PDCD1",
         "CD14", "FCGR3A", "LYZ", "ITGAX", "CLEC9A", "CLEC10A", "FCER1A", "CD68",
         "MS4A1", "CD79A", "NKG7", "KLRD1", "NCAM1",
         "EPCAM", "KRT18", "KRT19", "NKX2-1", "SFTPC", "SCGB1A1", "FOXJ1", "COL1A1",
         "PECAM1", "VWF", "TPSAB1"]

EPI_LIKE = ["Type I cells", "Type II cells", "Club cells", "Ciliated cells"]
NON_IMMUNE = {"Fibroblasts", "Endothelial cells", "Smooth muscle cells", "ND",
              "tRBC", "bRBC", "bPlatelets"}

# marker 方向检查（辅助）：(组名关键字, marker 基因, 比较范围关键字列表 or None)
#   None 表示与"其余所有组的最大值"比较；给出列表则只在所列组内比较
#   —— LYZ 在中性粒/DC 同样高表达，属髓系家族，不应算作"错误"
MARKER_CHECK = [
    ("T cells", "CD3E", None),
    ("Type II cells", "SFTPC", None),
    ("Endothelial cells", "PECAM1", None),
    ("Fibroblasts", "COL1A1", None),
    ("MoMacDC", "LYZ", ["T cells", "B cells", "NK cells", "Type II cells",
                        "Club cells", "Ciliated cells"]),
]


def load_gene_index():
    """返回 (panel_genes, j2k, mt_idx)。j2k: 基因列号(1 基) -> 面板内序号。"""
    names = [l.rstrip() for l in gzip.open(GENES, "rt")]
    sym2idx = {}
    for i, g in enumerate(names, start=1):
        if g not in sym2idx:
            sym2idx[g] = i
    panel_genes = [g for g in PANEL if g in sym2idx]
    # ⚠️ 修正：键必须是**基因列号**，不能是基因名
    j2k = {sym2idx[g]: k for k, g in enumerate(panel_genes)}
    mt_idx = {i for i, g in enumerate(names, start=1) if g.startswith("MT-")}
    miss = [g for g in PANEL if g not in sym2idx]
    print(f"[基因] 文件 {len(names)} 行；面板 {len(PANEL)} 个，命中 {len(panel_genes)}；"
          f"未命中: {miss or '无'}")
    print(f"[基因] MT- 基因列 {len(mt_idx)} 个")
    return panel_genes, j2k, mt_idx


def save_fig(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".tif"), dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"[图] results/figures/{name}.png / .tif")


def main():
    os.makedirs(FIG, exist_ok=True)
    print("=" * 70)
    print("独立数据集复核：GSE127465（Zilionis 2019, NSCLC）")
    print("=" * 70)

    meta = pd.read_csv(META, sep="\t")
    n_cells_meta = len(meta)
    print(f"[元数据] {meta.shape}；Major cell type {meta['Major cell type'].nunique()} 类")
    print("  Tissue:", meta["Tissue"].value_counts().to_dict())

    panel_genes, j2k, mt_idx = load_gene_index()
    kgene = len(panel_genes)

    major = meta["Major cell type"].astype(str).values
    minor = (meta["Minor subset"].astype(str).values if "Minor subset" in meta.columns
             else np.array(["NA"] * n_cells_meta))
    major_list = sorted(set(major))
    minor_list = sorted(set(minor))
    mpos = {c: i for i, c in enumerate(major_list)}
    spos = {c: i for i, c in enumerate(minor_list)}
    n_major, n_minor = len(major_list), len(minor_list)
    mcode = np.array([mpos[c] for c in major], dtype=np.int32)
    scode = np.array([spos[c] for c in minor], dtype=np.int32)
    cnt_maj = np.bincount(mcode, minlength=n_major).astype(float)
    cnt_min = np.bincount(scode, minlength=n_minor).astype(float)

    # 面板基因 × 组：累加 [表达和, 阳性数]
    acc_maj = np.zeros((kgene, n_major, 2))
    acc_min = np.zeros((kgene, n_minor, 2))
    # 逐细胞诊断量
    cell_sum = np.zeros(n_cells_meta)
    cell_mt = np.zeros(n_cells_meta)
    cell_det = np.zeros(n_cells_meta, dtype=np.int32)

    print("\n[流式解析 mtx] i=细胞，j=基因 ...")
    n_line = 0
    n_keep = 0
    nrow = ncol = nnz = 0
    header_read = False
    row_max = col_max = 0
    with gzip.open(MTX, "rt", errors="replace") as f:
        for line in f:
            if not header_read:
                if line.startswith("%"):
                    continue
                d = line.split()
                nrow, ncol, nnz = int(d[0]), int(d[1]), int(d[2])
                print(f"  维度: {nrow} 行(细胞) × {ncol} 列(基因)，非零 {nnz:,}")
                header_read = True
                continue
            n_line += 1
            a, b, c = line.split()
            i = int(a); j = int(b); v = float(c)
            if i > row_max:
                row_max = i
            if j > col_max:
                col_max = j
            ci = i - 1
            cell_sum[ci] += v
            cell_det[ci] += 1
            if j in mt_idx:
                cell_mt[ci] += v
            k = j2k.get(j)
            if k is not None:
                n_keep += 1
                cm = mcode[ci]
                acc_maj[k, cm, 0] += v
                acc_maj[k, cm, 1] += 1.0 if v > 0 else 0.0
                sm = scode[ci]
                acc_min[k, sm, 0] += v
                acc_min[k, sm, 1] += 1.0 if v > 0 else 0.0

    print(f"[完成] 读到 EOF；数据行 {n_line:,}；命中面板的行 {n_keep:,}")
    print(f"  索引范围: i ≤ {row_max}（应 = {nrow}），j ≤ {col_max}（应 = {ncol}）")
    if n_keep == 0:
        raise RuntimeError("面板命中仍为 0 —— 索引映射仍有问题，停止")

    # ---------------- 对齐校验 ----------------
    from scipy import stats as st

    tc = pd.to_numeric(meta["Total counts"], errors="coerce").values
    ok_tc = np.isfinite(tc) & (tc > 0) & (cell_sum > 0)
    rho_tc, p_tc = st.spearmanr(cell_sum[ok_tc], tc[ok_tc])

    pmt = pd.to_numeric(meta["Percent counts from mitochondrial genes"],
                        errors="coerce").values
    my_mt = np.divide(cell_mt, cell_sum, out=np.zeros_like(cell_mt), where=cell_sum > 0)
    ok_mt = np.isfinite(pmt) & (cell_sum > 0)
    rho_mt, p_mt = st.spearmanr(my_mt[ok_mt], pmt[ok_mt])

    print(f"\n[对齐诊断 A] 自算总量 vs 元数据 Total counts: Spearman rho = {rho_tc:.4f} "
          f"(P = {p_tc:.3g}, n = {ok_tc.sum()})  ← 该列语义存疑，不作判据")
    print(f"[对齐诊断 B] 自算 MT 比例 vs 元数据 MT 比例: Spearman rho = {rho_mt:.4f} "
          f"(P = {p_mt:.3g}, n = {ok_mt.sum()})  ← 硬门槛")
    print(f"  Total counts 列统计: 唯一值 {pd.Series(tc).nunique()}，"
          f"中位 {np.nanmedian(tc):.1f}，范围 {np.nanmin(tc):.1f}–{np.nanmax(tc):.1f}")

    # 硬门槛：MT 比例逐细胞吻合（行序一致的必要充分证据）
    if not (rho_mt > 0.99):
        raise RuntimeError(f"MT 比例对齐校验未通过（rho={rho_mt:.4f}）——"
                           "元数据行序与矩阵行序不一致，停止")

    print("\n[对齐辅助] marker 方向检查")
    maj_tab = pd.DataFrame({
        "gene": np.repeat(panel_genes, n_major),
        "group": np.tile(major_list, kgene),
        "n_cells": np.tile(cnt_maj, kgene),
        "mean_expr": acc_maj[:, :, 0].ravel() / np.tile(cnt_maj, kgene),
        "pct_positive": 100 * acc_maj[:, :, 1].ravel() / np.tile(cnt_maj, kgene),
    })
    maj_tab = maj_tab[maj_tab["n_cells"] > 0]

    marker_rows = []
    for key, gene, ck in MARKER_CHECK:
        sel = maj_tab["group"].str.contains(key, case=False, na=False)
        hit = maj_tab[(maj_tab["gene"] == gene) & sel]
        other = maj_tab[(maj_tab["gene"] == gene) & ~sel]
        if ck is not None:
            mask = pd.Series(False, index=other.index)
            for c in ck:
                mask |= other["group"].str.contains(c, case=False, na=False)
            other = other[mask]
        if len(hit) == 0 or len(other) == 0:
            print(f"  {gene}: 未匹配到组 —— 跳过")
            continue
        m_in = float(hit["mean_expr"].iloc[0])
        m_out = float(other["mean_expr"].max())
        ratio = m_in / m_out if m_out > 0 else np.inf
        good = m_in > m_out
        marker_rows.append({"gene": gene, "group_key": key,
                            "group": str(hit["group"].iloc[0]),
                            "mean_in": m_in, "max_other": m_out,
                            "ratio": ratio, "pass": bool(good),
                            "compare_scope": "all" if ck is None else ck})
        print(f"  {gene:8s} {str(hit['group'].iloc[0]):18s} 均值 {m_in:9.4f}；"
              f"比较组最大 {m_out:9.4f}；比值 {ratio:7.2f}  → {'通过' if good else '未通过'}")

    with open(os.path.join(RES, "m8_gse127465_alignment_check.json"), "w",
              encoding="utf-8") as f:
        json.dump({"nrow": nrow, "ncol": ncol, "nnz": nnz, "lines_read": int(n_line),
                   "panel_lines": int(n_keep),
                   "spearman_totalcounts": {"rho": float(rho_tc), "p": float(p_tc),
                                            "n": int(ok_tc.sum()),
                                            "used_as_gate": False},
                   "spearman_mt_fraction": {"rho": float(rho_mt), "p": float(p_mt),
                                            "n": int(ok_mt.sum()),
                                            "used_as_gate": True},
                   "marker_check": marker_rows},
                  f, ensure_ascii=False, indent=2)

    print("\n[对齐] MT 比例 rho = %.4f → 元数据行序与矩阵行序一致，可用于分组分析" % rho_mt)

    # ---------------- 汇总表 ----------------
    min_tab = pd.DataFrame({
        "gene": np.repeat(panel_genes, n_minor),
        "group": np.tile(minor_list, kgene),
        "n_cells": np.tile(cnt_min, kgene),
        "mean_expr": acc_min[:, :, 0].ravel() / np.tile(cnt_min, kgene),
        "pct_positive": 100 * acc_min[:, :, 1].ravel() / np.tile(cnt_min, kgene),
    })
    min_tab = min_tab[min_tab["n_cells"] > 0]
    maj_out = maj_tab.assign(level="Major").rename(columns={"group": "group"})
    min_out = min_tab.assign(level="Minor")
    out = pd.concat([maj_out, min_out], ignore_index=True)[
        ["gene", "level", "group", "n_cells", "mean_expr", "pct_positive"]]
    out.to_csv(os.path.join(RES, "m8_gse127465_panel_by_celltype.csv"), index=False)
    print(f"[输出] results/m8_gse127465_panel_by_celltype.csv  ({out.shape})")

    maj = out[out["level"] == "Major"]
    print("\n=== JAML 与 CXADR 在各 Major cell type ===")
    for gene in ["JAML", "CXADR"]:
        d = maj[maj["gene"] == gene].sort_values("pct_positive", ascending=False)
        print(f"\n--- {gene} ---")
        print(d[["group", "n_cells", "mean_expr", "pct_positive"]].head(16)
              .round(4).to_string(index=False))

    def grp_mean(gene, sel):
        d = maj[(maj["gene"] == gene) & (maj["group"].isin(sel))]
        w = d["n_cells"].values
        if w.sum() == 0:
            return np.nan, np.nan, 0
        return (float(np.average(d["mean_expr"], weights=w)),
                float(np.average(d["pct_positive"], weights=w)), int(w.sum()))

    immune = [c for c in major_list if not c.startswith("Patient")
              and c not in EPI_LIKE and c not in NON_IMMUNE]
    print(f"\n[分组] 免疫类 {len(immune)} 个；上皮类 {len(EPI_LIKE)} 个")
    cmp_rows = []
    for gene in ["JAML", "CXADR"]:
        for lab, sel in [("immune", immune), ("epithelial-like", EPI_LIKE)]:
            me, pe, ne = grp_mean(gene, sel)
            cmp_rows.append({"gene": gene, "compartment": lab, "n_cells": ne,
                             "mean_expr": me, "pct_positive": pe})
    cmp_df = pd.DataFrame(cmp_rows)
    cmp_df.to_csv(os.path.join(RES, "m8_gse127465_vs_gse131907.csv"), index=False)
    print("\n=== 免疫 vs 上皮样（GSE127465）===")
    print(cmp_df.round(4).to_string(index=False))

    print("\n=== 髓系/DC 与 T 细胞（JAML，Major 层）===")
    for key in ["MoMacDC", "T cells", "DC"]:
        d = maj[(maj["gene"] == "JAML") & maj["group"].str.contains(key, case=False, na=False)]
        for _, r in d.iterrows():
            print(f"  {r['group']:20s} n={int(r['n_cells']):6d}  "
                  f"阳性率 {r['pct_positive']:.1f}%  均值 {r['mean_expr']:.4f}")

    print("\n=== CD4 相关 Minor subset 的 JAML（n_cells ≥ 20）===")
    cd4 = out[(out["level"] == "Minor") & (out["gene"] == "JAML") &
              (out["group"].str.contains("CD4|Tfh", case=False, na=False))]
    cd4 = cd4[cd4["n_cells"] >= 20].sort_values("pct_positive", ascending=False)
    print(cd4[["group", "n_cells", "mean_expr", "pct_positive"]].round(4).to_string(index=False)
          if len(cd4) else "  未找到 CD4 亚群标签")

    # ---------------- 图 ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sel = ["tMoMacDC", "tpDC", "tNeutrophils", "bNeutrophils", "tT cells", "bT cells",
           "tNK cells", "bNK cells", "tB cells", "tPlasma cells", "tMast cells",
           "Fibroblasts", "Endothelial cells", "Type II cells", "Club cells",
           "Ciliated cells", "Type I cells"]
    sel = [c for c in sel if c in set(maj["group"])]
    piv = (maj[maj["group"].isin(sel)]
           .pivot_table(index="group", columns="gene", values="pct_positive")
           .reindex(sel))

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5),
                             gridspec_kw={"width_ratios": [1.1, 1]})
    ax = axes[0]
    y = np.arange(len(piv))
    h = 0.38
    ax.barh(y + h / 2, piv["CXADR"].values, height=h, color="#4575B4",
            edgecolor="#333", linewidth=0.5, label="CXADR (receptor)")
    ax.barh(y - h / 2, piv["JAML"].values, height=h, color="#D73027",
            edgecolor="#333", linewidth=0.5, label="JAML (ligand)")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c}  (n={int(cnt_maj[mpos[c]]):,})" for c in piv.index], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Positive cells (%)")
    ax.set_title("Independent replication: JAML\u2013CXADR compartmental separation\n"
                 "GSE127465 (Zilionis 2019, NSCLC, n = 54,773 cells)", fontsize=8.8)
    ax.legend(frameon=False, fontsize=7.5, loc="center right")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    spec = [("JAML", [("tMoMacDC", "Myeloid/DC (t)"), ("tT cells", "T cells (t)"),
                      ("tNeutrophils", "Neutrophils (t)"), ("tB cells", "B cells (t)"),
                      ("tNK cells", "NK cells (t)"),
                      ("Type II cells", "Type II epi."), ("Club cells", "Club")]),
            ("CXADR", [("Type II cells", "Type II epi."), ("Club cells", "Club"),
                       ("Ciliated cells", "Ciliated"), ("Type I cells", "Type I epi."),
                       ("tMoMacDC", "Myeloid/DC (t)"), ("tT cells", "T cells (t)")])]
    labels, vals, cols = [], [], []
    for gene, items in spec:
        for c, lab in items:
            d = maj[(maj["gene"] == gene) & (maj["group"] == c)]
            if len(d) == 0:
                continue
            labels.append(lab)
            vals.append(float(d["pct_positive"].iloc[0]))
            cols.append("#D73027" if gene == "JAML" else "#4575B4")
    yy = np.arange(len(labels))
    ax.barh(yy, vals, color=cols, edgecolor="#333", linewidth=0.5)
    for i, vv in enumerate(vals):
        ax.text(vv + 0.5, i, f"{vv:.1f}", va="center", fontsize=6.4)
    ax.set_yticks(yy)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Positive cells (%)")
    ax.set_title("Key compartments (red = JAML, blue = CXADR)", fontsize=8.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.98, 0.02, "few epithelial cells in this dataset",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.4, color="#777")
    fig.tight_layout()
    save_fig(fig, "m8_gse127465_jaml_cxadr")

    if len(cd4):
        fig, ax = plt.subplots(figsize=(4.8, 3.4))
        cols2 = plt.get_cmap("YlOrRd")(np.linspace(0.3, 0.95, len(cd4)))
        ax.bar(range(len(cd4)), cd4["pct_positive"].values, color=cols2,
               edgecolor="#333", linewidth=0.5)
        for i, (pv, nv) in enumerate(zip(cd4["pct_positive"].values, cd4["n_cells"].values)):
            ax.text(i, pv + 0.5, f"{pv:.1f}%\nn={int(nv)}", ha="center", va="bottom",
                    fontsize=6.4)
        ax.set_xticks(range(len(cd4)))
        ax.set_xticklabels(cd4["group"].values, rotation=25, ha="right", fontsize=7.5)
        ax.set_ylabel("JAML-positive cells (%)")
        ax.set_ylim(0, max(cd4["pct_positive"].max() * 1.4, 12))
        ax.set_title("JAML in CD4\u207a T subsets\nGSE127465 (independent dataset)", fontsize=8.8)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        save_fig(fig, "m8_gse127465_cd4_jaml")

    print("\nDONE")


if __name__ == "__main__":
    main()
