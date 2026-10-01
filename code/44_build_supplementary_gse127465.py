# -*- coding: utf-8 -*-
"""把 GSE127465 独立复核的证据整理成两个投稿包的补充材料。

产出（两个包统一用 "Additional file N" 命名——JTM 包的实际文件命名即如此，
其 README 里的 "SupplementaryTableS1–S8" 与实际不符，见本脚本末尾的提示）：
  Additional file 11 —— STROBE-MR 清单（JTM 包原本缺此文件，本次补齐以对齐编号）
  Additional file 12 —— 对齐校验证据（逐细胞 MT 比例 ρ = 1.0000 + 5 项 marker + 结构核验 + 独立重实现对拍）
  Additional file 13 —— 两数据集区室阳性率对照（核心结论表）
  Additional file 14 —— GSE127465 逐细胞类型表达面板（35 基因，Major + Minor）
  Additional file 15 —— 跨数据集区室对照图（PNG + 300 dpi LZW TIFF）
  JTM 另在 03_Figures/ 放 SupplementaryFigureS2_JAML_CXADR_compartments.{png,tiff}

用法：python 44_build_supplementary_gse127465.py
"""
import io
import json
import os
import shutil
import zipfile

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")
BM = os.path.join(ROOT, "BMCCancer投稿", "05_Supplementary")
JTM_SUP = os.path.join(ROOT, "JTM投稿", "05_Supplementary")
JTM_FIG = os.path.join(ROOT, "JTM投稿", "03_Figures")

CAMPAIGN = "GSE127465 (Zilionis 2019, NSCLC; 54,773 cells)"
NOTE_ALIGN = (f"Independent replication dataset: {CAMPAIGN}. Rows = cells, columns = genes in the "
              "source MatrixMarket file (verified empirically: in the first 300,000 data lines the "
              "column index spans only 1-766 while the row index already spans the full 1-54,773).")


def line_block(title, body):
    """生成带题注的多行文本，便于直接并入 CSV 的注释行。"""
    return f"# {title}\n# {body}\n"


def build_alignment():
    with io.open(os.path.join(RES, "m8_gse127465_alignment_check.json"), encoding="utf-8") as f:
        a = json.load(f)
    with io.open(os.path.join(RES, "verify_gse127465_targets_summary.json"), encoding="utf-8") as f:
        v = json.load(f)

    rows = []
    rows.append(("Dataset", CAMPAIGN, "source GEO accession GSE127465"))
    rows.append(("Matrix file", "GSE127465_human_counts_normalized_54773x41861.mtx.gz",
                 "528,303,938 bytes; MatrixMarket coordinate; header '54773 41861 44663765'"))
    rows.append(("Matrix orientation", "rows = cells (54,773), columns = genes (41,861)",
                 NOTE_ALIGN))
    rows.append(("Non-zero entries read", f"{a['lines_read']:,}",
                 f"equals the nnz declared in the MatrixMarket header ({a['nnz']:,})"))
    rows.append(("Panel lines matched", f"{a['panel_lines']:,}", "data lines belonging to the 35-gene panel"))
    rows.append(("Index range check", f"row <= {a['nrow']}, column <= {a['ncol']}", "both equal the declared dimensions"))
    tc = a["spearman_totalcounts"]
    rows.append(("Alignment test A (rejected)", f"Spearman rho = {tc['rho']:.4f} (n = {tc['n']:,})",
                 "self-computed per-cell total vs metadata 'Total counts' column; NOT used as a criterion "
                 "because the column's definition could not be established (see limitations)"))
    mt = a["spearman_mt_fraction"]
    rows.append(("Alignment test B (primary criterion)", f"Spearman rho = {mt['rho']:.4f} (n = {mt['n']:,})",
                 "self-computed per-cell mitochondrial fraction vs metadata 'Percent counts from "
                 "mitochondrial genes' (37 MT- genes); a per-cell ratio can only reach rho = 1.0000 if the "
                 "row order of the metadata and of the matrix are identical"))
    for m in a["marker_check"]:
        rows.append((f"Marker direction check: {m['gene']}",
                     f"{m['group']}: {m['mean_in']:.4f} vs max of comparison groups {m['max_other']:.4f} (ratio {m['ratio']:.2f})",
                     f"comparison scope = {m['compare_scope'] if m['compare_scope'] != 'all' else 'all other groups'}"))
    rows.append(("Independent re-implementation",
                 f"{v['n_compared']} cells compared, {v['n_mismatch']} mismatches",
                 "second implementation (column-index ordered scan with early stop, reporting absolute "
                 "positive-cell counts) reproduced the first implementation cell-for-cell"))
    rows.append(("Apparent coincidence resolved", "JAML and CXADR share identical positive-cell percentages in Endothelial cells and Fibroblasts",
                 "absolute counts show 4/535 and 3/585 for both genes; a genuine coincidence, not an index error"))

    df = pd.DataFrame(rows, columns=["Item", "Value", "Notes / source"])
    return df


def build_compartments():
    d = pd.read_csv(os.path.join(RES, "cross_dataset_compartment_pct.csv"))
    piv = d.pivot_table(index=["gene", "compartment"], columns="dataset",
                        values=["pct_positive", "n_cells"], observed=False).reset_index()
    piv.columns = ["gene", "compartment"] + [
        f"{a}_{b}" for a, b in piv.columns[2:]]
    order = ["Myeloid/DC", "Mast", "B lineage", "T cells", "NK cells",
             "Fibroblasts", "Endothelial", "Epithelial"]
    piv["compartment"] = pd.Categorical(piv["compartment"], categories=order, ordered=True)
    piv = piv.sort_values(["gene", "compartment"])
    # 增列两数据集比值
    piv["epithelial_to_myeloid_GSE131907"] = np.nan
    piv["epithelial_to_myeloid_GSE127465"] = np.nan
    for g in ["JAML", "CXADR"]:
        m = piv["gene"] == g
        for ds in ["GSE131907", "GSE127465"]:
            ep = piv.loc[m & (piv["compartment"] == "Epithelial"), f"pct_positive_{ds}"]
            my = piv.loc[m & (piv["compartment"] == "Myeloid/DC"), f"pct_positive_{ds}"]
            if len(ep) and len(my) and float(my.iloc[0]) > 0:
                piv.loc[m, f"epithelial_to_myeloid_{ds}"] = float(ep.iloc[0]) / float(my.iloc[0])
    cols = ["gene", "compartment", "pct_positive_GSE131907", "pct_positive_GSE127465",
            "n_cells_GSE131907", "n_cells_GSE127465",
            "epithelial_to_myeloid_GSE131907", "epithelial_to_myeloid_GSE127465"]
    return piv[cols].round(6)


def build_panel():
    d = pd.read_csv(os.path.join(RES, "m8_gse127465_panel_by_celltype.csv"))
    d = d[["gene", "level", "group", "n_cells", "mean_expr", "pct_positive"]]
    d = d.sort_values(["level", "group", "gene"])
    return d.round(6)


def write_af(dirpath, df, n, title, note=""):
    p = os.path.join(dirpath, f"Additional file {n} - {title}.csv")
    with io.open(p, "w", encoding="utf-8-sig", newline="") as f:
        f.write(f"# Additional file {n}. {title}.\n")
        if note:
            f.write(f"# {note}\n")
        df.to_csv(f, index=False)
    print(f"[写入] {os.path.relpath(p, ROOT)}  ({df.shape[0]} 行)")
    return p


def retire(dirpath, names):
    """把命名不符的文件移到 _obsolete/（沙箱禁用删除，改用移动）。"""
    ob = os.path.join(dirpath, "_obsolete")
    moved = []
    for n in names:
        p = os.path.join(dirpath, n)
        if os.path.isfile(p):
            os.makedirs(ob, exist_ok=True)
            shutil.move(p, os.path.join(ob, n))
            moved.append(n)
    if moved:
        print(f"[退役] {os.path.relpath(dirpath, ROOT)} 下移入 _obsolete/：{moved}")


def main():
    for d in (BM, JTM_SUP, JTM_FIG):
        os.makedirs(d, exist_ok=True)

    align = build_alignment()
    comp = build_compartments()
    panel = build_panel()
    print(f"[数据] 对齐证据 {align.shape}；区室对照 {comp.shape}；面板 {panel.shape}")

    T12 = "GSE127465 alignment verification"
    T13 = "JAML-CXADR compartment comparison across two single-cell datasets"
    T14 = "GSE127465 cell-type expression panel (35 genes)"
    for dirpath, label in [(BM, "BMCCancer"), (JTM_SUP, "JTM")]:
        write_af(dirpath, align, 12, T12)
        write_af(dirpath, comp, 13, T13)
        write_af(dirpath, panel, 14, T14)

    # JTM 包原本缺 STROBE-MR 清单（AF11）→ 从 BMCCancer 包补齐，使编号与母稿一致
    strobe = os.path.join(BM, "Additional file 11 - STROBE-MR checklist.csv")
    if os.path.isfile(strobe) and not os.path.isfile(
            os.path.join(JTM_SUP, "Additional file 11 - STROBE-MR checklist.csv")):
        shutil.copy2(strobe, os.path.join(JTM_SUP, "Additional file 11 - STROBE-MR checklist.csv"))
        print("[补齐] JTM 包新增 'Additional file 11 - STROBE-MR checklist.csv'（原缺失）")

    # 清理上一版误用的命名
    retire(JTM_SUP, ["SupplementaryTableS9.csv", "SupplementaryTableS10.csv",
                     "SupplementaryTableS11.csv", "SupplementaryFigureS2.png"])

    # ---------------- 图 ----------------
    src_png = os.path.join(FIG, "cross_dataset_compartment.png")
    src_tif = os.path.join(FIG, "cross_dataset_compartment.tif")
    copies = [
        (src_png, os.path.join(BM, "Additional file 15 - JAML-CXADR cross-dataset.png")),
        (src_tif, os.path.join(BM, "Additional file 15 - JAML-CXADR cross-dataset.tiff")),
        (src_png, os.path.join(JTM_SUP, "Additional file 15 - JAML-CXADR cross-dataset.png")),
        (src_tif, os.path.join(JTM_SUP, "Additional file 15 - JAML-CXADR cross-dataset.tiff")),
        (src_png, os.path.join(JTM_FIG, "SupplementaryFigureS2_JAML_CXADR_compartments.png")),
        (src_tif, os.path.join(JTM_FIG, "SupplementaryFigureS2_JAML_CXADR_compartments.tiff")),
    ]
    for s, d in copies:
        shutil.copy2(s, d)
        print(f"[复制] {os.path.relpath(d, ROOT)}  ({os.path.getsize(d)/1024:.1f} KB)")

    # ---------------- 重建 BMCCancer 的补充材料 zip ----------------
    zp = os.path.join(ROOT, "BMCCancer投稿", "AdditionalFiles.zip")
    names = sorted(os.listdir(BM))
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
        for n in names:
            z.write(os.path.join(BM, n), n)
    print(f"[打包] {os.path.relpath(zp, ROOT)}  ({len(names)} 个文件, {os.path.getsize(zp)/1048576:.2f} MB)")

    # ---------------- 备份到项目 results/ ----------------
    for src, dst in [(align, "supp_gse127465_alignment.csv"),
                     (comp, "supp_gse127465_compartments.csv"),
                     (panel, "supp_gse127465_panel.csv")]:
        src.to_csv(os.path.join(RES, dst), index=False)
    print("[备份] results/supp_gse127465_{alignment,compartments,panel}.csv")
    print("DONE")


if __name__ == "__main__":
    main()
