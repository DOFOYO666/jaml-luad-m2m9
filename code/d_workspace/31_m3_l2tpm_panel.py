# -*- coding: utf-8 -*-
"""
31_m3_l2tpm_panel.py — 用官方归一化 log2TPM 矩阵复核 M3 的基因面板（交叉验证）

目的：
  1. 验证 GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz 的完整性
     （流式读到 EOF；截断会抛 EOFError）
  2. 抽取同一基因面板在归一化空间的表达，与 raw UMI 版（scripts/22）对比
     —— 若两版结论方向一致，则 M3 的区室分离结论对归一化方式不敏感

输入：rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz
      rawdata/scRNA/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz
输出：results/m3_gse131907_l2tpm_panel_expression.csv
      results/m3_panel_raw_vs_l2tpm_consistency.csv

⚠️ 低内存策略：逐行流式，绝不整矩阵入内存。
⚠️ 符号：GSE131907 用 JAML（非 AMICA1）。
"""
import os
import gzip
import json
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "rawdata", "scRNA", "GSE131907")
RES = os.path.join(ROOT, "results")
ANN = os.path.join(DATA, "GSE131907_Lung_Cancer_cell_annotation.txt.gz")
L2 = os.path.join(DATA, "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz")
RAW_PANEL = os.path.join(RES, "m3_gse131907_panel_expression.csv")

PANEL = [
    # 核心配体-受体
    "JAML", "CXADR",
    # 免疫谱系 marker
    "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "TRDC", "FOXP3", "IL7R", "CCR7", "PDCD1",
    "CD14", "FCGR3A", "LYZ", "ITGAX", "CLEC9A", "CLEC10A", "FCER1A", "CD68", "CD163",
    "MS4A1", "CD79A", "NKG7", "KLRD1", "NCAM1",
    # 上皮 / 基质
    "EPCAM", "KRT18", "KRT19", "NKX2-1", "SFTPC", "SCGB1A1", "FOXJ1", "COL1A1",
    "PECAM1", "VWF", "TPSAB1",
]


def main():
    print("=" * 70)
    print("M3 交叉验证：官方 log2TPM 归一化矩阵复核")
    print("=" * 70)

    ann = pd.read_csv(ANN, sep="\t")
    ann.columns = [c.strip().lstrip("\ufeff") for c in ann.columns]
    # ⚠️ 对齐键必须是 Index 列（形如 AAACCTGCAAGGTGTG_LUNG_N01），
    #    与矩阵表头（AAACCTGAGAAACCGC_LN_05）同格式；Barcode 列缺样本后缀，对不上。
    key = "Index" if "Index" in ann.columns else ann.columns[0]
    bc = ann[key].astype(str).values
    cts = ann["Cell_type"].astype(str).values
    subs = ann["Cell_subtype"].astype(str).values
    print(f"[注释] {len(ann)} 细胞；Cell_type {len(set(cts))} 类；Cell_subtype {len(set(subs))} 类")

    with gzip.open(L2, "rt", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        cells = header[1:]
        print(f"[矩阵] {len(cells)} 个细胞列")

        # barcode 对齐（矩阵表头形如 AAACCTGAGAAACCGC_LN_05，注释 Barcode 同格式）
        pos = {b: i for i, b in enumerate(cells)}
        cell_idx = np.array([pos.get(b, -1) for b in bc])
        matched = int((cell_idx >= 0).sum())
        print(f"[对齐] 注释细胞在矩阵中命中 {matched} / {len(bc)}")
        if matched == 0:
            raise RuntimeError("barcode 对齐为 0 —— 检查对齐键（应为注释的 Index 列）")

        # 预先算好分组掩码
        ct_list = sorted(set(cts))
        ct_mask = {c: (cts == c) & (cell_idx >= 0) for c in ct_list}
        sub_list = sorted(set(subs))
        sub_mask = {s: (subs == s) & (cell_idx >= 0) for s in sub_list}

        acc_ct, acc_sub = {}, {}
        found = []
        nline = 0
        for line in f:
            nline += 1
            if nline % 5000 == 0:
                print(f"  ...已读 {nline} 行", flush=True)
            p = line.find("\t")
            name = line[:p]
            if name not in PANEL:
                continue
            vals = np.fromstring(line[p + 1:], dtype=np.float32, sep="\t")
            if vals.size != len(cells):
                print(f"  [警告] {name} 长度 {vals.size} != {len(cells)}，跳过")
                continue
            found.append(name)
            for c in ct_list:
                m = ct_mask[c]
                if m.sum() == 0:
                    continue
                v = vals[cell_idx[m]]
                acc_ct.setdefault(name, {})[c] = [float(v.mean()), int((v > 0).sum()), int(m.sum())]
            for s in sub_list:
                m = sub_mask[s]
                if m.sum() == 0:
                    continue
                v = vals[cell_idx[m]]
                acc_sub.setdefault(name, {})[s] = [float(v.mean()), int((v > 0).sum()), int(m.sum())]

    print(f"\n[完成] 读到 EOF（矩阵完整）；共 {nline} 行基因")
    print(f"[面板] 命中 {len(found)} / {len(PANEL)}：{sorted(found)}")
    missing = sorted(set(PANEL) - set(found))
    if missing:
        print(f"[面板] 未命中：{missing}")

    rows = []
    for g in sorted(acc_ct):
        for c, (mn, npz, n) in acc_ct[g].items():
            rows.append({"gene": g, "cell_type": c, "mean_log2tpm": mn,
                         "pos": npz, "n": n, "pos_frac": npz / n if n else np.nan})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "m3_gse131907_l2tpm_panel_expression.csv"), index=False)
    print(f"[输出] results/m3_gse131907_l2tpm_panel_expression.csv  ({df.shape})")

    rows2 = []
    for g in sorted(acc_sub):
        for c, (mn, npz, n) in acc_sub[g].items():
            rows2.append({"gene": g, "cell_subtype": c, "mean_log2tpm": mn,
                          "pos": npz, "n": n, "pos_frac": npz / n if n else np.nan})
    pd.DataFrame(rows2).to_csv(
        os.path.join(RES, "m3_gse131907_l2tpm_panel_subtype.csv"), index=False)
    print(f"[输出] results/m3_gse131907_l2tpm_panel_subtype.csv  ({len(rows2)} 行)")

    # ---------------- 与 raw UMI 版一致性 ----------------
    if os.path.exists(RAW_PANEL):
        raw = pd.read_csv(RAW_PANEL)
        print("\n[一致性] raw 面板列:", list(raw.columns))
        raw = raw[raw["level"] == "Cell_type"].copy()
        raw["pct_raw"] = pd.to_numeric(raw["pct_positive"], errors="coerce")
        l2 = df.copy()
        l2["pct_l2"] = l2["pos_frac"] * 100.0
        mg = raw[["gene", "group", "n_cells", "mean_counts", "pct_raw"]].merge(
            l2[["gene", "cell_type", "mean_log2tpm", "pct_l2"]],
            left_on=["gene", "group"], right_on=["gene", "cell_type"], how="inner")
        mg = mg.rename(columns={"group": "cell_type", "n_cells": "n_raw_cells"})
        if len(mg):
            a = mg["pct_raw"].values
            b = mg["pct_l2"].values
            ok = ~np.isnan(a) & ~np.isnan(b)
            r = np.corrcoef(a[ok], b[ok])[0, 1] if ok.sum() > 2 else np.nan
            print("\n=== raw UMI vs log2TPM 一致性（基因 × 细胞类型的阳性率 %）===")
            print(f"  配对点数 {ok.sum()}；Pearson r = {r:.4f}")
            print(f"  raw 均值 {a[ok].mean():.4f}% vs l2TPM 均值 {b[ok].mean():.4f}%")
            print(f"  最大绝对差 {np.abs(a[ok] - b[ok]).max():.4f} 个百分点")
            mg.to_csv(os.path.join(RES, "m3_panel_raw_vs_l2tpm_consistency.csv"), index=False)
            print(f"[输出] results/m3_panel_raw_vs_l2tpm_consistency.csv ({mg.shape})")
            # JAML / CXADR 两版并列
            print("\n=== JAML / CXADR 两版对照（Cell_type 层）===")
            sub = mg[mg["gene"].isin(["JAML", "CXADR"])]
            print(sub.sort_values(["gene", "pct_l2"], ascending=[True, False])
                     .to_string(index=False))
        else:
            print("\n[一致性] 无可配对记录")
    else:
        print(f"\n[一致性] 未找到 {RAW_PANEL}，跳过对比")

    print("\nDONE")


if __name__ == "__main__":
    main()
