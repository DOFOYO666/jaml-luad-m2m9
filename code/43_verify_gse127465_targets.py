# -*- coding: utf-8 -*-
"""对 GSE127465 复核结果做**独立重实现的定向核对**。

动机：主脚本（38_gse127465_replication.py）的输出里，JAML 与 CXADR 在
`Endothelial cells` 与 `Fibroblasts` 上的**阳性细胞百分比完全相同**
（0.747664%、0.512821%），而 mean_expr 不同。这可能是巧合（两者恰好各有 4 / 3 个阳性细胞），
也可能是索引串味的征兆。**不做独立核对就无法区分。**

本脚本用**不同写法**重算这 4 个组合：
  * 不用"按列累加 + j2k 字典"，改为**先把目标基因的列号排好，再按 i 是否属于目标细胞组**计数；
  * 只扫描到最大目标列号即停（节省时间）；
  * 输出阳性细胞的**绝对计数**（整数），便于人工核对 4/535、3/585 这种量级。

用法：python 43_verify_gse127465_targets.py
"""
import gzip
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "rawdata", "scRNA", "GSE127465")
RES = os.path.join(ROOT, "results")
MTX = os.path.join(DATA, "GSE127465_human_counts_normalized_54773x41861.mtx.gz")
GENES = os.path.join(DATA, "GSE127465_gene_names_human_41861.tsv.gz")
META = os.path.join(DATA, "GSE127465_human_cell_metadata_54773x25.tsv.gz")

TARGETS = ["JAML", "CXADR", "CD3E", "SFTPC", "EPCAM", "LYZ"]
# 需要核对的细胞组（含疑似巧合的四个组合，以及作为对照的高表达组）
GROUPS = ["Endothelial cells", "Fibroblasts", "Smooth muscle cells", "ND",
          "Type I cells", "Type II cells", "Club cells", "Ciliated cells",
          "tMoMacDC", "tpDC", "bMonocytes", "tT cells", "bT cells",
          "tNK cells", "bB cells", "tB cells", "tMast cells"]


def main():
    meta = pd.read_csv(META, sep="\t")
    major = meta["Major cell type"].astype(str).values
    n_cells = len(major)

    names = [l.rstrip() for l in gzip.open(GENES, "rt")]
    sym2idx = {}
    for i, g in enumerate(names, start=1):
        sym2idx.setdefault(g, i)
    tgt_j = {g: sym2idx[g] for g in TARGETS if g in sym2idx}
    miss = [g for g in TARGETS if g not in sym2idx]
    print(f"[基因] 目标 {len(tgt_j)} 个；未命中 {miss or '无'}")
    print("  列号:", tgt_j)
    jmax = max(tgt_j.values())
    inv = {v: g for g, v in tgt_j.items()}

    gsel = [g for g in GROUPS if g in set(major)]
    gcode = np.full(n_cells, -1, dtype=np.int32)
    for k, g in enumerate(gsel):
        gcode[major == g] = k
    n_g = len(gsel)
    print(f"[细胞组] {n_g} 个；细胞数:", {g: int((major == g).sum()) for g in gsel})
    print(f"[扫描] 只扫到列号 {jmax}（文件按列聚集，可提前停止）")

    acc = np.zeros((len(tgt_j), n_g, 2))     # [表达和, 阳性数]
    order = {j: k for k, j in enumerate(sorted(tgt_j.values()))}
    n_line = 0
    stopped_at = None
    with gzip.open(MTX, "rt", errors="replace") as f:
        hdr = False
        for line in f:
            if not hdr:
                if line.startswith("%"):
                    continue
                hdr = True
                continue
            a, b, c = line.split()
            j = int(b)
            if j > jmax:                    # 列号已越过全部目标 → 可安全停止
                stopped_at = j
                break
            k = order.get(j)
            if k is None:
                continue
            n_line += 1
            i = int(a) - 1
            ci = gcode[i]
            if ci < 0:
                continue
            v = float(c)
            acc[k, ci, 0] += v
            if v > 0:
                acc[k, ci, 1] += 1.0

    print(f"[完成] 目标列数据行 {n_line:,}；在列号 {stopped_at} 处停止（> {jmax}）")

    counts = np.array([int((major == g).sum()) for g in gsel], dtype=float)
    print("\n=== 定向核对：阳性细胞数（绝对计数）与百分比 ===")
    rows = []
    for j in sorted(tgt_j.values()):
        g = inv[j]
        k = order[j]
        for ci, grp in enumerate(gsel):
            pos = int(acc[k, ci, 1])
            rows.append({"gene": g, "group": grp, "n_cells": int(counts[ci]),
                         "n_positive": pos,
                         "pct_positive": 100.0 * pos / counts[ci],
                         "mean_expr": acc[k, ci, 0] / counts[ci]})
    df = pd.DataFrame(rows)

    focus = [("JAML", "Endothelial cells"), ("CXADR", "Endothelial cells"),
             ("JAML", "Fibroblasts"), ("CXADR", "Fibroblasts"),
             ("JAML", "tMoMacDC"), ("JAML", "tT cells"),
             ("CXADR", "Type II cells"), ("CXADR", "Club cells"),
             ("CD3E", "tT cells"), ("SFTPC", "Type II cells"),
             ("EPCAM", "Type II cells"), ("LYZ", "tMoMacDC")]
    print(f"{'gene':8s} {'group':22s} {'n':>6s} {'pos':>5s} {'pct':>9s} {'mean':>10s}")
    for g, grp in focus:
        r = df[(df["gene"] == g) & (df["group"] == grp)]
        if len(r) == 0:
            print(f"{g:8s} {grp:22s}   --")
            continue
        r = r.iloc[0]
        print(f"{g:8s} {grp:22s} {int(r['n_cells']):6d} {int(r['n_positive']):5d} "
              f"{r['pct_positive']:9.6f} {r['mean_expr']:10.6f}")
    df.to_csv(os.path.join(RES, "verify_gse127465_targets.csv"), index=False)
    print("\n[输出] results/verify_gse127465_targets.csv")

    # 与主脚本输出逐格比对
    main_p = os.path.join(RES, "m8_gse127465_panel_by_celltype.csv")
    if os.path.isfile(main_p):
        mp = pd.read_csv(main_p)
        mp = mp[mp["level"] == "Major"]
        comp = df.merge(mp, on=["gene", "group"], suffixes=("_verify", "_main"))
        comp["d_pct"] = (comp["pct_positive_verify"] - comp["pct_positive_main"]).abs()
        comp["d_n"] = (comp["n_cells_verify"] - comp["n_cells_main"]).abs()
        bad = comp[(comp["d_pct"] > 1e-9) | (comp["d_n"] > 0)]
        print(f"\n[比对] 与主脚本输出可比 {len(comp)} 格；不一致 {len(bad)} 格")
        if len(bad):
            print(bad[["gene", "group", "n_positive", "pct_positive_verify",
                       "pct_positive_main"]].to_string(index=False))
        else:
            print("        → 全部逐格一致（独立重实现通过）")
        with open(os.path.join(RES, "verify_gse127465_targets_summary.json"),
                  "w", encoding="utf-8") as f:
            json.dump({"n_compared": int(len(comp)), "n_mismatch": int(len(bad)),
                       "targets": TARGETS, "n_groups": n_g,
                       "note": "独立重实现（按列号排序 + 提前停止），用于排除索引串味"},
                      f, ensure_ascii=False, indent=2)
    print("\nDONE")


if __name__ == "__main__":
    main()
