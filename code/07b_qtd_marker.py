# -*- coding: utf-8 -*-
"""
07b_qtd_marker.py - 用细胞类型 marker 基因检测 QTD → OneK1K 细胞类型
统计每个 QTD 显著 eQTL 基因集合中各类别 marker 基因的命中比例，
命中率最高者即为该 QTD 的细胞类型。
"""
import os, glob, json
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EQTC_DIR = os.path.join(BASE, "data", "onek1k", "eqtl_catalogue")
ENSG2SYM = json.load(open(os.path.join(BASE, "data", "onek1k", "ensg2sym.json"), encoding="utf-8"))
RES = os.path.join(BASE, "results")

MARKERS = {
    "CD4_T": ["CD4", "IL7R", "LEF1", "TCF7", "CCR7", "KLRB1"],
    "CD8_T": ["CD8A", "CD8B", "LTB", "NKG7", "GZMK", "GZMH"],
    "NK": ["KLRD1", "KLRF1", "NCAM1", "GNLY", "XCL1", "XCL2", "FGFBP2"],
    "B_cell": ["CD19", "MS4A1", "CD79A", "CD79B", "TCL1A", "FCER2", "CD27", "CD24"],
    "Plasma": ["JCHAIN", "SDC1", "XBP1", "MZB1", "IGHG1", "IGHA1"],
    "Mono_C": ["CD14", "S100A8", "S100A9", "LYZ", "VCAN"],
    "Mono_NC": ["FCGR3A", "MS4A7", "MS4A4A", "LST1"],
    "DC": ["CD1C", "CLEC9A", "FLT3", "FCER1A", "ITGAX"],
    "SOX4_Treg": ["SOX4", "FOXP3", "IL2RA", "CTLA4"],
    "S100B_CD8": ["S100B"],
}


def load_symbols(qid):
    fp = os.path.join(EQTC_DIR, f"{qid}.permuted.tsv.gz")
    df = pd.read_csv(fp, sep="\t", compression="gzip", usecols=["molecular_trait_id"],
                     low_memory=False)
    syms = {ENSG2SYM.get(v) for v in df["molecular_trait_id"].astype(str)}
    syms.discard(None)
    return syms


def main():
    qtd_list = sorted(glob.glob(os.path.join(EQTC_DIR, "QTD*.permuted.tsv.gz")))
    results = []
    for fp in qtd_list:
        qid = os.path.basename(fp).split(".")[0]
        syms = load_symbols(qid)
        n = len(syms)
        row = {"QTD": qid, "n_genes": n}
        for mtype, mlist in MARKERS.items():
            hit = [g for g in mlist if g in syms]
            row[mtype] = ",".join(hit)
            row[f"{mtype}_n"] = len(hit)
        results.append(row)

    df = pd.DataFrame(results)
    # 决策：取命中数最多的 marker 组（>0 时），并列时取多
    mcols = [f"{m}_n" for m in MARKERS]
    df["best_type"] = df[mcols].idxmax(axis=1).str.replace("_n", "")
    df["best_n"] = df[mcols].max(axis=1)
    print(df[["QTD", "n_genes"] + [f"{m}_n" for m in MARKERS] + ["best_type", "best_n"]].to_string())
    df.to_csv(os.path.join(RES, "qtd_marker_annotation.csv"), index=False)


if __name__ == "__main__":
    main()
