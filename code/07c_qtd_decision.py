# -*- coding: utf-8 -*-
"""
07c_qtd_decision.py - QTD → 细胞类型最终判别
基于 OneK1K 官方细胞类型定义的 marker 基因（是否有显著 eQTL, p_perm<0.05）判定。
输出存在矩阵 + 建议标签。
"""
import os, glob, json
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EQTC_DIR = os.path.join(BASE, "data", "onek1k", "eqtl_catalogue")
ENSG2SYM = json.load(open(os.path.join(BASE, "data", "onek1k", "ensg2sym.json"), encoding="utf-8"))
RES = os.path.join(BASE, "results")

# OneK1K 官方 marker（Yazar 2022 注释）
MARKERS = {
    "CD4": ["CD4", "IL7R"],
    "CD8": ["CD8A", "CD8B"],
    "KLRB1(CD4ET)": ["KLRB1"],
    "LTB(CD8NC)": ["LTB"],
    "NKG7/GNLY(CD8ET/NK)": ["NKG7", "GNLY"],
    "S100B(CD8S100B)": ["S100B"],
    "SOX4(CD4SOX4)": ["SOX4"],
    "NK": ["KLRD1", "NCAM1", "KLRF1"],
    "XCL1(NKR)": ["XCL1"],
    "B": ["CD19", "MS4A1"],
    "TCL1A/FCER2(BIN)": ["TCL1A", "FCER2"],
    "CD27(BMEM)": ["CD27"],
    "Plasma": ["JCHAIN", "SDC1", "MZB1"],
    "MonoC": ["CD14"],
    "MonoNC": ["FCGR3A"],
    "DC": ["CD1C", "CLEC9A", "FCER1A"],
}


def load_qtd(qid):
    fp = os.path.join(EQTC_DIR, f"{qid}.permuted.tsv.gz")
    df = pd.read_csv(fp, sep="\t", compression="gzip", low_memory=False)
    df["gene_symbol"] = df["molecular_trait_id"].map(ENSG2SYM)
    sig = set(df.loc[df["p_perm"] < 0.05, "gene_symbol"].dropna())
    return sig, df


def main():
    qtd_list = sorted(glob.glob(os.path.join(EQTC_DIR, "QTD*.permuted.tsv.gz")))
    rows = []
    for fp in qtd_list:
        qid = os.path.basename(fp).split(".")[0]
        sig, _ = load_qtd(qid)
        row = {"QTD": qid, "n_sig_genes": len(sig)}
        for label, genes in MARKERS.items():
            row[label] = ",".join(g for g in genes if g in sig)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "qtd_marker_matrix.csv"), index=False)
    pd.set_option("display.max_colwidth", 40)
    print(df.to_string())


if __name__ == "__main__":
    main()
