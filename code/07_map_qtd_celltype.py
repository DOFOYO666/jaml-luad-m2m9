# -*- coding: utf-8 -*-
"""
07_map_qtd_celltype.py v2 - QTD → 细胞类型映射（基因集合重叠法）
用 gene symbol 集合重叠（Zenodo 版 vs eQTL Catalogue 版）确定 QTD 对应的 OneK1K 细胞类型，
并生成每基因 top cis-eQTL 工具变量表（含 variant, beta, p_perm, ENSG, symbol, cell_type）。
"""
import os, gzip, glob, json
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ONEX1K_DIR = os.path.join(BASE, "data", "onek1k", "extracted", "OneK1K_TensorQTL_top_eQTL_summary")
EQTC_DIR = os.path.join(BASE, "data", "onek1k", "eqtl_catalogue")
ENSG2SYM = json.load(open(os.path.join(BASE, "data", "onek1k", "ensg2sym.json"), encoding="utf-8"))
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

CELL_TYPES = ["B_IN", "B_MEM", "CD4_ET", "CD4_NC", "CD4_SOX4", "CD8_ET", "CD8_NC",
              "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK", "NK_R", "Plasma"]


def zenodo_symbols(ct):
    """Zenodo 版某细胞类型的 gene symbol 集合"""
    files = glob.glob(os.path.join(ONEX1K_DIR, f"OneK1K_{ct}.sig_cis_qtl_pairs.chr*.csv"))
    syms = set()
    for fp in files:
        df = pd.read_csv(fp, sep="\t", usecols=["phenotype_id"], low_memory=False)
        for v in df["phenotype_id"].astype(str):
            syms.add(v.split(".")[0])
    return syms


def qtd_symbols(qid):
    """eQTL Catalogue QTD 的 gene symbol 集合"""
    fp = os.path.join(EQTC_DIR, f"{qid}.permuted.tsv.gz")
    df = pd.read_csv(fp, sep="\t", compression="gzip", usecols=["molecular_trait_id"],
                     low_memory=False)
    syms = set()
    for v in df["molecular_trait_id"].astype(str):
        s = ENSG2SYM.get(v, None)
        if s:
            syms.add(s)
    return syms


def main():
    zen = {ct: zenodo_symbols(ct) for ct in CELL_TYPES}
    for ct, s in zen.items():
        print(f"[Zenodo] {ct}: {len(s)} symbols")

    qtd_list = sorted(glob.glob(os.path.join(EQTC_DIR, "QTD*.permuted.tsv.gz")))
    mapping = {}
    iv_parts = []
    for fp in qtd_list:
        qid = os.path.basename(fp).split(".")[0]
        qs = qtd_symbols(qid)
        best_ct, best_j = None, -1
        for ct, zs in zen.items():
            inter = len(qs & zs)
            union = len(qs | zs)
            j = inter / union if union else 0
            if j > best_j:
                best_j, best_ct = j, ct
        mapping[qid] = (best_ct, best_j, len(qs))
        print(f"[QTD] {qid}: -> {best_ct} (Jaccard {best_j:.3f}, n_sym {len(qs)})")
        df = pd.read_csv(fp, sep="\t", compression="gzip", low_memory=False)
        df["cell_type"] = best_ct
        df["gene_symbol"] = df["molecular_trait_id"].map(ENSG2SYM)
        iv_parts.append(df)

    iv = pd.concat(iv_parts, ignore_index=True)
    iv = iv.dropna(subset=["gene_symbol"])
    iv = iv.drop_duplicates(subset=["cell_type", "molecular_trait_id"], keep="first")
    # 仅保留 permutation 显著 eQTL
    iv = iv[iv["p_perm"] < 0.05]
    iv.to_csv(os.path.join(RES, "onek1k_top_eqtl_instruments.csv"), index=False)
    print(f"\n[完成] 工具变量表: {iv.shape}")
    print(f"  细胞类型分布:\n{iv['cell_type'].value_counts().to_string()}")
    with open(os.path.join(RES, "qtd_celltype_mapping.txt"), "w") as f:
        for qid, (ct, j, n) in mapping.items():
            f.write(f"{qid}\t{ct}\t{j:.3f}\t{n}\n")


if __name__ == "__main__":
    main()
