# -*- coding: utf-8 -*-
"""07d_qtd_official_validation.py - eQTL Catalogue 官方数据交叉复核
对 24 个 QTD（QTS000038/OneK1K）的显著 eQTL 基因集合，
与 OneK1K 官方 TensorQTL 输出（14 大类显著基因集合）做重叠评分，
推断每个 QTD 最可能的细胞类型（官方数据对官方数据的独立验证）。

输出: results/qtd_official_validation.csv + 控制台对照现有标注
"""
import gzip
import json
import os
import re
import zipfile

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
EQTC = os.path.join(ROOT, "data", "onek1k", "eqtl_catalogue")
ZIP = os.path.join(ROOT, "data", "onek1k", "OneK1K_TensorQTL_top_eQTL_summary.zip")
ENSMAP = os.path.join(ROOT, "data", "onek1k", "ensg2sym.json")

# 现用标注（02 脚本）
CURRENT = {
    "QTD000606": "B cell", "QTD000607": "B cell", "QTD000608": "B cell",
    "QTD000609": "CD4 T", "QTD000611": "CD4 T", "QTD000613": "CD4 T",
    "QTD000612": "CD8 T", "QTD000614": "CD8 T", "QTD000619": "CD8 T",
    "QTD000615": "CD8 S100B", "QTD000616": "NK", "QTD000617": "NK",
    "QTD000620": "NK", "QTD000621": "NK", "QTD000628": "NK",
    "QTD000623": "Plasma", "QTD000626": "DC", "QTD000629": "CD4 SOX4",
}


def load_tensorqtl_genes():
    """TensorQTL 官方: 14 大类 -> 显著基因集合(qval<0.05)"""
    ct_genes = {}
    with zipfile.ZipFile(ZIP) as z:
        for n in z.namelist():
            if not n.endswith(".csv") or "OneK1K_" not in n or "__MACOSX" in n or "/._" in n:
                continue
            m = re.search(r"OneK1K_(.+?)\.sig_cis_qtl", n)
            if not m:
                continue
            ct = m.group(1).split("/")[-1]
            genes = ct_genes.setdefault(ct, set())
            with z.open(n) as f:
                header = f.readline().decode().rstrip().split("\t")
                qi = header.index("qval")
                pi = header.index("phenotype_id")
                for line in f:
                    c = line.decode(errors="replace").rstrip().split("\t")
                    if len(c) <= pi:
                        continue
                    if len(c) > qi and c[qi] not in ("NA", ""):
                        try:
                            if float(c[qi]) < 0.05:
                                genes.add(c[pi])
                        except ValueError:
                            pass
    return ct_genes


def load_qtd_genes():
    """24 个 QTD permuted -> 显著基因集合(p_perm<0.05)"""
    ensg2sym = json.load(open(ENSMAP, encoding="utf-8"))
    qtd_genes = {}
    for fn in sorted(os.listdir(EQTC)):
        if not fn.endswith(".permuted.tsv.gz"):
            continue
        qtd = fn.split(".")[0]
        genes = set()
        with gzip.open(os.path.join(EQTC, fn), "rt", errors="replace") as f:
            header = f.readline().split("\t")
            gi = header.index("molecular_trait_id")
            pi = header.index("p_perm") if "p_perm" in header else header.index("pvalue")
            for line in f:
                c = line.split("\t")
                if len(c) <= max(gi, pi):
                    continue
                try:
                    if float(c[pi]) < 0.05:
                        sym = ensg2sym.get(c[gi], c[gi])
                        genes.add(sym)
                except ValueError:
                    continue
        qtd_genes[qtd] = genes
    return qtd_genes


def main():
    print("[1] 读取 TensorQTL 官方 14 类显著基因集 ...")
    ct_genes = load_tensorqtl_genes()
    print(f"    细胞类型: {sorted(ct_genes.keys())}")
    for ct, g in sorted(ct_genes.items()):
        print(f"    {ct}: {len(g)} 基因")
    print("\n[2] 读取 24 QTD 显著基因集 ...")
    qtd_genes = load_qtd_genes()
    print(f"    QTD 数: {len(qtd_genes)}")

    print("\n[3] 重叠评分（富集比：实际重叠/期望重叠，背景~1.2万基因）...")
    rows = []
    BACKGROUND = 12000  # OneK1K 测试基因数近似
    for qtd in sorted(qtd_genes):
        qg = qtd_genes[qtd]
        nq = len(qg)
        best = []
        for ct, cg in ct_genes.items():
            nc = len(cg)
            inter = len(qg & cg)
            expected = max(1.0, nq * nc / BACKGROUND)
            ratio = inter / expected
            jac = inter / max(1, len(qg | cg))
            best.append((ratio, inter, jac, ct))
        best.sort(key=lambda x: (-x[0], -x[1]))
        b = best[0]
        b2 = best[1] if len(best) > 1 else None
        rows.append({
            "QTD": qtd, "n_sig_genes": nq,
            "top_celltype": b[3], "enrichment": round(b[0], 2),
            "overlap": b[1], "jaccard": round(b[2], 3),
            "second": f"{b2[3]}({b2[1]})" if b2 else "",
            "现用标注": CURRENT.get(qtd, "排除"),
            "一致": "✓" if CURRENT.get(qtd) and CURRENT[qtd] in b[3].replace("_", " ") else "✗"
        })
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "qtd_official_validation.csv"), index=False)
    pd.set_option("display.width", 200)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
