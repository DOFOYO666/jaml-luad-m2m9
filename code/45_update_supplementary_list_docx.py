# -*- coding: utf-8 -*-
"""把新增的 Additional files 条目插入两个投稿包的正文 .docx。

做法：定位每条清单的最后一个条目，**deepcopy 该段落的 XML**（沿用其字体/缩进/bold 标签），
在其后依次插入新段落——保持与现有条目完全一致的排版，不做任何其他改动。

- BMCCancer：已有 AF1–11（含 STROBE-MR）→ 追加 **AF12–15**
- JTM：只有 AF1–10（原本缺 STROBE-MR）→ 追加 **AF11–15**（使编号与母稿一致）

**先备份**为 `01_Manuscript_main.docx.bak_YYYYMMDD`。

用法：python 45_update_supplementary_list_docx.py [--dry-run]
"""
import copy
import os
import shutil
import sys
from datetime import datetime

import docx
from docx.text.paragraph import Paragraph

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DRY = "--dry-run" in sys.argv

ENTRIES = [
    (11, "STROBE-MR checklist of recommended items to address in reports of Mendelian randomization "
         "studies, completed with the location of each item in this manuscript. "
         "(Additional file 11 - STROBE-MR checklist.csv)"),
    (12, "Alignment verification for the independent single-cell replication dataset (GSE127465; "
         "Zilionis et al. 2019, NSCLC, 54,773 cells). Includes the matrix-integrity check "
         "(44,663,765 non-zero entries, matching the MatrixMarket header), the primary row-order "
         "criterion (self-computed per-cell mitochondrial fraction versus the metadata mitochondrial "
         "percentage: Spearman \u03c1 = 1.0000 across 54,773 cells), five marker-direction checks, the "
         "reason the metadata \u201cTotal counts\u201d column was not used as a criterion, and the "
         "independent re-implementation comparison (102 cells, 0 mismatches). Comma-separated values "
         "(CSV). (results/supp_gse127465_alignment.csv)"),
    (13, "JAML and CXADR positivity across eight immune and non-immune compartments in the discovery "
         "dataset (GSE131907; LUAD tissue, 208,506 cells) and in the independent replication dataset "
         "(GSE127465; NSCLC, 54,773 cells), with epithelial-to-myeloid positivity ratios. "
         "Comma-separated values (CSV). (results/supp_gse127465_compartments.csv)"),
    (14, "GSE127465 cell-type expression panel: positivity and mean expression for the 35-gene marker "
         "panel across 35 major and 50 minor cell-type annotations (5,005 rows). Comma-separated "
         "values (CSV). (results/supp_gse127465_panel.csv)"),
    (15, "Compartmental separation of JAML (immune side) and CXADR (epithelial side) in the discovery "
         "and replication datasets. Portable Network Graphics (PNG) and Tagged Image File Format "
         "(TIFF). (results/figures/cross_dataset_compartment.png)"),
]


def find_para(doc, prefix):
    for i, p in enumerate(doc.paragraphs):
        if p.text.strip().startswith(prefix):
            return i, p
    return None, None


def clone_after(ref_para, model_para, n, body):
    """以 model_para 为模板，在 ref_para 之后插入一条新条目。返回新段落。"""
    new_p = copy.deepcopy(model_para._p)
    ref_para._p.addnext(new_p)
    para = Paragraph(new_p, ref_para._parent)
    runs = para.runs
    if not runs:
        raise RuntimeError("模板段落没有 run，无法复用其格式")
    runs[0].text = f"Additional file {n}."
    if len(runs) >= 2:
        runs[1].text = " " + body
        for r in runs[2:]:
            r._element.getparent().remove(r._element)
    else:
        runs[0].text = f"Additional file {n}. {body}"
    return para


def process(pkg, last_label, entries):
    path = os.path.join(ROOT, pkg, "01_Manuscript_main.docx")
    print(f"\n==== {pkg} ====")
    if not os.path.isfile(path):
        print("  [跳过] 找不到", path)
        return
    doc = docx.Document(path)
    idx, ref = find_para(doc, last_label)
    if ref is None:
        print(f"  [跳过] 未找到 '{last_label}'")
        return
    print(f"  末条 '{last_label}' 位于段落 {idx}；样式 style={ref.style.name if ref.style else None}，"
          f"run 数={len(ref.runs)}，run0.bold={ref.runs[0].bold if ref.runs else None}")

    n_existing = len(doc.paragraphs)
    if DRY:
        for n, body in entries:
            print(f"  [dry-run] 将插入：Additional file {n}. {body[:70]}…")
        return

    bak = path + f".bak_{datetime.now():%Y%m%d_%H%M%S}"
    shutil.copy2(path, bak)
    print(f"  [备份] {os.path.relpath(bak, ROOT)}")

    cur = ref
    for n, body in entries:
        cur = clone_after(cur, ref, n, body)
        print(f"  [插入] Additional file {n}.")
    doc.save(path)
    print(f"  [保存] 段落数 {n_existing} → {len(docx.Document(path).paragraphs)}")

    # 校验
    d2 = docx.Document(path)
    ok = True
    for n, _ in entries:
        i, _p = find_para(d2, f"Additional file {n}.")
        mark = "OK" if i is not None else "MISS"
        if i is None:
            ok = False
        print(f"  [校验] Additional file {n}. → {mark}")
    print("  [结果]", "全部插入成功" if ok else "有缺失，请人工检查")


def main():
    process("BMCCancer投稿", "Additional file 11.", ENTRIES[1:])          # 追加 12–15
    process("JTM投稿", "Additional file 10.", ENTRIES)                    # 追加 11–15
    print("\nDONE" + (" (dry-run)" if DRY else ""))


if __name__ == "__main__":
    main()
