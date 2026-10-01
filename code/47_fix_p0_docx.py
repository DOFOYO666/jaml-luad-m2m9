# -*- coding: utf-8 -*-
"""把母稿 md 的 P0 修正同步到两个投稿包的 01_Manuscript_main.docx。

改动四项（全部为 run 级精确替换，不用跨行正则）：
  A. Limitations：删除无数据支撑的 "OR up to ~3.0"，改为实测范围与方向（衰减而非膨胀）
  B. Results §3.7：TCGA 生存由"全阴性"改为如实报告（18/22 可评估；JAML 由 log-rank 显著、
     Cox 不显著；MAP4K4 相反）
  C. Methods §2.9：补上"只在矩阵内的基因被检验 + Cox 与 log-rank 两种检验 + FDR 分母"
  D. Table 1 汇总行：同步 TCGA 结果
  E. 新增 Methods §2.10（单细胞表达分析）——因正文已包含 GSE127465 新结果
     （见 scripts/46_insert_replication_sentence_docx.py），按用户规则须补 Methods 描述。

安全措施：改前做带时间戳备份；每处替换都断言"原文只落在单个 run 内"，不满足即报错退出，
不做任何跨 run 的正则拼接。改后逐项回读校验。

用法：python 47_fix_p0_docx.py [--dry-run]
"""
import copy
import io
import os
import shutil
import sys
import time

import docx
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKGS = ["BMCCancer投稿", "JTM投稿"]
DOC = "01_Manuscript_main.docx"

DRY = "--dry-run" in sys.argv

# ---------------------------------------------------------------- 替换内容

OLD_A = ("IVW estimates were inflated and unstable (OR up to ~3.0), reflecting "
         "weak-instrument bias (NOME violation) and winner's curse inherent to cis-MR;")
NEW_A = ("IVW estimates became less precise and were generally attenuated towards the null "
         "rather than inflated (median |log OR| 0.14 versus 0.20 with F > 10; all OR estimates "
         "between 0.48 and 1.16), as expected under weak-instrument bias (NOME violation and "
         "regression dilution);")

OLD_B = ("None of the 22 genes showed a significant TCGA-LUAD expression\u2013survival association "
         "after multiple-testing correction, consistent with the notion that germline-predicted "
         "risk effects need not translate into tumor-expression prognostic markers.")
NEW_B = [
    ("Eighteen of the 22 genes could be evaluated in TCGA-LUAD (RNA-seq, n = 563; 208 deaths); "
     "the remaining four loci (", None),
    ("ZNRD1ASP", True), (", ", None), ("CTC-490E21.14", True), (", ", None),
    ("RP1-167A14.3", True), (", ", None), ("RP11-514O12.4", True),
    (") were absent from the expression matrix and were therefore not tested. Using Cox "
     "proportional-hazards and log-rank tests with Benjamini\u2013Hochberg correction across the "
     "18 evaluable genes, no gene was significant by both tests: higher ", None),
    ("JAML", True),
    (" tumour expression was associated with longer overall survival by log-rank (HR per SD 0.85; "
     "log-rank P = 6.8\u00d710\u207b\u2074, FDR = 0.012) but not by Cox regression (P = 0.017, "
     "FDR = 0.11), whereas ", None),
    ("MAP4K4", True),
    (" showed the converse pattern (Cox P = 0.001, FDR = 0.018; log-rank FDR = 0.17). These "
     "observational associations are non-causal and need not mirror the direction of "
     "germline-predicted risk effects.", None),
]

OLD_C = ("LUAD tumor expression (RNA-seq, n = 563) and overall survival from The Cancer Genome "
         "Atlas (TCGA) were used in Cox regression to test whether tumor expression of the "
         "identified genes relates to prognosis (descriptive, non-causal).")
NEW_C = ("LUAD tumor expression (RNA-seq, n = 563) and overall survival (208 deaths) from The "
         "Cancer Genome Atlas (TCGA) were used to test whether tumor expression of the identified "
         "genes relates to prognosis (descriptive, non-causal). Only genes present in the TCGA "
         "expression matrix were tested; both Cox proportional-hazards regression (per-SD "
         "expression) and log-rank tests (median split) were performed, and "
         "Benjamini\u2013Hochberg correction was applied across the genes tested.")

OLD_D = "No significant association (after correction)"
NEW_D = ("JAML: better OS by log-rank (FDR = 0.012), not significant by Cox (FDR = 0.11); "
         "MAP4K4: significant by Cox (FDR = 0.018), not by log-rank")

HEAD_210 = "2.10 Tumour single-cell RNA-seq expression analyses"
BODY_210 = [
    ("To test whether the cell-type-restricted expression of the candidate genes is reproducible "
     "in independent tumour tissue, two publicly available single-cell RNA-seq datasets from "
     "treatment-na\u00efve NSCLC were analysed. GSE117570 (8 samples from 4 patients; 11,485 "
     "immune cells) was used to describe ", None),
    ("JAML", True),
    (" expression across tumour-infiltrating immune subsets. GSE127465 (Zilionis et al.; 54,773 "
     "cells from tumour and blood) was used as an independent replication dataset: cell-level "
     "expression was extracted from the GEO-supplied normalised count matrix and summarised as "
     "positivity and mean expression across the author-annotated major and minor cell types, "
     "which were then grouped into eight immune and non-immune compartments. Because the matrix "
     "and the cell metadata are supplied as separate files, row correspondence was verified "
     "before grouping by (i) correlating per-cell mitochondrial read fractions computed from the "
     "matrix against the corresponding metadata values (Spearman \u03c1 = 1.0000, n = 54,773) and "
     "(ii) confirming expected marker\u2013compartment directions (e.g. ", None),
    ("CD3E", True), (" in T cells, ", None), ("EPCAM", True),
    (" in epithelium); an independent re-implementation of the extraction reproduced all "
     "compartmental summaries exactly (Additional file 12). These expression analyses are "
     "descriptive and are not used for causal inference.", None),
]


# ---------------------------------------------------------------- 工具函数

def _run_span(full_text, old):
    idx = full_text.find(old)
    if idx < 0:
        return None
    return idx, idx + len(old)


def _set_run_text(run, text, italic=None):
    """清空 run 的子元素（保留 rPr），写入单段文本；italic=True/False 显式设置。"""
    for child in list(run._element):
        if child.tag != qn("w:rPr"):
            run._element.remove(child)
    t = OxmlElement("w:t")
    t.text = text
    t.set(qn("xml:space"), "preserve")
    run._element.append(t)
    if italic is not None:
        rpr = run._element.find(qn("w:rPr"))
        if rpr is None:
            rpr = OxmlElement("w:rPr")
            run._element.insert(0, rpr)
        for tag in ("w:i", "w:iCs"):
            el = rpr.find(qn(tag))
            if italic and el is None:
                rpr.append(OxmlElement(tag))
            elif not italic and el is not None:
                rpr.remove(el)


def replace_rich(par, old, new_parts, where):
    """在 par 内把 old 替换为 new_parts 描述的富文本；要求 old 完全落在单个 run 内。"""
    runs = par.runs
    texts = [r.text for r in runs]
    full = "".join(texts)
    sp = _run_span(full, old)
    if sp is None:
        return f"[跳过] {where}：未找到原文（可能已改过）"
    idx, end = sp
    pos, target = 0, None
    for i, r in enumerate(runs):
        s, e = pos, pos + len(r.text)
        pos = e
        if s <= idx and end <= e:
            target = (i, s, e)
            break
    if target is None:
        raise AssertionError(f"{where}：原文跨越多个 run，拒绝跨 run 拼接替换")
    i, s, e = target
    run = runs[i]
    pre = run.text[: idx - s]
    suf = run.text[end - s:]
    _set_run_text(run, pre if pre else "", italic=None)
    anchor = run._element
    last = anchor
    for txt, it in new_parts:
        new_el = copy.deepcopy(run._element)
        last.addnext(new_el)
        last = new_el
        # 写入文本
        for child in list(new_el):
            if child.tag != qn("w:rPr"):
                new_el.remove(child)
        t = OxmlElement("w:t")
        t.text = txt
        t.set(qn("xml:space"), "preserve")
        new_el.append(t)
        rpr = new_el.find(qn("w:rPr"))
        if rpr is None:
            rpr = OxmlElement("w:rPr")
            new_el.insert(0, rpr)
        for tag in ("w:i", "w:iCs"):
            el = rpr.find(qn(tag))
            if it is True and el is None:
                rpr.append(OxmlElement(tag))
            elif it is False and el is not None:
                rpr.remove(el)
    if suf:
        new_el = copy.deepcopy(run._element)
        last.addnext(new_el)
        for child in list(new_el):
            if child.tag != qn("w:rPr"):
                new_el.remove(child)
        t = OxmlElement("w:t")
        t.text = suf
        t.set(qn("xml:space"), "preserve")
        new_el.append(t)
    return f"[替换] {where}"


def insert_paragraph_after(par, template_par, parts):
    """在 par 之后插入一个段落（以 template_par 的段落格式为模板），内容由 parts 给出。"""
    new_el = copy.deepcopy(template_par._element)
    par._element.addnext(new_el)
    newpar = Paragraph(new_el, par._parent)
    runs = newpar.runs
    if not runs:
        r = newpar.add_run("")
        runs = [r]
    keep = runs[0]
    for r in runs[1:]:
        r._element.getparent().remove(r._element)
    _set_run_text(keep, parts[0][0] if parts else "", italic=parts[0][1] if parts else None)
    anchor = keep._element
    for txt, it in parts[1:]:
        e = copy.deepcopy(keep._element)
        anchor.addnext(e)
        anchor = e
        for child in list(e):
            if child.tag != qn("w:rPr"):
                e.remove(child)
        t = OxmlElement("w:t")
        t.text = txt
        t.set(qn("xml:space"), "preserve")
        e.append(t)
        rpr = e.find(qn("w:rPr"))
        if rpr is None:
            rpr = OxmlElement("w:rPr"); e.insert(0, rpr)
        for tag in ("w:i", "w:iCs"):
            el = rpr.find(qn(tag))
            if it is True and el is None:
                rpr.append(OxmlElement(tag))
            elif it is False and el is not None:
                rpr.remove(el)
    return newpar


# ---------------------------------------------------------------- 主流程

def process(pkg):
    path = os.path.join(ROOT, pkg, DOC)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    log = []
    if not DRY:
        bak = path + ".bak_" + time.strftime("%Y%m%d_%H%M%S")
        shutil.copy2(path, bak)
        log.append(f"[备份] {os.path.relpath(bak, ROOT)}")
    d = docx.Document(path)

    # --- 定位 ---
    paras = d.paragraphs
    i_lim = next(i for i, p in enumerate(paras) if OLD_A in p.text)
    i_res = next(i for i, p in enumerate(paras) if OLD_B in p.text)
    i_m29 = next(i for i, p in enumerate(paras) if OLD_C in p.text)
    i_h29 = next(i for i, p in enumerate(paras) if p.text.strip().startswith("2.9 "))
    i_res_h = next(i for i, p in enumerate(paras) if p.text.strip().startswith("3. Results"))

    log.append(f"[定位] Limitations=P{i_lim}  Results§3.7=P{i_res}  Methods§2.9=P{i_m29}  "
               f"§2.9标题=P{i_h29}  3.Results=P{i_res_h}")

    if DRY:
        for tag, p in [("A", paras[i_lim]), ("B", paras[i_res]), ("C", paras[i_m29])]:
            print(f"  --- 将替换 {tag} (P{paras.index(p)}) 片段：{p.text[:80]!r}")
        for ti, t in enumerate(d.tables):
            for ri, row in enumerate(t.rows):
                for ci, cell in enumerate(row.cells):
                    if cell.text.strip() == OLD_D:
                        print(f"  --- 将替换 D：table T{ti} r{ri} c{ci}")
        print("  --- 将插入 §2.10（标题 + 正文）于 P%d 之后" % i_m29)
        print("  --- dry-run：未写入")
        return log

    # --- A/B/C：顺序执行（后面的索引不受影响，因 P 索引在替换后不变） ---
    log.append(replace_rich(paras[i_lim], OLD_A, [(NEW_A, None)], "A Limitations"))
    log.append(replace_rich(paras[i_res], OLD_B, NEW_B, "B Results §3.7"))
    log.append(replace_rich(paras[i_m29], OLD_C, [(NEW_C, None)], "C Methods §2.9"))

    # --- D：表格单元格 ---
    hit_d = 0
    for t in d.tables:
        for row in t.rows:
            for cell in row.cells:
                if cell.text.strip() == OLD_D:
                    replace_rich(cell.paragraphs[0], OLD_D, [(NEW_D, None)], "D Table row")
                    hit_d += 1
    if hit_d != 1:
        raise AssertionError(f"D：表格单元格命中 {hit_d} 处，期望 1 处")
    log.append("[替换] D Table 1 汇总行")

    # --- E：插入 §2.10（先插标题，再插正文，避免顺序错乱） ---
    body_par = d.paragraphs[i_m29]
    head = insert_paragraph_after(body_par, d.paragraphs[i_h29], [(HEAD_210, None)])
    insert_paragraph_after(head, body_par, BODY_210)
    log.append("[插入] E Methods §2.10（标题 + 正文）")

    d.save(path)
    log.append(f"[保存] {os.path.relpath(path, ROOT)}")
    return log


def main():
    for pkg in PKGS:
        print(f"════════ {pkg} ════════")
        for line in process(pkg):
            print("  " + line)
        print()
    if not DRY:
        print("=== 回读校验 ===")
        for pkg in PKGS:
            d = docx.Document(os.path.join(ROOT, pkg, DOC))
            t = "\n".join(p.text for p in d.paragraphs)
            checks = {
                "含 'up to ~3.0'": "up to ~3.0" in t,
                "含旧 TCGA 全阴性句": "None of the 22 genes showed" in t,
                "含新 Limitations 句": "generally attenuated towards the null" in t,
                "含新 TCGA 句（JAML log-rank）": "log-rank P = 6.8" in t,
                "含 §2.10 标题": "2.10 Tumour single-cell RNA-seq" in t,
                "含 §2.10 正文（ρ = 1.0000）": "Spearman \u03c1 = 1.0000" in t,
                "含 'No significant association (after correction)'": "No significant association (after correction)" in t,
            }
            print(f"  {pkg}:")
            for k, v in checks.items():
                flag = "OK  " if ((k.startswith("含 'up to") or k.startswith("含旧") or
                                   "after correction" in k) and not v) or \
                              (not (k.startswith("含 'up to") or k.startswith("含旧") or
                                    "after correction" in k) and v) else "!!  "
                print(f"    {flag}{k}: {v}")
            ns = [i for i, p in enumerate(d.paragraphs) if p.text.strip().startswith("2.")]
            print(f"    Methods 段落索引（前 4）: {ns[:4]}")


if __name__ == "__main__":
    main()
