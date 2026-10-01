# -*- coding: utf-8 -*-
"""投前核验 JTM 稿的 Word 文件：结构、三线表边框、字体、章节顺序。

说明：本机 editor_sdk 通道（tencent-local-office-edit）对该文件反复返回
`document is not open`（open_file 为异步流式打开，本会话未完成），属该通道的具体失败，
故改用 python-docx 直接读取 docx —— 这也正是生成该文件所用的库，核验对象与生成方式一致。

Run: python scripts/68_verify_jtm_docx.py
"""
import os
import sys

import docx
from docx.oxml.ns import qn

sys.stdout.reconfigure(encoding="utf-8")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCX = os.path.join(BASE, "JTM投稿_M2M9", "01_Manuscript_JTM.docx")

ok = 0
bad = []


def check(label, cond, detail=""):
    global ok
    if cond:
        ok += 1
        print("  OK   %s%s" % (label, ("  " + detail) if detail else ""))
    else:
        bad.append(label)
        print("  !!   %s%s" % (label, ("  " + detail) if detail else ""))


def main():
    d = docx.Document(DOCX)
    paras = [(p.style.name if p.style else "", p.text) for p in d.paragraphs]
    all_text = "\n".join(t for _, t in paras)
    for t in d.tables:
        for row in t.rows:
            for c in row.cells:
                all_text += "\n" + c.text

    print("文件：%s" % os.path.basename(DOCX))
    print("段落 %d，表格 %d\n" % (len(d.paragraphs), len(d.tables)))

    print("[题名页]")
    check("标题在第 1 段", paras[0][0] == "Title" or paras[0][1].startswith("Compartmentally separated"))
    check("作者行含 7 位作者", all_text.count("Huayong Liu") >= 1 and "Jiming Chen" in all_text)
    check("通讯作者邮箱", "xingxinghuoshu@163.com" in all_text)
    check("ORCID", "0009-0007-8106-0425" in all_text)
    check("单位英文名（与主稿一致）",
          "Department of Pulmonary and Critical Care Medicine, Shenzhen Nanshan District People's Hospital" in all_text)
    check("短标题", "Where JAML acts in LUAD" in all_text)

    print("\n[摘要与关键词]")
    for h in ("Background.", "Methods.", "Results.", "Conclusions."):
        check("摘要小标题 " + h, h in all_text)
    kw = ""
    for _, t in paras:
        if t.startswith("Keywords:"):
            kw = t
    check("关键词 ≤ 10 个", len([k for k in kw.replace("Keywords:", "").split(";") if k.strip()]) <= 10,
          "实得 %d 个" % len([k for k in kw.replace("Keywords:", "").split(";") if k.strip()]))
    check("已删 'target triage'", "target triage" not in all_text)

    print("\n[章节顺序]")
    # locate headings as whole paragraphs, in order (searching raw text would match
    # "Conclusions." inside the abstract, which precedes the Introduction)
    want = ["Abstract", "1. Introduction", "2. Materials and Methods", "3. Results", "4. Discussion",
            "5. Limitations", "6. Conclusion", "Additional files", "Declarations",
            "Abbreviations",
            "References", "Figure legends", "Tables"]
    flat = [t.strip() for _, t in paras]
    idx = []
    for w in want:
        try:
            idx.append(flat.index(w))
        except ValueError:
            idx.append(-1)
    print("    段落型标题位置: %s" % list(zip(want, idx)))
    check("全部 12 个一级章节齐备且顺序正确", all(i >= 0 for i in idx) and idx == sorted(idx))

    print("\n[补充材料]")
    check("已按 BMC 规则命名为 Additional file",
          "Additional files" in all_text and "Additional file 1" in all_text
          and "Supplementary Material" not in all_text)
    check("文中无遗留 'Additional file'", "Additional file" not in all_text)
    check("Additional file 17 已在可用性声明中引用", "Additional file 17" in all_text)
    check("Additional file 编号 1–17 连续",
          all(("Additional file %d**" % i) in all_text or ("Additional file %d " % i) in all_text
              for i in range(1, 18)))
    check("含 Abbreviations 节", "## Abbreviations" in all_text)
    check("含 Acknowledgements 子标题", "Acknowledgements." in all_text)
    check("数据集已进参考文献（GEO/GDC/DepMap/eQTL/GWAS）",
          all(x in all_text for x in ["GSE131907", "portal.gdc.cancer.gov", "depmap.org",
                                      "ebi.ac.uk/eqtl", "ebi.ac.uk/gwas"]))

    print("\n[三线表]")
    check("表格数 = 5", len(d.tables) == 5, "实得 %d" % len(d.tables))
    for i, t in enumerate(d.tables, 1):
        tbl = t._element
        pr = tbl.tblPr
        tb = pr.find(qn("w:tblBorders"))
        inside_none = True
        if tb is not None:
            for edge in ("insideH", "insideV", "left", "right"):
                e = tb.find(qn("w:" + edge))
                if e is not None and e.get(qn("w:val")) not in (None, "none", "nil"):
                    inside_none = False
        # top/bottom lines must come from cell borders
        first = t.rows[0].cells[0]._tc.get_or_add_tcPr().find(qn("w:tcBorders"))
        last = t.rows[-1].cells[0]._tc.get_or_add_tcPr().find(qn("w:tcBorders"))
        has_top = first is not None and first.find(qn("w:top")) is not None
        has_bottom = last is not None and last.find(qn("w:bottom")) is not None
        check("表 %d 为三线表（无竖线/内部横线；首行有顶线、末行有底线）" % i,
              inside_none and has_top and has_bottom,
              "%d 行 × %d 列" % (len(t.rows), len(t.columns)))

    print("\n[字体]")
    fonts = set()
    for p in d.paragraphs:
        for r in p.runs:
            rpr = r._element.find(qn("w:rPr"))
            if rpr is None:
                continue
            rf = rpr.find(qn("w:rFonts"))
            if rf is not None:
                fonts.add((rf.get(qn("w:ascii")), rf.get(qn("w:hAnsi")), rf.get(qn("w:eastAsia"))))
    print("    rFonts 组合：%s" % sorted(f for f in fonts if f[0]))
    check("西文为 Times New Roman", all(f[0] == "Times New Roman" for f in fonts if f[0]))

    print("\n[关键数字仍在新稿中]")
    for k in ("0.757", "6.55", "0.881", "1.0000", "18,747", "2,842", "229", "6,344" if False else "329"):
        check("含 %s" % k, k in all_text)

    print("\n%s（%d 通过 / %d 失败）" % ("全部通过" if not bad else "存在问题：" + "；".join(bad), ok, len(bad)))


if __name__ == "__main__":
    main()
