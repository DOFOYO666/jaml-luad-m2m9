# -*- coding: utf-8 -*-
"""Finalise the second-round review paperwork and wire in a repository DOI.

Modes
  (default)               refresh the Supplementary Material 17 facts (file count, size,
                          SHA-256) in the three review documents, now that the archive also
                          carries the D-drive M4/M8 scripts.
  --doi X --url Y         replace the "being deposited" placeholder in both masters and in
                          the JTM manuscript with the real repository URL and DOI, then remind
                          you which builders to re-run.

Run: python scripts/67_finalise_review2.py
     python scripts/67_finalise_review2.py --doi 10.5281/zenodo.1234567 --url https://github.com/me/jaml
"""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(BASE, "M2-M9深化研究稿")
JTM = os.path.join(BASE, "JTM投稿_M2M9")
EN = os.path.join(PKG, "manuscript_EN.md")
CN = os.path.join(PKG, "manuscript_CN.md")
JTM_MD = os.path.join(JTM, "manuscript_JTM.md")

ARCHIVE = ("125 个文件：79 个 C 盘工作区脚本 + 19 个 D 盘 M4/M8 脚本 + 10 个派生结果表 "
           "+ 16 个补充数据表 + `MANIFEST.txt`；604,458 B；sha256[:16] = `ec1107c39f5daca2`"
           "（完整值见 `04_Supplementary/AdditionalFile11_code_and_results_SHA256.txt`）")

FACT_OPS = [
    (os.path.join(PKG, "README_投稿说明.md"),
     "（100 个文件：73 个 Python 脚本 + 16 个补充数据表 + 10 个派生结果表 + `MANIFEST.txt`；"
     "494,324 B；sha256[:16] = `bbf06e038548facf`）",
     "（" + ARCHIVE + "）"),
    (os.path.join(PKG, "审稿意见与逐条修改说明_第二轮.md"),
     "（代码与派生结果表压缩包；100 个文件 = 73 个 Python 脚本 + 16 个补充数据表 + "
     "10 个派生结果表 + 清单；494,324 字节；sha256[:16] = `bbf06e038548facf`）",
     "（代码与派生结果表压缩包；" + ARCHIVE + "）"),
    (os.path.join(PKG, "审稿意见与逐条修改说明_第二轮.md"),
     "AdditionalFile11_code_and_results.zip    （新增，100 个文件）",
     "AdditionalFile11_code_and_results.zip    （新增，125 个文件）"),
    (os.path.join(PKG, "投稿期刊推荐.md"),
     "（100 个文件，494,324 B，sha256[:16] = `bbf06e038548facf`）",
     "（" + ARCHIVE + "）"),
]

PLACEHOLDER_EN = ("An identical archive (release 1.0.0) is being deposited in a public repository; "
                  "its DOI will be quoted here once issued, and the archive remains available from "
                  "the corresponding author on reasonable request in the interim.")
PLACEHOLDER_CN = ("内容相同的归档（release 1.0.0）正存入公共仓库；其 DOI 一经取得即在此处引用，"
                  "在此之前可向通讯作者索取。")


def log(*a):
    print(*a, flush=True)


def sub_once(path, old, new, label):
    t = io.open(path, encoding="utf-8").read()
    if old not in t and new in t:
        log("  ..  already applied: " + label)
        return
    n = t.count(old)
    if n != 1:
        raise SystemExit("[%s] matched %d times (expected 1)" % (label, n))
    io.open(path, "w", encoding="utf-8").write(t.replace(old, new, 1))
    log("  OK  " + label)


def update_facts():
    log("[1] refresh the Supplementary Material 17 facts in the review documents")
    for path, old, new in FACT_OPS:
        sub_once(path, old, new, os.path.basename(path) + " :: archive facts")


def set_doi(doi, url):
    log("[2] wire the repository DOI into the manuscripts")
    cite_en = ("An identical archive (release 1.0.0) is openly available at %s (DOI: %s) and as "
               "Supplementary Material 17; the archive also remains available from the "
               "corresponding author on reasonable request." % (url, doi))
    cite_cn = ("内容相同的归档（release 1.0.0）已在 %s 公开发布（DOI：%s），并同时作为补充材料 17 提供；"
               "亦可向通讯作者索取。" % (url, doi))
    for path, old, new in ((EN, PLACEHOLDER_EN, cite_en), (CN, PLACEHOLDER_CN, cite_cn)):
        sub_once(path, old, new, os.path.basename(path) + " :: DOI")
    # the JTM manuscript is generated from EN, so regenerate it rather than patching it by hand
    if os.path.exists(JTM_MD):
        log("  note: manuscript_JTM.md is generated from the EN master - re-run scripts/64 then 65")
        log("        (and scripts/49 for the companion docx pair)")
    if not doi.startswith("10."):
        log("  WARNING: '%s' does not look like a DOI" % doi)


def main():
    if "--doi" in sys.argv:
        doi = sys.argv[sys.argv.index("--doi") + 1]
        url = sys.argv[sys.argv.index("--url") + 1] if "--url" in sys.argv else "(URL not given)"
        set_doi(doi, url)
    else:
        update_facts()


if __name__ == "__main__":
    main()
