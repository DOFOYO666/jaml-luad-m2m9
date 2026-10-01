# -*- coding: utf-8 -*-
"""Finalise the second-round review paperwork and wire in a repository DOI.

Modes
  (default)               recompute the Supplementary Material 17 facts (file count, size,
                          SHA-256) **from the archive itself** and refresh every place the
                          review documents quote them.
  --doi X --url Y         replace the "being deposited" placeholder in both masters with the
                          real repository URL and DOI, then remind you which builders to re-run.

Run: python scripts/67_finalise_review2.py
     python scripts/67_finalise_review2.py --doi 10.5281/zenodo.1234567 --url https://github.com/me/jaml
"""
import hashlib
import io
import os
import re
import sys
import zipfile

sys.stdout.reconfigure(encoding="utf-8")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(BASE, "M2-M9深化研究稿")
JTM = os.path.join(BASE, "JTM投稿_M2M9")
EN = os.path.join(PKG, "manuscript_EN.md")
CN = os.path.join(PKG, "manuscript_CN.md")
JTM_MD = os.path.join(JTM, "manuscript_JTM.md")

ARCHIVE = os.path.join(PKG, "04_Supplementary", "AdditionalFile11_code_and_results.zip")

DOC_README = os.path.join(PKG, "README_投稿说明.md")
DOC_REVIEW = os.path.join(PKG, "审稿意见与逐条修改说明_第二轮.md")
DOC_JOURNAL = os.path.join(PKG, "投稿期刊推荐.md")

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
        raise SystemExit("[%s] matched %d times (expected 1): %r" % (label, n, old[:70]))
    io.open(path, "w", encoding="utf-8").write(t.replace(old, new, 1))
    log("  OK  " + label)


def archive_facts():
    """Recompute everything the review documents quote about Supplementary Material 17."""
    z = zipfile.ZipFile(ARCHIVE)
    names = [n for n in z.namelist() if not n.endswith("/")]
    raw = open(ARCHIVE, "rb").read()
    f = {
        "n": len(names),
        "c": len([n for n in names if n.startswith("scripts/")]),
        "d": len([n for n in names if n.startswith("scripts_d_workspace/")]),
        "res": len([n for n in names if n.startswith("results/")]),
        "sup": len([n for n in names if n.startswith("04_Supplementary/")]),
        "size": os.path.getsize(ARCHIVE),
        "sha": hashlib.sha256(raw).hexdigest(),
    }
    f["short"] = f["sha"][:16]
    return f


def update_facts():
    """Refresh the archival facts in the review documents.

    Matching is **by shape, not by the previous literal value**: hardcoding "the string that is
    there now" is what let these numbers go stale twice (the refresh script only knew how to
    replace the wording it had written itself). Every pattern below must match exactly once.
    """
    log("[1] refresh the Supplementary Material 17 facts (recomputed from the archive)")
    f = archive_facts()
    core = (r"\d+ 个文件：\d+ 个 C 盘工作区脚本 \+ \d+ 个 D 盘 M4/M8 脚本 \+ \d+ 个派生结果表 "
            r"\+ \d+ 个补充数据表 \+ `MANIFEST\.txt`；[\d,]+ B；sha256\[:16\] = `[0-9a-f]{16}`")
    core_new = ("%d 个文件：%d 个 C 盘工作区脚本 + %d 个 D 盘 M4/M8 脚本 + %d 个派生结果表 "
                "+ %d 个补充数据表 + `MANIFEST.txt`；%d B；sha256[:16] = `%s`"
                % (f["n"], f["c"], f["d"], f["res"], f["sup"], f["size"], f["short"]))

    edits = [
        (DOC_README, core, core_new, "README :: archive facts"),
        (DOC_REVIEW, core, core_new, "review doc :: archive facts"),
        (DOC_JOURNAL, core, core_new, "journal doc :: archive facts"),
        (DOC_REVIEW, r"（新增，\d+ 个文件）", "（新增，%d 个文件）" % f["n"], "review doc :: tree listing"),
        (DOC_REVIEW, r"\| 文件数 \| \*\*\d+\*\*（[^|]*） \|",
         "| 文件数 | **%d**（%d C 盘脚本 + %d D 盘脚本 + %d 派生结果表 + %d 补充数据表 + `MANIFEST.txt`） |"
         % (f["n"], f["c"], f["d"], f["res"], f["sup"]), "review doc :: fact table row 1"),
        (DOC_REVIEW, r"\| 字节数 \| \*\*[\d,]+\*\* \|", "| 字节数 | **%d** |" % f["size"],
         "review doc :: fact table row 2"),
        (DOC_REVIEW, r"\| sha256 \| `[0-9a-f]{64}`（前 16 位 `[0-9a-f]{16}`） \|",
         "| sha256 | `%s`（前 16 位 `%s`） |" % (f["sha"], f["short"]),
         "review doc :: fact table row 3"),
    ]
    for path, pattern, repl, label in edits:
        t = io.open(path, encoding="utf-8").read()
        hits = re.findall(pattern, t)
        if len(hits) == 1 and hits[0] == repl:
            log("  ..  already current: " + label)
            continue
        if len(hits) != 1:
            raise SystemExit("[%s] pattern matched %d times (expected 1)" % (label, len(hits)))
        # lambda replacement: never let re.sub interpret backslashes in the new text
        io.open(path, "w", encoding="utf-8").write(re.sub(pattern, lambda m: repl, t, count=1))
        log("  OK  " + label)
    log("  facts now: %d files (%d C + %d D + %d results + %d supplementary + MANIFEST), "
        "%d B, sha256[:16] = %s" % (f["n"], f["c"], f["d"], f["res"], f["sup"], f["size"], f["short"]))


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
