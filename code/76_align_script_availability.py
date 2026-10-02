"""Align the "where are the scripts" statements with the fact of a public archive.

Two defects found on 2026-10-02, both consequences of the DOI write-back:

1. **The manuscripts contradict their own Availability section.** The Availability sections say
   that all analysis scripts, the derived result tables and the reference-verification record are
   deposited with the additional files and archived publicly at the Zenodo DOI -- yet
   Methods 2.13 (EN and CN) still ended with "scripts are available from the corresponding
   author", and the CN Availability section had a duplicated "provided as Supplementary
   material 17" plus a stale "and can also be requested from the corresponding author".
   A statement that the code is available on request, next to a statement that it is already
   public, is internally inconsistent and weaker than the journal's data-availability policy.

2. **scripts/67_finalise_review2.py --doi now aborts.** That script used to patch the two
   Availability sections from a stored placeholder. Round 5 (scripts/74) rewrote both sections
   into the journal's prescribed wording, so the placeholder no longer exists: `sub_once()` finds
   0 matches and raises SystemExit. scripts/70_github_publish.py --writeback calls that command,
   so the whole write-back path was broken. set_doi() is therefore turned into a *verifier*.

Replacement is exact-string, exactly-once, and aborts without writing if any anchor is missing.
"""

import io
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(BASE, "M2-M9深化研究稿")
EN = os.path.join(PKG, "manuscript_EN.md")
CN = os.path.join(PKG, "manuscript_CN.md")
S67 = os.path.join(BASE, "scripts", "67_finalise_review2.py")

EN_OLD = "Scripts are available from the corresponding author."
EN_NEW = ("The scripts that produced every figure and table, and the derived result tables behind "
          "them, are deposited with the additional files and publicly archived; see Availability "
          "of data and materials.")

CN_OLD_METHOD = "脚本可向通讯作者索取。"
CN_NEW_METHOD = "全部脚本与派生结果表均已随补充材料提交并公开存档，见“数据与材料可用性”声明。"

CN_OLD_AVAIL = ("均作为补充材料 17 提供。内容相同的归档（release 1.0.0）已在 "
                "https://github.com/DOFOYO666/jaml-luad-m2m9 公开发布"
                "（DOI：10.5281/zenodo.23076399），并同时作为补充材料 17 提供；"
                "亦可向通讯作者索取。")
CN_NEW_AVAIL = ("均作为补充材料 17 提供；内容相同的归档（release 1.0.0）已在 "
                "https://github.com/DOFOYO666/jaml-luad-m2m9 公开发布"
                "（DOI：10.5281/zenodo.23076399）。")

S67_SET_DOI = r'''def set_doi(doi, url):
    """Verify -- do NOT rewrite -- that the repository DOI is wired into both manuscripts.

    This function used to patch the two Availability sections from a stored placeholder. Round 5
    (scripts/74_apply_jtm_guidelines.py) rewrote them in the journal's prescribed wording, so the
    placeholder no longer exists and patching them from here would now either abort (0 matches)
    or re-insert the retired "available from the corresponding author" sentence. The write-back
    path in scripts/70 calls this command, so it must succeed and be idempotent when the DOI is
    already present. Patching is scripts/74's job; this checks its work.
    """
    log("[2] verify the repository DOI is wired into the manuscripts")
    if not doi.startswith("10."):
        log("  WARNING: '%s' does not look like a DOI" % doi)
    missing = []
    for path in (EN, CN):
        t = io.open(path, encoding="utf-8").read()
        has_doi, has_url = doi in t, url in t
        log("  %-16s DOI: %-5s URL: %-5s" % (os.path.basename(path), has_doi, has_url))
        if not (has_doi and has_url):
            missing.append(os.path.basename(path))
    if missing:
        raise SystemExit("DOI/URL not present in: %s\n"
                         "  -> rewrite the Availability sections with scripts/74_apply_jtm_guidelines.py\n"
                         "     and re-run the write-back" % ", ".join(missing))
    if os.path.exists(JTM_MD):
        log("  note: manuscript_JTM.md is generated from the EN master - re-run scripts/64 then 65")
        log("        (and scripts/49 for the companion docx pair)")
'''

S67_OLD_CONST_HEAD = "PLACEHOLDER_EN = ("
S67_OLD_CONST_TAIL = "\n\n\ndef log(*a):"
S67_NEW_CONST = '''# The Availability sections are written in full by scripts/74_apply_jtm_guidelines.py in the
# journal's prescribed wording. The placeholders that used to live here were retired on
# 2026-10-02: that is why set_doi() verifies instead of patching (see its docstring).'''


def sub_once(path, old, new, label):
    t = io.open(path, encoding="utf-8").read()
    n = t.count(old)
    if n != 1:
        raise SystemExit("[%s] matched %d times (expected 1): %r" % (label, n, old[:80]))
    io.open(path, "w", encoding="utf-8").write(t.replace(old, new, 1))
    print("  OK   %s" % label)


def main():
    print("[1] English master - Methods 2.13")
    sub_once(EN, EN_OLD, EN_NEW, "EN :: script availability")

    print("[2] Chinese master - Methods 2.13 and the Availability section")
    sub_once(CN, CN_OLD_METHOD, CN_NEW_METHOD, "CN :: script availability")
    sub_once(CN, CN_OLD_AVAIL, CN_NEW_AVAIL, "CN :: Availability redundancy + stale request clause")

    print("[3] scripts/67 - turn set_doi() into a verifier")
    t = io.open(S67, encoding="utf-8").read()
    a, b = t.index(S67_OLD_CONST_HEAD), t.index(S67_OLD_CONST_TAIL)
    i, j = t.index("def set_doi(doi, url):"), t.index("\n\ndef main():")
    assert a < b < i < j
    t2 = t[:a] + S67_NEW_CONST + t[b:i] + S67_SET_DOI + t[j:]
    io.open(S67, "w", encoding="utf-8").write(t2)
    print("  OK   67 :: set_doi -> verifier, placeholder constants retired")

    print()
    print("[4] checks")
    en = io.open(EN, encoding="utf-8").read()
    cn = io.open(CN, encoding="utf-8").read()
    s67 = io.open(S67, encoding="utf-8").read()
    checks = [
        ("EN no longer says 'available from the corresponding author'",
         "corresponding author" not in en.split("## 1. Introduction")[1].split("## Figure")[0]),
        ("EN points at the Availability section", "see Availability of data and materials." in en),
        ("CN no longer says 脚本可向通讯作者索取", "脚本可向通讯作者索取" not in cn),
        ("CN no longer says 亦可向通讯作者索取", "亦可向通讯作者索取" not in cn),
        ("CN Availability mentions 补充材料 17 only once",
         cn.count("均作为补充材料 17 提供") == 1),
        ("CN still carries the DOI", "10.5281/zenodo.23076399" in cn),
        ("EN still carries the DOI", "10.5281/zenodo.23076399" in en),
        ("67 has no leftover PLACEHOLDER constant", "PLACEHOLDER_EN" not in s67),
        ("67 no longer writes the corresponding-author sentence",
         "corresponding author on reasonable request" not in s67),
        ("67 verifies instead", "verify the repository DOI is wired" in s67),
    ]
    bad = 0
    for name, ok in checks:
        print(("  OK   " if ok else "  FAIL ") + name)
        bad += (not ok)
    if bad:
        raise SystemExit("%d check(s) failed" % bad)


if __name__ == "__main__":
    main()
