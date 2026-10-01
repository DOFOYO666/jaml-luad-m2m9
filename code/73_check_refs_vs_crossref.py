# -*- coding: utf-8 -*-
"""Compare every reference in the manuscript against the Crossref record produced by
scripts/56_verify_refs_m2m9.py, field by field.

Checks per entry: first-author family name, journal (abbreviation expanded word by word),
year (a +/-1 difference is tolerated and reported separately, because online-first and
issue years legitimately differ), volume, issue and pages.

Cross-project reusable: point DOC at any manuscript using the same reference style.
"""
import io
import json
import os
import re
import sys
import unicodedata

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "M2-M9深化研究稿", "manuscript_EN.md")
REC = os.path.join(ROOT, "results", "refcheck_m2m9.json")

# Abbreviations actually used in this reference list -> words of the container title.
ABBR = {
    "nat rev immunol": "nature reviews immunology",
    "nat rev cancer": "nature reviews cancer",
    "nat genet": "nature genetics",
    "nat commun": "nature communications",
    "nat methods": "nature methods",
    "nat immunol": "nature immunology",
    "nat biotechnol": "nature biotechnology",
    "j exp med": "journal of experimental medicine",
    "j clin invest": "journal of clinical investigation",
    "j cell biol": "journal of cell biology",
    "genome biol": "genome biology",
    "nucleic acids res": "nucleic acids research",
    "cancer cell int": "cancer cell international",
    "med sci monit": "medical science monitor",
    "am j epidemiol": "american journal of epidemiology",
    "int j epidemiol": "international journal of epidemiology",
    "j r stat soc series b": "journal of the royal statistical society series b",
    "bioinform adv": "bioinformatics advances",
    "lung cancer: targets and therapy": "lung cancer: targets and therapy",
}


def deacc(s):
    return "".join(c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c))


def norm(s):
    """Fold the differences that are typographic rather than substantive."""
    s = deacc(str(s or "")).lower()
    for a, b in (("–", "-"), ("—", "-"), ("−", "-"), ("‐", "-")):
        s = s.replace(a, b)
    s = re.sub(r"\s+", " ", s).strip().rstrip(".")
    if s.startswith("the "):
        s = s[4:]
    return s


t = io.open(DOC, encoding="utf-8").read()
block = t[t.index("## References"):t.index("## Figure legends")]
entries = {}
for m in re.finditer(r"(?m)^(\d+)\.\s+(.*)$", block):
    entries[int(m.group(1))] = m.group(2)

recs = {r["ref"]: r for r in json.load(io.open(REC, encoding="utf-8"))}

fail = 0
notes = []
print("%-4s %-8s %-30s %s" % ("ref", "fields", "journal", "issue"))
for n in sorted(entries):
    line = entries[n]
    r = recs.get(n)
    if not r or not r.get("ok"):
        print("%-4d %-8s NO CROSSREF RECORD" % (n, "-"))
        fail += 1
        continue
    problems = []

    # first author family name
    fam = deacc(line.split(",")[0].strip().split()[0]).lower()
    _cr = (r.get("first_author") or "").split()
    if not _cr:
        notes.append("ref %d: Crossref has no personal first author (consortium record); "
                     "author field not checked" % n)
    elif fam != deacc(_cr[0]).lower():
        problems.append("author %r vs %r" % (fam, deacc(_cr[0]).lower()))

    # journal
    jm = re.search(r"\*([^*]+)\*", line)
    if not jm:
        problems.append("no journal field")
    else:
        ab = norm(jm.group(1))
        full = norm(ABBR.get(ab, ab))
        crj = norm(r.get("journal"))
        if full != crj and not crj.startswith(full) and not full.startswith(crj):
            problems.append("journal %r vs %r" % (full, r.get("journal")))

    # year
    ym = re.search(r"(\d{4});", line)
    if not ym:
        problems.append("no year field")
    else:
        dy = int(ym.group(1)) - int(r.get("year") or 0)
        if dy != 0:
            (notes if abs(dy) == 1 else problems).append(
                "year %s vs Crossref %s (offline/online)" % (ym.group(1), r.get("year")))

    # volume / issue / pages sit in  YYYY;VOL(ISSUE):PAGES.
    vm = re.search(r";\s*([^();]+)\(([^()]*)\):([^.\s]+)", line)
    if vm:
        vol, iss, pgs = vm.group(1).strip(), vm.group(2).strip(), vm.group(3).strip()
        if vol and r.get("vol") and vol != str(r["vol"]):
            problems.append("volume %s vs %s" % (vol, r["vol"]))
        if iss and r.get("issue") and iss != str(r["issue"]):
            problems.append("issue %s vs %s" % (iss, r["issue"]))
        cp = norm(r.get("page"))
        pn = norm(pgs)
        if pn and cp and pn not in cp and cp not in pn:
            problems.append("pages %s vs %s" % (pn, cp))
    else:
        notes.append("no volume/issue/page pattern (article-number style)")

    status = "OK" if not problems else "FAIL"
    if problems:
        fail += 1
    print("%-4d %-8s %-30s %s" % (n, status, (r.get("journal") or "")[:28],
                                  "; ".join(problems) if problems else ""))

print()
print("entries         : %d" % len(entries))
print("crossref records: %d" % len(recs))
print("mismatches      : %d" % fail)
if notes:
    print("tolerated notes :")
    for x in notes:
        print("   - " + x)
sys.exit(1 if fail else 0)
