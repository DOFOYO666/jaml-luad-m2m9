# -*- coding: utf-8 -*-
"""Fetch a recent Journal of Translational Medicine open-access paper and dump its
section/back-matter skeleton plus the first few references, so that the manuscript
can be restyled to the journal's actual conventions (verified, not remembered).

Run: python scripts/61_jtm_style_probe.py
"""
import json
import os
import re
import sys
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")
UA = {"User-Agent": "fmt-check/1.0 (mailto:fmtcheck@example.org)"}
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "jtm_style")

QUERY = ('JOURNAL:"Journal of translational medicine" AND OPEN_ACCESS:Y AND HAS_FT:Y '
         'AND (PUB_YEAR:2024 OR PUB_YEAR:2025) AND HAS_SUPPL:Y')


def get(url, raw=False):
    r = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90)
    b = r.read()
    return b if raw else json.loads(b.decode("utf-8"))


def main():
    os.makedirs(OUT, exist_ok=True)
    u = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
         + urllib.parse.quote(QUERY) + "&format=json&pageSize=25&resultType=core")
    d = get(u)
    res = d["resultList"]["result"]
    print("hits:", d.get("hitCount"), "returned:", len(res))
    picked = None
    for r in res:
        pmcid = r.get("pmcid")
        title = (r.get("title") or "")[:100]
        print(" ", pmcid, r.get("pubYear"), "|", title)
        if pmcid and picked is None:
            picked = pmcid
    if not picked:
        sys.exit("no PMC full text found")

    xml = get("https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML" % picked, raw=True).decode("utf-8", "replace")
    open(os.path.join(OUT, picked + ".xml"), "w", encoding="utf-8").write(xml)
    print("\nsaved", picked + ".xml", len(xml), "chars")

    # section headings, in document order
    print("\n--- sec titles ---")
    for t in re.findall(r"<title>(.*?)</title>", xml):
        print("  ", re.sub(r"<[^>]+>", "", t)[:110])

    # abstract structure
    m = re.search(r"<abstract.*?</abstract>", xml, re.S)
    if m:
        print("\n--- abstract skeleton ---")
        for t in re.findall(r"<title>(.*?)</title>", m.group(0)):
            print("  *", re.sub(r"<[^>]+>", "", t))
        txt = re.sub(r"<[^>]+>", " ", m.group(0))
        print("  abstract word count:", len(re.findall(r"[A-Za-z][A-Za-z'\-]*", txt)))

    # back matter headings (BMC "Declarations")
    print("\n--- back matter headings (kw: Declarations/Ethics/Consent/Availability/Competing/Funding/Contributions/Acknowledg/Additional) ---")
    for t in re.findall(r"<title>(.*?)</title>", xml):
        s = re.sub(r"<[^>]+>", "", t)
        if re.search(r"Declar|Ethic|Consent|Availab|Competing|Funding|Contribut|Acknowledg|Additional|Abbreviat", s, re.I):
            print("  ", s[:120])

    # additional file naming
    print("\n--- additional-file strings ---")
    for s in set(re.findall(r"Additional file[^<]{0,90}", xml)):
        print("  ", s.strip()[:110])

    # first references, to read the style off real output
    refs = re.findall(r"<ref\b.*?</ref>", xml, re.S)
    print("\n--- n refs:", len(refs), "| first 4 rendered ---")
    for r in refs[:4]:
        s = re.sub(r"<[^>]+>", "", r)
        s = re.sub(r"\s+", " ", s).strip()
        print("  -", s[:230])


if __name__ == "__main__":
    main()
