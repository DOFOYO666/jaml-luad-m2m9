# -*- coding: utf-8 -*-
"""How do recent Journal of Translational Medicine papers name and cite their
supplementary files?  Evidence for the submission package's labelling.

Run: python scripts/62_jtm_supp_probe.py
"""
import json
import os
import re
import sys
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")
UA = {"User-Agent": "fmt-probe/1.0 (mailto:fmtprobe@example.org)"}
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "jtm_style")
os.makedirs(OUT, exist_ok=True)

QUERY = ('JOURNAL:"Journal of translational medicine" AND OPEN_ACCESS:Y AND HAS_FT:Y '
         'AND (PUB_YEAR:2024 OR PUB_YEAR:2025)')

u = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
     + urllib.parse.quote(QUERY) + "&format=json&pageSize=12")
ids = [r["pmcid"] for r in json.load(urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60))
       ["resultList"]["result"] if r.get("pmcid")]
print("papers:", ids)

tot = {"Supplementary Material": 0, "Additional file": 0, "Supplementary Table": 0,
       "Supplementary Fig": 0, "Supplementary Data": 0, "MOESM": 0}
for pid in ids:
    try:
        x = urllib.request.urlopen(urllib.request.Request(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML" % pid,
            headers=UA), timeout=90).read().decode("utf-8", "replace")
    except Exception as e:
        print(pid, "ERR", repr(e)[:70])
        continue
    row = {k: x.count(k) for k in tot}
    for k, v in row.items():
        tot[k] += v
    caps = re.findall(r"<caption><p>(.*?)</p></caption>", x)
    intext = set(re.findall(r"[Ss]upplementary (?:Material|Table|Fig(?:ure)?|Data|File)[^<.,;)]{0,18}", x))
    print("%s | %s" % (pid, row))
    print("    caps:", caps[:3])
    print("    in-text:", sorted(intext)[:6])
print("\nTOTALS:", tot)
