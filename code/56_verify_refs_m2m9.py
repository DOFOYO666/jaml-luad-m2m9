# -*- coding: utf-8 -*-
"""Verify every reference DOI in M2-M9深化研究稿/manuscript_EN.md against Crossref.
Outputs a machine-readable JSON and a human-readable report.
Cross-project reusable: point REFS at any {n: doi} mapping.
"""
import json, sys, time, urllib.request, urllib.parse, os

sys.stdout.reconfigure(encoding="utf-8")
UA = {"User-Agent": "refcheck/1.0 (mailto:refcheck@example.org)"}

REFS = {
    1: "10.1126/science.abf3041",
    2: "10.1038/ng.3892",
    3: "10.1038/s41586-022-04713-1",
    4: "10.1084/jem.20202644",
    5: "10.1186/s12935-022-02517-x",
    6: "10.12659/MSM.933503",
    7: "10.2147/LCTT.S616193",
    8: "10.1182/blood-2002-11-3462",
    9: "10.1083/jcb.200805061",
    10: "10.3390/cells7040025",
    11: "10.1038/nri2096",
    12: "10.1126/science.275.5304.1320",
    13: "10.1126/science.1187996",
    14: "10.1126/science.1192698",
    15: "10.1038/s41568-019-0153-5",
    16: "10.1038/s41467-020-16164-1",
    17: "10.1016/j.immuni.2019.03.009",
    18: "10.1016/j.cell.2017.06.010",
    19: "10.1038/s41467-019-12159-9",
    20: "10.1038/s41588-021-00924-w",
    21: "10.1038/s41588-021-00913-z",
    22: "10.1186/s13059-016-1070-5",
    23: "10.7554/eLife.26476",
    24: "10.1038/nmeth.3337",
    25: "10.1038/nbt.2203",
    26: "10.1093/aje/kws412",
    27: "10.1093/ije/dyx206",
    28: "10.1038/s41467-022-30755-0",
    29: "10.1038/nmeth.3971",
    30: "10.1186/s13059-017-1382-0",
    31: "10.1093/bioadv/vbac016",
    32: "10.1093/nar/gkad841",
    33: "10.1038/s41586-022-05688-9",
    34: "10.1186/s13059-021-02540-7",
    35: "10.1111/j.2517-6161.1995.tb02031.x",
    36: "10.1001/jama.2021.18236",
    37: "10.1016/j.cell.2021.01.010",
    38: "10.1038/s41586-020-2134-y",
    39: "10.1038/s41586-021-03651-8",
    40: "10.1038/s41577-019-0127-6",
    41: "10.1038/s41467-020-19015-1",
    42: "10.1093/bioinformatics/btz363",
    43: "10.1038/ncomms9971",
    44: "10.1038/nri2507",
    45: "10.1038/ni1403",
    46: "10.1038/ni903",
    47: "10.1038/s41592-019-0690-6",
    48: "10.1172/JCI91190",
    49: "10.1016/j.immuni.2018.09.024",
}

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
os.makedirs(OUTDIR, exist_ok=True)

out = []
for n, doi in REFS.items():
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    rec = {"ref": n, "doi": doi}
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=45) as r:
                m = json.load(r)["message"]
            rec["title"] = (m.get("title") or [""])[0]
            ct = m.get("container-title") or [""]
            rec["journal"] = ct[0] if ct else ""
            rec["year"] = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
            rec["vol"] = m.get("volume", "")
            rec["issue"] = m.get("issue", "")
            rec["page"] = m.get("page", "")
            au = m.get("author") or []
            rec["first_author"] = (au[0].get("family", "") + " " + au[0].get("given", "")) if au else ""
            rec["n_auth"] = len(au)
            rec["ok"] = True
            break
        except Exception as e:
            rec["ok"] = False
            rec["err"] = repr(e)[:200]
            time.sleep(2)
    out.append(rec)
    print(json.dumps(rec, ensure_ascii=False), flush=True)
    time.sleep(0.4)

with open(os.path.join(OUTDIR, "refcheck_m2m9.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

lines = ["# Reference DOI verification (Crossref)", ""]
lines.append("| # | DOI | Crossref title | Journal | Year | Vol(Issue):Pages | First author |")
lines.append("|---|---|---|---|---|---|---|")
for r in out:
    if r.get("ok"):
        lines.append("| {ref} | {doi} | {title} | {journal} | {year} | {vol}({issue}):{page} | {first_author} |".format(**r))
    else:
        lines.append("| {ref} | {doi} | **FAILED** {err} | | | | |".format(**r))
with open(os.path.join(OUTDIR, "refcheck_m2m9.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print("WROTE", os.path.join(OUTDIR, "refcheck_m2m9.md"))
