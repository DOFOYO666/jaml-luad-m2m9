# -*- coding: utf-8 -*-
"""Regenerate the REFS map inside scripts/56_verify_refs_m2m9.py from the manuscript's
current reference list, so the Crossref verification record matches the new numbering."""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EN = os.path.join(ROOT, "M2-M9深化研究稿", "manuscript_EN.md")
SCRIPT = os.path.join(ROOT, "scripts", "56_verify_refs_m2m9.py")

t = io.open(EN, encoding="utf-8").read()
block = t[t.index("## References"):t.index("## Figure legends")]
rows = re.findall(r"(?m)^(\d+)\.\s+(.*)$", block)
print("reference entries found:", len(rows))

pairs = []
for num, line in rows:
    m = re.search(r"doi:(10\.[^\s]+)", line)
    if not m:
        raise SystemExit("entry %s has no DOI: %s" % (num, line[:100]))
    pairs.append((int(num), m.group(1)))

assert [n for n, _ in pairs] == list(range(1, len(pairs) + 1)), "numbering not 1..N"
assert len({d for _, d in pairs}) == len(pairs), "duplicate DOI in the reference list"

s = io.open(SCRIPT, encoding="utf-8").read()
new_map = "REFS = {\n" + "".join('    %d: "%s",\n' % (n, d) for n, d in pairs) + "}\n"
s2, n = re.subn(r"REFS = \{.*?\n\}\n", new_map, s, count=1, flags=re.S)
assert n == 1, "REFS block not replaced"
io.open(SCRIPT, "w", encoding="utf-8").write(s2)
print("patched 56: %d DOIs written" % len(pairs))
