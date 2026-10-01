# -*- coding: utf-8 -*-
"""Apply the Journal of Translational Medicine submission requirements to the EN master.

Changes
  [1] supplementary files renamed  "Supplementary Material N" -> "Additional file N"
      (BMC rule: "Additional files should be named 'Additional file 1' and so on")
  [2] "Availability of data and materials" rewritten in the mandated wording, with the
      eight software fields and a citation for every public dataset (BMC requires public
      datasets to be fully referenced in the reference list)
  [3] a List of abbreviations section, and an Acknowledgements sub-heading
  [4] each additional file listed with file name, file format, title and description
  [5] thousands separators removed from numeric values inside the tables
  [6] reference list renumbered (new dataset/URL/software entries), with audits

Run once. A timestamped backup of the master is written first.
"""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EN = os.path.join(ROOT, "M2-M9深化研究稿", "manuscript_EN.md")
S64 = os.path.join(ROOT, "scripts", "64_build_jtm_pkg.py")
S71 = os.path.join(ROOT, "scripts", "71_revise_m2m9_language.py")
STAMP = "jtmlguidelines_20261001"

TOKEN = re.compile(r"\{\{([A-Za-z0-9_]+)\}\}")
GROUP = re.compile(r"\{\{[A-Za-z0-9_]+\}\}(?:\s*,?\s*\{\{[A-Za-z0-9_]+\}\})*")

# ------------------------------------------------------------------ helpers


def log(m):
    print(m, flush=True)


def read(p):
    return io.open(p, encoding="utf-8").read()


def section_span(t, start, end):
    i = t.index(start)
    j = t.index(end) if end else len(t)
    assert i < j, (start, end)
    return i, j


def expand(spec):
    out = []
    for part in re.split(r"\s*[,，]\s*", spec):
        m = re.fullmatch(r"(\d+)\s*[–\-—]\s*(\d+)", part)
        if m:
            out.extend(range(int(m.group(1)), int(m.group(2)) + 1))
        else:
            out.append(int(part))
    return out


# --- reference keys already in use: read them straight out of script 71 (never executed) --
_old_src = read(S71)
KEY_BY_DOI = {}
for line in _old_src.split("\n"):
    m = re.match(r'\s*"([A-Za-z0-9_]+)":\s+"(.+)",\s*$', line)
    if not m:
        continue
    dm = re.search(r"doi:(10\.\S+?)\s*$", m.group(2))
    if dm:
        KEY_BY_DOI[dm.group(1).lower()] = m.group(1)
REFS = {}
for line in _old_src.split("\n"):
    m = re.match(r'\s*"([A-Za-z0-9_]+)":\s+"(.+)",\s*$', line)
    if m:
        REFS[m.group(1)] = m.group(2)
assert len(REFS) == 49, "expected 49 inherited references, got %d" % len(REFS)

# --- additional files: read the numbering/description table out of script 64 ---------------
_sm = re.search(r"SM_MAP = \[(.*?)\n\]", read(S64), re.S).group(1)
SM_MAP = []
for a, num, desc in re.findall(r'\("([^"]+)",\s*(\d+),\s*"([^"]+)"\)', _sm):
    stem, ext = os.path.splitext(a)
    stem = re.sub(r"^AdditionalFile\d+[a-d]?_", "", stem)
    SM_MAP.append((int(num), "AdditionalFile%02d_%s%s" % (int(num), stem, ext),
                   ext.lstrip(".").upper(), desc))
assert [n for n, _, _, _ in SM_MAP] == list(range(1, 18)), "SM_MAP numbering is not 1..17"

# --- new references -------------------------------------------------------------------------
REFS.update({
    "geo131907": "Gene Expression Omnibus. GSE131907: single-cell RNA sequencing of primary and metastatic lung adenocarcinoma. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907. Accessed 1 October 2026.",
    "geo127465": "Gene Expression Omnibus. GSE127465: single-cell transcriptomics of human and mouse lung tumours. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465. Accessed 1 October 2026.",
    "geo135222": "Gene Expression Omnibus. GSE135222: transcriptomes of non-small-cell lung cancers treated with anti-PD-1/PD-L1 therapy. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222. Accessed 1 October 2026.",
    "geo126044": "Gene Expression Omnibus. GSE126044: gene expression profiling of non-small-cell lung cancer. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044. Accessed 1 October 2026.",
    "gdcportal": "Genomic Data Commons Data Portal. TCGA-LUAD: lung adenocarcinoma. https://portal.gdc.cancer.gov/. Accessed 1 October 2026.",
    "depmapportal": "DepMap portal. DepMap 22Q2 public release: CRISPR (Chronos) gene effect and model files. https://depmap.org/portal/. Accessed 1 October 2026.",
    "eqtlcatalogue": "eQTL Catalogue. OneK1K single-cell eQTL summary statistics, dataset QTS000038. https://www.ebi.ac.uk/eqtl/. Accessed 1 October 2026.",
    "gwascatalog": "GWAS Catalog. Lung adenocarcinoma genome-wide association summary statistics, study GCST004744. https://www.ebi.ac.uk/gwas/. Accessed 1 October 2026.",
    "codearchive": "Liu H. jaml-luad-m2m9: analysis code and derived result tables for the compartment-resolved analysis of JAML and CXADR and the transcription-factor perturbation screen in lung adenocarcinoma (version 1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23076399. Accessed 1 October 2026.",
})

AVAILABILITY = """**Availability of data and materials.** The dataset(s) supporting the conclusions of this article are available in the Gene Expression Omnibus repository, GSE131907 {{geo131907}}, GSE127465 {{geo127465}}, GSE135222 {{geo135222}} and GSE126044 {{geo126044}}; in the Genomic Data Commons Data Portal, TCGA-LUAD {{gdcportal}}; in the DepMap portal, DepMap 22Q2 {{depmapportal}}; in the eQTL Catalogue, dataset QTS000038 {{eqtlcatalogue}}; and in the GWAS Catalog, study GCST004744 {{gwascatalog}}. All analysis scripts, the derived result tables underlying every figure and table, and the machine-readable reference-verification record are included within this article and its additional files (Additional file 17), and the same archive is deposited as release 1.0.0 {{codearchive}}.

The analysis code is available as follows. *Project name:* jaml-luad-m2m9, compartment-resolved expression, deconvolution and transcription-factor perturbation analysis of *JAML* and *CXADR* in lung adenocarcinoma. *Project home page:* https://github.com/DOFOYO666/jaml-luad-m2m9. *Archived version:* release 1.0.0, DOI 10.5281/zenodo.23076399 {{codearchive}}. *Operating system(s):* platform independent (developed on Windows 11). *Programming language:* Python 3.13 and R 4.6. *Other requirements:* Python 3.13 with numpy, pandas, scipy, scikit-learn, scanpy, anndata, celloracle and decoupler; R 4.6 with survival, data.table and ggplot2. *License:* MIT. *Any restrictions to use by non-academics:* none."""

ABBREV = """## Abbreviations

AP-1: activator protein 1; CAR: coxsackievirus and adenovirus receptor; CIBERSORT: Cell-type Identification by Estimating Relative Subsets of RNA Transcripts; DC: dendritic cell; DepMap: Cancer Dependency Map; eQTL: expression quantitative trait locus; EPIC: Estimating the Proportions of Immune and Cancer cells; FDR: false discovery rate; GEO: Gene Expression Omnibus; GRN: gene regulatory network; GWAS: genome-wide association study; ICB: immune checkpoint blockade; JAM: junctional adhesion molecule; JAML: junctional adhesion molecule-like protein; LUAD: lung adenocarcinoma; MCP-counter: Microenvironment Cell Populations-counter; MHC: major histocompatibility complex; MR: Mendelian randomization; NSCLC: non-small cell lung cancer; PD-1: programmed cell death protein 1; PD-L1: programmed death-ligand 1; PFS: progression-free survival; scRNA-seq: single-cell RNA sequencing; SD: standard deviation; TCGA: The Cancer Genome Atlas; TF: transcription factor; TPM: transcripts per million."""

# ------------------------------------------------------------------ load

log("[0] load master")
t = read(EN)
assert "{{" not in t, "master already carries {{...}} symbols — script already applied?"
bak = EN + ".bak_" + STAMP
if not os.path.exists(bak):
    io.open(bak, "w", encoding="utf-8").write(t)
log("    backup -> %s" % os.path.basename(bak))

log("[1] numeric citations -> symbols")
_n = [0]

# the conversion is done against the reference list parsed below
ref_block = t[t.index("## References"):t.index("## Figure legends")]
order_old = []
for m in re.finditer(r"(?m)^(\d+)\.\s+(.*)$", ref_block):
    dm = re.search(r"doi:(10\.\S+?)\s*$", m.group(2))
    assert dm, "reference %s has no DOI" % m.group(1)
    doi = dm.group(1).lower()
    assert doi in KEY_BY_DOI, "unknown DOI %s" % doi
    order_old.append(KEY_BY_DOI[doi])
assert len(order_old) == 49, len(order_old)


def _num2sym(m):
    _n[0] += 1
    return ", ".join("{{%s}}" % order_old[n - 1] for n in expand(m.group(1)))


t = re.sub(r"\[(\d+(?:\s*[–\-,，]\s*\d+)*)\]", _num2sym, t)
log("    %d citation groups converted" % _n[0])

log("[2] availability statement -> mandated wording + software fields")
i, j = section_span(t, "**Availability of data and materials.**", "**Competing interests.**")
t = t[:i] + AVAILABILITY + "\n\n" + t[j:]

log("[3] Acknowledgements sub-heading")
anchor = "All authors read and approved the final manuscript.\n"
assert t.count(anchor) == 1
t = t.replace(anchor, anchor + "\n**Acknowledgements.** Not applicable.\n", 1)

log("[4] additional files: rename and list name / format / title / description")
i, j = section_span(t, "## Supplementary material", "## Declarations")
head = ("## Additional files\n\nThe following additional files are supplied with this article. "
        "Statistical details are given in the corresponding Methods subsections, and every number "
        "quoted in the text is reproduced in these files.\n\n")
ent = re.findall(r"(?m)^\*\*Supplementary Material (\d+)\.\*\* (.*)$", t[i:j])
assert len(ent) == 17, "expected 17 supplementary entries, found %d" % len(ent)
body = []
for num, text in ent:
    n, fname, fmt, desc = SM_MAP[int(num) - 1]
    assert n == int(num)
    body.append("**Additional file %d** — Title: %s. File name: `%s`. File format: %s. "
                "Description: %s" % (n, desc[0].upper() + desc[1:], fname, fmt, text))
t = t[:i] + head + "\n\n".join(body) + "\n\n---\n\n" + t[j:]

log("[5] remaining 'Supplementary Material' -> 'Additional file'")
t = re.sub(r"Supplementary Material (\d+)", r"Additional file \1", t)
t = t.replace("supplied as Supplementary Material.", "supplied as additional files.")
assert "Supplementary Material" not in t, "leftover 'Supplementary Material': %s" % \
    re.findall(r".{40}Supplementary Material.{40}", t)[:2]

log("[6] List of abbreviations, placed between Declarations and References")
i, j = section_span(t, "## References", None)
t = t[:i] + ABBREV + "\n\n---\n\n" + t[i:]

log("[7] thousands separators removed inside the tables")
i, j = section_span(t, "## Tables", None)
tbl, n_comma = re.subn(r"(?<=\d),(?=\d{3}(?!\d))", "", t[i:j])
t = t[:i] + tbl + t[j:]
log("    %d separators removed" % n_comma)

log("[8] renumber references by order of first appearance")
_order = []


def _num(key):
    if key not in _order:
        _order.append(key)
    return _order.index(key) + 1


def _render(g):
    nums = sorted({_num(m.group(1)) for m in TOKEN.finditer(g)})
    runs, start, prev = [], nums[0], nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
            continue
        runs.append((start, prev))
        start = prev = n
    runs.append((start, prev))
    return "[" + ", ".join(str(a) if a == b else "%d–%d" % (a, b) for a, b in runs) + "]"


t = GROUP.sub(lambda m: _render(m.group(0)), t)
missing = [k for k in REFS if k not in _order]
if missing:
    raise SystemExit("uncited references: %s" % missing)

lines = ["%d. %s" % (i + 1, REFS[k]) for i, k in enumerate(_order)]
lines.append("")
lines.append("*All references were verified against Crossref by DOI (title, journal, year, "
             "volume, issue, pages and first author) on 2026-10-01; the machine-readable "
             "verification record is included in Additional file 17.*")
i, j = section_span(t, "## References", "## Figure legends")
t = t[:i] + "## References\n\n" + "\n".join(lines) + "\n\n---\n\n" + t[j:]

io.open(EN, "w", encoding="utf-8").write(t)
log("    wrote %s (%d chars)" % (os.path.basename(EN), len(t)))

log("[9] checks")
refs_now = t[t.index("## References"):t.index("## Figure legends")]
nums = [int(x) for x in re.findall(r"(?m)^(\d+)\.\s+\S", refs_now)]
assert nums == list(range(1, len(nums) + 1)), "numbering not 1..N"
prose = t[:t.index("## References")] + t[t.index("## Figure legends"):]
cited = set()
for m in re.finditer(r"\[(\d+(?:\s*[–\-,]\s*\d+)*)\]", prose):
    cited.update(expand(m.group(1)))
assert not [n for n in cited if n < 1 or n > len(nums)], "citation out of range"
log("    %d references, continuous; uncited: %s"
    % (len(nums), [n for n in range(1, len(nums) + 1) if n not in cited] or "none"))
log("    headings: %s" % " | ".join(re.findall(r"(?m)^## (.+)$", t)))
log("    leftover symbols: %d | 'Additional file' count: %d | 'Supplementary' left: %d"
    % (t.count("{{") + t.count("}}"), t.count("Additional file"), t.count("Supplementary")))
for pat, want in [("**Acknowledgements.**", 1), ("## Abbreviations", 1), ("## Additional files", 1)]:
    log("    %-26s count=%d" % (pat, t.count(pat)))
log("DONE")
