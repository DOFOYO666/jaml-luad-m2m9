# -*- coding: utf-8 -*-
"""Second-round supplementary artefacts for the M2-M9 manuscript.

 1. AdditionalFile10d - gene evaluability of GSE126044 (backs the Limitations statement).
 2. Integrity note appended to AdditionalFile06 (degenerate within-positive p/FDR columns).
 3. AdditionalFile11 - archive of analysis code + derived result tables + reference-verification record.

Run:  python scripts/59_build_review2_supplementary.py
"""
import csv
import gzip
import hashlib
import os
import shutil
import sys
import zipfile

sys.stdout.reconfigure(encoding="utf-8")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUP = os.path.join(BASE, "M2-M9深化研究稿", "04_Supplementary")
SCRIPTS = os.path.join(BASE, "scripts")
RESULTS = os.path.join(BASE, "results")
os.makedirs(SUP, exist_ok=True)

CAND22 = ["JAML", "FUBP1", "MAP4K4", "NUMBL", "PARVA", "IREB2", "RPS6KA2", "EID1", "HYKK",
          "STMN3", "GALK2", "RNASET2", "DNAJA4", "SECISBP2L", "RAB31", "CMIP", "HLA-C",
          "RAB4B", "ZNRD1ASP", "CTC-490E21.14", "RP1-167A14.3", "RP11-514O12.4"]
ALIASES = {"JAML": ["JAML", "AMICA1"], "HYKK": ["HYKK", "AGPHD1"]}

# ---------------------------------------------------------------- AF10d ----
def build_af10d():
    path = os.path.join(BASE, "rawdata", "icb", "GSE126044", "GSE126044_counts.txt.gz")
    genes = set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        f.readline()
        for line in f:
            genes.add(line.split("\t", 1)[0].strip())
    out = os.path.join(SUP, "AdditionalFile10d_GSE126044_gene_evaluability.csv")
    with open(out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["# Additional file 10d. Gene evaluability of GSE126044 (NSCLC, anti-PD-1, 16 samples). "
                    "The matrix carries %d gene symbols; GSE126044 uses legacy annotation, in which JAML is AMICA1."
                    % len(genes)])
        w.writerow(["# source: rawdata/icb/GSE126044/GSE126044_counts.txt.gz ; script: scripts/57_check_gse126044_genes.py"])
        w.writerow(["candidate_locus", "symbols_tried", "present_in_GSE126044"])
        n_present = 0
        for g in CAND22:
            tried = ALIASES.get(g, [g])
            hit = [s for s in tried if s in genes]
            present = "yes" if hit else "no"
            if hit:
                n_present += 1
            w.writerow([g, "|".join(tried), present])
        w.writerow([])
        w.writerow(["# summary", "candidate loci present: %d of %d" % (n_present, len(CAND22)), ""])
        w.writerow(["# note", "JAML/AMICA1 is absent, so GSE126044 could not be analysed for the gene of interest", ""])
    print("wrote", out, "| loci present =", n_present, "of", len(CAND22))


# ---------------------------------------------------------------- AF06 -----
NOTE = ("# INTEGRITY NOTE (2026-10-01): the p and fdr columns for split=within_positive_median are degenerate "
        "(identical to 17 significant digits for all 498 TFs) and must not be interpreted; the corresponding test "
        "statistic was not usable in that run. Only the effect sizes (diff) from that contrast are reported in the "
        "manuscript. The p/fdr columns for split=positive_vs_negative are valid.")


def patch_af06():
    p = os.path.join(SUP, "AdditionalFile06_TF_activity_group_differences.csv")
    with open(p, encoding="utf-8-sig") as f:
        lines = f.read().split("\n")
    if any(l.startswith("# INTEGRITY NOTE") for l in lines[:4]):
        print("AF06 already annotated")
        return
    shutil.copy2(p, p + ".bak_review2")
    lines.insert(2, NOTE)
    with open(p, "w", encoding="utf-8-sig", newline="") as f:
        f.write("\n".join(lines))
    print("annotated", p)


# ---------------------------------------------------------------- AF11 -----
def build_af11():
    out = os.path.join(SUP, "AdditionalFile11_code_and_results.zip")
    members = []
    for root, dirs, files in os.walk(SCRIPTS):
        dirs[:] = [d for d in dirs if d != "__pycache__"]
        for fn in sorted(files):
            if fn.endswith(".py"):
                members.append(os.path.join(root, fn))
    derived = [
        ("m2_tcga_jaml_partial_corr.csv", "M2 partial correlations"),
        ("m2_marker_specificity.csv", "M2 adjustment-validity control"),
        ("m2_marker_profile_similarity.csv", "M2 profile similarity"),
        ("m2_tcga_stage_jaml.csv", "M2 stage-wise expression"),
        ("cross_dataset_compartment_pct.csv", "M3/M9 two-dataset compartments"),
        ("reviewer_R1_ID2_relative_effect.csv", "R1 relative effect"),
        ("reviewer_R2_TCGA_parameterisation.csv", "R2 survival parameterisations"),
        ("reviewer_R3_M5_power.json", "R3 power analysis"),
        ("refcheck_m2m9.json", "Crossref DOI verification (machine readable)"),
        ("refcheck_m2m9.md", "Crossref DOI verification (human readable)"),
    ]
    manifest = ["Additional file 11 - analysis code and derived result tables",
                "assembled 2026-10-01", "", "scripts/"]
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for m in members:
            z.write(m, "scripts/" + os.path.basename(m))
            manifest.append("  scripts/" + os.path.basename(m))
        manifest += ["", "results/"]
        for fn, desc in derived:
            src = os.path.join(RESULTS, fn)
            if os.path.exists(src):
                z.write(src, "results/" + fn)
                manifest.append("  results/%s   # %s" % (fn, desc))
            else:
                manifest.append("  (missing) results/%s   # %s" % (fn, desc))
        # the supplementary CSVs as shipped
        manifest += ["", "04_Supplementary/"]
        for fn in sorted(os.listdir(SUP)):
            if fn.endswith(".csv"):
                z.write(os.path.join(SUP, fn), "04_Supplementary/" + fn)
                manifest.append("  04_Supplementary/" + fn)
        z.writestr("MANIFEST.txt", "\n".join(manifest) + "\n")
    size = os.path.getsize(out)
    h = hashlib.sha256(open(out, "rb").read()).hexdigest()[:16]
    print("wrote %s (%d bytes, sha256[:16]=%s, %d members)"
          % (out, size, h, len(members) + len(derived) + 1))


if __name__ == "__main__":
    build_af10d()
    patch_af06()
    build_af11()
