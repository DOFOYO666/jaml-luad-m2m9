# -*- coding: utf-8 -*-
"""Build the Journal of Translational Medicine submission package for the M2-M9 manuscript.

Steps
  1. Normalise three reference entries in both masters to the JTM house author rule
     (verified on 49 element-citations of a 2025 JTM paper: every entry with <etal/>
     lists exactly 3 names) and align the department's English/Chinese name with the
     one already used on the companion submission's title page.
  2. Convert manuscript_EN.md -> JTM投稿_M2M9/manuscript_JTM.md: BMC title page,
     keywords trimmed to 10, "Additional file N" -> "Supplementary Material M",
     section order ... Conclusions -> Supplementary material -> Declarations ->
     References -> Figure legends -> Tables, and Figure 1 / Figure 2 legends relabelled
     for the merged two-panel figures that are actually supplied.
     The same conversion (with Chinese headings) is applied to manuscript_CN.md.
  3. Merge the 15 figure files into the 9 figures the legends describe (PNG + LZW TIFF, 300 dpi).
  4. Copy the 17 supplementary files as SupplementaryMaterial01..17 and write the mapping table.
  5. Write cover letter, README, build log.

Then: python scripts/65_build_jtm_docx.py
"""
import io
import os
import re
import shutil
import sys

sys.stdout.reconfigure(encoding="utf-8")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(BASE, "M2-M9深化研究稿")
EN = os.path.join(PKG, "manuscript_EN.md")
CN = os.path.join(PKG, "manuscript_CN.md")
FIG_SRC = os.path.join(PKG, "03_Figures")
SUP_SRC = os.path.join(PKG, "04_Supplementary")
OUT = os.path.join(BASE, "JTM投稿_M2M9")

AFFIL_EN = ("Department of Pulmonary and Critical Care Medicine, Shenzhen Nanshan District People's "
            "Hospital (Shenzhen University Affiliated Nanshan Hospital), Shenzhen, Guangdong, China")
AFFIL_CN = ("\u6df1\u5733\u5e02\u5357\u5c71\u533a\u4eba\u6c11\u533b\u9662\uff08\u6df1\u5733\u5927\u5b66"
            "\u9644\u5c5e\u5357\u5c71\u533b\u9662\uff09\u547c\u5438\u4e0e\u5371\u91cd\u75c7\u533b\u5b66\u79d1"
            "\uff0c\u5e7f\u4e1c \u6df1\u5733")
CORRESP = ("Huayong Liu, Department of Pulmonary and Critical Care Medicine, Shenzhen Nanshan District "
           "People's Hospital (Shenzhen University Affiliated Nanshan Hospital), No. 89 Taoyuan Road, "
           "Nanshan District, Shenzhen 518000, Guangdong, China. Email: xingxinghuoshu@163.com. "
           "ORCID: 0009-0007-8106-0425. Tel: +86-15989240414.")
AUTHORS_EN = ("Huayong Liu\u00b9\\*, Xuefang Gong\u00b9, Bo Xiong\u00b9, Jiadun Wang\u00b9, "
              "Chunyan Tao\u00b9, Zhuoming Peng\u00b9, Jiming Chen\u00b9")
AUTHORS_CN = ("\u5218\u534e\u52c7\u00b9\\*\uff0c\u9f9a\u96ea\u82b3\u00b9\uff0c\u718a\u6ce2\u00b9\uff0c"
              "\u738b\u5bb6\u76fe\u00b9\uff0c\u9676\u6625\u71d5\u00b9\uff0c\u5f6d\u5353\u660e\u00b9\uff0c"
              "\u9648\u7ee7\u660e\u00b9")

SM_MAP = [
    ("AdditionalFile01_GSE127465_alignment_verification.csv", 1, "GSE127465 alignment verification"),
    ("AdditionalFile02_compartment_positivity_two_datasets.csv", 2, "compartmental positivity, both datasets"),
    ("AdditionalFile03_partial_correlation_all_populations.csv", 3, "partial correlations, all populations"),
    ("AdditionalFile04a_adjustment_control_marker_genes.csv", 4, "marker-gene control for the adjustment"),
    ("AdditionalFile04b_profile_similarity_39populations.csv", 5, "profile similarity, 39 populations"),
    ("AdditionalFile05a_TCGA_survival_all_genes.csv", 6, "TCGA survival, all evaluable genes"),
    ("AdditionalFile05b_TCGA_stage_JAML.csv", 7, "stage-wise JAML expression"),
    ("AdditionalFile06_TF_activity_group_differences.csv", 8, "TF activity group differences"),
    ("AdditionalFile07a_TF_perturbation_three_configurations.csv", 9, "perturbation, three configurations"),
    ("AdditionalFile07b_paired_randomized_control.csv", 10, "paired randomized control"),
    ("AdditionalFile08_cross_module_concordance.csv", 11, "cross-module concordance"),
    ("AdditionalFile09_DepMap_JAML_summary.csv", 12, "DepMap gene-effect summary"),
    ("AdditionalFile10a_ID2_relative_effect.csv", 13, "relative effect of ID2 perturbation"),
    ("AdditionalFile10b_TCGA_survival_parameterisations.csv", 14, "TCGA survival, four parameterizations"),
    ("AdditionalFile10c_M5_power_analysis.csv", 15, "power analysis, checkpoint cohort"),
    ("AdditionalFile10d_GSE126044_gene_evaluability.csv", 16, "GSE126044 gene evaluability"),
    ("AdditionalFile11_code_and_results.zip", 17, "analysis code and derived result tables"),
]

SM_LIST_EN = """## Supplementary material

The following files are supplied as Supplementary Material. Statistical details are given in the corresponding Methods subsections, and every number quoted in the text is reproduced in these files.

**Supplementary Material 1.** Alignment verification for the independent replication dataset (GSE127465): per-cell mitochondrial read fraction versus metadata (rho = 1.0000, n = 54,773), marker-compartment checks, matrix-structure verification, the note on why the `Total counts` field was not used as a criterion, and the independent reimplementation that reproduced 102/102 compartment counts.

**Supplementary Material 2.** Compartmental positivity and mean expression of *JAML* and *CXADR* in both datasets.

**Supplementary Material 3.** Correlation of *JAML* expression with every deconvolved population: raw, adjusted for global immune content, and adjusted for tumour purity.

**Supplementary Material 4.** Marker-gene control for the adjustment procedure, under each candidate adjustment.

**Supplementary Material 5.** Correlation-profile similarity between *JAML* and 15 marker genes across 39 populations.

**Supplementary Material 6.** TCGA-LUAD expression-survival results for all evaluable genes.

**Supplementary Material 7.** Stage-wise *JAML* expression in TCGA-LUAD.

**Supplementary Material 8.** TF activity in *JAML*-positive versus *JAML*-negative CD4+ T cells and within positive cells, with the permutation negative control and the *PDCD1* positive control. The p and FDR columns of the within-positive median split are degenerate and must not be interpreted; only that contrast's effect sizes are used.

**Supplementary Material 9.** Perturbation results for all candidate transcription factors across the three network configurations.

**Supplementary Material 10.** Paired randomized-network control and the two-criterion classification of candidate regulators.

**Supplementary Material 11.** Cross-module concordance between TF activity and perturbation effects.

**Supplementary Material 12.** DepMap 22Q2 gene-effect summary for *JAML*.

**Supplementary Material 13.** Relative effect size of *ID2* perturbation on *JAML* (per-cell distribution, fraction of mean expression, rank among 3,074 genes).

**Supplementary Material 14.** TCGA-LUAD overall survival for all 18 evaluable genes under four parameterizations.

**Supplementary Material 15.** Power analysis for the checkpoint-blockade cohort (minimum detectable hazard ratio).

**Supplementary Material 16.** Gene evaluability of GSE126044: 17 of the 22 candidate loci are present, while *JAML*/*AMICA1* is absent.

**Supplementary Material 17.** Analysis code and derived result tables (archive), including the machine-readable reference-verification record.
"""

SM_LIST_CN = """## \u8865\u5145\u6750\u6599

\u4ee5\u4e0b\u6587\u4ef6\u4f5c\u4e3a\u8865\u5145\u6750\u6599\u968f\u7a3f\u63d0\u4ea4\u3002\u7edf\u8ba1\u7ec6\u8282\u89c1\u76f8\u5e94\u7684 Methods \u5c0f\u8282\uff1b\u6b63\u6587\u5f15\u7528\u7684\u6bcf\u4e00\u4e2a\u6570\u5b57\u5747\u53ef\u5728\u8fd9\u4e9b\u6587\u4ef6\u4e2d\u590d\u73b0\u3002

**\u8865\u5145\u6750\u6599 1.** \u72ec\u7acb\u590d\u6838\u6570\u636e\u96c6\uff08GSE127465\uff09\u7684\u5bf9\u9f50\u6838\u9a8c\uff1a\u9010\u7ec6\u80de\u7ebf\u7c92\u4f53\u8bfb\u6bb5\u6bd4\u4f8b vs \u5143\u6570\u636e\uff08\u03c1 = 1.0000\uff0cn = 54,773\uff09\u3001marker\u2013\u533a\u5ba4\u68c0\u67e5\u3001\u77e9\u9635\u7ed3\u6784\u6838\u9a8c\u3001`Total counts` \u4e0d\u4f5c\u4e3a\u5224\u636e\u7684\u8bf4\u660e\uff0c\u4ee5\u53ca\u590d\u73b0 102/102 \u4e2a\u533a\u5ba4\u8ba1\u6570\u7684\u72ec\u7acb\u5b9e\u73b0\u3002

**\u8865\u5145\u6750\u6599 2.** \u4e24\u5957\u6570\u636e\u4e2d *JAML* \u4e0e *CXADR* \u7684\u533a\u5ba4\u9633\u6027\u7387\u4e0e\u5e73\u5747\u8868\u8fbe\u3002

**\u8865\u5145\u6750\u6599 3.** *JAML* \u8868\u8fbe\u4e0e\u5168\u90e8\u53cd\u5377\u79ef\u7fa4\u4f53\u7684\u76f8\u5173\uff1a\u539f\u59cb\u3001\u6821\u6b63\u6574\u4f53\u514d\u75ab\u542b\u91cf\u3001\u6821\u6b63\u80bf\u7624\u7eaf\u5ea6\u3002

**\u8865\u5145\u6750\u6599 4.** \u6821\u6b63\u6d41\u7a0b\u7684 marker \u57fa\u56e0\u5bf9\u7167\uff0c\u9010\u79cd\u5019\u9009\u6821\u6b63\u91cf\u3002

**\u8865\u5145\u6750\u6599 5.** *JAML* \u4e0e 15 \u4e2a marker \u57fa\u56e0\u5728 39 \u4e2a\u7fa4\u4f53\u4e0a\u7684\u76f8\u5173\u8c31\u76f8\u4f3c\u5ea6\u3002

**\u8865\u5145\u6750\u6599 6.** TCGA-LUAD \u5168\u90e8\u53ef\u8bc4\u4f30\u57fa\u56e0\u7684\u8868\u8fbe\u2013\u751f\u5b58\u7ed3\u679c\u3002

**\u8865\u5145\u6750\u6599 7.** TCGA-LUAD \u4e2d *JAML* \u7684\u5206\u671f\u8868\u8fbe\u3002

**\u8865\u5145\u6750\u6599 8.** *JAML* \u9633\u6027 vs \u9634\u6027 CD4\u207aT \u7ec6\u80de\u53ca\u9633\u6027\u7ec6\u80de\u5185\u90e8\u7684 TF \u6d3b\u6027\uff0c\u542b\u7f6e\u6362\u9634\u6027\u5bf9\u7167\u4e0e *PDCD1* \u9633\u6027\u5bf9\u7167\u3002\u9633\u6027\u7ec6\u80de\u5185\u4e2d\u4f4d\u6570\u5207\u5206\u7684 p \u4e0e FDR \u4e24\u5217\u5df2\u9000\u5316\u3001\u4e0d\u53ef\u89e3\u8bfb\uff1b\u8be5\u5bf9\u6bd4\u4ec5\u91c7\u7528\u5176\u6548\u5e94\u91cf\u3002

**\u8865\u5145\u6750\u6599 9.** \u5168\u90e8\u5019\u9009 TF \u5728\u4e09\u5957\u7f51\u7edc\u914d\u7f6e\u4e2d\u7684\u6270\u52a8\u7ed3\u679c\u3002

**\u8865\u5145\u6750\u6599 10.** \u914d\u5bf9\u6253\u4e71\u7f51\u7edc\u5bf9\u7167\u4e0e\u5019\u9009\u8c03\u63a7\u56e0\u5b50\u7684\u53cc\u5224\u636e\u5206\u7c7b\u3002

**\u8865\u5145\u6750\u6599 11.** TF \u6d3b\u6027\u4e0e\u6270\u52a8\u6548\u5e94\u4e4b\u95f4\u7684\u6a21\u5757\u95f4\u4e00\u81f4\u6027\u3002

**\u8865\u5145\u6750\u6599 12.** DepMap 22Q2 \u4e2d *JAML* \u7684\u57fa\u56e0\u6548\u5e94\u6c47\u603b\u3002

**\u8865\u5145\u6750\u6599 13.** *ID2* \u6270\u52a8\u5bf9 *JAML* \u7684\u76f8\u5bf9\u6548\u5e94\u91cf\uff08\u9010\u7ec6\u80de\u5206\u5e03\u3001\u5360\u5e73\u5747\u8868\u8fbe\u6bd4\u4f8b\u3001\u5728 3,074 \u4e2a\u57fa\u56e0\u4e2d\u7684\u6392\u540d\uff09\u3002

**\u8865\u5145\u6750\u6599 14.** TCGA-LUAD \u5168\u90e8 18 \u4e2a\u53ef\u8bc4\u4f30\u57fa\u56e0\u5728\u56db\u79cd\u53c2\u6570\u5316\u4e0b\u7684\u603b\u751f\u5b58\u7ed3\u679c\u3002

**\u8865\u5145\u6750\u6599 15.** \u68c0\u67e5\u70b9\u6291\u5236\u5242\u961f\u5217\u7684\u529f\u6548\u5206\u6790\uff08\u6700\u5c0f\u53ef\u68c0\u51fa HR\uff09\u3002

**\u8865\u5145\u6750\u6599 16.** GSE126044 \u7684\u57fa\u56e0\u53ef\u8bc4\u4f30\u6027\uff1a22 \u4e2a\u5019\u9009\u4f4d\u70b9\u4e2d 17 \u4e2a\u5728\u4f4d\uff0c*JAML*/*AMICA1* \u7f3a\u5931\u3002

**\u8865\u5145\u6750\u6599 17.** \u5206\u6790\u4ee3\u7801\u4e0e\u6d3e\u751f\u7ed3\u679c\u8868\uff08\u538b\u7f29\u5305\uff09\uff0c\u542b\u673a\u8bfb\u7684\u53c2\u8003\u6587\u732e\u6838\u9a8c\u8bb0\u5f55\u3002
"""

FIG_PLAN = [
    # layout is "v" (stacked) throughout: the source panels are already wide, and placing
    # two of them side by side would push the merged file to 14-19 inches, so that at the
    # journal's column width the axis text would fall below ~4 pt.
    ("Figure1", [("Figure1_compartment_discovery", "a"), ("Figure1b_CD4_subtypes_discovery", "b")], "v", None),
    ("Figure2", [("Figure2_compartment_replication", "a"), ("Figure2b_replication_panel", "b")], "v", None),
    ("Figure3", [("Figure3_TCGA_raw_vs_adjusted", "a")], "v", None),
    ("Figure4", [("Figure4a_adjustment_validity", "a"), ("Figure4b_profile_similarity", "b")], "v", None),
    ("Figure5", [("Figure5_pseudotime_CD4", "a")], "v", None),
    ("Figure6", [("Figure6a_TF_activity", "a"), ("Figure6b_TF_activity_negative_control", "b")], "v", 0.50),
    ("Figure7", [("Figure7a_TF_perturbation", "a"), ("Figure7b_paired_randomized_control", "b")], "v", None),
    ("Figure8", [("Figure8a_activity_vs_perturbation", "a"), ("Figure8b_cross_configuration_reproducibility", "b")], "v", None),
    ("Figure9", [("Figure9_DepMap", "a")], "v", None),
]

REF_FIXES = [
    ("Verdino P, Witherden DA, Havran WL, Wilson IA.",
     "Verdino P, Witherden DA, Havran WL, et al."),
    ("Racle J, de Jonge K, Baumgaertner P, Speiser DE, Gfeller D.",
     "Racle J, de Jonge K, Baumgaertner P, et al."),
    ("Haghverdi L, B\u00fcttner M, Wolf FA, Buettner F, Theis FJ.",
     "Haghverdi L, B\u00fcttner M, Wolf FA, et al."),
]


def log(*a):
    print(*a, flush=True)


def exactly_once(path, old, new, label):
    """Apply old->new exactly once; if it is already applied, say so and move on
    (so the script can be re-run after a later step fails)."""
    t = io.open(path, encoding="utf-8").read()
    if old not in t and new in t:
        log("  ..  already applied: " + label)
        return
    n = t.count(old)
    if n != 1:
        raise SystemExit("[%s] matched %d times (expected 1): %r" % (label, n, old[:70]))
    io.open(path, "w", encoding="utf-8").write(t.replace(old, new, 1))
    log("  OK  " + label)


# ------------------------------------------------------------------ step 1 --
def fix_refs_and_affiliation():
    log("[1] reference house style + affiliation alignment")
    for k, (old, new) in enumerate(REF_FIXES):
        for f in (EN, CN):
            exactly_once(f, old, new, "ref fix %d in %s" % (k + 1, os.path.basename(f)))
    exactly_once(EN,
                 "\u00b9 Department of Respiratory Medicine, Shenzhen Nanshan People's Hospital "
                 "(Shenzhen University Affiliated Nanshan Hospital), Shenzhen, Guangdong, China",
                 "\u00b9 " + AFFIL_EN, "EN affiliation")
    exactly_once(CN,
                 "\u00b9 \u6df1\u5733\u5e02\u5357\u5c71\u533a\u4eba\u6c11\u533b\u9662\uff08\u6df1\u5733\u5927"
                 "\u5b66\u9644\u5c5e\u5357\u5c71\u533b\u9662\uff09\u547c\u5438\u5185\u79d1\uff0c\u5e7f\u4e1c "
                 "\u6df1\u5733",
                 "\u00b9 " + AFFIL_CN, "CN affiliation")


# ----------------------------------------------------------- shared md work --
def _heading(text, name):
    k = text.index(name)
    nxt = text.find("\n## ", k + 1)
    return k, (nxt if nxt > 0 else len(text))


def _reorder(text, tables_h, figs_h, legends_h):
    """Move the Tables and Figures blocks to the very end, after References."""
    kt, nt = _heading(text, tables_h)
    tables = text[kt:nt]
    text = text[:kt] + text[nt:]
    kf, nf = _heading(text, figs_h)
    figs = text[kf:nf].replace(figs_h, legends_h, 1)
    text = text[:kf] + text[nf:]
    return text.rstrip() + "\n\n---\n\n" + figs.strip() + "\n\n---\n\n" + tables.strip() + "\n"


def soft_replace(text, old, new, expected=1, label=""):
    """Replace old->new; return (text, applied). Skip silently if already applied."""
    if old not in text and new in text:
        return text, False
    n = text.count(old)
    if n != expected:
        raise SystemExit("[%s] matched %d times (expected %d): %r" % (label, n, expected, old[:70]))
    return text.replace(old, new), True


AVAIL_OLD_EN = ("are provided as Additional file 11; that archive will be deposited in a public "
                "repository with a DOI before publication and remains available from the corresponding "
                "author on reasonable request in the interim.")
AVAIL_NEW_EN = ("are provided as Supplementary Material 17. An identical archive (release 1.0.0) is being "
                "deposited in a public repository; its DOI will be quoted here once issued, and the archive "
                "remains available from the corresponding author on reasonable request in the interim.")
AVAIL_OLD_CN = ("\u5747\u4f5c\u4e3a\u9644\u52a0\u6587\u4ef6 11 \u63d0\u4f9b\uff1b\u540c\u4e00\u538b\u7f29"
                "\u5305\u5c06\u5728\u6295\u7a3f\u524d\u5b58\u5165\u5e26 DOI \u7684\u516c\u5171\u4ed3\u5e93")
AVAIL_NEW_CN = ("\u5747\u4f5c\u4e3a\u8865\u5145\u6750\u6599 17 \u63d0\u4f9b\u3002\u5185\u5bb9\u76f8\u540c"
                "\u7684\u5f52\u6863\uff08release 1.0.0\uff09\u6b63\u5b58\u5165\u516c\u5171\u4ed3\u5e93")


def normalize_master_en():
    """Bring the English master to the same conventions as the Chinese one, in place."""
    log("[2a-0] normalise manuscript_EN.md (supplementary naming, order, availability)")
    t = io.open(EN, encoding="utf-8").read()
    if "## Supplementary material" in t:
        log("  ..  already normalised")
        return
    t, _ = soft_replace(t, "; DepMap; target triage", "; DepMap", 1, "EN keywords")
    t, _ = soft_replace(t, "Additional file 1)", "Supplementary Material 1)", 3, "EN citations")
    t, _ = soft_replace(t, AVAIL_OLD_EN, AVAIL_NEW_EN, 1, "EN availability")
    i = t.index("## Additional files")
    j = t.index("## Declarations")
    t = t[:i] + SM_LIST_EN + "\n---\n\n" + t[j:]
    t = _reorder(t, "## Tables", "## Figures", "## Figure legends")
    io.open(EN, "w", encoding="utf-8").write(t)
    log("  EN master normalised; headings: " + " | ".join(re.findall(r"(?m)^## (.+)$", t)))


CN_AVAIL_FINAL = ("\u5747\u4f5c\u4e3a\u8865\u5145\u6750\u6599 17 \u63d0\u4f9b\u3002\u5185\u5bb9\u76f8\u540c\u7684"
                  "\u5f52\u6863\uff08release 1.0.0\uff09\u6b63\u5b58\u5165\u516c\u5171\u4ed3\u5e93\uff1b\u5176 DOI "
                  "\u4e00\u7ecf\u53d6\u5f97\u5373\u5728\u6b64\u5904\u5f15\u7528\uff0c\u5728\u6b64\u4e4b\u524d\u53ef\u5411"
                  "\u901a\u8baf\u4f5c\u8005\u7d22\u53d6\u3002")
CN_AVAIL_CANDIDATES = [
    "\u5747\u4f5c\u4e3a\u9644\u52a0\u6587\u4ef6 11 \u63d0\u4f9b\uff1b\u540c\u4e00\u538b\u7f29\u5305\u5c06\u5728"
    "\u6295\u7a3f\u524d\u5b58\u5165\u5e26 DOI \u7684\u516c\u5171\u4ed3\u5e93\uff0c\u5728\u6b64\u4e4b\u524d\u53ef"
    "\u5411\u901a\u8baf\u4f5c\u8005\u7d22\u53d6\u3002",
    "\u5747\u4f5c\u4e3a\u8865\u5145\u6750\u6599 17 \u63d0\u4f9b\uff1b\u540c\u4e00\u538b\u7f29\u5305\u5c06\u5728"
    "\u6295\u7a3f\u524d\u5b58\u5165\u5e26 DOI \u7684\u516c\u5171\u4ed3\u5e93\uff0c\u5728\u6b64\u4e4b\u524d\u53ef"
    "\u5411\u901a\u8baf\u4f5c\u8005\u7d22\u53d6\u3002",
]


def normalize_master_cn(t):
    """Return the CN master text with the final availability wording."""
    # 若文中已出现真实 DOI（`67 --doi` 的成果），说明可用性声明已定稿，不要再改措辞 ——
    # 否则 `--writeback` 跑到这一步会因为"认不出新措辞"而中止（实测踩过）。
    if CN_AVAIL_FINAL in t or "10.5281/zenodo." in t:
        return t, False
    for old in CN_AVAIL_CANDIDATES:
        if t.count(old) == 1:
            return t.replace(old, CN_AVAIL_FINAL, 1), True
    raise SystemExit("CN availability statement: no known source wording matched")


def build_en_md():
    log("[2a] manuscript_JTM.md (EN)")
    t = io.open(EN, encoding="utf-8").read()
    head, body = t.split("\n---\n", 1)
    title = head.split("\n", 1)[0].strip()
    running = re.search(r"\*\*Running title:\*\*\s*(.+)", head).group(1).strip()
    t = (title + "\n\n" + AUTHORS_EN + "\n\n" + "\u00b9 " + AFFIL_EN + "\n\n"
         + "**Corresponding author:** " + CORRESP + "\n\n"
         + "**Running title:** " + running + "\n\n---\n" + body)

    # Figure legends must describe the merged panels that are actually supplied, and the
    # numbers printed in Figure 7a must be attributed to the run they came from.  Verified
    # against the source tables: Figure 1 = m3_jaml_cxadr_normalized + m3_jaml_cd4_subtypes;
    # Figure 2 = cross_dataset_compartment + m8_gse127465_jaml_cxadr; Figure 7a =
    # m4_tf_perturbation_summary.csv (3,013 genes, bagging = 20, 34 perturbable factors,
    # randomized-control max 1.18e-3) while Table 4 / panel (b) use the final configuration
    # and the paired control (global ceiling 9.59e-3).
    ops = [
        # Figure 1
        ("Mean normalized expression (log2(TPM+1)) and positivity for *JAML* and *CXADR* across ten author-annotated cell types in GSE131907 (208,506 cells; LUAD). Points are cell-type means; the x-axis is the fraction of positive cells.",
         "(a) Mean normalized expression (log2(TPM+1)) and positivity for *JAML* and *CXADR* across ten author-annotated cell types in GSE131907 (208,506 cells; LUAD); points are cell-type means and the x-axis is the fraction of positive cells."),
        ("Panel b: *JAML* positivity across CD4-lineage subsets annotated by the authors (`Cell_subtype`).",
         "(b) *JAML* mean expression (bars), annotated with the fraction of positive cells, across the CD4-lineage subsets defined by the authors' `Cell_subtype` labels."),
        # Figure 2 — panel (a) is handled below by a regex, because its text quotes a
        # citation number that changes whenever the reference list is renumbered.
        ("Horizontal bars are positivity percentages; the discovery and replication values for the same compartment are shown side by side.",
         "(b) The same replication data at the authors' own annotation granularity, separately resolving the key epithelial and myeloid compartments; the note 'few epithelial cells in this dataset' flags the small epithelial compartment discussed in the Limitations. Horizontal bars are positivity percentages."),
        # Figure 7a — the master legend now describes the final network configuration, so the only
        # JTM-specific change is the name of the supplementary table it points to (the master calls
        # it "Additional file 7a", JTM calls it "Supplementary Material 9").
        ("is in Additional file 7a.", "is in Supplementary Material 9."),
    ]
    for old, new in ops:
        t, applied = soft_replace(t, old, new, 1, "legend op")

    # Figure 2, panel (a). Regex rather than a literal string: the op used to hard-code
    # "Zilionis et al. [8]", which silently stopped matching once the reference list was
    # renumbered (the script aborted and the JTM file was left stale). Match the number-free
    # part and carry whatever citation number the master happens to use.
    fig2 = re.compile(
        r"Positivity of \*JAML\* and \*CXADR\* across eight compartments in GSE127465 "
        r"\(54,773 cells; treatment-naïve NSCLC; (Zilionis et al\. \[\d+(?:[–\-]\d+)*\])\) "
        r"compared with GSE131907\.")
    if "(a) Positivity of *JAML* and *CXADR* across the eight compartments" in t:
        log("  ..  already applied: Figure 2 legend (a)")
    else:
        t, n_fig2 = fig2.subn(
            r"(a) Positivity of *JAML* and *CXADR* across the eight compartments, with the discovery "
            r"(GSE131907) and replication (GSE127465; 54,773 cells; treatment-naïve NSCLC; \1) values "
            r"for the same compartment shown side by side.",
            t, count=1)
        if n_fig2 != 1:
            raise SystemExit("[Figure 2 legend] regex matched %d times (expected 1)" % n_fig2)
        log("  OK  Figure 2 legend (a) split out")

    assert "Additional file" not in t, "EN -> JTM still contains 'Additional file'"
    log("  headings: " + " | ".join(re.findall(r"(?m)^## (.+)$", t)))

    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, "manuscript_JTM.md")
    io.open(p, "w", encoding="utf-8").write(t)
    log("  wrote %s (%d chars)" % (os.path.basename(p), len(t)))


def build_cn_md():
    log("[2b] manuscript_CN.md companion (Chinese headings)")
    t = io.open(CN, encoding="utf-8").read()
    if "## \u8865\u5145\u6750\u6599" in t:
        log("  ..  structure already converted")
    else:
        t = t.replace("\uff1b\u9776\u70b9\u7b5b\u9009", "", 1)          # ；靶点筛选
        n = t.count("\u9644\u52a0\u6587\u4ef6 1\uff09")            # 附加文件 1）
        assert n == 3, "CN in-text citations: found %d (expected 3)" % n
        t = t.replace("\u9644\u52a0\u6587\u4ef6 1\uff09", "\u8865\u5145\u6750\u6599 1\uff09")
        i = t.index("## \u9644\u52a0\u6587\u4ef6")
        j = t.index("## \u58f0\u660e")
        t = t[:i] + SM_LIST_CN + "\n---\n\n" + t[j:]
        t = _reorder(t, "## \u8868\u683c", "## \u56fe", "## \u56fe\u6ce8")

    t, applied = normalize_master_cn(t)
    io.open(CN, "w", encoding="utf-8").write(t)
    log("  availability wording %s" % ("updated" if applied else "already final"))
    log("  headings: " + " | ".join(re.findall(r"(?m)^## (.+)$", t)))


def build_cn_md_unused():
    pass


# ------------------------------------------------------------------ step 3 --
def merge_figures():
    log("[3] figures -> 9 merged files")
    from PIL import Image, ImageDraw, ImageFont
    rs = getattr(Image, "Resampling", Image)
    dst = os.path.join(OUT, "02_Figures")
    os.makedirs(dst, exist_ok=True)

    for name, panels, layout, scale_b in FIG_PLAN:
        ims = [Image.open(os.path.join(FIG_SRC, s + ".png")).convert("RGB") for s, _ in panels]
        marks = []                                   # (x, y, letter)
        if len(ims) == 1:
            canvas = ims[0]
        elif layout == "h":
            h = max(im.height for im in ims)
            sc = [im.resize((round(im.width * h / im.height), h), rs.LANCZOS) for im in ims]
            gap = 40
            W = sum(im.width for im in sc) + gap * (len(sc) - 1)
            canvas = Image.new("RGB", (W, h), "white")
            x = 0
            for im, (_, letter) in zip(sc, panels):
                canvas.paste(im, (x, 0))
                marks.append((x + 14, 10, letter))
                x += im.width + gap
        else:
            W = max(im.width for im in ims)
            sc = []
            for k, im in enumerate(ims):
                w = round(W * scale_b) if (k == 1 and scale_b) else W
                sc.append(im.resize((w, round(im.height * w / im.width)), rs.LANCZOS))
            gap = 40
            H = sum(im.height for im in sc) + gap * (len(sc) - 1)
            canvas = Image.new("RGB", (W, H), "white")
            y = 0
            for im, (_, letter) in zip(sc, panels):
                canvas.paste(im, ((W - im.width) // 2, y))
                marks.append((14, y + 8, letter))
                y += im.height + gap

        if len(panels) > 1:
            fsize = max(36, round(canvas.width * 0.018))
            font = None
            for cand in (r"C:\Windows\Fonts\arialbd.ttf", r"C:\Windows\Fonts\arial.ttf"):
                if os.path.exists(cand):
                    font = ImageFont.truetype(cand, fsize)
                    break
            if font is not None:
                dr = ImageDraw.Draw(canvas)
                for x, y, letter in marks:
                    bb = dr.textbbox((x, y), letter, font=font)
                    dr.rectangle([bb[0] - 6, bb[1] - 4, bb[2] + 6, bb[3] + 4], fill="white")
                    dr.text((x, y), letter, fill="black", font=font)

        canvas.save(os.path.join(dst, name + ".png"), dpi=(300, 300))
        canvas.save(os.path.join(dst, name + ".tiff"), dpi=(300, 300),
                    **{"compression": "tiff_lzw"})
        log("  %s  <- %s   %dx%d" % (name + ".png", " + ".join(s for s, _ in panels),
                                    canvas.width, canvas.height))


# ------------------------------------------------------------------ step 4 --
def copy_supplementary():
    log("[4] supplementary files -> Supplementary Material 1..17")
    dst = os.path.join(OUT, "03_Supplementary")
    os.makedirs(dst, exist_ok=True)
    rows = []
    for src_name, num, desc in SM_MAP:
        src = os.path.join(SUP_SRC, src_name)
        if not os.path.exists(src):
            raise SystemExit("missing supplementary file: " + src_name)
        stem, ext = os.path.splitext(src_name)
        stem = re.sub(r"^AdditionalFile\d+[a-d]?_", "", stem)      # drop the legacy prefix
        new = "SupplementaryMaterial%02d_%s%s" % (num, stem, ext)
        shutil.copy2(src, os.path.join(dst, new))
        rows.append((num, new, src_name, desc, os.path.getsize(src)))
        log("  %2d  %s" % (num, new))
    return rows


# ------------------------------------------------------------------ step 5 --
COVER = """To the Editors
*Journal of Translational Medicine*

Dear Editors,

We submit for your consideration our manuscript, "{title}".

**The question.** Germline variation at *JAML* (*AMICA1*) has been linked to lung adenocarcinoma risk through CD4+ T-cell expression, and the molecule has since been proposed as an immunotherapeutic target. Before such a target is pursued, three things need settling: where the molecule is expressed in tumour tissue, which cells determine its tissue abundance, and whether its expression is under identifiable transcriptional control. We set out to answer those three questions and to state plainly what the resulting evidence does and does not support.

**What we found.** Using public data only, we interrogated this single gene with eight complementary analyses, each paired with an explicit negative or positive control. Two findings carry the paper. First, *JAML* and its ligand *CXADR* occupy different compartments, and the separation replicates in an independently processed NSCLC cohort (*JAML* positivity 47.0% versus 1.7% in myeloid versus epithelium in the discovery set, 30.1% versus 1.5% in replication, epithelial-to-myeloid ratios 0.04 and 0.05), which reframes the therapeutic question as one about the myeloid-epithelial interface rather than about the tumour cell. Second, of 39 perturbed transcription factors only *ID2* survives a paired randomized-network control together with a global noise ceiling, and it does so in all three network configurations. *JAML* is not a dependency of LUAD cell lines, and we report that and two further null results as null results.

**Why we think it fits this journal.** The paper is a target-triage study: it does not nominate another gene, it narrows where a candidate should be tested. It is also, deliberately, a study about controls. Two of our control designs failed at first pass, for reasons we now quantify — a randomized coefficient matrix that the software caches, and a ratio whose denominator collapses for a low-detection factor — as did a composite adjustment variable that erased the very lineage signal it was meant to preserve. Each failure was caught only by inspecting control values rather than test statistics. We describe these in Methods and Discussion, and we flag the degrees of freedom a reader should be sceptical about, including that our paired control rests on five randomized replicates, so that no empirical P value below 0.17 is attainable and the quantity of interest is the size of the excess rather than a significance test.

**What we do not claim.** *ID2* is a candidate, not a mechanism: the supporting edge is of moderate strength (6.3% of mean *JAML* expression; rank 33 of 3,074 genes by relative effect). The 229 transcription factors flagged by activity inference are not regulators of *JAML*. *JAML* is not a T-cell infiltration marker. No therapeutic claim is made: no protein-level measurement was performed and the interface has not been tested physically. The experimental follow-up we consider warranted — perturbing *ID2* in primary human CD4+ T cells or monocyte-derived cells, and blocking the *JAML*-*CXADR* interface in co-culture — is set out in the Discussion.

All data are public, and all analysis code together with the derived result tables underlying every figure and table is supplied as Supplementary Material 17. No new human participant data were collected. The manuscript is original, is not under consideration elsewhere, and all authors have approved its submission.

We hope the work suits *Journal of Translational Medicine* and look forward to your assessment.

Yours sincerely,

Huayong Liu, on behalf of all authors
Department of Pulmonary and Critical Care Medicine
Shenzhen Nanshan District People's Hospital (Shenzhen University Affiliated Nanshan Hospital)
No. 89 Taoyuan Road, Nanshan District, Shenzhen 518000, Guangdong, China
Email: xingxinghuoshu@163.com
"""


def write_docs(sm_rows):
    log("[5] cover letter / README / build log")
    title = io.open(os.path.join(OUT, "manuscript_JTM.md"), encoding="utf-8").read().split("\n", 1)[0]
    title = title.replace("# ", "").strip()
    io.open(os.path.join(OUT, "Cover_Letter_JTM.md"), "w", encoding="utf-8").write(COVER.format(title=title))

    tbl = "\n".join("| %d | `%s` | `%s` | %s | %s |" % (n, new, old, d, "{:,}".format(sz))
                    for n, new, old, d, sz in sm_rows)
    readme = README_TMPL.format(title=title, tbl=tbl)
    io.open(os.path.join(OUT, "README_JTM_投稿说明.md"), "w", encoding="utf-8").write(readme)

    lines = ["# JTM 投稿包构建记录", "",
             "- 来源：`M2-M9深化研究稿/manuscript_EN.md`（EN 权威源）；CN 母稿同步改名与重排",
             "- 脚本：`scripts/64_build_jtm_pkg.py`（本包）、`scripts/65_build_jtm_docx.py`（docx）",
             "- 格式依据：`scripts/61_jtm_style_probe.py`、`scripts/62_jtm_supp_probe.py`（12 篇 JTM 全文 XML）",
             "- 图：9 张（PNG + 300 dpi LZW TIFF）；补充材料：17 个", ""]
    for n, new, old, d, sz in sm_rows:
        lines.append("- [补充] SupplementaryMaterial%02d  <-  %s  (%.1f KB)" % (n, old, sz / 1024))
    for name, panels, layout, sb in FIG_PLAN:
        lines.append("- [图] %s  <-  %s  (%s)" % (name, " + ".join(p[0] for p in panels),
                                                 "side-by-side" if layout == "h" else "stacked"))
    io.open(os.path.join(OUT, "BUILD_LOG_JTM.md"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
    log("  cover letter / README / BUILD_LOG written")


README_TMPL = """# 《M2–M9 深化研究》Journal of Translational Medicine 投稿包

**稿件：** {title}
**目标期刊：** Journal of Translational Medicine（BMC / Springer Nature）
**编制日期：** 2026-10-01
**正文来源：** `../M2-M9深化研究稿/manuscript_EN.md`（由 `scripts/64_build_jtm_pkg.py` 转换，科学内容未改动）

---

## 一、本包按 JTM/BMC 要求做了什么（规则均有来源，非凭记忆）

格式依据：抓取 **12 篇 2024–2025 年 JTM 开放获取论文的全文 XML** 逐条核对（`scripts/61`、`scripts/62`）。

| 项目 | 核到的事实 | 本包的处理 |
|---|---|---|
| 摘要 | 结构式四标题 Background / Methods / Results / Conclusions | 一致（329 词，≤ 350 上限） |
| 关键词 | 3–10 个 | 11 → 10（删 "target triage"） |
| 正文章节 | Introduction → Methods → Results → Discussion → Conclusions | 一致；Limitations 作为独立小节保留 |
| 参考文献作者 | **> 3 位时列前 3 位 + et al.**（1 篇 JTM 论文 49 条 `element-citation` 中，39 条带 `<etal/>` 者**全部**为 3 名） | 全部条目均按此规则：> 3 位作者列前 3 位 + et al.，作者数 ≤ 3 时完整列出 |
| 参考文献编号 | Vancouver，按正文首次出现顺序 | 本次修订新增 26 条文献，编号已整体重排；正文引用与文献表**双向对账**通过（无未引用条目、无悬空引用），并逐条经 Crossref DOI 反查 |
| 补充材料命名 | 一律 **"Supplementary Material N"**（12 篇抽样中 "Additional file" 出现 **0** 次） | "Additional file N" → "Supplementary Material 1–17"，文件按 `SupplementaryMaterialNN_*.csv` 重命名，正文引用同步 |
| 后置章节顺序 | Conclusions → Electronic supplementary material → Acknowledgements → Author contributions → Funding → Data availability → Declarations → Abbreviations → References | Declarations 四小标题采用 JTM 原文措辞；References 置于声明之后 |
| 图 | 图件单独上传、图注随正文 | 15 个图件**合并为图注所述的 9 张图**；图注集中于 References 之后的 "Figure legends" |
| 表 | 可编辑格式，非图片 | 5 张均为真 Word 三线表 |
| 报告规范 | 观察性 + 计算分析，不强制 MR 清单 | Methods §2.14 已声明不报告 MR 估计，故不适用 STROBE-MR |

---

## 二、目录结构

```
JTM投稿_M2M9/
├── manuscript_JTM.md        权威排版源（md）
├── 01_Manuscript_JTM.docx   Word 稿（题名页 → 摘要 → 正文 → 声明 → 参考文献 → 图注 → 表）
├── Cover_Letter_JTM.md      投稿信（To the Editors，不具名）
├── 02_Figures/              Figure1–Figure9（.png + 300 dpi LZW .tiff）
├── 03_Supplementary/        SupplementaryMaterial01–17
├── README_JTM_投稿说明.md    本文件
└── BUILD_LOG_JTM.md         构建记录
```

---

## 三、图的合并方式（正文图号 → 源文件）

| 正文图 | 源文件 | 排布 |
|---|---|---|
| Figure 1 | Figure1_compartment_discovery + Figure1b_CD4_subtypes_discovery | 上下堆叠（a, b） |
| Figure 2 | Figure2_compartment_replication + Figure2b_replication_panel | 上下堆叠（a, b） |
| Figure 3 | Figure3_TCGA_raw_vs_adjusted | 单幅 |
| Figure 4 | Figure4a_adjustment_validity + Figure4b_profile_similarity | 上下堆叠（a, b） |
| Figure 5 | Figure5_pseudotime_CD4 | 单幅（内含 a/b/c） |
| Figure 6 | Figure6a_TF_activity + Figure6b_TF_activity_negative_control | 上下堆叠（a, b；b 缩至 50% 宽） |
| Figure 7 | Figure7a_TF_perturbation + Figure7b_paired_randomized_control | 上下堆叠（a, b） |
| Figure 8 | Figure8a_activity_vs_perturbation + Figure8b_cross_configuration_reproducibility | 上下堆叠（a, b） |
| Figure 9 | Figure9_DepMap | 单幅 |

> **为什么统一上下堆叠**：源面板本身已偏宽（例如 Figure 8b 单张即 10.8 英寸宽）。若左右并排，
> 合并图会达到 14–19 英寸；按期刊栏宽缩排后坐标轴文字将降到约 4 pt。上下堆叠后每张图宽度控制在
> 6.8–12.5 英寸，缩放后文字约 5–9 pt。这一点已用像素数与字号反算核对（脚本内注释）。

Figure 1 与 Figure 2 的图注已改写为 (a)/(b) 结构，与合并后的版式一致。

---

## 四、补充材料编号对照（旧 → 新）

| 新编号 | 文件 | 原编号与文件名 | 内容 | 字节 |
|---|---|---|---|---|
{tbl}

---

## 五、投稿前仍需确认

1. **通讯作者信息**按主稿（BMC Cancer 投稿包）填写：`xingxinghuoshu@163.com`、ORCID `0009-0007-8106-0425`、Tel `+86-15989240414`、`No. 89 Taoyuan Road, Nanshan District, Shenzhen 518000`。如需修改，改 `manuscript_JTM.md` 题名页后重跑 `scripts/65`。
2. **科室英文名已统一并经作者确认**（2026-10-01）：全文一律用 **"Department of Pulmonary and Critical Care Medicine, Shenzhen Nanshan District People's Hospital (Shenzhen University Affiliated Nanshan Hospital)"**，中文稿为"深圳市南山区人民医院（深圳大学附属南山医院）呼吸与危重症医学科"。原 M2–M9 稿的 "Department of Respiratory Medicine, Shenzhen Nanshan People's Hospital" 与"呼吸内科"已全部替换（含本包 Cover letter、代码发布树的 `CITATION.cff` 与 `.zenodo.json`）。
3. **Cover letter 不写日期**（投稿系统自带）；如需可在首行补。

---

## 六、代码发布物与 DOI（投稿系统会校验这一项）

代码归档已备好两步之外的**全部内容**：

```
JAML_M2M9_code_release/            发布树（README / DATA_SOURCES / LICENSE(MIT，可改) /
                                   CITATION.cff / .zenodo.json / SHA256SUMS.txt）
   ├── code/                       82 个 C 盘工作区脚本
   ├── code/d_workspace/           20 个 D 盘 M4/M8 脚本（含 m4_fig7_replot.py 等；以 D 盘版本为准）
   ├── results/                    10 个派生结果表 + Crossref 核验记录
   └── supplementary/              Supplementary Material 1–17
JAML_M2M9_code_release.tar.gz      约 1.1 MB。**校验值每次重建都会变，请以发布树内的
                                   `SHA256SUMS.txt` 与重建脚本打印的 sha256 为准，勿从本文件转抄。**
```

> 建仓、取 DOI、许可确认的**逐步操作**见 `M2-M9深化研究稿/代码发布与DOI操作手册.md`
> 与 `M2-M9深化研究稿/MIT许可证确认流程.md`。

**剩余两步需要你的账号**（本机无 GitHub 凭据、无 `gh` CLI，代码不作伪 DOI）：

```
cd "…/JAML_M2M9_code_release"
git remote add origin https://github.com/<你的账号>/<仓库名>.git
git branch -M main && git push -u origin main
# Zenodo：用 GitHub 登录 → Settings → GitHub → 打开该仓库开关；
#         回 GitHub 建 Release（tag v1.0.0）→ Zenodo 自动颁发 DOI。
python scripts/67_finalise_review2.py --doi 10.5281/zenodo.XXXXXXX --url https://github.com/<你>/<仓库>
```

第 3 条命令会把仓库地址与 DOI 写回英文母稿、中文母稿的可用性声明，随后重跑
`scripts/49`（中英母稿 docx）、`scripts/64`（本包）、`scripts/65`（本包 docx）即全部同步。
在此之前，稿件中的表述为"内容相同的归档（release 1.0.0）正存入公共仓库；其 DOI 一经取得即在此处引用"
——**如实描述当前状态，不预留假 DOI**。
"""


def main():
    fix_refs_and_affiliation()
    normalize_master_en()
    build_en_md()
    build_cn_md()
    merge_figures()
    rows = copy_supplementary()
    write_docs(rows)
    print("\nDONE -> " + OUT)


if __name__ == "__main__":
    main()
