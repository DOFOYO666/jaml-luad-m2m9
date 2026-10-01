# -*- coding: utf-8 -*-
"""Revise the M2-M9 English master: richer Introduction/Discussion, de-AI prose,
native polish, tighter Methods/Results, and a renumbered reference list.

Pipeline
  [0] sanity-check the section markers
  [1] convert the existing numeric citations       [7] / [1-5]  ->  {{key}}
  [2] splice the rewritten Abstract / Introduction / Methods / Results / Discussion
  [3] assign new numbers by order of first appearance across the whole document
  [4] rebuild the reference list in that order
  [5] write, then run structural checks (count, continuity, both-way citation audit)

Every rewrite partner file lives in  _rewrite/  and is kept so the change is auditable.
"""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EN = os.path.join(ROOT, "M2-M9深化研究稿", "manuscript_EN.md")
NEW = os.path.join(ROOT, "_rewrite")
DATE = "2026-10-01"

# ---------------------------------------------------------------- reference data
# old number -> key (the 23 references already in the manuscript)
OLD = {
    1: "weber", 2: "verdino", 3: "mcgraw", 4: "feng", 5: "yazar", 6: "kerimov",
    7: "kim", 8: "zilionis", 9: "becht", 10: "racle", 11: "newman", 12: "carter",
    13: "dimitrov", 14: "decoupler", 15: "collecTRI", 16: "celloracle", 17: "tsherniak",
    18: "haghverdi", 19: "scanpy", 20: "bh", 21: "vosa", 22: "jung", 23: "strobe",
}

REFS = {
    # ---- already cited in the previous version (wording unchanged) ----
    "weber": "Weber C, Fraemohs L, Dejana E. The role of junctional adhesion molecules in vascular inflammation. *Nat Rev Immunol.* 2007;7(6):467–477. doi:10.1038/nri2096",
    "verdino": "Verdino P, Witherden DA, Havran WL, et al. The molecular interaction of CAR and JAML recruits the central cell signal transducer PI3K. *Science.* 2010;329(5996):1210–1214. doi:10.1126/science.1187996",
    "mcgraw": "McGraw JM, Witherden DA, Havran WL, et al. JAML promotes CD8 and γδ T cell antitumor immunity and is a novel target for cancer immunotherapy. *J Exp Med.* 2021;218(10):e20202644. doi:10.1084/jem.20202644",
    "feng": "Feng Z, Zhang Y, He M, et al. AMICA1 is a diagnostic and prognostic biomarker and induces immune cells infiltration by activating cGAS-STING signaling in lung adenocarcinoma. *Cancer Cell Int.* 2022;22(1):111. doi:10.1186/s12935-022-02517-x",
    "yazar": "Yazar S, Alquicira-Hernandez J, Wing K, et al. Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease. *Science.* 2022;376(6589):eabf3041. doi:10.1126/science.abf3041",
    "kerimov": "Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human gene expression and splicing quantitative trait loci. *Nat Genet.* 2021;53(9):1290–1299. doi:10.1038/s41588-021-00924-w",
    "kim": "Kim N, Kim HK, Lee K, et al. Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. *Nat Commun.* 2020;11(1):2285. doi:10.1038/s41467-020-16164-1",
    "zilionis": "Zilionis R, Engblom C, Pfirschke C, et al. Single-cell transcriptomics of human and mouse lung cancers reveals conserved myeloid populations across individuals and species. *Immunity.* 2019;50(5):1317–1334.e10. doi:10.1016/j.immuni.2019.03.009",
    "becht": "Becht E, Giraldo NA, Lacroix L, et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. *Genome Biol.* 2016;17(1):218. doi:10.1186/s13059-016-1070-5",
    "racle": "Racle J, de Jonge K, Baumgaertner P, et al. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. *eLife.* 2017;6:e26476. doi:10.7554/eLife.26476",
    "newman": "Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell subsets from tissue expression profiles. *Nat Methods.* 2015;12(5):453–457. doi:10.1038/nmeth.3337",
    "carter": "Carter SL, Cibulskis K, Helman E, et al. Absolute quantification of somatic DNA alterations in human cancer. *Nat Biotechnol.* 2012;30(5):413–421. doi:10.1038/nbt.2203",
    "dimitrov": "Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat Commun.* 2022;13(1):3224. doi:10.1038/s41467-022-30755-0",
    "decoupler": "Badia-i-Mompel P, Vélez Santiago J, Braunger J, et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinform Adv.* 2022;2(1):vbac016. doi:10.1093/bioadv/vbac016",
    "collecTRI": "Müller-Dott S, Tsirvouli E, Vazquez M, et al. Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities. *Nucleic Acids Res.* 2023;51(20):10934–10949. doi:10.1093/nar/gkad841",
    "celloracle": "Kamimoto K, Stringa B, Hoffmann CM, et al. Dissecting cell identity via network inference and in silico gene perturbation. *Nature.* 2023;614(7949):742–751. doi:10.1038/s41586-022-05688-9",
    "tsherniak": "Tsherniak A, Vazquez F, Montgomery PG, et al. Defining a cancer dependency map. *Cell.* 2017;170(3):564–576.e16. doi:10.1016/j.cell.2017.06.010",
    "haghverdi": "Haghverdi L, Büttner M, Wolf FA, et al. Diffusion pseudotime robustly reconstructs lineage branching. *Nat Methods.* 2016;13(10):845–848. doi:10.1038/nmeth.3971",
    "scanpy": "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* 2018;19(1):15. doi:10.1186/s13059-017-1382-0",
    "bh": "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B.* 1995;57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x",
    "vosa": "Võsa U, Claringbould A, Westra HJ, et al. Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression. *Nat Genet.* 2021;53(9):1300–1310. doi:10.1038/s41588-021-00913-z",
    "jung": "Jung H, Kim HS, Kim JY, et al. DNA methylation loss promotes immune evasion of tumours with high mutation and copy number load. *Nat Commun.* 2019;10(1):4278. doi:10.1038/s41467-019-12159-9",
    "strobe": "Skrivankova VW, Richmond RC, Woolf BAR, et al. Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization: the STROBE-MR statement. *JAMA.* 2021;326(16):1614–1621. doi:10.1001/jama.2021.18236",

    # ---- added in this revision; every entry verified against Crossref by DOI on 2026-10-01 ----
    "mooglutz": "Moog-Lutz C, Cavé-Riant F, Guibal FC, et al. JAML, a novel protein with characteristics of a junctional adhesion molecule, is induced during differentiation of myeloid leukemia cells. *Blood.* 2003;102(9):3371–3378. doi:10.1182/blood-2002-11-3462",
    "witherden": "Witherden DA, Verdino P, Rieder SE, et al. The junctional adhesion molecule JAML is a costimulatory receptor for epithelial γδ T cell activation. *Science.* 2010;329(5996):1205–1210. doi:10.1126/science.1192698",
    "kummer": "Kummer D, Ebnet K. Junctional adhesion molecules (JAMs): the JAM–integrin connection. *Cells.* 2018;7(4):25. doi:10.3390/cells7040025",
    "luissint": "Luissint AC, Lutz PG, Calderwood DA, et al. JAM-L-mediated leukocyte adhesion to endothelial cells is regulated in cis by α4β1 integrin activation. *J Cell Biol.* 2008;183(6):1159–1173. doi:10.1083/jcb.200805061",
    "silvasantos": "Silva-Santos B, Mensurado S, Coffelt SB. γδ T cells: pleiotropic immune effectors with therapeutic potential in cancer. *Nat Rev Cancer.* 2019;19(7):392–404. doi:10.1038/s41568-019-0153-5",
    "bergelson": "Bergelson JM, Cunningham JA, Droguett G, et al. Isolation of a common receptor for coxsackie B viruses and adenoviruses 2 and 5. *Science.* 1997;275(5304):1320–1323. doi:10.1126/science.275.5304.1320",
    "mckay": "McKay JD, Hung RJ, Han Y, et al. Large-scale association analysis identifies new lung cancer susceptibility loci and heterogeneity in genetic susceptibility across histological subtypes. *Nat Genet.* 2017;49(7):1126–1132. doi:10.1038/ng.3892",
    "nathan": "Nathan A, Asgari S, Ishigaki K, et al. Single-cell eQTL models reveal dynamic T cell state dependence of disease loci. *Nature.* 2022;606(7912):120–128. doi:10.1038/s41586-022-04713-1",
    "fang": "Fang L, Yu W, Yu G, et al. Junctional adhesion molecule-like protein (JAML) is correlated with prognosis and immune infiltrates in lung adenocarcinoma. *Med Sci Monit.* 2021;27:e933503. doi:10.12659/MSM.933503",
    "tian": "Tian Q, Zhang F, Han C. Identification of JAML as an immune-associated prognostic marker in non-small cell lung cancer. *Lung Cancer: Targets and Therapy.* 2026;17:1–16. doi:10.2147/LCTT.S616193",
    "maier": "Maier B, Leader AM, Chen ST, et al. A conserved dendritic-cell regulatory program limits antitumour immunity. *Nature.* 2020;580(7802):257–262. doi:10.1038/s41586-020-2134-y",
    "cheng": "Cheng S, Li Z, Gao R, et al. A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells. *Cell.* 2021;184(3):792–809.e23. doi:10.1016/j.cell.2021.01.010",
    "garris": "Garris CS, Arlauckas SP, Kohler RH, et al. Successful anti-PD-1 cancer immunotherapy requires T cell–dendritic cell crosstalk involving the cytokines IFN-γ and IL-12. *Immunity.* 2018;49(6):1148–1161.e7. doi:10.1016/j.immuni.2018.09.024",
    "casanova": "Casanova-Acebes M, Dalla E, Leader AM, et al. Tissue-resident macrophages provide a pro-tumorigenic niche to early NSCLC cells. *Nature.* 2021;595(7868):578–584. doi:10.1038/s41586-021-03651-8",
    "denardo": "DeNardo DG, Ruffell B. Macrophages as regulators of tumour immunity and immunotherapy. *Nat Rev Immunol.* 2019;19(6):369–382. doi:10.1038/s41577-019-0127-6",
    "kee": "Kee BL. E and ID proteins branch out. *Nat Rev Immunol.* 2009;9(3):175–184. doi:10.1038/nri2507",
    "cannarile": "Cannarile MA, Lind NA, Rivera R, et al. Transcriptional regulator Id2 mediates CD8+ T cell immunity. *Nat Immunol.* 2006;7(12):1317–1325. doi:10.1038/ni1403",
    "hacker": "Hacker C, Kirsch RD, Ju X, et al. Transcriptional profiling identifies Id2 function in dendritic cell development. *Nat Immunol.* 2003;4(4):380–386. doi:10.1038/ni903",
    "westreich": "Westreich D, Greenland S. The Table 2 fallacy: presenting and interpreting confounder and modifier coefficients. *Am J Epidemiol.* 2013;177(4):292–298. doi:10.1093/aje/kws412",
    "munafo": "Munafò MR, Tilling K, Taylor AE, et al. Collider scope: when selection bias can substantially influence observed associations. *Int J Epidemiol.* 2018;47(1):226–235. doi:10.1093/ije/dyx206",
    "cobos": "Avila Cobos F, Alquicira-Hernandez J, Powell JE, et al. Benchmarking of cell type deconvolution pipelines for transcriptomics data. *Nat Commun.* 2020;11(1):5650. doi:10.1038/s41467-020-19015-1",
    "sturm": "Sturm G, Finotello F, Petitprez F, et al. Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. *Bioinformatics.* 2019;35(14):i436–i445. doi:10.1093/bioinformatics/btz363",
    "aran": "Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. *Nat Commun.* 2015;6:8971. doi:10.1038/ncomms9971",
    "dempster": "Dempster JM, Boyle I, Vazquez F, et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biol.* 2021;22(1):343. doi:10.1186/s13059-021-02540-7",
    "pratapa": "Pratapa A, Jalihal AP, Law JN, et al. Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data. *Nat Methods.* 2020;17(2):147–154. doi:10.1038/s41592-019-0690-6",
    "ayers": "Ayers M, Lunceford J, Nebozhyn M, et al. IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade. *J Clin Invest.* 2017;127(8):2930–2940. doi:10.1172/JCI91190",
}

TOKEN = re.compile(r"\{\{([A-Za-z0-9_]+)\}\}")
GROUP = re.compile(r"\{\{[A-Za-z0-9_]+\}\}(?:\s*,?\s*\{\{[A-Za-z0-9_]+\}\})*")

_order = []


def log(m):
    print(m, flush=True)


def read(p):
    return io.open(p, encoding="utf-8").read()


def load(name):
    """Read a rewrite partner file and normalise its trailing separator."""
    return read(os.path.join(NEW, name)).rstrip("\n") + "\n\n"


def section_span(t, start, end):
    """Return (i, j) with t[i:j] the region from `start` up to (not incl.) `end`."""
    i = t.index(start)
    j = t.index(end) if end else len(t)
    assert i < j, (start, end)
    return i, j


def expand(spec):
    """'1-5' / '18,19' -> [1,2,3,4,5] / [18,19]"""
    out = []
    for part in re.split(r"\s*[,，]\s*", spec):
        m = re.fullmatch(r"(\d+)\s*[–\-—]\s*(\d+)", part)
        if m:
            out.extend(range(int(m.group(1)), int(m.group(2)) + 1))
        else:
            out.append(int(part))
    return out


def _num(key):
    """Reserve a number for `key`, in the order keys are first encountered."""
    if key not in _order:
        _order.append(key)
    return _order.index(key) + 1


def _render_group(g):
    """Turn one citation group into '[3,7]' / '[12–15]' with ranges compressed."""
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


# ---------------------------------------------------------------- steps

log("[0] markers")
t = read(EN)
for marker in ["## Abstract", "## 1. Introduction", "## 2. Materials and Methods",
               "## Supplementary material", "## References", "## Figure legends"]:
    n = t.count(marker)
    assert n == 1, "marker %r found %d times" % (marker, n)
log("    all section markers unique")

if "{{" in t:
    raise SystemExit("文件里已存在 {{...}} 记号，本脚本只应运行一次（请先恢复 .bak）")

log("[1] numeric citations -> symbols")
seen = []


def _to_sym(m):
    seen.append(m.group(0))
    return ", ".join("{{%s}}" % OLD[n] for n in expand(m.group(1)))


t = re.sub(r"\[(\d+(?:\s*[–\-,，]\s*\d+)*)\]", _to_sym, t)
log("    %d citation groups converted: %s" % (len(seen), " ".join(seen)))

log("[2] splice rewritten sections")
i, j = section_span(t, "## Abstract", "## 1. Introduction")
t = t[:i] + load("abstract.md") + t[j:]
i, j = section_span(t, "## 1. Introduction", "## 2. Materials and Methods")
t = t[:i] + load("intro.md") + t[j:]
i, j = section_span(t, "## 2. Materials and Methods", "## Supplementary material")
t = t[:i] + load("methods.md") + load("results.md") + load("discussion.md") + t[j:]
log("    Abstract / Introduction / Methods / Results / Discussion spliced")

log("[3] renumber by order of first appearance")
t = GROUP.sub(lambda m: _render_group(m.group(0)), t)
log("    %d references cited" % len(_order))
missing = [k for k in REFS if k not in _order]
if missing:
    raise SystemExit("these references ended up uncited: %s" % missing)

log("[4] rebuild the reference list")
lines = ["%d. %s" % (i + 1, REFS[k]) for i, k in enumerate(_order)]
lines.append("")
lines.append("*All references were verified against Crossref by DOI (title, journal, year, "
             "volume, issue, pages and first author) on %s; the machine-readable verification "
             "record is included in Supplementary Material 17.*" % DATE)
new_refs = "## References\n\n" + "\n".join(lines) + "\n\n---\n\n"
i, j = section_span(t, "## References", "## Figure legends")
t = t[:i] + new_refs + t[j:]

bak = EN + ".bak_langrev_20261001"
if not os.path.exists(bak):
    io.open(bak, "w", encoding="utf-8").write(read(EN))
io.open(EN, "w", encoding="utf-8").write(t)
log("    wrote %s (%d chars); backup %s" % (os.path.basename(EN), len(t), os.path.basename(bak)))

log("[5] structural checks")
ref_block = t[t.index("## References"):t.index("## Figure legends")]
nums = [int(x) for x in re.findall(r"(?m)^(\d+)\.\s+\S", ref_block)]
assert nums == list(range(1, len(nums) + 1)), "reference numbering is not 1..N"
prose = t[:t.index("## References")] + t[t.index("## Figure legends"):]
cited = set()
forms = set()
for m in re.finditer(r"\[(\d+(?:\s*[–\-,]\s*\d+)*)\]", prose):
    forms.add(m.group(0))
    cited.update(expand(m.group(1)))
assert not [n for n in cited if n < 1 or n > len(nums)], "citation out of range"
uncited = [n for n in range(1, len(nums) + 1) if n not in cited]
log("    references: %d entries, numbering 1..%d continuous" % (len(nums), len(nums)))
log("    citation forms now in use: %s" % " ".join(sorted(forms, key=len)))
log("    uncited entries: %s" % (uncited if uncited else "none"))
log("    leftover symbols: %d" % (t.count("{{") + t.count("}}")))

L = t.split("\n")


def words(a, b):
    i = next(k for k, l in enumerate(L) if l.startswith(a))
    j = next(k for k, l in enumerate(L) if l.startswith(b))
    return len("\n".join(L[i:j]).split())


parts = [("Abstract", "## Abstract", "## 1. Introduction"),
         ("Introduction", "## 1. Introduction", "## 2. Materials"),
         ("Methods", "## 2. Materials", "## 3. Results"),
         ("Results", "## 3. Results", "## 4. Discussion"),
         ("Discussion", "## 4. Discussion", "## 5. Limitations"),
         ("Limitations", "## 5. Limitations", "## 6. Conclusion"),
         ("Conclusion", "## 6. Conclusion", "## Supplementary material")]
tot = 0
for name, a, b in parts:
    w = words(a, b)
    log("    %-13s %5d words" % (name, w))
    if name != "Abstract":
        tot += w
log("    BODY (Introduction -> Conclusion) %d words" % tot)
log("DONE")
