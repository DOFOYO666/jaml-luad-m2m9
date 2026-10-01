# Data sources

All input data are public. Nothing in this release redistributes raw data.

| Dataset | Accession / URL | Used for |
|---|---|---|
| GSE131907 | GEO | discovery single-cell cohort (LUAD, 208,506 cells) |
| GSE127465 | GEO | independent replication (treatment-naive NSCLC, 54,773 cells) |
| GSE135222 | GEO | checkpoint-blockade cohort (n = 27, PFS) |
| GSE126044 | GEO | evaluated and rejected: the gene of interest is absent from its matrix |
| TCGA-LUAD | https://portal.gdc.cancer.gov/ (expression via UCSC Xena) | tissue-level deconvolution and survival |
| DepMap 22Q2 | https://depmap.org/ | CRISPR (Chronos) gene effect |
| OneK1K single-cell eQTL | eQTL Catalogue, study QTS000038 | genetic anchor (not re-estimated here) |
| LUAD GWAS (McKay et al.) | GWAS Catalog, GCST004744 | genetic anchor (not re-estimated here) |

Reference resources pulled at run time: CollecTRI (via decoupleR), the liana consensus
ligand-receptor resource, the CellOracle human promoter base GRN (hg38), EPIC / MCP-counter /
CIBERSORT-LM22 signatures, ABSOLUTE tumour purity (TCGA clinical matrix).
