# -*- coding: utf-8 -*-
"""Count how many of the 22 LUAD candidate genes are evaluable in GSE126044.
Resolves an EN/CN inconsistency in the M2-M9 manuscript Limitations.
"""
import gzip, os, sys, io
sys.stdout.reconfigure(encoding="utf-8")

BASE = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化"
PATH = os.path.join(BASE, "rawdata", "icb", "GSE126044", "GSE126044_counts.txt.gz")

CAND = ["JAML", "AMICA1", "FUBP1", "MAP4K4", "NUMBL", "PARVA", "IREB2", "RPS6KA2",
        "EID1", "HYKK", "AGPHD1", "STMN3", "GALK2", "RNASET2", "DNAJA4", "SECISBP2L",
        "RAB31", "CMIP", "HLA-C", "RAB4B",
        "ZNRD1ASP", "CTC-490E21.14", "RP1-167A14.3", "RP11-514O12.4"]

print("exists:", os.path.exists(PATH), "size:", os.path.getsize(PATH) if os.path.exists(PATH) else None)

hit = set()
n_genes = 0
with gzip.open(PATH, "rt", encoding="utf-8", errors="replace") as f:
    header = f.readline().rstrip("\n").split("\t")
    print("header cols:", len(header), "->", header[:4], "...")
    for line in f:
        g = line.split("\t", 1)[0].strip()
        n_genes += 1
        if g in CAND:
            hit.add(g)

print("genes in matrix:", n_genes)
print("hits:", sorted(hit))
present_canonical = [g for g in ["JAML","FUBP1","MAP4K4","NUMBL","PARVA","IREB2","RPS6KA2","EID1",
                                 "HYKK","STMN3","GALK2","RNASET2","DNAJA4","SECISBP2L","RAB31",
                                 "CMIP","HLA-C","RAB4B"] if g in hit or (g=="JAML" and "AMICA1" in hit)
                                                       or (g=="HYKK" and "AGPHD1" in hit)]
print("canonical present:", len(present_canonical), present_canonical)
print("of 22 candidate loci:", len([g for g in hit]))
