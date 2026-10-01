# -*- coding: utf-8 -*-
"""
02_scmr_luad_risk.py v4-fine-celltype - 细胞类型特异 sc-eQTL MR 主分析（流式内存优化）
"""
import os, glob, json, sys
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

EQTC_DIR = os.path.join(DATA, "onek1k", "eqtl_catalogue")
GWAS = os.path.join(DATA, "luad_gwas", "GCST004744.h.tsv.gz")
ENSG2SYM = json.load(open(os.path.join(DATA, "onek1k", "ensg2sym.json"), encoding="utf-8"))

# QTD → 细胞类型大类（v6：依据 eQTL Catalogue 官方 dataset_id_map.tsv + dataset_metadata_r8.tsv 复核，2026-08-21）
# 官方 sample_group: 606=B_intermediate 607=B_memory 608=B_naive 609=CD14_Mono 610=CD16_Mono
#   611=CD4_CTL 612=CD4_Naive 613=CD4_TCM 614=CD4_TEM 615=CD8_Naive 616=CD8_TCM 617=CD8_TEM
#   618=HSPC 619=MAIT 620=NK 621=NK_CD56bright 622=NK_Proliferating 623=Plasmablast 624=Platelet
#   625=Treg 626=cDC2 627=dnT 628=gdT 629=pDC
# 大类映射：MAIT/gdT→CD8 T（T 细胞系）；Plasmablast→Plasma；cDC2/pDC→DC；
#   HSPC/Platelet/dnT 非 8 大类免疫细胞，排除
QTD_MAP = {
    "QTD000606": "B_intermediate", "QTD000607": "B_memory", "QTD000608": "B_naive",
    "QTD000609": "CD14_Mono", "QTD000610": "CD16_Mono",
    "QTD000611": "CD4_CTL", "QTD000612": "CD4_Naive", "QTD000613": "CD4_TCM", "QTD000614": "CD4_TEM",
    "QTD000615": "CD8_Naive", "QTD000616": "CD8_TCM", "QTD000617": "CD8_TEM",
    "QTD000619": "MAIT", "QTD000620": "NK", "QTD000621": "NK_CD56bright", "QTD000622": "NK_Proliferating",
    "QTD000623": "Plasmablast", "QTD000625": "Treg", "QTD000626": "cDC2", "QTD000628": "gdT", "QTD000629": "pDC",
}
EXCLUDED = {"QTD000618", "QTD000624", "QTD000627"}  # HSPC / Platelet / dnT
QTD_FILES = [f for f in sorted(glob.glob(os.path.join(EQTC_DIR, "QTD*.permuted.tsv.gz")))
             if os.path.basename(f).split(".")[0] not in EXCLUDED
             and os.path.basename(f).split(".")[0] in QTD_MAP]


def parse_variant(v):
    s = str(v).replace("chr", "")
    parts = s.split("_")
    if len(parts) < 4:
        return None
    return parts[0], parts[1], parts[2], parts[3]


def load_eqtl_need():
    """读取所有 QTD 的 eQTL，返回：需要 vid 集合 + (chr,ref,alt) 组合"""
    need_vid = set()
    need_cra = set()
    for fp in QTD_FILES:
        df = pd.read_csv(fp, sep="\t", compression="gzip", usecols=["variant"], low_memory=False)
        for v in df["variant"].astype(str):
            p = parse_variant(v)
            if p:
                need_vid.add(f"{p[0]}_{p[1]}_{p[2]}_{p[3]}")
                need_cra.add((p[0], p[2], p[3]))
    print(f"[eQTL] 需要精确 vid: {len(need_vid)}, chr_ref_alt 组合: {len(need_cra)}", flush=True)
    return need_vid, need_cra


def stream_gwas(need_vid, need_cra):
    """流式读取 GWAS，只保留需要的行；返回 (df_sub, g_cra)"""
    cols = ["hm_variant_id", "hm_effect_allele", "hm_other_allele", "hm_beta",
            "standard_error", "p_value"]
    keep = []
    g_cra = {}
    for ch in pd.read_csv(GWAS, sep="\t", compression="gzip", chunksize=1_000_000,
                          low_memory=False, usecols=cols):
        ch = ch.rename(columns={"hm_beta": "beta", "standard_error": "se", "p_value": "p"})
        ch["vid"] = ch["hm_variant_id"].astype(str).str.replace("chr", "", regex=False)
        ch["ea"] = ch["hm_effect_allele"].astype(str).str.upper()
        ch["oa"] = ch["hm_other_allele"].astype(str).str.upper()
        for c in ["beta", "se", "p"]:
            ch[c] = pd.to_numeric(ch[c], errors="coerce")
        hit = ch[ch["vid"].isin(need_vid)]
        keep.append(hit)
    # 合并命中行
    g = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
    if not g.empty:
        g = g.drop_duplicates(subset=["vid"], keep="first")
        parts = g["vid"].str.split("_", expand=True)
        if parts.shape[1] >= 4:
            cra_keys = list(zip(parts[0], parts[2], parts[3]))
            for i, key in enumerate(cra_keys):
                if key in need_cra and key not in g_cra:
                    g_cra[key] = g["vid"].iloc[i]
    print(f"[GWAS] 命中行: {len(g)}", flush=True)
    return g, g_cra


def main():
    need_vid, need_cra = load_eqtl_need()
    gwas, g_cra = stream_gwas(need_vid, need_cra)
    if gwas.empty:
        print("GWAS 无命中"); return
    gwas_map = gwas.set_index("vid")
    g_vid = set(gwas["vid"])

    results = []
    for fp in QTD_FILES:
        qid = os.path.basename(fp).split(".")[0]
        ct = QTD_MAP[qid]
        df = pd.read_csv(fp, sep="\t", compression="gzip", low_memory=False)
        df["gene_symbol"] = df["molecular_trait_id"].map(ENSG2SYM)
        df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
        df["p_perm"] = pd.to_numeric(df["p_perm"], errors="coerce")
        df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
        df = df[(df["gene_symbol"].notna()) & (df["p_perm"] < 0.05) & (df["beta"] != 0)
                & (df["pvalue"] > 0) & (df["pvalue"] < 1)]
        if df.empty:
            continue
        z_e = np.abs(stats.norm.ppf(df["pvalue"] / 2))
        df["se"] = np.abs(df["beta"] / z_e)
        df["F"] = z_e ** 2
        parsed = df["variant"].apply(parse_variant)
        df["chr"] = parsed.apply(lambda x: x[0] if x else None)
        df["pos"] = parsed.apply(lambda x: x[1] if x else None)
        df["ref"] = parsed.apply(lambda x: x[2] if x else None)
        df["alt"] = parsed.apply(lambda x: x[3] if x else None)
        df["vid"] = df["chr"] + "_" + df["pos"] + "_" + df["ref"] + "_" + df["alt"]
        df = df.dropna(subset=["vid"])

        m1 = df[df["vid"].isin(g_vid)].copy()
        m1["match_type"] = "exact"
        rest = df[~df["vid"].isin(g_vid)].copy()
        m2 = rest[rest.apply(lambda r: (r["chr"], r["ref"], r["alt"]) in g_cra, axis=1)].copy()
        if not m2.empty:
            m2["vid"] = m2.apply(lambda r: g_cra[(r["chr"], r["ref"], r["alt"])], axis=1)
            m2["match_type"] = "ref_alt"
        merged = pd.concat([m1, m2], ignore_index=True) if not m2.empty else m1
        if merged.empty:
            print(f"  [{qid}->{ct}] 无匹配", flush=True)
            continue
        merged = merged.reset_index(drop=True).merge(
            gwas_map.reset_index(), on="vid", how="left", suffixes=("", "_g"))
        merged = merged.dropna(subset=["beta_g", "se_g"])
        ea = merged["alt"].str.upper()
        oa = merged["ref"].str.upper()
        same = (ea == merged["ea"]) & (oa == merged["oa"])
        flip = (ea == merged["oa"]) & (oa == merged["ea"])
        merged["beta_g_aligned"] = np.where(same, merged["beta_g"],
                                            np.where(flip, -merged["beta_g"], np.nan))
        merged = merged.dropna(subset=["beta_g_aligned"])
        if merged.empty:
            continue
        merged = merged.sort_values("p_perm").drop_duplicates(subset=["molecular_trait_id"], keep="first")
        merged["b_mr"] = merged["beta_g_aligned"] / merged["beta"]
        merged["se_mr"] = np.abs(merged["se_g"] / merged["beta"])
        merged["or"] = np.exp(merged["b_mr"])
        merged["p_mr"] = 2 * (1 - stats.norm.cdf(np.abs(merged["b_mr"] / merged["se_mr"])))
        merged["cell_type"] = ct
        merged["qtd"] = qid
        results.append(merged[["qtd", "cell_type", "molecular_trait_id", "gene_symbol", "variant",
                               "beta", "se", "F", "b_mr", "se_mr", "or", "p_mr", "match_type"]])
        print(f"  [{qid}->{ct}] 匹配基因 {len(merged)}", flush=True)

    if not results:
        print("无结果"); return
    res = pd.concat(results, ignore_index=True)
    res = res[res["F"] >= 10]
    res["fdr"] = stats.false_discovery_control(res["p_mr"], method="bh")
    res = res.sort_values("p_mr")
    res.to_csv(os.path.join(RES, "scmr_luad_risk_all.csv"), index=False)
    sig = res[res["fdr"] < 0.05]
    sig.to_csv(os.path.join(RES, "scmr_luad_risk_sig_fine.csv"), index=False)

    print(f"\n[完成] 总基因-细胞对: {len(res)} | FDR<0.05: {len(sig)}")
    print("各细胞类型基因数:")
    print(res.groupby("cell_type")["gene_symbol"].count().to_string())
    if not sig.empty:
        print("\n显著结果按细胞类型:")
        print(sig.groupby("cell_type")["gene_symbol"].count().to_string())
        print("\nTop 20 显著结果:")
        print(sig.head(20).to_string())


if __name__ == "__main__":
    main()
