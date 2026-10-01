# -*- coding: utf-8 -*-
"""
06_drug_annotation.py - 药物靶点注释
对 scMR 显著基因查询 Open Targets Platform 已知药物（若 API 可达），
否则用内置 curated 免疫/肿瘤靶点表。
输出：results/drug_annotation.csv
"""
import os, json, urllib.request
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
SCMR_SIG = os.path.join(RES, "scmr_luad_risk_sig.csv")

OT_URL = "https://platform.opentargets.org/api/v4/graphql"
OT_QUERY = """
query target($ensemblIds: [String!]!) {
  targets(ensemblIds: $ensemblIds) {
    id approvedSymbol approvedName
    knownDrugs(size: 5) {
      count
      rows { drug { id name } mechanismOfAction }
    }
  }
}
"""

# ENSG -> symbol 映射（复用 gencode）
ENSG2SYM = json.load(open(os.path.join(BASE, "data", "onek1k", "ensg2sym.json"), encoding="utf-8"))
SYM2ENSG = {}
for e, s in ENSG2SYM.items():
    SYM2ENSG.setdefault(s, e)

# curated 表（免疫/肿瘤相关可成药靶点标注）
CURATED = {
    "PDCD1": ("PD-1", "免疫检查点抑制剂（已批准 NSCLC）"),
    "CD274": ("PD-L1", "免疫检查点抑制剂（已批准 NSCLC）"),
    "CTLA4": ("CTLA-4", "免疫检查点抑制剂"),
    "EGFR": ("EGFR", "TKI 已批准 LUAD"),
    "ALK": ("ALK", "TKI 已批准 LUAD"),
    "KRAS": ("KRAS", "KRAS G12C 抑制剂已批准"),
    "MET": ("MET", "c-MET 抑制剂在研"),
    "VEGFA": ("VEGF", "抗血管生成"),
    "TNFRSF9": ("4-1BB", "共刺激靶点在研"),
    "JAML": ("JAML", "连接粘附分子样蛋白，无获批药物"),
    "RNASET2": ("RNASET2", "核糖核酸酶，无获批药物"),
    "IREB2": ("IREB2", "铁反应元件结合蛋白2，无获批药物"),
}


def query_ot(genes):
    """Open Targets GraphQL 查询，返回 gene -> known drug 数"""
    ensembl_ids = [SYM2ENSG.get(g) for g in genes if SYM2ENSG.get(g)]
    if not ensembl_ids:
        return {}
    body = json.dumps({"query": OT_QUERY, "variables": {"ensemblIds": ensembl_ids}}).encode()
    req = urllib.request.Request(OT_URL, data=body,
                                 headers={"Content-Type": "application/json",
                                          "User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            data = json.loads(r.read().decode())
        out = {}
        for t in data.get("data", {}).get("targets", []):
            sym = t.get("approvedSymbol")
            kd = t.get("knownDrugs") or {}
            drugs = [(d.get("drug", {}).get("name"), d.get("mechanismOfAction"))
                     for d in (kd.get("rows") or [])]
            out[sym] = {"drug_count": kd.get("count", 0), "drugs": drugs}
        return out
    except Exception as e:
        print(f"[Open Targets API 不可达: {e}] 使用 curated 表")
        return {}


def main():
    if not os.path.exists(SCMR_SIG):
        print("[跳过] 无 scMR 显著结果"); return
    sig = pd.read_csv(SCMR_SIG)
    genes = sorted(set(sig["gene_symbol"]))
    print(f"[输入] {len(genes)} 个显著基因")

    ot = query_ot(genes)
    rows = []
    for g in genes:
        r = {"gene": g}
        if g in ot and ot[g]["drug_count"] > 0:
            r["druggable"] = "yes"
            r["n_drugs"] = ot[g]["drug_count"]
            r["drug_names"] = "; ".join(sorted(set(d[0] for d in ot[g]["drugs"] if d[0])))[:200]
            r["note"] = "Open Targets known drugs"
        elif g in CURATED:
            r["druggable"] = "yes" if CURATED[g][1] and "无获批" not in CURATED[g][1] else "potential"
            r["n_drugs"] = 0
            r["drug_names"] = ""
            r["note"] = CURATED[g][1]
        else:
            r["druggable"] = "no"
            r["n_drugs"] = 0
            r["drug_names"] = ""
            r["note"] = ""
        rows.append(r)
    res = pd.DataFrame(rows)
    res = res.merge(sig.drop_duplicates("gene_symbol")[["gene_symbol", "cell_type", "or", "p_mr", "fdr"]],
                    left_on="gene", right_on="gene_symbol", how="left")
    res = res.drop(columns=["gene_symbol"])
    res.to_csv(os.path.join(RES, "drug_annotation.csv"), index=False)
    print(f"[完成] 注释 {len(res)} 基因")
    print(res[["gene", "cell_type", "or", "druggable", "n_drugs", "drug_names", "note"]].to_string(index=False))


if __name__ == "__main__":
    main()
