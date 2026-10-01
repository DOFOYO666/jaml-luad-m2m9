"""M1 步骤 2：F 值多口径对照表 + winner's curse 专项表。

产出：
  results/m1_f_reconciliation.csv   —— JAML 三种 F 口径对照
  results/m1_winners_curse.csv      —— 全部 22 基因 ALL vs F>10 的 IVW 对照
  results/m1_supplementary.md       —— 可直接进 Supplementary 的表格
"""
import numpy as np
import pandas as pd
from scipy.stats import norm


def md_table(df):
    """极简 markdown 表格（避免依赖 tabulate）。"""
    cols = [str(c) for c in df.columns]
    out = ["| " + " | ".join(cols) + " |",
           "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]) for c in df.columns) + " |")
    return "\n".join(out)

# ---------------- 1. JAML 三种 F 口径 ----------------
qc = pd.read_csv("results/m1_jaml_instrument_qc.csv")

# 口径 A：OneK1K 官方 top-eQTL summary（本机可核），F = (slope/slope_se)^2
A = qc[["cell_class", "celltype_file", "variant_id", "MAF", "slope", "slope_se", "F",
        "pval_nominal", "pval_perm", "qval"]].copy()

# 口径 B：反推 F —— QTD000613(eQTL Catalogue) 的 JAML top cis-eQTL
#   beta 来自 permuted 文件（官方β），SE 由 nominal p 反推
BETA = 0.853937
P_NOM = 2.691150e-102
se_rev = abs(BETA) / norm.isf(P_NOM / 2)
F_rev = (BETA / se_rev) ** 2

# 口径 C：手稿 Table 1/3 报告的官方 F —— 取自 OneK1K all-variant 文件（本机已缺失）
F_ms = 590.5136331948422

tcga_pairs = []
for _, r in A.iterrows():
    tcga_pairs.append({
        "F口径": "A. OneK1K 官方 top-eQTL summary（本机可核）",
        "细胞类型": r["celltype_file"],
        "variant": r["variant_id"],
        "MAF": round(r["MAF"], 4),
        "beta/slope": round(r["slope"], 6),
        "SE": round(r["slope_se"], 6),
        "F": round(r["F"], 2),
        "置换P": r["pval_perm"],
        "qval": r["qval"],
    })
recon = pd.DataFrame(tcga_pairs)

extra = pd.DataFrame([
    {"F口径": "B. 反推 F（eQTL Catalogue QTD000613，SE 由 nominal P 反推）",
     "细胞类型": "CD4_TCM", "variant": "chr11_118234507_C_T", "MAF": np.nan,
     "beta/slope": BETA, "SE": round(se_rev, 6), "F": round(F_rev, 2),
     "置换P": 0.000999, "qval": np.nan},
    {"F口径": "C. 手稿报告官方 F（OneK1K all-variant 文件；本机已缺失）",
     "细胞类型": "CD4_TCM", "variant": "chr11:118234507 C>T", "MAF": np.nan,
     "beta/slope": np.nan, "SE": round(BETA / np.sqrt(F_ms), 6), "F": round(F_ms, 2),
     "置换P": np.nan, "qval": np.nan},
])
recon = pd.concat([recon, extra], ignore_index=True)
recon.to_csv("results/m1_f_reconciliation.csv", index=False)

# ---------------- 2. winner's curse 对照（22 基因） ----------------
ms = pd.read_csv("results/multisnp_mr_results_ld.csv")
piv = ms.pivot_table(index="gene", columns="tag",
                     values=["n_snp", "or_ivw", "p_ivw", "p_intercept"], aggfunc="first")
wc = pd.DataFrame({
    "n_snp_all": piv[("n_snp", "all")],
    "n_snp_f10": piv[("n_snp", "f10")],
    "OR_all": piv[("or_ivw", "all")],
    "OR_f10": piv[("or_ivw", "f10")],
    "p_intercept_f10": piv[("p_intercept", "f10")],
}).reset_index()
wc["OR膨胀比_all_over_f10"] = wc["OR_all"] / wc["OR_f10"]
wc["方向是否翻转"] = np.where(
    (wc["OR_all"] > 1) == (wc["OR_f10"] > 1), "否", "是")
wc = wc.sort_values("OR膨胀比_all_over_f10", ascending=False, na_position="last")
wc.to_csv("results/m1_winners_curse.csv", index=False)

# ---------------- 3. Summary markdown ----------------
jaml = wc[wc["gene"] == "JAML"].iloc[0]
lines = []
lines.append("# M1 补充材料：工具变量强度与 F 值规范化报告\n")
lines.append("## 表 S-M1-1　JAML 工具变量 MAF 与强度（OneK1K 官方 top-eQTL summary，14 种细胞类型）\n")
tbl = A.copy()
tbl["MAF"] = tbl["MAF"].round(4)
tbl["slope"] = tbl["slope"].round(6)
tbl["slope_se"] = tbl["slope_se"].round(6)
tbl["F"] = tbl["F"].round(2)
tbl["pval_perm"] = tbl["pval_perm"].map(lambda x: f"{x:.4f}")
tbl["qval"] = tbl["qval"].map(lambda x: f"{x:.3g}")
tbl = tbl.rename(columns={
    "cell_class": "细胞大类", "celltype_file": "细胞类型", "variant_id": "变异",
    "MAF": "MAF", "slope": "slope(β)", "slope_se": "slope_se(SE)", "F": "F=(β/SE)²",
    "pval_nominal": "nominal P", "pval_perm": "置换 P", "qval": "qval"})
lines.append(md_table(tbl))
lines.append("")
lines.append(f"- **MAF 范围：{A['MAF'].min():.4f} – {A['MAF'].max():.4f}；MAF < 0.05 的条目数：{(A['MAF'] < 0.05).sum()}**")
lines.append(f"- F 范围：{A['F'].min():.2f} – {A['F'].max():.2f}；F > 10 的细胞类型数：{(A['F'] > 10).sum()} / {len(A)}")
lines.append("")
lines.append("## 表 S-M1-2　JAML 的 F 值三口径对照\n")
lines.append(md_table(recon))
lines.append("")
lines.append("## 表 S-M1-3　winner's curse 专项：全部 LD-pruned 变体（ALL）vs 仅 F>10 变体（F10）\n")
w2 = wc.copy()
for c in ["OR_all", "OR_f10", "OR膨胀比_all_over_f10"]:
    w2[c] = w2[c].round(3)
w2["p_intercept_f10"] = w2["p_intercept_f10"].map(lambda x: "nan" if pd.isna(x) else f"{x:.3g}")
lines.append(md_table(w2))
lines.append("")
lines.append(f"- JAML：ALL（{int(jaml['n_snp_all'])} SNP）IVW OR = {jaml['OR_all']:.3f} → "
             f"F>10（{int(jaml['n_snp_f10'])} SNP）IVW OR = {jaml['OR_f10']:.3f}；"
             f"Egger 截距 P = {jaml['p_intercept_f10']:.3f}（无多效性证据）")
lines.append(f"- 方向翻转的基因数：{(wc['方向是否翻转'] == '是').sum()} / {len(wc)}")

with open("results/m1_supplementary.md", "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

print("=" * 80)
print("JAML F 值三口径对照（反推 SE 复算验证）")
print("=" * 80)
print(f"  反推 SE = |β| / z,  z = z(1 - P/2),  P = {P_NOM:.3g}")
print(f"          = {BETA} / {se_rev:.6f}")
print(f"  → 反推 F = (β/SE)² = {F_rev:.2f}   （sig 表记录值 461.16 → 复现{'成功' if abs(F_rev-461.16)<0.5 else '失败'}）")
print(f"  手稿官方 F = {F_ms:.2f}  → 隐含官方 SE = {BETA/np.sqrt(F_ms):.6f}")
print(f"  反推 SE / 官方 SE = {se_rev/(BETA/np.sqrt(F_ms)):.3f}")
print()
print("=" * 80)
print("winner's curse 对照（按 OR 膨胀比降序，前 8 行）")
print("=" * 80)
print(w2.head(8).to_string(index=False))
print()
print("已写出: results/m1_f_reconciliation.csv")
print("        results/m1_winners_curse.csv")
print("        results/m1_supplementary.md")
