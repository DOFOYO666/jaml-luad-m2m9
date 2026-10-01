# -*- coding: utf-8 -*-
"""构建《M2-M9 深化研究》投稿包：补充材料、图件与清单。

输入：主工作区 results/（M2/M3/M6/M7/M8-M9）与 D 盘工作区 results/（M4/M8）。
输出：M2-M9深化研究稿/04_Supplementary/*.csv、03_Figures/*.{png,tif}、MANIFEST.md

用法：python 48_build_new_pkg.py
"""
import io
import os
import shutil

ROOT = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化"
D = r"D:\workbuddy工作空间\JAML深度研究"
PKG = os.path.join(ROOT, "M2-M9深化研究稿")
SUPP = os.path.join(PKG, "04_Supplementary")
FIG = os.path.join(PKG, "03_Figures")
RES = os.path.join(ROOT, "results")
DRES = os.path.join(D, "results")
DFIG = os.path.join(D, "results", "figures")
RFIG = os.path.join(RES, "figures")
BMCSUPP = os.path.join(ROOT, "BMCCancer投稿", "05_Supplementary")

log = []


def mkd():
    for d in (SUPP, FIG):
        os.makedirs(d, exist_ok=True)


def copy_csv(src, dst_name, title, note=""):
    if not os.path.isfile(src):
        log.append(f"[缺失] {src}")
        return False
    with io.open(src, encoding="utf-8-sig") as f:
        body = f.read()
    with io.open(os.path.join(SUPP, dst_name), "w", encoding="utf-8-sig", newline="") as f:
        f.write(f"# {title}\n")
        if note:
            f.write(f"# {note}\n")
        f.write(f"# source: {os.path.relpath(src, ROOT) if src.startswith(ROOT) else src}\n")
        f.write(body)
    rows = max(0, body.count("\n") - 1)
    log.append(f"[补充] {dst_name}  ({rows} 行)")
    return True


def copy_fig(stem, out_name, search_dirs, exts=("png", "tif")):
    ok = []
    for e in exts:
        src = None
        for d in search_dirs:
            p = os.path.join(d, f"{stem}.{e}")
            if os.path.isfile(p):
                src = p
                break
        if src:
            ext = "tiff" if e == "tif" else e
            dst = os.path.join(FIG, f"{out_name}.{ext}")
            shutil.copy2(src, dst)
            ok.append(f"{out_name}.{ext} ({os.path.getsize(dst)/1024:.0f} KB)")
    if ok:
        log.append("[图] " + " / ".join(ok))
    else:
        log.append(f"[图缺失] {stem}")

    return bool(ok)


def main():
    mkd()
    print("=== 补充材料 ===")
    # AF1 对齐核验（复用已生成版本）
    copy_csv(os.path.join(BMCSUPP, "Additional file 12 - GSE127465 alignment verification.csv"),
             "AdditionalFile01_GSE127465_alignment_verification.csv",
             "Additional file 1. Alignment verification for the independent replication dataset (GSE127465).",
             "Per-cell mitochondrial read fraction vs metadata (Spearman rho), marker-compartment checks, matrix structure, Total counts caveat, independent re-implementation check.")
    # AF2 区室阳性率
    copy_csv(os.path.join(RES, "cross_dataset_compartment_pct.csv"),
             "AdditionalFile02_compartment_positivity_two_datasets.csv",
             "Additional file 2. Compartmental positivity and mean expression of JAML and CXADR in the discovery (GSE131907) and replication (GSE127465) datasets.")
    # AF3 偏相关
    copy_csv(os.path.join(RES, "m2_tcga_jaml_partial_corr.csv"),
             "AdditionalFile03_partial_correlation_all_populations.csv",
             "Additional file 3. Correlation of JAML expression with all deconvolved populations in TCGA-LUAD: raw, adjusted for global immune content (G_mc, G_epic) and for tumour purity (G_pur).")
    # AF4 marker 对照 + 相关谱相似度
    copy_csv(os.path.join(RES, "m2_marker_specificity.csv"),
             "AdditionalFile04a_adjustment_control_marker_genes.csv",
             "Additional file 4a. Marker-gene control for the adjustment procedure (rho_raw vs rho_adj under G_epic, rho_adj_mc under G_mc).")
    copy_csv(os.path.join(RES, "m2_marker_profile_similarity.csv"),
             "AdditionalFile04b_profile_similarity_39populations.csv",
             "Additional file 4b. Pearson similarity of 39-population correlation profiles between JAML and 15 marker genes.")
    # AF5 临床
    copy_csv(os.path.join(RES, "m2_tcga_survival_proper.csv"),
             "AdditionalFile05a_TCGA_survival_all_genes.csv",
             "Additional file 5a. TCGA-LUAD expression-survival analysis (Cox and log-rank) for all evaluable candidate genes; n = 563, 208 deaths.")
    copy_csv(os.path.join(RES, "m2_tcga_stage_jaml.csv"),
             "AdditionalFile05b_TCGA_stage_JAML.csv",
             "Additional file 5b. JAML expression by pathological stage in TCGA-LUAD.")
    # AF6 M8
    copy_csv(os.path.join(DRES, "m8_tf_group_diff.csv"),
             "AdditionalFile06_TF_activity_group_differences.csv",
             "Additional file 6. Transcription-factor activity in JAML-positive versus JAML-negative CD4+ T cells and within positive cells (decoupleR + CollecTRI; 498 TFs). Negative control: 200 label permutations (mean 0.04, max 2 significant TFs). Positive control: PDCD1, 5/9 expected TFs.")
    # AF7 M4
    copy_csv(os.path.join(DRES, "m4_round_compare.csv"),
             "AdditionalFile07a_TF_perturbation_three_configurations.csv",
             "Additional file 7a. Perturbation effect on JAML for all candidate TFs across three network configurations, with paired noise calibration.")
    copy_csv(os.path.join(DRES, "m4_noise_control_classify.csv"),
             "AdditionalFile07b_paired_randomized_control.csv",
             "Additional file 7b. Paired randomized-network control: real versus 5 freshly seeded randomized replicates per TF, with degenerate-control flag and global noise ceiling (9.59e-03).")
    # AF8 交叉
    copy_csv(os.path.join(DRES, "m4_m8_cross_tf.csv"),
             "AdditionalFile08_cross_module_concordance.csv",
             "Additional file 8. Concordance between TF activity differences (observational) and perturbation effects (interventional) across 32 shared TFs.")
    # AF9 DepMap
    p = os.path.join(RES, "m6_depmap_jaml_summary.json")
    if os.path.isfile(p):
        import json
        j = json.load(open(p, encoding="utf-8"))
        lines = ["key,value"]
        for k, v in j.items():
            if isinstance(v, dict):
                for k2, v2 in v.items():
                    lines.append(f'"{k}.{k2}","{v2}"')
            else:
                lines.append(f'"{k}","{v}"')
        with io.open(os.path.join(SUPP, "AdditionalFile09_DepMap_JAML_summary.csv"), "w",
                     encoding="utf-8-sig", newline="") as f:
            f.write("# Additional file 9. DepMap 22Q2 CRISPR gene-effect summary for JAML.\n")
            f.write("# source: results/m6_depmap_jaml_summary.json\n")
            f.write("\n".join(lines) + "\n")
        log.append("[补充] AdditionalFile09_DepMap_JAML_summary.csv")
    else:
        log.append("[缺失] m6_depmap_jaml_summary.json")

    # AF10 审稿人要求补充的三项分析
    copy_csv(os.path.join(RES, "reviewer_R1_ID2_relative_effect.csv"),
             "AdditionalFile10a_ID2_relative_effect.csv",
             "Additional file 10a. Relative effect size of ID2 perturbation on JAML (reviewer-requested).",
             "mean |delta JAML| as a fraction of mean JAML expression, per-cell distribution, and rank of JAML among 3,074 genes.")
    copy_csv(os.path.join(RES, "reviewer_R2_TCGA_parameterisation.csv"),
             "AdditionalFile10b_TCGA_survival_parameterisations.csv",
             "Additional file 10b. TCGA-LUAD overall survival for all 18 evaluable genes under four parameterizations (reviewer-requested).",
             "Pre-specified primary model = continuous per-SD Cox; secondary = per-log2-unit Cox, categorical Cox, and log-rank.")
    copy_csv(os.path.join(RES, "reviewer_R3_M5_power.json"),
             "AdditionalFile10c_M5_power_analysis.csv",
             "Additional file 10c. Power analysis for the checkpoint-blockade cohort (GSE135222).",
             "n = 27, 21 events; minimum detectable hazard ratio at 80% power, alpha = 0.05 two-sided.")

    print("\n=== 图件 ===")
    copy_fig("m3_jaml_cxadr_normalized", "Figure1_compartment_discovery", [RFIG])
    copy_fig("m3_jaml_cd4_subtypes", "Figure1b_CD4_subtypes_discovery", [RFIG])
    copy_fig("cross_dataset_compartment", "Figure2_compartment_replication", [RFIG])
    copy_fig("m8_gse127465_jaml_cxadr", "Figure2b_replication_panel", [RFIG])
    copy_fig("m2_tcga_jaml_raw_vs_partial", "Figure3_TCGA_raw_vs_adjusted", [RFIG])
    copy_fig("m2_marker_validity_Tcell", "Figure4a_adjustment_validity", [RFIG])
    copy_fig("m2_marker_specificity", "Figure4b_profile_similarity", [RFIG])
    copy_fig("m7_pseudotime_jaml", "Figure5_pseudotime_CD4", [RFIG])
    copy_fig("m8_tf_jaml", "Figure6a_TF_activity", [DFIG])
    copy_fig("m8_negative_control", "Figure6b_TF_activity_negative_control", [DFIG])
    copy_fig("m4_tf_ko_effect_on_jaml", "Figure7a_TF_perturbation", [DFIG])
    copy_fig("m4_noise_control", "Figure7b_paired_randomized_control", [DFIG])
    copy_fig("m4_m8_cross_tf", "Figure8a_activity_vs_perturbation", [DFIG])
    copy_fig("m4_round_compare", "Figure8b_cross_configuration_reproducibility", [DFIG])
    copy_fig("m6_depmap_jaml_effect", "Figure9_DepMap", [RFIG])

    for line in log:
        print("  " + line)

    n_supp = len([f for f in os.listdir(SUPP) if f.endswith(".csv")])
    n_tif = len([f for f in os.listdir(FIG) if f.endswith(".tiff")])
    n_png = len([f for f in os.listdir(FIG) if f.endswith(".png")])
    print(f"\n合计：补充材料 {n_supp} 个 CSV；图 {n_png} 个 PNG + {n_tif} 个 TIFF")
    with io.open(os.path.join(PKG, "BUILD_LOG.md"), "w", encoding="utf-8", newline="") as f:
        f.write("# 投稿包构建记录\n\n")
        f.write(f"- 补充材料：{n_supp} 个 CSV\n- 图：{n_png} PNG / {n_tif} TIFF\n\n## 明细\n\n")
        for line in log:
            f.write(f"- {line}\n")
    print("DONE")


if __name__ == "__main__":
    main()
