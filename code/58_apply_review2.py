# -*- coding: utf-8 -*-
"""Second-round reviewer revisions for the M2-M9 manuscript (EN + CN).

Discipline (project convention):
  * whole-string, exactly-once replacements; the script aborts if any fragment does not match once;
  * no cross-line regex, no line-number addressing;
  * long paragraphs are replaced whole (anchor -> next blank line) rather than by fragile substrings;
  * a residual scan and a re-count of uncited references run at the end.

Run:  python scripts/58_apply_review2.py
"""
import collections
import os
import re
import shutil
import sys
import time

sys.stdout.reconfigure(encoding="utf-8")

PKG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "M2-M9深化研究稿")
EN = os.path.join(PKG, "manuscript_EN.md")
CN = os.path.join(PKG, "manuscript_CN.md")

OPS = []          # ('sub', path, old, new, label) | ('para', path, anchor, new, label)


def sub(path, old, new, label):
    OPS.append(("sub", path, old, new, label))


def para(path, anchor, new, label):
    OPS.append(("para", path, anchor, new, label))


E, C = EN, CN

# ================================================================ ENGLISH ==
sub(E, "one gene with eight orthogonal analyses, each with an explicit control",
       "one gene with eight complementary analyses, each with an explicit control",
       "EN-abstract: orthogonal -> complementary")

sub(E, "with eight orthogonal analyses, each with an explicit negative or positive control",
       "with eight complementary analyses, each with an explicit negative or positive control",
       "EN-intro: orthogonal -> complementary")

# ---- 7 uncited references: cite the tools where they are used -------------
sub(E, "The OneK1K single-cell eQTL resource [5] and the eQTL Catalogue harmonised release [6] provided the genetic anchor",
       "The OneK1K single-cell eQTL resource [5], the eQTL Catalogue harmonised release [6] and the corresponding bulk blood eQTL compendium [21] provided the genetic anchor",
       "EN-2.1: cite ref 21")

sub(E, "were deconvolved with MCP-counter, EPIC (TRef reference, `otherCells` retained) and CIBERSORT (LM22 signature, `perm = 0`, no quantile normalisation, with EPIC and MCP-counter inputs expressed on the linear scale).",
       "were deconvolved with MCP-counter [9], EPIC (TRef reference, `otherCells` retained) [10] and CIBERSORT (LM22 signature, `perm = 0`, no quantile normalisation, with EPIC and MCP-counter inputs expressed on the linear scale) [11].",
       "EN-2.4: cite refs 9,10,11")

sub(E, "and `G_pur` (ABSOLUTE tumour purity, available for 170 of 514 tumours",
       "and `G_pur` (ABSOLUTE tumour purity [12], available for 170 of 514 tumours",
       "EN-2.4: cite ref 12")

sub(E, "using the consensus subset of the liana resource (1,747 ligand–receptor pairs) with the liana implementation (`rank_aggregate`)",
       "using the consensus subset of the liana resource [13] (1,747 ligand–receptor pairs) with the liana implementation (`rank_aggregate`)",
       "EN-2.6: cite ref 13")

sub(E, "No Mendelian randomization estimate is presented, so no STROBE-MR checklist applies;",
       "No Mendelian randomization estimate is presented, so the STROBE-MR checklist [23] does not apply;",
       "EN-2.14: cite ref 23")

# ---- reconcile 2,442 vs 2,842 cells --------------------------------------
sub(E, "with the naive population used to define the root.",
       "with the naive population used to define the root. This subsample (2,442 cells) is smaller than the CD4+ subset used below for TF activity and for perturbation (2,842 cells): the pseudotime analysis capped each subtype at 500 cells so that the diffusion map stayed tractable, whereas the TF analyses used every CD4+ T cell in the same annotation.",
       "EN-2.7: reconcile cell counts")

# ---- label the two TF-activity contrasts ---------------------------------
sub(E, "Group comparisons used *JAML*-positive versus *JAML*-negative cells and, within positive cells, a median split; the latter was used because zero-inflation makes an all-cell median split degenerate.",
       "Two contrasts were used and both are reported: *JAML*-positive versus *JAML*-negative cells, and a median split restricted to *JAML*-positive cells. The split is restricted to positive cells because zero-inflation (over 70% of these cells score zero) makes an all-cell median split degenerate, returning a single group. Because the two contrasts sit on different scales — the split is taken inside a smaller, more homogeneous group — every value quoted below is labelled with the contrast it comes from.",
       "EN-2.8: label contrasts")

# ---- the noise ceiling is a single order statistic ------------------------
sub(E, "the global noise ceiling defined by the largest randomized reading among all non-degenerate TFs.",
       "the global noise ceiling defined by the largest randomized reading among all non-degenerate TFs. Two properties of that ceiling are stated plainly. It is a single order statistic — here it is set by the largest of five randomized readings for one factor — so it is a conservative screen rather than a fitted null distribution. And with five randomized replicates per factor the smallest attainable one-sided empirical P value is 1/6 ≈ 0.17, so the procedure cannot in principle yield a significant P value: it is reported as a calibrated screen, not as a hypothesis test.",
       "EN-2.9: noise-ceiling caveat")

sub(E, "The affected TFs were nevertheless biologically coherent: within *JAML*-positive cells, higher *JAML* was accompanied by higher activity of the RFX complex (RFXAP +3.36, RFXANK +3.32, RFX5 +3.05 activity units), which drives MHC class II expression, together with *JUN*, *HIF1A*, *IRF6*, *RELA* and *IRF1* — an activation and effector programme rather than a resting one.",
       "The affected TFs were nevertheless biologically coherent, and the same programme appears in both contrasts. Between *JAML*-positive and *JAML*-negative cells the largest activity differences were *RFXANK* (+0.48), *RFXAP* (+0.47), *HIF1A* (+0.44), *JUN* (+0.42), *RFX5* (+0.41) and *IRF6* (+0.41 activity units). Within *JAML*-positive cells, splitting at the median, the same factors separated the two halves, headed by the RFX complex — which drives MHC class II expression — (*RFXAP* +3.36, *RFXANK* +3.32, *RFX5* +3.05). The two contrasts sit on different scales because the split is taken inside a smaller and more homogeneous group, so their values are not comparable with one another.",
       "EN-3.6: name both contrasts")

sub(E, "recovered 5 of the 9 pre-specified activation/exhaustion TFs (PRDM1, NR4A1, BATF, MAF, TBX21), confirming that the pipeline detects a known programme when one exists.",
       "recovered 5 of the 9 pre-specified activation/exhaustion TFs (PRDM1, NR4A1, BATF, MAF, TBX21), confirming that the pipeline detects a known programme when one exists. Two of the nine, *TOX* and *IKZF2*, are absent from the CollecTRI collection used here and could not be scored, so 7 of the 9 were testable; that denominator, not 9, is the one printed on the control panel of the figure.",
       "EN-3.6: positive-control denominator")

sub(E, "and 1.8× the global noise ceiling (9.59 × 10⁻³, defined as the largest randomized reading among non-degenerate TFs).",
       "and 1.8× the global noise ceiling (9.59 × 10⁻³, the largest randomized reading among non-degenerate TFs). That ceiling is a single order statistic, set by one factor's replicate, so it is a conservative screen rather than a fitted null, and with five replicates per factor no empirical P value below 0.17 is attainable.",
       "EN-3.7: noise-ceiling caveat")

sub(E, "The TFs with the largest activity differences in *JAML*-positive cells were the ones with negligible perturbation effects:",
       "Among those 32, the TFs with the largest activity differences between *JAML*-positive and *JAML*-negative cells were the ones with negligible perturbation effects:",
       "EN-3.8: name the contrast")

# ---- "eight analyses" vs 11 table rows -----------------------------------
sub(E, "### Table 1. Summary of the eight analyses and their controls",
       "### Table 1. Summary of the analyses and their controls",
       "EN-Table 1: title")

sub(E, "| Checkpoint blockade | GSE135222; n = 27 | log-rank P = 0.881 | GSE126044 not evaluable |",
       "| Checkpoint blockade | GSE135222; n = 27 | log-rank P = 0.881 | GSE126044 not evaluable |\n\nRows one to eight follow the eight questions listed at the end of the Introduction, in order; the last three rows (independent replication, cross-module concordance, checkpoint blockade) are reported in addition, which is why the table has eleven rows.",
       "EN-Table 1: footnote")

# ---- GSE126044, verified by direct re-count ------------------------------
sub(E, "and only one of 22 genes was evaluable in GSE126044.",
       "and *JAML* itself is absent from the GSE126044 expression matrix, which carries 18,747 genes and does include 17 of the other 21 candidate loci, so that cohort could not be analysed for the gene of interest.",
       "EN-Limitations: GSE126044")

# ---- Figure 6 legend, written to match the panels that exist -------------
para(E, "**Figure 6. Transcription-factor activity in *JAML*-positive CD4+ T cells",
     "**Figure 6. Transcription-factor activity in CD4+ T cells and its two controls.** (a, left) The 20 transcription factors whose inferred activity tracks *JAML* expression most closely, ranked by Spearman ρ between TF activity and *JAML* expression (498 CollecTRI factors scored with decoupleR, univariate linear model; 2,842 CD4+ T cells); red bars are positive ρ and blue bars negative ρ. Activity scores are strongly intercorrelated, so individual identities are not interpretable in isolation. (a, middle) z-scored activity of the 12 factors with the largest activity difference between *JAML*-negative (n = 2,231) and *JAML*-positive (n = 611) cells; these 12 are the top of the 229 factors reaching FDR < 0.05, not the whole significant set. (a, right) Positive control: Δ activity (*PDCD1*-positive minus negative) for the 12 highest-ranked factors of that comparison; the pipeline recovered 5 of the 9 pre-specified activation/exhaustion factors (PRDM1, NR4A1, BATF, MAF, TBX21). *TOX* and *IKZF2* are absent from the CollecTRI collection used, so 7 of the 9 were testable and the panel prints the denominator 7. (b) Number of transcription factors reaching FDR < 0.05 in 200 random permutations of the *JAML* label; observed 229 (red line), permutation mean 0.04 and maximum 2, so the empirical P lies below the resolution of 200 permutations (P < 0.005), and the panel annotation prints only the rounded value 0.0000.",
     "EN-Figure 6: legend")

# ---- Additional files: AF10 restated, AF11 added -------------------------
para(E, "**Additional file 10.** Reviewer-requested supplementary analyses:",
     "**Additional file 10.** Reviewer-requested supplementary analyses: relative effect size of *ID2* perturbation on *JAML* (mean |Δ| as a fraction of mean expression, per-cell distribution, and rank among 3,074 genes); survival results for all 18 evaluable genes in TCGA-LUAD under four parameterizations; the power analysis for the checkpoint-blockade cohort (minimum detectable hazard ratio); and the gene-evaluability check for GSE126044, which records that 17 of the 22 candidate loci are present in that matrix while *JAML*/*AMICA1* is not. (AdditionalFile10a–10d)\n\n**Additional file 11.** Analysis code and derived result tables (archive): the Python and R scripts that generate every number, figure and table in this manuscript, the machine-readable reference-verification record produced by Crossref DOI lookup on 2026-10-01, and the integrity notes on the supplementary data files. (AdditionalFile11_code_and_results.zip)",
     "EN: additional files")

sub(E, "Analysis code is available from the corresponding author on reasonable request.",
       "All analysis scripts, the derived result tables underlying every figure and table, and the machine-readable reference-verification record are provided as Additional file 11; that archive will be deposited in a public repository with a DOI before publication and remains available from the corresponding author on reasonable request in the interim.",
       "EN-availability: code")

# ================================================================ CHINESE ==
sub(C, "围绕单一基因完成八项正交分析，每项均配显式对照：",
       "围绕单一基因完成八项互补分析，每项均配显式对照：",
       "CN-摘要: 正交->互补")

sub(C, "而是用八项正交分析追问单一候选，每项分析均配明确的阴性或阳性对照，且仅使用公共数据。",
       "而是用八项互补分析追问单一候选，每项分析均配明确的阴性或阳性对照，且仅使用公共数据。",
       "CN-引言: 正交->互补")

sub(C, "OneK1K 单细胞 eQTL 资源[5]与 eQTL Catalogue 统一化发布[6]构成本研究的遗传学出发点；",
       "OneK1K 单细胞 eQTL 资源[5]、eQTL Catalogue 统一化发布[6]，以及与之对应的 bulk 血液 eQTL 汇总[21]，构成本研究的遗传学出发点；",
       "CN-2.1: 引文献 21")

sub(C, "分别用 MCP-counter、EPIC（TRef 参考、保留 `otherCells`）与 CIBERSORT（LM22 签名，`perm = 0`、不做分位数归一化；EPIC 与 MCP-counter 输入取线性尺度）反卷积。",
       "分别用 MCP-counter[9]、EPIC（TRef 参考、保留 `otherCells`）[10] 与 CIBERSORT（LM22 签名，`perm = 0`、不做分位数归一化；EPIC 与 MCP-counter 输入取线性尺度）[11] 反卷积。",
       "CN-2.4: 引文献 9,10,11")

sub(C, "`G_pur`（ABSOLUTE 肿瘤纯度，514 例中 170 例有值",
       "`G_pur`（ABSOLUTE 肿瘤纯度[12]，514 例中 170 例有值",
       "CN-2.4: 引文献 12")

sub(C, "在发现集中用 liana 资源共识子集（1,747 对配体–受体）",
       "在发现集中用 liana 资源[13]共识子集（1,747 对配体–受体）",
       "CN-2.6: 引文献 13")

sub(C, "本文不报告孟德尔随机化估计，故不适用 STROBE-MR 清单；",
       "本文不报告孟德尔随机化估计，故不适用 STROBE-MR 清单[23]；",
       "CN-2.14: 引文献 23")

sub(C, "拟合扩散图与扩散拟时序[18,19]，以初始型群体定根。",
       "拟合扩散图与扩散拟时序[18,19]，以初始型群体定根。本项子抽样（2,442 细胞）小于下文 TF 活性与扰动所用的 CD4⁺子集（2,842 细胞）：拟时序分析为控制扩散图的计算量，将每个亚型抽样上限设为 500 细胞，而 TF 分析使用同一注释下的全部 CD4⁺T 细胞。",
       "CN-2.7: 细胞数口径")

sub(C, '分组比较采用"*JAML* 阳性 vs 阴性"；在阳性细胞内部另做中位数切分——后者是必要的，因为零膨胀会使全细胞中位数切分退化为"全部为高"。',
       "共使用两种对比，且均予报告：*JAML* 阳性 vs 阴性；以及在 *JAML* 阳性细胞内部的中位数切分。之所以把切分限制在阳性细胞内，是因为零膨胀（该群 >70% 细胞为 0）会使全细胞中位数切分退化为单一分组。两种对比所用量纲不同——后者在更小、更同质的亚群内比较——故下文凡引用其数值处，均标明所属对比。",
       "CN-2.8: 对比口径")

sub(C, "与**全局噪声上界**（全部非退化 TF 中最大的随机读数）时，才被接受。",
       "与**全局噪声上界**（全部非退化 TF 中最大的随机读数）时，才被接受。有关该上界有两点须明说：它是一个**单一顺序统计量**——此处由某一个 TF 五次随机读数中的最大值决定——因此是保守的筛选门槛，而非拟合出的零分布；又，每个 TF 仅有 5 次随机重复，单侧经验 P 值的最小可能取值为 1/6 ≈ 0.17，故本流程在原理上无法给出显著的 P 值，只能作为校准过的筛选，而非假设检验。",
       "CN-2.9: 噪声上界说明")

sub(C, "受影响的 TF 在生物学上自洽：在 *JAML* 阳性细胞内部，*JAML* 越高，RFX 复合体活性越高（RFXAP +3.36、RFXANK +3.32、RFX5 +3.05 活性单位）——该复合体驱动 MHC II 类表达——并伴随 *JUN*、*HIF1A*、*IRF6*、*RELA*、*IRF1*，即一套活化/效应程序而非静息程序。",
       '受影响的 TF 在生物学上自洽，且两种对比下出现的是同一套程序。在"*JAML* 阳性 vs 阴性"对比中，活性差最大者为 *RFXANK*（+0.48）、*RFXAP*（+0.47）、*HIF1A*（+0.44）、*JUN*（+0.42）、*RFX5*（+0.41）、*IRF6*（+0.41 活性单位）。在 *JAML* 阳性细胞内部按中位数切分，分离两半的仍是同一批因子，其中驱动 MHC II 类表达的 RFX 复合体居首（*RFXAP* +3.36、*RFXANK* +3.32、*RFX5* +3.05）。两种对比量纲不同（后者在更小、更同质的亚群内比较），其数值彼此不可比。',
       "CN-3.6: 两种对比")

sub(C, "同一流程用于 *PDCD1* 时命中了 9 个预设活化/耗竭 TF 中的 5 个（PRDM1、NR4A1、BATF、MAF、TBX21），说明该流程在存在已知程序时能够检出。",
       "同一流程用于 *PDCD1* 时命中了 9 个预设活化/耗竭 TF 中的 5 个（PRDM1、NR4A1、BATF、MAF、TBX21），说明该流程在存在已知程序时能够检出。9 个预设 TF 中有 2 个（*TOX*、*IKZF2*）不在本研究所用的 CollecTRI 集合内、无法评分，故实际可检验者为 7 个；图中阳性对照面板标注的分母是 7 而非 9。",
       "CN-3.6: 阳性对照分母")

sub(C, "且为全局噪声上界的 1.8 倍（9.59×10⁻³，取全部非退化 TF 中最大的随机读数）。",
       "且为全局噪声上界的 1.8 倍（9.59×10⁻³，即全部非退化 TF 中最大的随机读数）。该上界是单一顺序统计量、由某一个 TF 的随机重复决定，属保守筛选门槛而非拟合零分布；每 TF 仅 5 次重复，无法得到低于 0.17 的经验 P 值。",
       "CN-3.7: 噪声上界说明")

sub(C, "活性差最大的 TF 恰恰是扰动效应可忽略的那些：",
       '在这 32 个 TF 中，"阳性 vs 阴性"活性差最大的那些，恰恰是扰动效应可忽略的：',
       "CN-3.8: 对比口径")

sub(C, "### 表 1. 八项分析及其对照一览",
       "### 表 1. 分析及其对照一览",
       "CN-表 1: 标题")

sub(C, "| 检查点抑制剂 | GSE135222；n = 27 | log-rank P = 0.881 | GSE126044 不可评估 |",
       "| 检查点抑制剂 | GSE135222；n = 27 | log-rank P = 0.881 | GSE126044 不可评估 |\n\n前 8 行依次对应引言末尾列出的 8 个问题；后 3 行（独立复核、模块间一致性、检查点抑制剂）为额外报告，故本表共 11 行。",
       "CN-表 1: 脚注")

sub(C, "且 GSE126044 中 22 个基因仅有 0 个可评估。",
       "且 *JAML* 本身不在 GSE126044 的表达矩阵中（该矩阵含 18,747 个基因，22 个候选位点中另外 21 个里有 17 个在位），故该队列无法就目标基因进行分析。",
       "CN-局限: GSE126044")

para(C, "**图 6. *JAML* 阳性 CD4⁺T 细胞的转录因子活性",
     "**图 6. CD4⁺T 细胞的转录因子活性及其两项对照。**（a 左）与 *JAML* 表达关联最密切的 20 个转录因子，按 TF 活性与 *JAML* 表达的 Spearman ρ 排序（498 个 CollecTRI 因子，decoupleR 单变量线性模型，2,842 个 CD4⁺T 细胞）；红柱为 ρ > 0，蓝柱为 ρ < 0。活性分数彼此高度相关，单个因子的身份不可孤立解读。（a 中）在 *JAML* 阴性（n = 2,231）与阳性（n = 611）细胞之间活性差最大的 12 个因子的 z 标准化活性；这 12 个是 229 个达 FDR < 0.05 因子中的前列，而非全部显著集合。（a 右）阳性对照：该对比中排名前 12 的因子的 Δ 活性（*PDCD1* 阳性减阴性）。流程命中了 9 个预设活化/耗竭因子中的 5 个（PRDM1、NR4A1、BATF、MAF、TBX21）；*TOX* 与 *IKZF2* 不在所用 CollecTRI 集合内，故 9 个中有 7 个可检验，面板标注的分母为 7。（b）对 *JAML* 标签做 200 次随机置换后达到 FDR < 0.05 的转录因子数目分布；实测 229（红线），置换均值 0.04、最大 2，故经验 P 低于 200 次置换的分辨率（P < 0.005）；面板标注为取整后的数值 0.0000。",
     "CN-图 6: 图注")

para(C, "**附加文件 10.** 应审稿人要求补充的分析：",
     "**附加文件 10.** 应审稿人要求补充的分析：*ID2* 扰动对 *JAML* 的相对效应量（平均 |Δ| 占该基因平均表达的比例、逐细胞分布、在 3,074 个基因中的排名）；TCGA-LUAD 全部 18 个可评估基因在四种参数化下的生存结果；检查点抑制剂队列的功效分析（最小可检出 HR）；以及 GSE126044 的基因可评估性核查（该矩阵含 22 个候选位点中的 17 个，但不含 *JAML*/*AMICA1*）。（AdditionalFile10a–10d）\n\n**附加文件 11.** 分析代码与派生结果表（压缩包）：生成本文每一个数字、每一张图与每一张表的 Python 与 R 脚本，2026-10-01 通过 Crossref DOI 反查产生的机读核验记录，以及补充数据文件的完整性说明。（AdditionalFile11_code_and_results.zip）",
     "CN: 附加文件")

sub(C, "分析代码可向通讯作者索取。",
       "全部分析脚本、支撑每张图与每张表的派生结果表，以及机读的参考文献核验记录，均作为附加文件 11 提供；同一压缩包将在投稿前存入带 DOI 的公共仓库，在此之前可向通讯作者索取。",
       "CN-数据可用性: 代码")

# ---- CN Discussion: repair the find/replace artefact ---------------------
sub(C, '仍属可行，而配体–受体分析与"髓系和树突状细胞是该界面的主要发送方。',
       "仍属可行；配体–受体分析也提示髓系与树突状细胞是该界面的主要发送方。",
       "CN-讨论: 修复破损句(前)")

sub(C, '这一外部注释。"这一判断一致。',
       "这一外部注释。”",
       "CN-讨论: 修复破损句(尾)")


def main():
    stamp = time.strftime("%Y%m%d_%H%M%S")
    for p in (EN, CN):
        shutil.copy2(p, p + ".bak_review2_" + stamp)
    print("backups: *.bak_review2_" + stamp)

    texts = {EN: open(EN, encoding="utf-8").read(),
             CN: open(CN, encoding="utf-8").read()}

    problems, ok = [], 0
    for kind, path, a, b, label in OPS:
        t = texts[path]
        if kind == "sub":
            n = t.count(a)
            if n != 1:
                problems.append("[%s] substring hit %d time(s), expected 1 :: %r" % (label, n, a[:70]))
                continue
            texts[path] = t.replace(a, b, 1)
        else:
            n = t.count(a)
            if n != 1:
                problems.append("[%s] paragraph anchor hit %d time(s), expected 1 :: %r" % (label, n, a[:70]))
                continue
            i = t.index(a)
            j = t.find("\n\n", i)
            if j < 0:
                problems.append("[%s] no paragraph terminator after anchor" % label)
                continue
            texts[path] = t[:i] + b + t[j:]
        ok += 1
        print("  OK  " + label)

    if problems:
        print("\nABORTED, nothing written. Unmatched edits:")
        for p in problems:
            print("   " + p)
        sys.exit(1)

    for p, t in texts.items():
        with open(p, "w", encoding="utf-8") as f:
            f.write(t)
    print("\napplied %d edits." % ok)

    for p, t in texts.items():
        body = t.split("## References")[0].split("## 参考文献")[0]
        cited = collections.Counter()
        for m in re.finditer(r"\[([0-9][0-9,\u2013\- ]*)\]", body):
            for part in m.group(1).split(","):
                part = part.strip().replace("\u2013", "-")
                if "-" in part:
                    try:
                        lo, hi = part.split("-")
                        for i in range(int(lo), int(hi) + 1):
                            cited[i] += 1
                    except Exception:
                        pass
                elif part.isdigit():
                    cited[int(part)] += 1
        print("  %s -> uncited refs: %s" % (os.path.basename(p),
                                            [i for i in range(1, 24) if i not in cited]))
        for bad in ("八项正交", "eight orthogonal", "这一判断一致", "而配体–受体分析与"):
            if bad in t:
                print("  !! residual fragment still present: %r" % bad)


if __name__ == "__main__":
    main()
