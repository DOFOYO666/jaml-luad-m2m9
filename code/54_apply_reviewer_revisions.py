# -*- coding: utf-8 -*-
"""按审稿意见修改 M2-M9 中英双稿（可重跑、逐条校验）。

支持四种操作：
  line  —— 替换"以 anchor 开头"的整行
  sub   —— 精确子串替换
  block —— 从以 anchor 开头的行起，替换连续的若干非空行（整块）
  after / before —— 在 anchor 行之后/之前插入若干行

任一操作未命中即报 MISS 并计入失败数（退出码非 0），避免静默漏改。

用法：python 54_apply_reviewer_revisions.py
"""
import io
import os
import sys

PKG = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化\M2-M9深化研究稿"
EN = "manuscript_EN.md"
CN = "manuscript_CN.md"

OPS = []          # (file, mode, anchor, new)


def op(f, mode, anchor, new):
    OPS.append((f, mode, anchor, new))


# =====================================================================
# 一、摘要（M5：压缩至 ≤350 词，按"一问一答"重构）
# =====================================================================
op(EN, "line", "**Background.** Germline variation at",
   "**Background.** Germline variation at *JAML* has been linked to lung adenocarcinoma (LUAD) risk "
   "through CD4+ T-cell expression, yet where the molecule acts, which cells determine its tissue "
   "abundance, and whether its expression is transcriptionally controlled remain unknown — the "
   "minimum needed before it can be treated as a target.")

op(EN, "line", "**Methods.** We assembled eight",
   "**Methods.** In a secondary analysis of public data we interrogated one gene with eight orthogonal "
   "analyses, each with an explicit control: single-cell expression in LUAD (GSE131907; 208,506 cells) "
   "replicated in NSCLC (GSE127465; 54,773 cells) with verified matrix–metadata alignment; three "
   "deconvolution methods with partial-correlation adjustment validated against lineage markers; "
   "TCGA-LUAD survival under four parameterizations; ligand–receptor inference with permutation; "
   "CD4-lineage pseudotime; transcription-factor (TF) activity with permutation and positive controls; "
   "*in silico* TF perturbation with a paired randomized-network control calibrated to a global noise "
   "ceiling; and DepMap dependency data.")

op(EN, "line", "**Results.** *JAML* and its receptor",
   "**Results.** *JAML* and *CXADR* occupied separate compartments in both datasets: *JAML* positivity "
   "was highest in myeloid/dendritic cells (47.0%, 30.1%) and rare in epithelium (1.7%, 1.5%), the "
   "epithelial-to-myeloid ratio being 0.04 and 0.05, whereas *CXADR* was epithelial (40.5%, 19.3%). "
   "Alignment was confirmed independently (per-cell mitochondrial fraction ρ = 1.0000; 5/5 marker "
   "checks; 102/102 counts reproduced exactly). *JAML* abundance tracked myeloid/DC content "
   "(ρ = 0.757) — the only association robust to adjustment for immune content and tumour purity "
   "(0.750, 0.738) — and its correlation profile resembled myeloid/DC genes rather than T-cell genes "
   "(*CLEC10A* 0.973 versus *CD8A* 0.623). Among 39 perturbed TFs, only *ID2* exceeded both its own "
   "randomized control (6.6×) and the global noise ceiling (1.8×), reproduced across three network "
   "configurations, and ranked 33rd of 3,074 genes by relative effect (6.3% of mean *JAML* expression); "
   "*GATA3* was suggestive, while *MAF*, *STAT1*, *IRF1* and *STAT3* were indistinguishable from noise. "
   "TF-activity inference flagged 229 TFs but showed no concordance with perturbation (r = −0.002). "
   "*JAML* was neither a LUAD cell-line dependency (0/50 below threshold) nor a consistent prognostic "
   "marker (categorical Cox FDR = 0.014; continuous Cox, the pre-specified model, FDR = 0.11; no gene "
   "significant under both), and was unrelated to progression-free survival under checkpoint blockade "
   "(n = 27; minimum detectable hazard ratio 3.4).")

op(EN, "line", "**Conclusions.** *JAML* is an immune-side molecule",
   "**Conclusions.** *JAML* is an immune-side molecule whose tissue abundance is dominated by myeloid "
   "and dendritic cells, whose counterpart in this setting is the epithelial protein *CXADR*, and whose "
   "expression is at best weakly controlled, with *ID2* the only candidate regulator surviving noise "
   "calibration. These results define the compartment and experimental system in which *JAML* biology "
   "should be tested, and they argue against pursuing it as a tumour-cell-intrinsic dependency.")

op(CN, "line", "**背景。**",
   "**背景。** *JAML* 位点的胚系变异通过 CD4⁺T 细胞表达与肺腺癌（LUAD）风险相关联，但该分子在何处发挥作用、"
   "哪些细胞决定其在肿瘤组织中的丰度、其表达是否受转录调控，三者均未确定——而这正是把它当作靶点之前"
   "至少需要回答的问题。")

op(CN, "line", "**方法。**",
   "**方法。** 本文是对公共数据的二次分析，围绕单一基因完成八项正交分析，每项均配显式对照："
   "LUAD 单细胞表达（GSE131907，208,506 细胞）及其在 NSCLC 中的独立复核（GSE127465，54,773 细胞，"
   "矩阵–元数据对齐经核验）；三种反卷积方法结合偏相关校正（校正流程以谱系 marker 基因验证）；"
   "TCGA-LUAD 生存分析在四种参数化下报告；配体–受体推断（含置换对照）；CD4 谱系拟时序；"
   "转录因子（TF）活性（含置换与阳性对照）；*in silico* TF 扰动（配对打乱网络对照 + 全局噪声上界校准）；"
   "以及 DepMap 依赖性数据。")

op(CN, "line", "**结果。**",
   "**结果。** 两套数据中 *JAML* 与 *CXADR* 均处于不同区室：*JAML* 阳性率在髓系/树突状细胞最高"
   "（47.0%、30.1%），上皮中罕见（1.7%、1.5%），上皮/髓系比值 0.04 与 0.05；*CXADR* 位于上皮"
   "（40.5%、19.3%）。对齐经独立核验（逐细胞线粒体读段比例 ρ = 1.0000；5/5 项 marker 检查；"
   "102/102 个计数逐格复现）。*JAML* 丰度跟随髓系/DC（ρ = 0.757），是唯一对整体免疫含量与肿瘤纯度"
   "校正均稳健的关联（0.750、0.738）；其相关谱更像髓系/DC 基因（*CLEC10A* 0.973）而非 T 细胞基因"
   "（*CD8A* 0.623）。在 39 个被扰动的 TF 中，只有 *ID2* 同时超过其自身随机对照（6.6 倍）与全局噪声"
   "上界（1.8 倍），并在三套网络配置中复现；其相对效应为 *JAML* 平均表达的 6.3%，在 3,074 个基因中"
   "排第 33 位。*GATA3* 仅提示性；*MAF*、*STAT1*、*IRF1*、*STAT3* 与噪声不可区分。TF 活性推断标出"
   "229 个 TF，但与扰动无一致性（r = −0.002）。*JAML* 既非 LUAD 细胞系依赖基因（0/50 低于阈值），"
   "也不是稳健的预后标志物（分类参数化 Cox FDR = 0.014；预设的连续 Cox FDR = 0.11；无基因在两种参数化下"
   "均显著），且与检查点抑制剂治疗下的无进展生存无关（n = 27；最小可检出 HR 3.4）。")

op(CN, "line", "**结论。**",
   "**结论。** *JAML* 是免疫侧分子，其组织丰度由髓系与树突状细胞主导，其在本文情境下的对应分子是上皮蛋白"
   "*CXADR*，其表达至多受微弱调控，*ID2* 是唯一通过噪声校准的候选调控因子。这些结果界定了检验 *JAML* "
   "生物学作用时应选用的区室与实验体系，也说明不应把 *JAML* 作为肿瘤细胞内在依赖来推进。")

# =====================================================================
# 二、Methods
# =====================================================================
op(EN, "sub", "do not re-estimate causal effects.",
   "do not re-estimate causal effects. This is a secondary analysis of published data; no new causal "
   "estimate is reported here.")
op(CN, "sub", "不重新估计因果效应。",
   "不重新估计因果效应。**本文为已发表数据的二次分析，不报告新的因果估计。**")

op(EN, "line", "Overall survival was analysed in TCGA-LUAD",
   "Overall survival was analysed in TCGA-LUAD (n = 563; 208 deaths). The pre-specified primary model "
   "was Cox proportional-hazards regression on continuous expression standardized per SD, with "
   "Benjamini–Hochberg correction across the genes actually tested; three secondary parameterizations "
   "were also run on the identical sample and gene set (Cox per 1 log2 unit; Cox on a median split, "
   "i.e. categorical; and the Mantel–Cox log-rank test on the same median split), because "
   "inconsistency between parameterizations would preclude any prognostic claim. Stage-wise expression "
   "was summarized by Kruskal–Wallis test and a rank trend test. Analyses are descriptive and "
   "non-causal.")

op(CN, "line", "在 TCGA-LUAD（n = 563",
   "在 TCGA-LUAD（n = 563；208 例死亡）分析总生存。**预设的主要模型**为连续表达（每 SD 标准化）的 "
   "Cox 比例风险回归，Benjamini–Hochberg 校正施于实际检验的基因数；另在同一批样本与同一基因集上运行"
   "三种次要参数化（连续每 1 个 log2 单位的 Cox；中位数分组即分类变量的 Cox；以及同一分组下的 "
   "Mantel–Cox log-rank 检验）——因为若不同参数化之间不一致，就不应作出任何预后主张。"
   "分期表达以 Kruskal–Wallis 检验与秩趋势检验汇总。分析为描述性、非因果。")

op(EN, "line", "Because *JAML* is an adhesion molecule",
   "Because *JAML* is an adhesion molecule and not a transcription factor, it cannot itself be "
   "perturbed in a TF-based gene regulatory network model. We therefore inverted the question and "
   "perturbed candidate TFs, reading out the change in *JAML*. CellOracle (v0.20.0) was run [16] on "
   "the CD4+ subset (2,842 cells × 3,013 genes, extended to 3,074 genes in a final configuration) "
   "using the human promoter base GRN (hg38, gimmemotifs v5, FPR 2%, score threshold 10; 1,094 TF "
   "columns), α = 10, k = 30 imputation neighbours, 50 principal components, and the five CD4+ "
   "subtypes as network units with `bagging = 20`. Because *JAML* has no regulatory connection "
   "reaching the tutorial threshold (its 125 candidate edges had a median P of 0.28 and none below "
   "P = 0.001), the network filter was relaxed to P < 0.05 with the top 20,000 edges retained per "
   "unit; this relaxation is stated explicitly because it was necessary to make *JAML* evaluable at "
   "all, and it is important that it was not the criterion of evidence: after relaxation each unit "
   "retained roughly 12,000–20,000 edges (about 15–24% of the 83,545 candidate edges per unit), and "
   "the final judgement rested on the paired randomized control and the global noise ceiling rather "
   "than on any P value. Effect sizes are reported both on the count scale and relative to the "
   "target's own expression, as the ratio of mean |Δ*JAML*| to mean *JAML* expression in the same "
   "cells, together with the rank of *JAML* among all 3,074 genes by that ratio, so that the "
   "magnitude can be judged independently of the model's units.")

op(CN, "line", "由于 *JAML* 是黏附分子而非转录因子",
   "由于 *JAML* 是黏附分子而非转录因子，它本身无法在基于 TF 的基因调控网络模型中被敲除。因此我们反转问题："
   "扰动候选 TF，读出 *JAML* 的变化。用 CellOracle（v0.20.0）[16]在 CD4⁺子集（2,842 细胞 × 3,013 基因，"
   "最终配置扩展至 3,074 基因）上运行，使用人 promoter base GRN（hg38、gimmemotifs v5、FPR 2%、"
   "分数阈值 10；1,094 个 TF 列），α = 10，k = 30 个插补近邻，50 个主成分，以 5 个 CD4⁺亚型为网络单元，"
   "`bagging = 20`。由于 *JAML* 没有任何调控连接达到教程默认阈值（其 125 条候选边的 P 值中位数为 0.28，"
   "无一条低于 P = 0.001），网络过滤放宽至 P < 0.05 且每单元保留前 20,000 条边；**这一放宽被显式声明**，"
   "因为它是使 *JAML* 可被评估的前提；更重要的是，它**不是**证据判据：放宽后每单元保留约 12,000–20,000 条边"
   "（约占每单元 83,545 条候选边的 15–24%），最终判断依据的是配对随机对照与全局噪声上界，而非任何 P 值。"
   "效应量同时以 count 尺度与相对目标基因自身表达的比值报告（即平均 |Δ*JAML*| 与同一批细胞 *JAML* 平均表达之比），"
   "并给出 *JAML* 在该比值上于全部 3,074 个基因中的排名，使量级可不依赖模型单位来判断。")

op(EN, "line", "In GSE135222 (n = 27 NSCLC patients",
   "In GSE135222 (n = 27 NSCLC patients treated with anti-PD-1/PD-L1 therapy; 21 progression events), "
   "*JAML* expression above the median was not associated with progression-free survival (log-rank "
   "χ² = 0.022; P = 0.881; median PFS 65 versus 51 days; Cox HR = 0.94, 95% CI 0.40–2.22). This "
   "cohort is small, and its negative result must be read with its power: at 80% power and α = 0.05 "
   "(two-sided) with 21 events, the minimum detectable hazard ratio is approximately 3.4, so the "
   "analysis excludes only a strong association and has limited power for hazard ratios between 1.0 "
   "and 3.4. A second cohort, GSE126044, was not evaluable because *JAML* is absent from its "
   "expression matrix.")

op(CN, "line", "GSE135222（27 例",
   "GSE135222（27 例接受抗 PD-1/PD-L1 治疗的 NSCLC 患者；21 例进展事件）中，*JAML* 表达高于中位数与"
   "无进展生存无关（log-rank χ² = 0.022；P = 0.881；中位 PFS 65 天 vs 51 天；Cox HR = 0.94，"
   "95% CI 0.40–2.22）。该队列样本量小，其阴性结果必须结合功效解读：在 21 个事件、80% 功效、"
   "α = 0.05（双侧）下，**最小可检出 HR 约为 3.4**，故本分析只能排除强关联，对 HR 在 1.0–3.4 之间的"
   "效应功效不足。另一队列 GSE126044 因表达矩阵中不含 *JAML* 而不可评估。")

op(EN, "after", "All P values are two-sided",
   "\n### 2.14 Reporting and reproducibility\n\n"
   "Because several analyses here rest on statistical models rather than on direct measurement, three "
   "reporting rules were fixed in advance. First, one primary parameterization was pre-specified for "
   "each family of tests (continuous per-SD Cox for survival, Spearman correlation for composition "
   "analyses, the paired randomized control for perturbation) and all secondary parameterizations are "
   "reported alongside it rather than selected from. Second, every effect is reported together with "
   "its control value on the same scale; a P value without its control is not treated as evidence. "
   "Third, since several analyses have sample sizes in the hundreds to thousands, effect sizes rather "
   "than P values carry the interpretation, and no claim is made from a P value alone. Multiplicity "
   "was controlled by Benjamini–Hochberg [20] within each analysis family, and the nominal families "
   "are stated in the corresponding Methods subsection. No Mendelian randomization estimate is "
   "presented, so no STROBE-MR checklist applies; the clinical and expression analyses are reported "
   "with sample sizes, event counts, exclusions and the number of genes tested.")

op(CN, "after", "所有 P 值为双侧",
   "\n### 2.14 报告规范与可复现性\n\n"
   "由于本文多项分析依赖统计模型而非直接测量，我们预先固定三条报告规则。第一，每一族检验**预设一个主要参数化**"
   "（生存分析为连续每 SD 的 Cox，构成分析为 Spearman 相关，扰动分析为配对随机对照），其余参数化与主要参数化"
   "一并报告，而非从中挑选。第二，每个效应量必须与同一尺度上的对照值同报；没有对照值的 P 值不作为证据。"
   "第三，本文多项分析的样本量在数百至数千，因此由**效应量**承担解读，不得仅凭 P 值作结论。"
   "多重比较在每个分析族内用 Benjamini–Hochberg 校正[20]，各族的具体范围在相应方法小节中写明。"
   "本文不报告孟德尔随机化估计，故不适用 STROBE-MR 清单；临床与表达分析均报告样本量、事件数、排除标准"
   "与受检基因数。")

# =====================================================================
# 三、Results
# =====================================================================
op(EN, "line", "Eighteen of the 22 candidate genes were present",
   "Eighteen of the 22 candidate genes were present in the TCGA expression matrix. Because the choice "
   "of survival parameterization turned out to determine the answer, we report four parameterizations "
   "on the identical sample and gene set (Table 5). Under the pre-specified continuous model, a higher "
   "*JAML* level was associated with longer overall survival (HR per SD 0.85; P = 0.017), but the "
   "association did not survive correction (FDR = 0.11); the only gene passing FDR under this "
   "parameterization was *MAP4K4* (FDR = 0.018). Under a categorical parameterization (median split) "
   "the picture inverted: *JAML* was significant (HR = 0.62; P = 0.0008; FDR = 0.014), matching the "
   "log-rank result (χ² = 11.55; P = 6.8 × 10⁻⁴; FDR = 0.012), while *MAP4K4* was not (FDR = 0.17). "
   "**No gene was significant under both Cox parameterizations.** The discrepancy is therefore a "
   "property of the parameterization rather than of the gene: *JAML*'s association with survival is "
   "confined to the lower half of its expression range and is not log-linear, so a per-SD continuous "
   "term dilutes it. The defensible conclusion is that no candidate here is a robust prognostic "
   "marker. *JAML* expression also decreased monotonically across stage (means 8.56, 8.32, 8.18, 7.75 "
   "for stages I–IV; Kruskal–Wallis P = 0.0020; trend ρ = −0.171). All of these associations are "
   "observational and non-causal, and their direction need not match that of a germline instrument.")

op(CN, "line", "22 个候选基因中有 18 个",
   "22 个候选基因中有 18 个存在于 TCGA 表达矩阵。由于生存分析的参数化选择会直接决定答案，我们在同一批样本"
   "与同一基因集上报告四种参数化（表 5）。在预设的连续模型下，*JAML* 高表达与更长总生存相关"
   "（每 SD 的 HR 0.85；P = 0.017），但未通过校正（FDR = 0.11）；该参数化下唯一通过 FDR 的基因是 "
   "*MAP4K4*（FDR = 0.018）。而在分类参数化（中位数分组）下结论反转：*JAML* 显著（HR = 0.62；"
   "P = 0.0008；FDR = 0.014），与 log-rank 结果一致（χ² = 11.55；P = 6.8×10⁻⁴；FDR = 0.012），"
   "而 *MAP4K4* 不显著（FDR = 0.17）。**没有任何基因在两种 Cox 参数化下均显著。**因此这一分歧属于"
   "**参数化的属性而非基因的属性**：*JAML* 与生存的关联集中于其表达量的下半区、并非对数线性，"
   "以每 SD 连续项建模会稀释该效应。可辩护的结论是：本文没有哪个候选是稳健的预后标志物。"
   "*JAML* 表达还随分期单调下降（I–IV 期均值 8.56、8.32、8.18、7.75；Kruskal–Wallis P = 0.0020；"
   "趋势 ρ = −0.171）。以上关联均为观察性、非因果，其方向亦不必与胚系工具变量的方向一致。")

op(EN, "sub", "Two further checks supported *ID2*.",
   "The magnitude should be read in relative rather than absolute terms. Mean |Δ*JAML*| equalled 6.31% "
   "of mean *JAML* expression in the same cells (1.36% of the mean of *JAML*-positive cells); the "
   "change was non-zero in 918 of 2,842 cells (32%), with a cell-level 95th percentile of 0.086 and a "
   "maximum of 0.242, so the effect is concentrated rather than uniform. Ranked by relative effect "
   "across all 3,074 genes, *JAML* placed 33rd — notable, but far behind *ID2* itself, whose "
   "self-perturbation moved its own expression by 45.6%. The *ID2*→*JAML* edge is therefore a "
   "moderate-strength edge in this network, not a dominant one. Two further checks supported *ID2*.")
op(CN, "sub", "两项进一步检验支持 *ID2*。",
   "量级应以相对而非绝对值解读。平均 |Δ*JAML*| 等于同一批细胞 *JAML* 平均表达的 6.31%"
   "（为 *JAML* 阳性细胞均值的 1.36%）；变化在 2,842 个细胞中的 918 个（32%）非零，逐细胞第 95 百分位为 "
   "0.086、最大 0.242，说明效应是**集中**而非均匀的。在全部 3,074 个基因中按相对效应排序，*JAML* 位列第 33 —— "
   "值得注意，但远低于 *ID2* 扰动自身时的效应（45.6%）。因此 *ID2*→*JAML* 在该网络中属**中等强度边**，"
   "而非主导边。两项进一步检验支持 *ID2*。")

op(EN, "line", "### 3.9 *JAML* is not a dependency of LUAD cells, and its blockade biomarker value is unproven",
   "### 3.9 *JAML* is not a dependency of LUAD cells")
op(CN, "line", "### 3.9 *JAML* 不是 LUAD 细胞的依赖基因，其作为阻断疗效标志物的价值亦未确立",
   "### 3.9 *JAML* 不是 LUAD 细胞的依赖基因")

op(EN, "before", "In GSE135222 (n = 27 NSCLC patients",
   "### 3.10 Clinical context: two negative results\n")
op(CN, "before", "GSE135222（27 例",
   "### 3.10 临床语境：两项阴性结果\n")

# =====================================================================
# 四、图注（M6：补完整图注）
# =====================================================================
op(EN, "block", "**Figure 1.**",
   "**Figure 1. Compartmental expression of *JAML* and *CXADR* in LUAD tissue.** Mean normalized "
   "expression (log2(TPM+1)) and positivity for *JAML* and *CXADR* across ten author-annotated cell "
   "types in GSE131907 (208,506 cells; LUAD). Points are cell-type means; the x-axis is the fraction "
   "of positive cells. Oligodendrocytes (n = 716) originate from brain metastases and are shown but "
   "excluded from primary-tumour statements. Panel b: *JAML* positivity across CD4-lineage subsets "
   "annotated by the authors (`Cell_subtype`). No statistical test is applied; these are descriptive "
   "summaries of a complete dataset.\n\n"
   "**Figure 2. Independent replication of the compartmental separation (GSE127465).** Positivity of "
   "*JAML* and *CXADR* across eight compartments in GSE127465 (54,773 cells; treatment-naïve NSCLC; "
   "Zilionis et al. [8]) compared with GSE131907. Cells were grouped from the author-annotated major "
   "and minor cell types; positivity is expression > 0 in the GEO-supplied normalised count matrix. "
   "Row correspondence between matrix and metadata was verified before grouping (Additional file 1). "
   "Horizontal bars are positivity percentages; the discovery and replication values for the same "
   "compartment are shown side by side.\n\n"
   "**Figure 3. *JAML* expression versus immune populations in TCGA-LUAD: raw and adjusted "
   "correlations.** Spearman ρ between *JAML* expression and each deconvolved population in n = 514 "
   "primary tumours, computed with MCP-counter, EPIC and CIBERSORT. Filled circles are raw ρ; open "
   "markers are partial ρ controlling for global immune content (1 − EPIC `otherCells`) or, where "
   "available, tumour purity (ABSOLUTE; n = 170). Bars connect the raw and immune-adjusted values. "
   "Populations are ordered by raw ρ. All correlations shown had raw FDR < 0.05.\n\n"
   "**Figure 4. Validity of the adjustment and similarity of correlation profiles.** (a) Mean "
   "Spearman ρ between true lineage marker genes (*CD3D*, *CD3E*, *CD8A* for T cells; *LYZ* for "
   "myeloid cells) and T-cell populations, before adjustment (grey bars) and after adjustment with "
   "the structurally independent leukocyte-fraction proxy (red diamonds) or with the composite "
   "MCP-counter score (blue). The composite control removes the lineage signal of genuine T-cell "
   "markers and is therefore not a valid adjustment; the independent proxy preserves it. (b) Pairwise "
   "Pearson similarity of 39-population correlation profiles between *JAML* and 15 marker genes; "
   "*JAML* most resembles myeloid/DC genes (*CLEC10A*, *LYZ*) and is opposite in sign to epithelial "
   "markers (*EPCAM*, *KRT18*).\n\n"
   "**Figure 5. Pseudotime within the CD4+ lineage and the *JAML* gradient.** Diffusion pseudotime "
   "(50 principal components on 2,006 genes after filtering; 2,442 cells sampled to at most 500 per "
   "subtype) with the naive population as root. (a) Pseudotime by subtype (median and interquartile "
   "range). (b) *JAML* positivity by subtype. (c) Scaled *JAML* expression against pseudotime with "
   "the fraction of zero values indicated, because more than 75% of CD4+ T cells were *JAML*-negative "
   "and the correlation (ρ = 0.128, P = 2.6 × 10⁻¹⁰) is driven by the positive minority. The gradient "
   "is supported by this dataset only.\n\n"
   "**Figure 6. Transcription-factor activity in *JAML*-positive CD4+ T cells and its negative "
   "control.** (a) Activity of 498 CollecTRI transcription factors (decoupleR, univariate linear "
   "model) compared between *JAML*-positive and *JAML*-negative cells; the top-ranked factors by "
   "activity difference are shown. Because activity scores are intercorrelated, individual identities "
   "are not interpretable and the figure is annotated accordingly. (b) Distribution of the number of "
   "transcription factors reaching FDR < 0.05 in 200 random permutations of the *JAML* label (observed "
   "value 229 shown as a dashed line); the permutation mean is 0.04 and the maximum is 2. The positive "
   "control (*PDCD1*, identical pipeline) recovered 5 of 9 pre-specified activation/exhaustion factors.\n\n"
   "**Figure 7. Perturbation of candidate transcription factors and readout of *JAML*.** (a) Mean "
   "|Δ*JAML*| after knocking out each of 39 candidate transcription factors in the CD4+ subset "
   "(2,842 cells) using CellOracle with bagging = 20; the dashed line is the global noise ceiling "
   "(9.59 × 10⁻³), defined as the largest randomized reading among non-degenerate factors. (b) Paired "
   "randomized-network control: for each factor, the real effect (red bars) against the mean and "
   "standard deviation of five freshly seeded randomized replicates (open markers). Factors whose "
   "randomized baseline collapsed below 1 × 10⁻⁴ are flagged as degenerate (shaded region); ratios "
   "derived from a degenerate denominator are not interpreted. Only *ID2* exceeded both its own "
   "randomized maximum and the global ceiling.\n\n"
   "**Figure 8. Reproducibility across network configurations and absence of concordance between "
   "association and intervention.** (a) Perturbation effect on *JAML* for the 32 transcription factors "
   "evaluable in both the perturbation model and the activity analysis, plotted against the activity "
   "difference; the dashed line is y = 0. No relationship is present (Pearson r = −0.002; sign "
   "agreement 12/32). (b) Effect sizes for candidate factors across three network configurations "
   "(bagging 3 with 3,013 genes; bagging 20 with 3,013 genes; bagging 20 with 3,074 genes), showing "
   "the rank stability of *ID2* and the position of the 60 transcription factors newly included in the "
   "final configuration relative to the noise ceiling.\n\n"
   "**Figure 9. *JAML* gene effect across DepMap 22Q2 cell lines.** Violin and box plots of CRISPR "
   "(Chronos) gene effect for *JAML* in LUAD lines (`lineage_sub_subtype == NSCLC_adenocarcinoma`; "
   "n = 50) versus all other lineages (n = 1,036); points are individual cell lines. The dashed line "
   "marks the conventional dependency threshold (−0.5); no LUAD line falls below it. RPL5 (−2.25) and "
   "CDK1 (−2.04) are shown in the text as positive controls for the measurement range, and KRAS "
   "(−0.40), a driver rather than a common dependency in this lineage, as a negative control.")

op(CN, "block", "**图 1.**",
   "**图 1. LUAD 组织中 *JAML* 与 *CXADR* 的区室表达。** GSE131907（208,506 细胞；LUAD）中 *JAML* 与 "
   "*CXADR* 在 10 个作者注释细胞类型上的平均归一化表达（log2(TPM+1)）与阳性率。点为细胞类型均值，横轴为"
   "阳性细胞比例。少突胶质细胞（n = 716）来自脑转移灶，图中显示但在原发灶相关表述中排除。b 图：按作者 "
   "`Cell_subtype` 注释的 CD4 谱系亚型中 *JAML* 阳性率。本图不做统计检验，是全数据集的描述性汇总。\n\n"
   "**图 2. 区室分离的独立复核（GSE127465）。** GSE127465（54,773 细胞；未经治疗的 NSCLC；"
   "Zilionis 等[8]）中 *JAML* 与 *CXADR* 在 8 个区室的阳性率，与 GSE131907 对照。细胞由作者注释的"
   "主要与次要细胞类型归组；阳性定义为 GEO 提供的归一化计数矩阵中表达 > 0。分组前已核验矩阵与元数据的"
   "行对应（附加文件 1）。横条为阳性率百分比，同一区室的发现集与复发集数值并列显示。\n\n"
   "**图 3. TCGA-LUAD 中 *JAML* 表达与免疫群体的相关：原始与校正后。** n = 514 例原发瘤中 *JAML* 表达与"
   "各反卷积群体的 Spearman ρ，分别用 MCP-counter、EPIC 与 CIBERSORT 计算。实心圆为原始 ρ；空心标记为"
   "控制整体免疫含量（1 − EPIC `otherCells`）或（有值时）肿瘤纯度（ABSOLUTE；n = 170）的偏 ρ。"
   "连线连接原始值与免疫含量校正值。群体按原始 ρ 排序。图中所有相关的原始 FDR 均 < 0.05。\n\n"
   "**图 4. 校正的有效性与相关谱相似度。**（a）真谱系 marker 基因（T 细胞：*CD3D*、*CD3E*、*CD8A*；"
   "髓系：*LYZ*）与 T 细胞群体的平均 Spearman ρ：未校正（灰柱）、以结构独立的白细胞比例代理校正"
   "（红色菱形）、以 MCP-counter 复合分数校正（蓝色）。复合校正会抹去真 T 细胞 marker 的谱系信号，"
   "因而不是有效校正；独立代理则完整保留。（b）*JAML* 与 15 个 marker 基因在 39 个群体上的相关谱两两 "
   "Pearson 相似度；*JAML* 最像髓系/DC 基因（*CLEC10A*、*LYZ*），与上皮 marker（*EPCAM*、*KRT18*）符号相反。\n\n"
   "**图 5. CD4 谱系拟时序与 *JAML* 梯度。** 扩散拟时序（过滤后 2,006 基因上计算 50 个主成分；2,442 细胞，"
   "每亚型抽样上限 500），以初始型群体定根。（a）各亚型拟时序（中位数与四分位距）。（b）各亚型 *JAML* 阳性率。"
   "（c）scaled *JAML* 表达对拟时序的散点，并标出零值比例——因为 >75% 的 CD4⁺T 细胞 *JAML* 阴性，"
   "该相关（ρ = 0.128，P = 2.6×10⁻¹⁰）由阳性少数细胞驱动。该梯度仅由本数据集支持。\n\n"
   "**图 6. *JAML* 阳性 CD4⁺T 细胞的转录因子活性及其阴性对照。**（a）498 个 CollecTRI 转录因子"
   "（decoupleR，单变量线性模型）在 *JAML* 阳性与阴性细胞间的活性比较，按活性差展示前列因子。"
   "由于活性分数彼此相关，单个因子的身份不可解读，图中已据此标注。（b）对 *JAML* 标签做 200 次随机置换后，"
   "达到 FDR < 0.05 的转录因子数目分布（实测值 229 以虚线标出）；置换均值 0.04、最大 2。阳性对照"
   "（*PDCD1*，同一流程）命中了 9 个预设活化/耗竭因子中的 5 个。\n\n"
   "**图 7. 候选转录因子扰动与 *JAML* 读出。**（a）用 CellOracle（bagging = 20）在 CD4⁺子集（2,842 细胞）中"
   "逐个敲除 39 个候选转录因子后的平均 |Δ*JAML*|；虚线为全局噪声上界（9.59×10⁻³），即全部非退化因子中最大的"
   "随机读数。（b）配对打乱网络对照：每个因子的真实效应（红柱）对 5 次独立换种子打乱重复的均值与标准差"
   "（空心标记）。随机基线低于 1×10⁻⁴ 的因子被标为对照退化（阴影区），分母退化的比值不作解读。"
   "只有 *ID2* 同时超过其自身随机最大值与全局上界。\n\n"
   "**图 8. 跨网络配置的可复现性，以及「关联」与「干预」之间缺乏一致性。**（a）在扰动模型与活性分析中均可评估的 "
   "32 个转录因子，其 *JAML* 扰动效应对活性差作图；虚线为 y = 0。二者无关系（Pearson r = −0.002；"
   "符号一致 12/32）。（b）候选因子在三套网络配置（bagging 3 / 3,013 基因；bagging 20 / 3,013 基因；"
   "bagging 20 / 3,074 基因）中的效应量，显示 *ID2* 的排名稳定性，以及最终配置中新纳入的 60 个转录因子"
   "相对于噪声上界的位置。\n\n"
   "**图 9. DepMap 22Q2 细胞系中 *JAML* 的基因效应。** LUAD 细胞系"
   "（`lineage_sub_subtype == NSCLC_adenocarcinoma`；n = 50）与其他谱系（n = 1,036）中 *JAML* 的 "
   "CRISPR（Chronos）基因效应小提琴/箱线图；点为单个细胞系。虚线为常规依赖阈值（−0.5），无 LUAD 细胞系"
   "低于该阈值。正文以 RPL5（−2.25）与 CDK1（−2.04）作为量程阳性对照，以 KRAS（−0.40，该谱系中的驱动基因"
   "而非普遍依赖基因）作为阴性对照。")

# =====================================================================
# 五、Discussion（m4：HPA 旁证；M8：不主张清单）
# =====================================================================
op(EN, "sub", "myeloid and dendritic cells being the dominant senders in that interface.",
   "myeloid and dendritic cells being the dominant senders in that interface. An independent source "
   "agrees on the compartment: the Human Protein Atlas classifies *JAML* as cell-type enhanced in "
   "neutrophils, monocytes, dendritic cells and T cells (Tau specificity score 0.86, RNA–protein "
   "concordance 0.87), i.e. myeloid-leaning but not myeloid-restricted, which is what the present "
   "data show. That concordance is at the level of mRNA and, for the atlas, of protein staining; "
   "we did not measure *JAML* protein ourselves, and the compartmental claim therefore rests on "
   "transcript-level evidence from two datasets plus this external annotation.")
op(CN, "sub", "髓系和树突状细胞是该界面的主要发送方",
   "髓系和树突状细胞是该界面的主要发送方。一个独立来源与本文的区室判断一致：Human Protein Atlas 将 "
   "*JAML* 归为在中性粒细胞、单核细胞、树突状细胞与 T 细胞中增强表达（Tau 特异度 0.86，"
   "RNA–蛋白一致性 0.87），即**偏向髓系但非髓系专有**——这正是本文数据所显示的。该一致性建立在 mRNA 层面，"
   "对图谱而言还包括蛋白染色；但本文自身未检测 *JAML* 蛋白，因此区室结论依据的是两套数据的转录本证据"
   "加上这一外部注释。")

op(EN, "before", "Three negative results have direct translational bearing.",
   "**What this study does not claim.** *ID2* is a candidate regulator, not a demonstrated one: the "
   "evidence is a computational perturbation model, and the edge it rests on is of moderate strength. "
   "The 229 transcription factors flagged by activity inference are not regulators of *JAML*; they "
   "describe the activation state of cells that express it, and their lack of concordance with "
   "perturbation (r = −0.002) is the reason we say so. *JAML* is not a T-cell infiltration marker, "
   "and we do not claim that it drives CD4+ T-cell infiltration; its correlation profile resembles "
   "myeloid/DC genes and its T-cell association is weaker than that of genuine T-cell markers. The "
   "CD4+ T-cell gradient is single-dataset evidence driven by a minority of cells. The compartmental "
   "and tissue-level analyses are observational, so no causal statement follows from them. And no "
   "therapeutic claim is made: *JAML* is not required by LUAD cells, has no approved drug or "
   "registered trial, and the interface described here has not been tested physically.\n")
op(CN, "before", "三项阴性结果具有直接的转化意义。",
   "**本文不主张什么。***ID2* 是候选调控因子而非已证实的调控因子：证据来自计算机扰动模型，"
   "且其依据的边属中等强度。活性推断标出的 229 个转录因子**不是** *JAML* 的调控因子——它们描述的是"
   "表达 *JAML* 的细胞所处的活化状态，而它们与扰动效应缺乏一致性（r = −0.002）正是我们如此表述的原因。"
   "*JAML* 不是 T 细胞浸润标志物，我们也不主张它驱动 CD4⁺T 浸润：其相关谱更像髓系/DC 基因，"
   "与 T 细胞的关联弱于真正的 T 细胞 marker。CD4⁺T 梯度仅为单一数据集证据，且由少数阳性细胞驱动。"
   "区室与组织层面分析属观察性，不能据此推出因果。也不作任何治疗学主张：*JAML* 并非 LUAD 细胞所必需，"
   "无已批准药物或注册临床试验，本文所述的界面亦未经物理验证。\n")

# =====================================================================
# 六、Limitations 补充
# =====================================================================
op(EN, "sub", "*JAML*–*CXADR* binding would require biochemical confirmation.",
   "*JAML*–*CXADR* binding would require biochemical confirmation. Eleventh, the survival result "
   "depends on the parameterization: *JAML* was significant under a categorical model and not under "
   "the pre-specified continuous model, and no gene was significant under both, so none should be "
   "presented as a prognostic marker; the same caution applies to *MAP4K4*. Twelfth, the "
   "checkpoint-blockade cohort (n = 27; 21 events) can detect only a hazard ratio of about 3.4 or "
   "larger at 80% power, so its null result is uninformative about moderate effects. Thirteenth, the "
   "perturbation analysis is restricted to transcription factors, is built on a promoter-based network "
   "in which the filter threshold had to be relaxed for *JAML* to be evaluable, and its effect sizes "
   "are modest in relative terms (6.3% of mean *JAML* expression; *JAML* ranks 33rd of 3,074 genes), "
   "so *ID2* should be treated as a candidate to be tested rather than a mechanism. Fourteenth, no "
   "protein-level measurement of *JAML* was made in this study.")
op(CN, "sub", "*JAML*–*CXADR* 的结合需生化实验确认。",
   "*JAML*–*CXADR* 的结合需生化实验确认。第十一，生存结论依赖参数化：*JAML* 仅在分类模型下显著、"
   "在预设的连续模型下不显著，且没有任何基因在两种参数化下均显著，因此均不应作为预后标志物呈现；"
   "对 *MAP4K4* 同样适用。第十二，检查点抑制剂队列（n = 27；21 例事件）在 80% 功效下只能检出约 3.4 "
   "或更大的 HR，故其阴性结果对中等效应不具信息量。第十三，扰动分析局限于转录因子、建立在 promoter "
   "网络上、且为使 *JAML* 可评估而放宽了过滤阈值，其效应量在相对意义上并不大（为 *JAML* 平均表达的 "
   "6.3%；*JAML* 在 3,074 个基因中排第 33 位），因此 *ID2* 应视为**待检验的候选**而非机制。"
   "第十四，本文未做 *JAML* 的蛋白层面检测。")

# =====================================================================
# 七、表格与附加文件
# =====================================================================
op(EN, "line", "| Clinical association | TCGA-LUAD; n = 563",
   "| Clinical association | TCGA-LUAD; n = 563; 208 deaths | Parameterization-dependent: continuous "
   "Cox FDR = 0.11; categorical Cox FDR = 0.014; log-rank FDR = 0.012; no gene significant under both "
   "(Table 5) | Four parameterizations on identical sample; stage trend |")
op(CN, "line", "| 临床关联 | TCGA-LUAD；n = 563",
   "| 临床关联 | TCGA-LUAD；n = 563；208 例死亡 | 依赖参数化：连续 Cox FDR = 0.11；分类 Cox FDR = 0.014；"
   "log-rank FDR = 0.012；无基因在两种参数化下均显著（表 5） | 同一批样本上的四种参数化；分期趋势 |")

op(EN, "after", "Global noise ceiling = 9.59 × 10⁻³ (largest randomized reading among non-degenerate TFs).",
   "\n### Table 5. TCGA-LUAD overall survival for *JAML* and *MAP4K4* under four parameterizations\n\n"
   "| Gene | Cox, continuous per SD | Cox, per 1 log2 unit | Cox, median split | log-rank (median split) |\n"
   "|---|---|---|---|---|\n"
   "| *JAML* | HR 0.85; P = 0.017; FDR = 0.11 | HR 0.89; P = 0.017; FDR = 0.11 | HR 0.62; "
   "P = 0.0008; **FDR = 0.014** | χ² = 11.55; P = 6.8 × 10⁻⁴; **FDR = 0.012** |\n"
   "| *MAP4K4* | HR 1.26; P = 0.001; **FDR = 0.018** | HR 1.26; P = 0.001; FDR = 0.018 | HR 1.35; "
   "P = 0.031; FDR = 0.17 | χ² = 4.72; P = 0.030; FDR = 0.17 |\n\n"
   "n = 563; 208 deaths; 18 genes tested, Benjamini–Hochberg correction within each parameterization. "
   "The pre-specified primary model was the continuous per-SD Cox. No gene was significant under both "
   "Cox parameterizations, and the gene passing FDR differs between them.")

op(CN, "after", "全局噪声上界 = 9.59×10⁻³（全部非退化 TF 中最大的随机读数）。",
   "\n### 表 5. TCGA-LUAD 总生存：*JAML* 与 *MAP4K4* 在四种参数化下的结果\n\n"
   "| 基因 | Cox 连续（每 SD） | Cox 连续（每 1 log2 单位） | Cox 中位数分组 | log-rank（中位数分组） |\n"
   "|---|---|---|---|---|\n"
   "| *JAML* | HR 0.85；P = 0.017；FDR = 0.11 | HR 0.89；P = 0.017；FDR = 0.11 | HR 0.62；"
   "P = 0.0008；**FDR = 0.014** | χ² = 11.55；P = 6.8×10⁻⁴；**FDR = 0.012** |\n"
   "| *MAP4K4* | HR 1.26；P = 0.001；**FDR = 0.018** | HR 1.26；P = 0.001；FDR = 0.018 | HR 1.35；"
   "P = 0.031；FDR = 0.17 | χ² = 4.72；P = 0.030；FDR = 0.17 |\n\n"
   "n = 563；208 例死亡；共检验 18 个基因，Benjamini–Hochberg 校正施于各参数化内部。"
   "预设的主要模型为连续每 SD 的 Cox。**没有任何基因在两种 Cox 参数化下均显著**，且通过 FDR 的基因"
   "在两种参数化之间并不相同。")

op(EN, "after", "**Additional file 9.** DepMap 22Q2 gene-effect summary for *JAML*.",
   "\n**Additional file 10.** Reviewer-requested supplementary analyses: relative effect size of *ID2* "
   "perturbation on *JAML* (mean |Δ| as a fraction of mean expression, per-cell distribution, and rank "
   "among 3,074 genes); six-parameterization survival results for all 18 evaluable genes in TCGA-LUAD; "
   "and the power analysis for the checkpoint-blockade cohort (minimum detectable hazard ratio). "
   "(AdditionalFile10_reviewer_analyses.csv)")
op(CN, "after", "**附加文件 9.** DepMap 22Q2 中 *JAML* 基因效应汇总。",
   "\n**附加文件 10.** 应审稿人要求补充的分析：*ID2* 扰动对 *JAML* 的相对效应量（平均 |Δ| 占该基因"
   "平均表达的比例、逐细胞分布、在 3,074 个基因中的排名）；TCGA-LUAD 全部 18 个可评估基因在多种参数化下的"
   "生存结果；以及检查点抑制剂队列的功效分析（最小可检出 HR）。（AdditionalFile10_reviewer_analyses.csv）")


# =====================================================================
def run():
    files = {}
    for f in (EN, CN):
        p = os.path.join(PKG, f)
        files[f] = io.open(p, encoding="utf-8").read().split("\n")
    ok = miss = 0
    for f, mode, anchor, new in OPS:
        lines = files[f]
        idx = [i for i, ln in enumerate(lines) if ln.startswith(anchor)]
        if mode in ("line", "block", "after", "before") and not idx:
            idx = [i for i, ln in enumerate(lines) if anchor in ln]
        if not idx and mode != "sub":
            print(f"  [MISS] {f} [{mode}] {anchor[:60]}")
            miss += 1
            continue
        if mode == "sub":
            joined = "\n".join(lines)
            if anchor not in joined:
                print(f"  [MISS] {f} [sub] {anchor[:60]}")
                miss += 1
                continue
            files[f] = joined.replace(anchor, new, 1).split("\n")
            ok += 1
            continue
        i = idx[0]
        if mode == "line":
            lines[i] = new
        elif mode == "after":
            lines[i + 1:i + 1] = new.split("\n")
        elif mode == "before":
            lines[i:i] = new.split("\n")
        elif mode == "block":
            j = i
            while j < len(lines) and lines[j].strip():
                j += 1
            lines[i:j] = new.split("\n")
        ok += 1
    for f, content in files.items():
        io.open(os.path.join(PKG, f), "w", encoding="utf-8", newline="\n").write("\n".join(content))
    print(f"\n应用 {ok} 处修改；未命中 {miss} 处")
    return miss


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
