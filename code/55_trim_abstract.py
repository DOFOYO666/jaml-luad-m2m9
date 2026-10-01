# -*- coding: utf-8 -*-
"""把英文摘要压缩到 ≤350 词（审稿意见 M5），并同步压缩中文摘要的对应段落。"""
import io
import os
import re

PKG = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化\M2-M9深化研究稿"

EN_NEW = {
    "**Background.**":
        "**Background.** Germline variation at *JAML* has been linked to lung adenocarcinoma (LUAD) "
        "risk through CD4+ T-cell expression, but where the molecule acts, which cells determine its "
        "tissue abundance, and whether it is transcriptionally controlled remain unknown.",
    "**Methods.**":
        "**Methods.** In a secondary analysis of public data we interrogated one gene with eight "
        "orthogonal analyses, each with an explicit control: single-cell expression in LUAD (GSE131907; "
        "208,506 cells) and NSCLC (GSE127465; 54,773 cells) with verified alignment; three "
        "deconvolution methods with validated partial-correlation adjustment; TCGA-LUAD survival under "
        "four parameterizations; ligand\u2013receptor inference and CD4-lineage pseudotime; "
        "transcription-factor (TF) activity with permutation and positive controls; *in silico* TF "
        "perturbation against a paired randomized-network control and a global noise ceiling; and "
        "DepMap dependency data.",
    "**Results.**":
        "**Results.** *JAML* and *CXADR* were compartmentally separated in both datasets: *JAML* "
        "positivity was highest in myeloid/dendritic cells (47.0%, 30.1%) and rare in epithelium "
        "(1.7%, 1.5%), the epithelial-to-myeloid ratio being 0.04 and 0.05, whereas *CXADR* was "
        "epithelial (40.5%, 19.3%) (alignment verified; Additional file 1). *JAML* abundance tracked "
        "myeloid/DC content (\u03c1 = 0.757), the only association robust to adjustment for immune "
        "content and purity (0.750, 0.738), and its correlation profile resembled myeloid/DC genes "
        "(*CLEC10A* 0.973 versus *CD8A* 0.623). Of 39 perturbed TFs, only *ID2* exceeded both its own "
        "randomized control (6.6\u00d7) and the global noise ceiling (1.8\u00d7), reproduced in three "
        "network configurations, and ranked 33rd of 3,074 genes by relative effect (6.3% of mean "
        "expression); *GATA3* was suggestive, while *MAF*, *STAT1*, *IRF1* and *STAT3* were "
        "indistinguishable from noise. TF-activity inference flagged 229 TFs without concordance with "
        "perturbation (r = \u22120.002). *JAML* was not a LUAD cell-line dependency (0/50) and not a "
        "consistent prognostic marker (categorical Cox FDR = 0.014; continuous Cox, pre-specified, "
        "FDR = 0.11; none significant under both), and was unrelated to progression-free survival on "
        "checkpoint blockade (n = 27; minimum detectable hazard ratio 3.4).",
    "**Conclusions.**":
        "**Conclusions.** *JAML* is an immune-side molecule whose tissue abundance is dominated by "
        "myeloid and dendritic cells, whose counterpart here is epithelial *CXADR*, and whose "
        "expression is at best weakly controlled, with *ID2* the only candidate regulator surviving "
        "noise calibration. These results define where *JAML* biology should be tested, and argue "
        "against treating it as a tumour-cell-intrinsic dependency.",
}

CN_NEW = {
    "**结果。**":
        "**结果。** 两套数据中 *JAML* 与 *CXADR* 均处于不同区室：*JAML* 阳性率在髓系/树突状细胞最高"
        "（47.0%、30.1%），上皮中罕见（1.7%、1.5%），上皮/髓系比值 0.04 与 0.05；*CXADR* 位于上皮"
        "（40.5%、19.3%）（对齐经核验，附加文件 1）。*JAML* 丰度跟随髓系/DC（ρ = 0.757），是唯一对整体"
        "免疫含量与肿瘤纯度校正均稳健的关联（0.750、0.738）；其相关谱更像髓系/DC 基因（*CLEC10A* 0.973）"
        "而非 T 细胞基因（*CD8A* 0.623）。在 39 个被扰动的 TF 中，只有 *ID2* 同时超过其自身随机对照"
        "（6.6 倍）与全局噪声上界（1.8 倍），并在三套网络配置中复现；其相对效应为 *JAML* 平均表达的 6.3%，"
        "在 3,074 个基因中排第 33 位。*GATA3* 仅提示性；*MAF*、*STAT1*、*IRF1*、*STAT3* 与噪声不可区分。"
        "TF 活性推断标出 229 个 TF，但与扰动无一致性（r = −0.002）。*JAML* 既非 LUAD 细胞系依赖基因"
        "（0/50 低于阈值），也不是稳健的预后标志物（分类参数化 Cox FDR = 0.014；预设的连续 Cox FDR = 0.11；"
        "无基因在两种参数化下均显著），且与检查点抑制剂治疗下的无进展生存无关（n = 27；最小可检出 HR 3.4）。",
}


def wc(txt):
    return len(re.findall(r"[A-Za-z][A-Za-z\-']+", txt))


def apply(fname, table):
    p = os.path.join(PKG, fname)
    lines = io.open(p, encoding="utf-8").read().split("\n")
    done = 0
    for tag, new in table.items():
        hit = [i for i, ln in enumerate(lines) if ln.startswith(tag)]
        if not hit:
            print(f"  [MISS] {fname} {tag}")
            continue
        lines[hit[0]] = new
        done += 1
    io.open(p, "w", encoding="utf-8", newline="\n").write("\n".join(lines))
    print(f"  {fname}: 替换 {done}/{len(table)} 段")
    return done


if __name__ == "__main__":
    apply("manuscript_EN.md", EN_NEW)
    apply("manuscript_CN.md", CN_NEW)
    tot = sum(wc(v) for v in EN_NEW.values())
    print(f"\n英文摘要词数 = {tot}（目标 ≤ 350）")
    for k, v in EN_NEW.items():
        print(f"   {k:18s} {wc(v):4d} 词")
    cn = io.open(os.path.join(PKG, "manuscript_CN.md"), encoding="utf-8").read().split("\n")
    ab = []
    ins = False
    for ln in cn:
        if ln.startswith("## 摘要"):
            ins = True
            continue
        if ins and ln.startswith("---"):
            break
        if ins and ln.strip():
            ab.append(ln)
    print(f"中文摘要非空白字符数 = {len(re.sub(r'\\s', '', ''.join(ab)))}")
