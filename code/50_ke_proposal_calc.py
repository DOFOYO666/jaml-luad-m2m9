# -*- coding: utf-8 -*-
"""课题样本量估算与前期关键数值的复算脚本。

用途：为《JAML机制课题研究》立项书 §6.1 的样本量估算提供可复现来源，
并复算正文引用的若干关键统计量（DepMap、ICB log-rank、对齐 ρ 等），
确保立项书中的每个数字都能指到文件或脚本。

用法：python 50_ke_proposal_calc.py
"""
import json
import math
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化"
OUT = os.path.join(ROOT, "JAML机制课题研究", "prelim_calc_output.txt")


def sample_sizes():
    lines = ["=== 样本量估算（双侧 alpha = 0.05，power = 0.80）===",
             "方法：statsmodels TTestIndPower / NormalIndPower；相关系数用 Fisher z 近似 +3 修正",
             "说明：效应量假设来自前期公共数据的观察值（阳性率差约 27 个百分点）与常规经验值，",
             "      不是文献推测的精确值；立项后应以预实验的实际离散度重算。"]
    try:
        from statsmodels.stats.power import TTestIndPower, NormalIndPower
        from statsmodels.stats.proportion import proportion_effectsize
    except Exception as e:                                   # pragma: no cover
        lines.append(f"[警告] statsmodels 不可用（{type(e).__name__}），改用正态近似公式")
        TTestIndPower = None

    def n_t(d):
        if TTestIndPower is not None:
            return math.ceil(TTestIndPower().solve_power(effect_size=d, alpha=0.05,
                                                         power=0.80, ratio=1.0,
                                                         alternative="two-sided"))
        z = 1.959964 + 0.841621
        return math.ceil(2 * (z / d) ** 2)

    lines.append("")
    lines.append("-- 两组比较（t 检验）--")
    for d in (0.6, 0.8, 1.0):
        n = n_t(d)
        lines.append(f"   Cohen's d = {d:.1f}  →  每组 n = {n}（总 {2*n}）")

    lines.append("")
    lines.append("-- 配对比较（同一患者肿瘤 vs 癌旁）--")
    for d in (0.6, 0.8, 1.0):
        z = 1.959964 + 0.841621
        lines.append(f"   d = {d:.1f}  →  n = {math.ceil((z/d)**2)} 例")

    lines.append("")
    lines.append("-- 相关分析（Fisher z 近似）--")
    for r in (0.4, 0.5, 0.6):
        zz = 0.5 * math.log((1 + r) / (1 - r))
        n = ((stats.norm.ppf(0.975) + stats.norm.ppf(0.80)) / zz) ** 2 + 3
        lines.append(f"   r = {r:.1f}  →  n = {math.ceil(n)}")

    lines.append("")
    lines.append("-- 阳性率比较（正态近似）--")
    pairs = [(0.47, 0.20, "JAML+ 髓系占比：肿瘤 vs 癌旁（前期实测 47.0%）"),
             (0.30, 0.10, "JAML+ 髓系占比（复发集 30.1% 情境）")]
    for p1, p2, tag in pairs:
        try:
            from statsmodels.stats.proportion import proportion_effectsize
            from statsmodels.stats.power import NormalIndPower
            es = proportion_effectsize(p1, p2)
            n = math.ceil(NormalIndPower().solve_power(effect_size=es, alpha=0.05,
                                                       power=0.80, ratio=1.0,
                                                       alternative="two-sided"))
        except Exception:
            h = 2 * math.asin(math.sqrt(p1)) - 2 * math.asin(math.sqrt(p2))
            n = math.ceil(2 * ((1.959964 + 0.841621) / h) ** 2)
        lines.append(f"   {p1:.2f} vs {p2:.2f}  →  每组 n = {n}（总 {2*n}）   [{tag}]")

    lines.append("")
    lines.append("结论：课题 1 拟定 60 例配对标本（覆盖最保守的 47 例/组，并预留约 20% 脱落）。")
    return lines


def recompute_key_numbers():
    lines = ["", "=== 关键数值复算（与立项书/稿件正文对照）==="]
    # DepMap
    p = os.path.join(ROOT, "results", "m6_depmap_jaml_summary.json")
    if os.path.isfile(p):
        j = json.load(open(p, encoding="utf-8"))
        lines.append(f"  [M6] DepMap {j.get('release')}: LUAD n = {j.get('n_LUAD')}, "
                     f"mean = {j.get('LUAD_mean_gene_effect')}, "
                     f"n< -0.5 = {j.get('LUAD_n_below_neg0.5')}, Wilcoxon P = {j.get('wilcox_p')}")
        pc = j.get("positive_control_median_gene_effect", {})
        lines.append(f"       阳性对照: RPL5 {pc.get('RPL5')}, CDK1 {pc.get('CDK1')}, KRAS {pc.get('KRAS')}")
    else:
        lines.append("  [M6] 摘要文件缺失")

    # ICB log-rank 独立复算
    p = os.path.join(ROOT, "results", "m5_gse135222_jaml_pfs.csv")
    if os.path.isfile(p):
        d = pd.read_csv(p)
        hi = (d["expr"] > d["expr"].median()).values
        t = d["pfs_time"].values.astype(float)
        e = d["pfs"].values.astype(int)
        O1 = E1 = V = 0.0
        for tt in np.sort(np.unique(t[e == 1])):
            at = t >= tt
            n, n1 = int(at.sum()), int((hi & at).sum())
            d_ = int(((t == tt) & (e == 1)).sum())
            d1 = int(((t == tt) & (e == 1) & hi).sum())
            if n < 2:
                continue
            O1 += d1
            E1 += d_ * n1 / n
            V += d_ * (n1 / n) * (1 - n1 / n) * (n - d_) / (n - 1)
        chi2 = (O1 - E1) ** 2 / V
        lines.append(f"  [M5] GSE135222: n = {len(d)}, 事件 {int(e.sum())}, "
                     f"log-rank chi2 = {chi2:.4f}, P = {1-stats.chi2.cdf(chi2,1):.4f}")

    # M4 配对对照
    p = os.path.join(r"D:\workbuddy工作空间\JAML深度研究\results", "m4_noise_control_classify.json")
    if os.path.isfile(p):
        j = json.load(open(p, encoding="utf-8"))
        lines.append(f"  [M4] 全局噪声上界 = {j.get('global_noise_ceiling_nondegenerate'):.4g}; "
                     f"通过配对的 TF = {j.get('pass_paired_nondegenerate')}; "
                     f"通过全局上界的 TF = {j.get('pass_global_ceiling')}")

    # M9 对齐
    p = os.path.join(ROOT, "results", "m8_gse127465_alignment_check.json")
    if os.path.isfile(p):
        j = json.load(open(p, encoding="utf-8"))
        mt = j["spearman_mt_fraction"]
        mc = j.get("marker_check", [])
        n_pass = sum(1 for m in mc if m.get("pass"))
        lines.append(f"  [M9] 对齐 MT 比例 rho = {mt['rho']:.7f} (n = {mt['n']}); "
                     f"marker 检查通过 {n_pass}/{len(mc)}; nnz = {j['nnz']} "
                     f"(读取行数 = {j.get('lines_read')})")
        tc = j.get("spearman_totalcounts", {})
        lines.append(f"       Total counts 列 rho = {tc.get('rho'):.4f}（语义未明，全程未用作判据）")
    return lines


if __name__ == "__main__":
    out = sample_sizes() + recompute_key_numbers()
    txt = "\n".join(out)
    print(txt)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(txt + "\n")
    print(f"\n[输出] {OUT}")
