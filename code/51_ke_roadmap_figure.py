# -*- coding: utf-8 -*-
"""生成《JAML机制课题研究》技术路线图（PNG + 300 dpi LZW TIFF）。

用法：python 51_ke_roadmap_figure.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化\JAML机制课题研究"

# 中文字体（Windows 常见可用字体，逐个尝试）
for fam in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "SimSun"]:
    try:
        matplotlib.font_manager.findfont(fam, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [fam]
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False

C_PREM = "#E8EEF6"; C_EDGE = "#4A6FA5"
C_1 = "#FBE9E7"; C_2 = "#FFF6E5"; C_3 = "#E8F5E9"; C_4 = "#F3E5F5"


def box(ax, x, y, w, h, text, fc, ec=C_EDGE, fs=8.4, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                linewidth=1.0, edgecolor=ec, facecolor=fc))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", linespacing=1.5)


def arrow(ax, x1, y1, x2, y2, style="-|>", lw=1.1):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=11, linewidth=lw, color="#37474F"))


fig, ax = plt.subplots(figsize=(7.4, 9.0))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# ---- 顶部：前期工作（已完成）
box(ax, 0.06, 0.895, 0.88, 0.085,
    "公共数据深化分析（已完成 M2–M9，八项正交分析，全部含对照）\n"
    "产出三条硬约束：① 干预位点在免疫侧（上皮阳性率 <2%、非肿瘤细胞依赖基因）\n"
    "② 读出以髓系/DC 为主（CD4+T 梯度仅单数据集支持）③ ID2 为唯一通过噪声校准的候选 TF",
    C_PREM, fs=8.0)
arrow(ax, 0.5, 0.893, 0.5, 0.862)

# ---- 课题 1
box(ax, 0.06, 0.735, 0.88, 0.12,
    "课题 1　临床标本区室与空间定位验证\n"
    "手术标本 60 例（肿瘤 + 配对癌旁）→ 单细胞悬液流式 + 多色免疫荧光（mIF）\n"
    "读出：JAML+ 髓系占比；JAML+ 髓系 – CXADR+ 上皮的空间邻接富集倍数\n"
    "作用：确立“在哪里”，并为课题 3/4 提供取材与判读标准",
    C_1, fs=8.0)
arrow(ax, 0.5, 0.733, 0.5, 0.702)

# ---- 课题 2
box(ax, 0.06, 0.545, 0.88, 0.15,
    "课题 2　ID2 → JAML 转录调控验证（核心机制）\n"
    "双路线敲低：siRNA + 慢病毒 CRISPRi（dCas9-KRAB）→ JAML mRNA / 表面蛋白\n"
    "过表达验证方向性；ChIP-qPCR / CUT&RUN 检测 ID2 在 JAML 调控区的结合\n"
    "报告基因：启动子野生型 vs 结合位点突变型；siRNA 抗性 ID2 回补\n"
    "强制对照：阴性（非靶基因）+ 阳性（已知 ID2 靶基因）",
    C_2, fs=8.0)
arrow(ax, 0.5, 0.543, 0.5, 0.512)

# ---- 课题 3
box(ax, 0.06, 0.355, 0.88, 0.15,
    "课题 3　JAML–CXADR 界面的功能后果与阻断\n"
    "共培养：上皮（CXADR 高/低、敲低/过表达）× 髓系（THP-1 来源巨噬、mo-DC）× 原代 CD4+T\n"
    "干预：抗 JAML 阻断抗体 / 重组 JAML-Fc / 抗 CXADR；对照：同型 IgG、抗 PD-L1\n"
    "读数：CD69 / CD25 / IFN-γ / GZMB；上皮 pAKT 与 PI3K 通路\n"
    "反向验证：髓系 JAML 敲低是否复现抗体效应；Fc 沉默突变体排除 Fc 依赖效应",
    C_3, fs=8.0)
arrow(ax, 0.5, 0.353, 0.5, 0.322)

# ---- 课题 4
box(ax, 0.06, 0.175, 0.88, 0.14,
    "课题 4　患者来源模型中的治疗相关性探索（探索性）\n"
    "LUAD 患者来源类器官（PDO）+ 自体 PBMC 来源巨噬 / T 细胞共培养\n"
    "重复课题 3 关键干预；与标本 JAML+ 髓系占比、邻接比例做相关\n"
    "预案：PDO 建成 < 10 例则降级为方法学报告，不影响课题 1–3 结题",
    C_4, fs=8.0)
arrow(ax, 0.5, 0.173, 0.5, 0.142)

# ---- 交付
box(ax, 0.06, 0.035, 0.88, 0.10,
    "交付物　机制结论（ID2→JAML 成立 / 不成立，二者均可结题，阴性结果照实报告并给出可检出效应量）\n"
    "SCI 论文 2–3 篇；发明专利 1 项；60 例 LUAD 配对标本 mIF/流式数据库；脚本与结果文件归档",
    "#ECEFF1", fs=8.0, bold=True)

# ---- 右侧纵向：贯穿全程的质量控制
ax.add_patch(FancyBboxPatch((0.955, 0.035), 0.038, 0.94,
                            boxstyle="round,pad=0.008,rounding_size=0.015",
                            linewidth=1.0, edgecolor="#795548", facecolor="#FBE9E7"))
ax.text(0.974, 0.505, "质量控制贯穿全程\n样本量与功效复核\n抗体特异度验证\n对照体系前置\n原始数据可追溯",
        ha="center", va="center", rotation=90, fontsize=7.6, linespacing=1.9)

fig.tight_layout()
fig.savefig(os.path.join(OUT, "技术路线图.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "技术路线图.tiff"), dpi=300, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
print("[图] 技术路线图.png / .tiff 已生成 ->", OUT)
