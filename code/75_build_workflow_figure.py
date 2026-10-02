"""Build the study-design / analysis-workflow flowchart for the M2-M9 (JAML/LUAD) paper.

Format follows the conventions already fixed for this project:
  * drawn at the JTM/BMC full-page width (170 mm), so every font size in the file
    IS its printed size; the minimum here is 6 pt, which meets the journal's
    "legible at these dimensions" rule (>= ~0.897 pt per inch of figure width);
  * 300 dpi raster, TIFF saved with LZW compression;
  * light background, flat fills, nothing but the frame inside the boxes.

Line breaks are found by MEASURING the rendered text (matplotlib renderer), not by
counting characters: character-count estimates came out ~10 % too generous for
mixed-case text and let strings run past the frame edge. Chinese has no word
spaces, so the wrapper tokenises per character for ideographs.

Outputs (Workflow_figure/):
    workflow_study_design.{pdf,svg,png,tiff}       English (figure for the paper)
    workflow_study_design_zh.{pdf,svg,png,tiff}    Chinese (for internal use)
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "Workflow_figure")

W = 170.0                      # mm; JTM/BMC full page width
MAX_H = 225.0                  # mm; JTM/BMC maximum height for figure + legend
PAD = 8.0
GAP = 12.0
COLW = (W - 2 * PAD - 2 * GAP) / 3.0
COLX = [PAD, PAD + COLW + GAP, PAD + 2 * (COLW + GAP)]
COLC = [x + COLW / 2.0 for x in COLX]
INNER = COLW - 4.4             # text width available inside a module box

INK = "#1F2933"
MUTED = "#5C6B7A"
RULE = "#B8C4CF"

F_DATA, S_DATA = "#F2F4F6", "#9AA7B4"
F_TRK, S_TRK, T_TRK = "#E8F1FB", "#4A7FA8", "#1F4E6E"
F_GATE, S_GATE, T_GATE = "#FCECEC", "#BC6363", "#8A3232"
F_CONC, S_CONC, T_CONC = "#EBF3E2", "#6E9440", "#3F5C1F"
BADGE = "#3E6F96"

FS_TITLE, FS_SUB, FS_SEC = 8.6, 6.0, 7.0
FS_BOX, FS_DET, FS_TRK = 6.3, 6.0, 6.6
LH = 2.8                       # mm per line at 6 pt

EN = dict(
    name="workflow_study_design",
    font="DejaVu Sans",
    title="Study design and analysis workflow",
    subtitle="Public data only; no new human participant data were generated, "
             "and no causal effect was re-estimated.",
    data_head="Data sources",
    data=[
        ("GSE131907 \u2014 discovery", "LUAD scRNA-seq; 208,506 cells, 29,634 genes"),
        ("GSE127465 \u2014 replication", "treatment-na\u00efve NSCLC scRNA-seq; 54,773 cells"),
        ("TCGA-LUAD (GDC / UCSC Xena)", "bulk RNA + clinical; n = 514 / 563"),
        ("DepMap 22Q2 (Chronos)", "CRISPR gene effect; 50 LUAD cell lines"),
        ("GSE135222", "anti-PD-(L)1 cohort; n = 27, progression-free survival"),
        ("OneK1K / eQTL Catalogue", "genetic anchor; no causal effect re-estimated"),
    ],
    tracks=[
        ("Track 1 \u00b7 where is JAML?", [
            ("Compartmental expression",
             "66-gene panel, 10 cell types; epithelial : myeloid ratio", "Fig 1"),
            ("Independent replication",
             "matrix alignment verified two ways; same ratio reproduced", "Fig 2"),
            ("Tissue deconvolution",
             "MCP-counter / EPIC / CIBERSORT \u00d7 3 adjustments", "Fig 3"),
            ("Is the adjustment valid?",
             "calibrated against lineage marker genes", "Fig 4"),
            ("Ligand\u2013receptor inference",
             "liana + independent re-implementation; JAML \u2192 CXADR", ""),
        ]),
        ("Track 2 \u00b7 clinical & dependency", [
            ("Survival: pre-specified model",
             "continuous per-SD Cox, BH-corrected", ""),
            ("Survival: secondary models",
             "per log2 unit, median split, log-rank", ""),
            ("Stage association",
             "Kruskal\u2013Wallis + rank trend", ""),
            ("Immune checkpoint blockade",
             "median split + log-rank; n = 27, underpowered", ""),
            ("CRISPR dependency",
             "JAML not required by LUAD cell lines", "Fig 9"),
        ]),
        ("Track 3 \u00b7 regulatory layer", [
            ("CD4+ diffusion pseudotime",
             "naive population as root; gradient in a minority", "Fig 5"),
            ("TF activity inference",
             "decoupleR / CollecTRI, 498 TFs, two contrasts", "Fig 6"),
            ("Null calibration and control",
             "200 permutations; PDCD1 expected set recovered", ""),
            ("TF perturbation",
             "CellOracle; three network configurations", "Fig 7, 8b"),
            ("Association vs intervention",
             "no concordance (r = \u22120.002; 12/32)", "Fig 8a"),
        ]),
    ],
    gate_title="Calibration and reporting rules, applied across every analysis",
    gate_lines=[
        "Paired randomised control, fresh seed each run \u2192 global noise ceiling "
        "9.59 \u00d7 10\u207b\u00b3, an empirical maximum over a finite number of "
        "draws and therefore a liberal benchmark.",
        "Five replicates per factor \u21d2 smallest attainable empirical P = "
        "1/6 \u2248 0.17: a calibrated screen, not a hypothesis test.",
        "One pre-specified parameterisation per family; every effect reported "
        "together with its control; eight complementary analyses \u2014 not "
        "independent, not orthogonal.",
    ],
    conc_title="Where this leaves JAML",
    conc_bullets=[
        "JAML abundance in tumour tissue is carried by myeloid and dendritic "
        "cells; its epithelial counterpart is CXADR.",
        "It is not a T-cell infiltration marker, not a LUAD dependency, and not "
        "a robust prognostic marker.",
        "ID2 is the only candidate regulator to survive noise calibration \u2014 "
        "a model-generated hypothesis, not a mechanism.",
        "Next test: perturb ID2 with JAML readout in primary myeloid / CD4+ "
        "cells, and block the JAML\u2013CXADR interface.",
    ],
)

ZH = dict(
    name="workflow_study_design_zh",
    font="Microsoft YaHei",
    title="\u7814\u7a76\u8bbe\u8ba1\u4e0e\u5206\u6790\u6d41\u7a0b",
    subtitle="\u4ec5\u7528\u516c\u5171\u6570\u636e\uff1b\u672a\u4ea7\u751f\u65b0\u7684"
             "\u4eba\u4f53\u53d7\u8bd5\u8005\u6570\u636e\uff0c\u4e5f\u672a\u91cd"
             "\u4f30\u56e0\u679c\u6548\u5e94\u3002",
    data_head="\u6570\u636e\u6765\u6e90",
    data=[
        ("GSE131907\u00a0\u00b7 \u53d1\u73b0\u96c6", "LUAD \u5355\u7ec6\u80de\uff1b208,506 \u7ec6\u80de\u300129,634 \u57fa\u56e0"),
        ("GSE127465\u00a0\u00b7 \u9a8c\u8bc1\u96c6", "\u672a\u6cbb\u7597 NSCLC \u5355\u7ec6\u80de\uff1b54,773 \u7ec6\u80de"),
        ("TCGA-LUAD\uff08GDC / Xena\uff09", "bulk \u8868\u8fbe + \u4e34\u5e8a\uff1bn = 514 / 563"),
        ("DepMap 22Q2\uff08Chronos\uff09", "CRISPR \u57fa\u56e0\u6548\u5e94\uff1b50 \u682a LUAD \u7ec6\u80de\u7cfb"),
        ("GSE135222", "\u6297 PD-(L)1 \u961f\u5217\uff1bn = 27\uff0c\u65e0\u8fdb\u5c55\u751f\u5b58"),
        ("OneK1K / eQTL Catalogue", "\u9057\u4f20\u5b66\u7ebf\u7d22\uff1b\u4e0d\u91cd\u4f30\u56e0\u679c\u6548\u5e94"),
    ],
    tracks=[
        ("\u901a\u9053 1 \u00b7 JAML \u5728\u54ea\u91cc\uff1f", [
            ("\u5206\u533a\u95f4\u8868\u8fbe",
             "66 \u57fa\u56e0\u9762\u677f\u300110 \u79cd\u7ec6\u80de\uff1b\u4e0a\u76ae : \u9ad3\u7cfb\u6bd4", "\u56fe 1"),
            ("\u72ec\u7acb\u9a8c\u8bc1",
             "\u77e9\u9635\u5bf9\u9f50\u53cc\u91cd\u6821\u9a8c\uff1b\u6bd4\u503c\u53ef\u590d\u73b0", "\u56fe 2"),
            ("\u7ec4\u7ec7\u53bb\u5377\u79ef",
             "MCP-counter / EPIC / CIBERSORT \u00d7 3 \u79cd\u6821\u6b63", "\u56fe 3"),
            ("\u6821\u6b63\u662f\u5426\u6709\u6548\uff1f",
             "\u4ee5\u8c31\u7cfb\u6807\u5fd7\u57fa\u56e0\u6807\u5b9a", "\u56fe 4"),
            ("\u914d\u4f53\u2013\u53d7\u4f53\u63a8\u65ad",
             "liana + \u72ec\u7acb\u590d\u73b0\uff1bJAML \u2192 CXADR", ""),
        ]),
        ("\u901a\u9053 2 \u00b7 \u4e34\u5e8a\u4e0e\u4f9d\u8d56\u6027", [
            ("\u751f\u5b58\uff1a\u9884\u8bbe\u6a21\u578b",
             "\u8fde\u7eed per-SD Cox\uff0cBH \u6821\u6b63", ""),
            ("\u751f\u5b58\uff1a\u6b21\u8981\u6a21\u578b",
             "per log2\u3001\u4e2d\u4f4d\u6570\u5206\u7ec4\u3001log-rank", ""),
            ("\u5206\u671f\u5173\u8054",
             "Kruskal\u2013Wallis + \u8d8b\u52bf\u68c0\u9a8c", ""),
            ("\u514d\u75ab\u68c0\u67e5\u70b9\u6cbb\u7597",
             "\u4e2d\u4f4d\u6570\u5206\u7ec4 + log-rank\uff1bn = 27\uff0c\u529f\u6548\u4e0d\u8db3", ""),
            ("CRISPR \u4f9d\u8d56\u6027",
             "LUAD \u7ec6\u80de\u7cfb\u4e0d\u4f9d\u8d56 JAML", "\u56fe 9"),
        ]),
        ("\u901a\u9053 3 \u00b7 \u8c03\u63a7\u5c42", [
            ("CD4+ \u6269\u6563\u62df\u65f6\u5e8f",
             "\u4ee5 naive \u4e3a\u6839\uff1b\u4ec5\u5c11\u6570\u7ec6\u80de\u6709\u68af\u5ea6", "\u56fe 5"),
            ("\u8f6c\u5f55\u56e0\u5b50\u6d3b\u6027",
             "decoupleR / CollecTRI\uff0c498 \u4e2a TF\uff0c\u4e24\u79cd\u5bf9\u6bd4", "\u56fe 6"),
            ("\u96f6\u5206\u5e03\u6807\u5b9a\u4e0e\u5bf9\u7167",
             "200 \u6b21\u7f6e\u6362\uff1bPDCD1 \u9884\u671f\u96c6\u53ef\u590d\u73b0", ""),
            ("\u865a\u62df\u6270\u52a8",
             "CellOracle\uff1b\u4e09\u5957\u7f51\u7edc\u914d\u7f6e", "\u56fe 7, 8b"),
            ("\u5173\u8054 vs \u5e72\u9884",
             "\u65e0\u4e00\u81f4\u6027\uff08r = \u22120.002\uff1b12/32\uff09", "\u56fe 8a"),
        ]),
    ],
    gate_title="\u6807\u5b9a\u4e0e\u62a5\u544a\u89c4\u5219\uff08\u9002\u7528\u4e8e\u5168\u90e8\u5206\u6790\uff09",
    gate_lines=[
        "\u914d\u5bf9\u968f\u673a\u5bf9\u7167\u3001\u6bcf\u6b21\u91cd\u62bd\u79cd\u5b50 "
        "\u2192 \u5168\u5c40\u566a\u58f0\u4e0a\u754c $9.59\\times10^{-3}$\uff1a"
        "\u6709\u9650\u6b21\u62bd\u53d6\u7684\u7ecf\u9a8c\u6700\u5927\u503c\uff0c"
        "\u5c5e\u504f\u5bbd\u677e\u7684\u5224\u636e\u3002",
        "\u6bcf\u4e2a\u56e0\u5b50 5 \u6b21\u91cd\u590d \u2192 \u6700\u5c0f\u7ecf\u9a8c "
        "P = 1/6 \u2248 0.17\uff1a\u8fd9\u662f\u6807\u5b9a\u8fc7\u7684\u7b5b\u9009\uff0c"
        "\u4e0d\u662f\u5047\u8bbe\u68c0\u9a8c\u3002",
        "\u6bcf\u4e2a\u5206\u6790\u65cf\u53ea\u9884\u8bbe\u4e00\u4e2a\u4e3b\u53c2\u6570"
        "\u5316\uff1b\u6bcf\u4e2a\u6548\u5e94\u90fd\u4e0e\u540c\u671f\u5bf9\u7167\u4e00"
        "\u5e76\u62a5\u544a\uff1b\u516b\u9879\u4e92\u8865\u5206\u6790\u2014\u2014\u5e76"
        "\u975e\u72ec\u7acb\uff0c\u4e5f\u975e\u6b63\u4ea4\u3002",
    ],
    conc_title="\u672c\u6587\u5bf9 JAML \u7684\u5b9a\u4f4d",
    conc_bullets=[
        "\u80bf\u7624\u7ec4\u7ec7\u4e2d JAML \u4e3b\u8981\u7531\u9ad3\u7cfb\u4e0e"
        "\u6811\u7a81\u72b6\u7ec6\u80de\u627f\u8f7d\uff1b\u5176\u4e0a\u76ae\u4fa7\u5bf9"
        "\u5e94\u5206\u5b50\u662f CXADR\u3002",
        "\u5b83\u4e0d\u662f T \u7ec6\u80de\u6d78\u6da6\u6807\u5fd7\u7269\uff0c\u4e0d"
        "\u662f LUAD \u7ec6\u80de\u4f9d\u8d56\u57fa\u56e0\uff0c\u4e5f\u4e0d\u662f"
        "\u7a33\u5065\u7684\u9884\u540e\u6807\u5fd7\u7269\u3002",
        "ID2 \u662f\u552f\u4e00\u901a\u8fc7\u566a\u58f0\u6807\u5b9a\u7684\u5019\u9009"
        "\u8c03\u63a7\u56e0\u5b50\u2014\u2014\u5c5e\u6a21\u578b\u751f\u6210\u7684\u5047"
        "\u8bbe\uff0c\u4e0d\u662f\u673a\u5236\u3002",
        "\u4e0b\u4e00\u6b65\uff1a\u5728\u539f\u4ee3\u9ad3\u7cfb / CD4+ \u7ec6\u80de\u4e2d"
        "\u6270\u52a8 ID2 \u5e76\u8bfb\u51fa JAML\uff0c\u5171\u57f9\u517b\u963b\u65ad"
        " JAML\u2013CXADR \u754c\u9762\u3002",
    ],
)

_probe = None


def _cls(ch):
    o = ord(ch)
    if 0x4E00 <= o <= 0x9FFF or 0x3400 <= o <= 0x4DBF:
        return "ideo"
    if 0x3000 <= o <= 0x303F or 0xFF00 <= o <= 0xFFEF or o in (0x2018, 0x2019,
                                                               0x201C, 0x201D,
                                                               0x2026):
        return "punct"
    return "latin"


def tokens(text):
    """Split into unbreakable tokens: Latin words, or one ideograph each."""
    out, buf = [], ""
    for ch in text:
        if ch.isspace():
            if buf:
                out.append(buf)
                buf = ""
        elif _cls(ch) == "ideo":
            if buf:
                out.append(buf)
                buf = ""
            out.append(ch)
        else:
            buf += ch
    if buf:
        out.append(buf)
    return out


def text_mm(ax, s, size):
    """Width of `s` in mm at `size` pt, measured by the renderer."""
    global _probe
    if _probe is None:
        _probe = ax.text(0, 0, "", fontsize=size)
    _probe.set_fontsize(size)
    _probe.set_text(s)
    bb = _probe.get_window_extent(renderer=ax.figure.canvas.get_renderer())
    return bb.width / ax.figure.dpi * 25.4


def wrap_mm(ax, text, size, width_mm):
    """Greedy wrap to a measured width (mm), handling CJK per character."""
    lines, cur = [], ""
    for tok in tokens(text):
        if not cur:
            cur = tok
            continue
        a, b = _cls(cur[-1]), _cls(tok[0])
        # no space inside an ideograph run, and punctuation always attaches
        gap = "" if a == "punct" or b == "punct" or (a == b == "ideo") else " "
        if text_mm(ax, cur + gap + tok, size) <= width_mm:
            cur += gap + tok
        else:
            lines.append(cur)
            cur = tok
    if cur:
        lines.append(cur)
    return lines


def box(ax, x, y, w, h, fill, stroke):
    ax.add_patch(FancyBboxPatch((x, y + h), w, -h,
                                boxstyle="round,pad=0,rounding_size=1.1",
                                linewidth=0.5, edgecolor=stroke, facecolor=fill,
                                zorder=2))


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=4.0, linewidth=0.6, color=RULE,
                                 shrinkA=0.0, shrinkB=0.0, zorder=1))


def line(ax, x0, y0, x1, y1):
    ax.plot([x0, x1], [y0, y1], color=RULE, linewidth=0.6, zorder=1,
            solid_capstyle="butt")


def plan_content(ax, C):
    data = [(t, wrap_mm(ax, s, FS_DET, INNER)) for t, s in C["data"]]
    tracks = []
    for header, mods in C["tracks"]:
        out = []
        for title, det, badge in mods:
            bw = text_mm(ax, badge, FS_DET) + 2.0 if badge else 0.0
            tl = wrap_mm(ax, title, FS_BOX, INNER - bw)
            dl = wrap_mm(ax, det, FS_DET, INNER)
            out.append((badge, tl, dl))
        tracks.append((header, out))
    gate_w = W - 2 * PAD - 6.0
    gate = [t for ln in C["gate_lines"] for t in wrap_mm(ax, ln, FS_DET, gate_w)]
    conc = []
    for b in C["conc_bullets"]:
        for i, t in enumerate(wrap_mm(ax, b, FS_DET, gate_w - 3.6)):
            conc.append((t, i == 0))
    return data, tracks, gate, conc


def draw(ax, C, plan, height):
    data, tracks, gate, conc = plan
    ax.set_xlim(0, W)
    ax.set_ylim(height, 0)
    ax.set_axis_off()

    ax.text(PAD, 7.6, C["title"], fontsize=FS_TITLE, color=INK, va="baseline")
    ax.text(PAD, 13.0, C["subtitle"], fontsize=FS_SUB, color=MUTED,
            va="baseline")
    ax.text(PAD, 20.8, C["data_head"], fontsize=FS_SEC, color=INK,
            va="baseline")

    dh = 2.6 + (1 + max(len(d) for _, d in data)) * LH
    dgap, dtop = 2.0, 23.6
    for i, (title, dl) in enumerate(data):
        r, c = divmod(i, 3)
        x, y = COLX[c], dtop + r * (dh + dgap)
        box(ax, x, y, COLW, dh, F_DATA, S_DATA)
        yy = y + (dh - (1 + len(dl)) * LH) / 2.0 + LH * 0.95
        ax.text(x + 2.2, yy, title, fontsize=FS_BOX, color=INK, va="baseline")
        yy += LH
        for ln in dl:
            ax.text(x + 2.2, yy, ln, fontsize=FS_DET, color=MUTED,
                    va="baseline")
            yy += LH
    dbot = dtop + dh + dgap + dh

    bus = dbot + 5.5
    line(ax, COLC[0], dbot, COLC[0], bus)
    line(ax, COLC[2], dbot, COLC[2], bus)
    line(ax, COLC[0], bus, COLC[2], bus)
    line(ax, COLC[1], dbot, COLC[1], bus)

    hh, htop = 7.4, bus + 3.0
    btop = htop + hh + 4.5
    bh = 2.6 + max(len(t) + len(d) for _, m in tracks for _, t, d in m) * LH
    bgap = 3.0
    for c, (header, mods) in enumerate(tracks):
        arrow(ax, COLC[c], bus, COLC[c], htop)
        box(ax, COLX[c], htop, COLW, hh, F_TRK, S_TRK)
        ax.text(COLC[c], htop + 4.9, header, fontsize=FS_TRK, color=T_TRK,
                va="baseline", ha="center")
        for i, (badge, tl, dl) in enumerate(mods):
            y = btop + i * (bh + bgap)
            if i:
                arrow(ax, COLC[c], y - bgap, COLC[c], y)
            box(ax, COLX[c], y, COLW, bh, "#FFFFFF", S_TRK)
            yy = y + (bh - (len(tl) + len(dl)) * LH) / 2.0 + LH * 0.95
            first = yy
            for ln in tl:
                ax.text(COLX[c] + 2.2, yy, ln, fontsize=FS_BOX, color=INK,
                        va="baseline")
                yy += LH
            for ln in dl:
                ax.text(COLX[c] + 2.2, yy, ln, fontsize=FS_DET, color=MUTED,
                        va="baseline")
                yy += LH
            if badge:
                ax.text(COLX[c] + COLW - 2.2, first, badge, fontsize=FS_DET,
                        color=BADGE, va="baseline", ha="right")
    colbot = btop + 4 * (bh + bgap) + bh

    gtop = colbot + 10.0
    gh = 3.0 + (1 + len(gate)) * LH
    for c in range(3):
        arrow(ax, COLC[c], colbot, COLC[c], gtop)
    box(ax, PAD, gtop, W - 2 * PAD, gh, F_GATE, S_GATE)
    ax.text(PAD + 3.0, gtop + 2.6 + LH * 0.95, C["gate_title"],
            fontsize=FS_TRK, color=T_GATE, va="baseline")
    yy = gtop + 2.6 + LH * 1.95
    for ln in gate:
        ax.text(PAD + 3.0, yy, ln, fontsize=FS_DET, color=INK, va="baseline")
        yy += LH

    ctop = gtop + gh + 9.0
    ch = 3.0 + (1 + len(conc)) * LH
    arrow(ax, W / 2.0, gtop + gh, W / 2.0, ctop)
    box(ax, PAD, ctop, W - 2 * PAD, ch, F_CONC, S_CONC)
    ax.text(PAD + 3.0, ctop + 2.6 + LH * 0.95, C["conc_title"],
            fontsize=FS_TRK, color=T_CONC, va="baseline")
    yy = ctop + 2.6 + LH * 1.95
    for ln, is_first in conc:
        if is_first:
            ax.text(PAD + 3.0, yy, "\u2022", fontsize=FS_DET, color=T_CONC,
                    va="baseline")
        ax.text(PAD + 5.6, yy, ln, fontsize=FS_DET, color=INK, va="baseline")
        yy += LH

    return ctop + ch + 7.0


def build(C):
    plt.rcParams["font.family"] = [C["font"]]
    global _probe
    _probe = None
    fig = plt.figure(figsize=(W / 25.4, MAX_H / 25.4))
    fig.canvas.draw()
    ax = fig.add_axes([0, 0, 1, 1])

    plan = plan_content(ax, C)
    height = round(draw(ax, C, plan, MAX_H) + 1.0, 1)
    assert height <= MAX_H, "%s: %.1f mm exceeds the %.0f mm cap" % (
        C["name"], height, MAX_H)

    fig.set_size_inches(W / 25.4, height / 25.4, forward=True)
    draw(ax, C, plan, height)

    stem = os.path.join(OUT, C["name"])
    fig.savefig(stem + ".pdf")
    fig.savefig(stem + ".svg")
    fig.savefig(stem + ".png", dpi=300)
    fig.savefig(stem + ".tiff", dpi=300,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    print("%-28s %.0f x %.1f mm (cap %.0f)" % (C["name"], W, height, MAX_H))
    for ext in (".pdf", ".svg", ".png", ".tiff"):
        print("    %-6s %8.1f kB" % (ext, os.path.getsize(stem + ext) / 1e3))
    return height


def main():
    os.makedirs(OUT, exist_ok=True)
    want = sys.argv[1:] or ["en", "zh"]
    for key in want:
        build({"en": EN, "zh": ZH}[key])


if __name__ == "__main__":
    main()
