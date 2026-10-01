# -*- coding: utf-8 -*-
"""Figure 6 has been re-drawn with corrected in-figure labels, so the two sentences
that described the *old* labels must go.  Exactly-once replacement, abort on mismatch.

Run: python scripts/63_fix_fig6_legend.py
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
PKG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "M2-M9深化研究稿")
EN = os.path.join(PKG, "manuscript_EN.md")
CN = os.path.join(PKG, "manuscript_CN.md")

OPS = [
    (EN,
     "so the empirical P lies below the resolution of 200 permutations (P < 0.005), and the panel annotation prints only the rounded value 0.0000.",
     "so the empirical P lies below the resolution of 200 permutations; the panel is annotated accordingly as P < 0.005.",
     "EN-Fig6: drop obsolete 'prints 0.0000'"),
    (CN,
     "故经验 P 低于 200 次置换的分辨率（P < 0.005）；面板标注为取整后的数值 0.0000。",
     "故经验 P 低于 200 次置换的分辨率；面板已据此标注为 P < 0.005。",
     "CN-图6: 删去已不适用的 0.0000 说明"),
    # the panels now carry the denominators, so keep the legend aligned with them
    (EN,
     "the pipeline recovered 5 of the 9 pre-specified activation/exhaustion factors (PRDM1, NR4A1, BATF, MAF, TBX21). *TOX* and *IKZF2* are absent from the CollecTRI collection used, so 7 of the 9 were testable and the panel prints the denominator 7.",
     "the pipeline recovered 5 of the 9 pre-specified activation/exhaustion factors (PRDM1, NR4A1, BATF, MAF, TBX21). *TOX* and *IKZF2* are absent from the CollecTRI collection used, so 7 of the 9 were testable, and the panel states 5/7 testable (5/9 pre-specified).",
     "EN-Fig6: align panel c wording"),
    (CN,
     "流程命中了 9 个预设活化/耗竭因子中的 5 个（PRDM1、NR4A1、BATF、MAF、TBX21）；*TOX* 与 *IKZF2* 不在所用 CollecTRI 集合内，故 9 个中有 7 个可检验，面板标注的分母为 7。",
     "流程命中了 9 个预设活化/耗竭因子中的 5 个（PRDM1、NR4A1、BATF、MAF、TBX21）；*TOX* 与 *IKZF2* 不在所用 CollecTRI 集合内，故 9 个中有 7 个可检验；面板标注为 5/7（可检验），并附 5/9（预设）。",
     "CN-图6: 与面板 c 对齐"),
]


def main():
    bad = []
    for path, old, new, label in OPS:
        t = open(path, encoding="utf-8").read()
        n = t.count(old)
        if n != 1:
            bad.append("%s -> hit %d" % (label, n))
            continue
        open(path, "w", encoding="utf-8").write(t.replace(old, new, 1))
        print("OK  " + label)
    if bad:
        print("FAILED:", bad)
        sys.exit(1)


if __name__ == "__main__":
    main()
