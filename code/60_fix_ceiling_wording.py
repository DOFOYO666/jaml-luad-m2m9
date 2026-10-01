# -*- coding: utf-8 -*-
"""Correction pass: the wording introduced in 58_apply_review2.py described the
global noise ceiling as 'conservative'.  That is wrong: the ceiling is an empirical
maximum over a finite number of randomized draws, so it sits BELOW the true null
maximum and the screen is therefore liberal (anti-conservative).  Fix in EN + CN.

Run: python scripts/60_fix_ceiling_wording.py
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")

PKG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "M2-M9深化研究稿")
EN = os.path.join(PKG, "manuscript_EN.md")
CN = os.path.join(PKG, "manuscript_CN.md")

OPS = [
    (EN,
     " Two properties of that ceiling are stated plainly. It is a single order statistic — here it is set by the largest of five randomized readings for one factor — so it is a conservative screen rather than a fitted null distribution. And with five randomized replicates per factor the smallest attainable one-sided empirical P value is 1/6 ≈ 0.17, so the procedure cannot in principle yield a significant P value: it is reported as a calibrated screen, not as a hypothesis test.",
     " Two properties of that ceiling are stated plainly. It is the largest randomized reading observed across all non-degenerate factors — an empirical maximum over a finite number of draws, not a fitted null distribution — and because a maximum taken over a finite sample sits below the true null maximum, the benchmark is biased low, which makes this screen liberal rather than conservative; the margin by which a factor clears it, rather than the bare fact of clearing it, is therefore what carries weight. Second, with five randomized replicates per factor the smallest attainable one-sided empirical P value is 1/6 ≈ 0.17, so the procedure cannot in principle yield a significant P value: it is reported as a calibrated screen, not as a hypothesis test.",
     "EN-2.9: ceiling is liberal, not conservative"),

    (EN,
     " That ceiling is a single order statistic, set by one factor's replicate, so it is a conservative screen rather than a fitted null, and with five replicates per factor no empirical P value below 0.17 is attainable.",
     " That ceiling is an empirical maximum over a finite number of randomized draws rather than a fitted null, so it is biased low and the screen is liberal rather than conservative; with five replicates per factor no empirical P value below 0.17 is attainable.",
     "EN-3.7: ceiling is liberal, not conservative"),

    (CN,
     "有关该上界有两点须明说：它是一个**单一顺序统计量**——此处由某一个 TF 五次随机读数中的最大值决定——因此是保守的筛选门槛，而非拟合出的零分布；又，每个 TF 仅有 5 次随机重复，单侧经验 P 值的最小可能取值为 1/6 ≈ 0.17，故本流程在原理上无法给出显著的 P 值，只能作为校准过的筛选，而非假设检验。",
     "有关该上界有两点须明说：它是**全部非退化因子中观测到的最大随机读数**，即在有限次抽取上取得的经验最大值，而非拟合出的零分布；由于有限样本的最大值低于真实的零分布最大值，该门槛被**低估**，即**偏宽松（liberal）而非保守**；因此真正有分量的是某个因子**超出该门槛的倍数**（*ID2* 为 1.8 倍），而不是仅仅看是否超过。第二，每个 TF 仅有 5 次随机重复，单侧经验 P 值的最小可能取值为 1/6 ≈ 0.17，故本流程在原理上无法给出显著的 P 值，只能作为校准过的筛选，而非假设检验。",
     "CN-2.9: 上界偏宽松而非保守"),

    (CN,
     "该上界是单一顺序统计量、由某一个 TF 的随机重复决定，属保守筛选门槛而非拟合零分布；每 TF 仅 5 次重复，无法得到低于 0.17 的经验 P 值。",
     "该上界是有限次随机抽取上的经验最大值、而非拟合零分布，故被低估、**偏宽松（liberal）而非保守**；每 TF 仅 5 次重复，无法得到低于 0.17 的经验 P 值。",
     "CN-3.7: 上界偏宽松而非保守"),
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
    for path in (EN, CN):
        t = open(path, encoding="utf-8").read()
        assert "conservative screen" not in t and "conservative" not in t.split("## L")[0].split("## 局")[0], path
        print("verified: %s no longer calls the ceiling conservative" % os.path.basename(path))


if __name__ == "__main__":
    main()
