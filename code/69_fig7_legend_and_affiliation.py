# -*- coding: utf-8 -*-
"""收尾改稿（第三轮）：

1. Figure 7a 换用最终网络配置（ext3074 / 3,074 基因 / 70 个可扰动 TF）后，重写中英母稿的图 7 图注。
   —— 同时纠正上一轮发现的第三处图注/图不符：图注写"全部可扰动因子 ranked"，而旧图只画了前 16 根柱。
2. 科室英文名统一为用户确认的写法（原 M2–M9 稿的两封投稿信仍是旧写法）。
3. 代码发布树 .zenodo.json 的 affiliation 与 CITATION.cff 对齐（补上科室）。

所有替换都是"整串精确 + 恰好命中一次"，否则中止且不写盘；已应用过的替换在重跑时跳过（幂等）。
改前对每个文件留 .bak_fig7_<时间戳> 备份。
"""
import io
import os
import re
import shutil
import sys
import time

sys.stdout.reconfigure(encoding="utf-8")

BASE = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化"
PKG = os.path.join(BASE, "M2-M9深化研究稿")
REL = os.path.join(BASE, "JAML_M2M9_code_release")
STAMP = time.strftime("%Y%m%d_%H%M%S")

LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def sub_once(path, old, new, label, expected=1):
    """old -> new，恰好 expected 次；若已应用（old 不在而 new 在）则跳过。

    ⚠️ 陷阱（本脚本踩过一次）：当 old 是 new 的子串时（如"医院名"⊂"科室+医院名"），
    幂等判据 `old not in t` 永远为假 —— 第二次运行会把 new 里的 old 再替换一遍，
    结果重复插入。凡属这种情况，调用方必须**自己判已应用**（见下面 .zenodo.json 那一步）。
    """
    t = io.open(path, encoding="utf-8").read()
    if old not in t and new in t:
        log("  ..  already applied: " + label)
        return
    n = t.count(old)
    if n != expected:
        raise SystemExit("[%s] 命中 %d 次（应为 %d）：%r" % (label, n, expected, old[:80]))
    bak = path + ".bak_fig7_" + STAMP
    if not os.path.exists(bak):
        shutil.copy2(path, bak)
    io.open(path, "w", encoding="utf-8").write(t.replace(old, new))
    log("  OK  %s  （命中 %d 次）" % (label, n))


# ---------------------------------------------------------------- 1. 图 7 图注 --
EN = os.path.join(PKG, "manuscript_EN.md")
CN = os.path.join(PKG, "manuscript_CN.md")

EN_FIG7_OLD = ("(a) Mean |Δ*JAML*| after knocking out each of 39 candidate transcription factors "
               "in the CD4+ subset (2,842 cells) using CellOracle with bagging = 20; the dashed "
               "line is the global noise ceiling (9.59 × 10⁻³), defined as the largest randomized "
               "reading among non-degenerate factors.")
EN_FIG7_NEW = ("(a) Mean |Δ*JAML*| after knocking out each of the 70 transcription factors that "
               "could be perturbed in the final network configuration (3,074 genes; bagging = 20; "
               "CD4+ subset, 2,842 cells; CellOracle 0.20.0), ranked; the dashed line is the global "
               "noise ceiling (9.59 × 10⁻³) derived from the paired randomized-network control run "
               "on that same configuration, and *ID2* (1.75 × 10⁻²) is the only factor above it, at "
               "5.3× the second-ranked factor. Factors whose knockout left *JAML* unchanged (Δ = 0) "
               "are drawn as open markers at the lower axis limit, and dark markers denote the "
               "factors for which a paired randomized control was run. The complete "
               "three-configuration table, including the factors assessed only in the intermediate "
               "configuration, is in Additional file 7a.")

CN_FIG7_OLD = ("（a）用 CellOracle（bagging = 20）在 CD4⁺子集（2,842 细胞）中逐个敲除 39 个候选转录因子后的"
               "平均 |Δ*JAML*|；虚线为全局噪声上界（9.59×10⁻³），即全部非退化因子中最大的随机读数。")
CN_FIG7_NEW = ("（a）在最终网络配置（3,074 基因；bagging = 20；CD4⁺子集 2,842 细胞；CellOracle 0.20.0）中"
               "可扰动的 70 个转录因子，按敲除后的平均 |Δ*JAML*| 排序；虚线为在同一配置上做的配对打乱网络"
               "对照所得的全局噪声上界（9.59×10⁻³），只有 *ID2*（1.75×10⁻²）越线，为第二名的 5.3 倍。"
               "敲除后 *JAML* 无变化者（Δ = 0）画为轴内下限处的空心标记；深色标记表示该因子做过配对打乱"
               "对照。三套网络配置的完整对照表（含只在中间配置中评估的因子）见附加文件 7a。")

# ---------------------------------------------------------------- 2. 科室英文名 --
CLE_EN = os.path.join(PKG, "Cover_Letter_EN.md")
CLE_CN = os.path.join(PKG, "Cover_Letter_CN.md")
AFFIL_OLD_EN = ("Department of Respiratory Medicine\n"
                "Shenzhen Nanshan People's Hospital (Shenzhen University Affiliated Nanshan Hospital)")
AFFIL_NEW_EN = ("Department of Pulmonary and Critical Care Medicine\n"
                "Shenzhen Nanshan District People's Hospital (Shenzhen University Affiliated "
                "Nanshan Hospital)")
AFFIL_OLD_CN = "深圳市南山区人民医院（深圳大学附属南山医院）呼吸内科"
AFFIL_NEW_CN = "深圳市南山区人民医院（深圳大学附属南山医院）呼吸与危重症医学科"

# ---------------------------------------------------------------- 3. Zenodo 元数据 --
ZEN = os.path.join(REL, ".zenodo.json")
ZEN_OLD = "Shenzhen Nanshan District People's Hospital (Shenzhen University Affiliated Nanshan Hospital)"
ZEN_NEW = ("Department of Pulmonary and Critical Care Medicine, Shenzhen Nanshan District People's "
           "Hospital (Shenzhen University Affiliated Nanshan Hospital)")


def main():
    log("[1] 图 7 图注（换用最终配置）")
    sub_once(EN, EN_FIG7_OLD, EN_FIG7_NEW, "EN Figure 7 panel (a)")
    sub_once(CN, CN_FIG7_OLD, CN_FIG7_NEW, "CN 图 7 面板 (a)")

    log("[2] 科室英文名统一")
    sub_once(CLE_EN, AFFIL_OLD_EN, AFFIL_NEW_EN, "Cover_Letter_EN.md 单位行")
    sub_once(CLE_CN, AFFIL_OLD_CN, AFFIL_NEW_CN, "Cover_Letter_CN.md 单位行")

    log("[3] Zenodo 元数据单位与 CITATION.cff 对齐")
    # ZEN_OLD ⊂ ZEN_NEW，故不能交给 sub_once 判幂等：这里显式以 ZEN_NEW 的出现次数为准。
    _zen = io.open(ZEN, encoding="utf-8").read()
    if _zen.count("Department of Pulmonary and Critical Care Medicine, Department of Pulmonary"):
        raise SystemExit(".zenodo.json 存在重复的科室前缀，请先修复")
    if _zen.count(ZEN_NEW) == 7:
        log("  ..  already applied: .zenodo.json affiliation")
    else:
        sub_once(ZEN, ZEN_OLD, ZEN_NEW, ".zenodo.json affiliation", expected=7)

    log("[4] 回读校验")
    en = io.open(EN, encoding="utf-8").read()
    cn = io.open(CN, encoding="utf-8").read()
    zen = io.open(ZEN, encoding="utf-8").read()
    checks = [
        ("EN 图注已写 70 个可扰动因子", "each of the 70 transcription factors" in en),
        ("EN 图注已写上界同源", "derived from the paired randomized-network control run" in en),
        ("EN 图注已写 5.3 倍", "at 5.3× the second-ranked factor" in en),
        ("EN 图注已指向 Additional file 7a", "is in Additional file 7a." in en),
        ("EN 旧图注（39 个 / 无 panel b 说明）已消失", "each of 39 candidate transcription factors" not in en),
        ("CN 图注已写 70 个可扰动因子", "可扰动的 70 个转录因子" in cn),
        ("CN 图注已写 5.3 倍", "为第二名的 5.3 倍" in cn),
        ("CN 图注已指向附加文件 7a", "见附加文件 7a" in cn),
        ("CN 旧图注已消失", "逐个敲除 39 个候选转录因子" not in cn),
        ("EN 投稿信旧科室名已消失", "Department of Respiratory Medicine" not in
         io.open(CLE_EN, encoding="utf-8").read()),
        ("CN 投稿信旧科室名已消失", "呼吸内科" not in io.open(CLE_CN, encoding="utf-8").read()),
        ("两封投稿信均为确认后的写法",
         "Department of Pulmonary and Critical Care Medicine" in io.open(CLE_EN, encoding="utf-8").read()
         and "呼吸与危重症医学科" in io.open(CLE_CN, encoding="utf-8").read()),
        # 注意：ZEN_NEW 里包含 ZEN_OLD 作为子串，故不能用 "ZEN_OLD 不在文中" 作判据；
        # 正确判据是带科室前缀的次数 = 不带前缀的医院串次数 = 作者数（7）。
        (".zenodo.json 的 7 条 affiliation 均带科室前缀",
         zen.count("Department of Pulmonary and Critical Care Medicine, Shenzhen Nanshan District") == 7
         and zen.count(ZEN_OLD) == 7),
    ]
    bad = [k for k, v in checks if not v]
    for k, v in checks:
        log(("  OK   " if v else "  FAIL ") + k)
    if bad:
        raise SystemExit("校验未过：%s" % bad)
    log("\nDONE")


if __name__ == "__main__":
    main()
