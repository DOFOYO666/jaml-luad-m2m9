# -*- coding: utf-8 -*-
"""
03_tcga_survival.py - TCGA-LUAD 表达-生存预后验证

⚠️ 已弃用（superseded）—— 请改用 scripts/18_m2_tcga_clinical.py
   本脚本存在两处方法学缺陷，2026-09-29 审查发现：
   (1) `if gene not in expr.index: continue` 会**静默丢弃**符号不匹配的基因。
       TCGA/Xena HiSeqV2 使用较老注释：JAML→AMICA1、HYKK→AGPHD1，
       导致核心基因 JAML **从未被实际分析**（旧结果 tcga_survival_validation.csv
       只有 16 个基因，不含 JAML）。
   (2) `cox_hr_quick` 用"事件率之比 + 二项 z 检验"近似，并把结果写入名为
       `logrank_p` 的列 —— 既不是 log-rank 检验，也不是 Cox 回归。
   18_m2_tcga_clinical.py 已用正规 Mantel-Cox log-rank + 一维 Cox PH 重做，
   并加入基因别名解析。本文件仅保留以追溯历史结果。

原始功能：对 scMR 显著（基因 × 细胞类型），在 TCGA-LUAD 中验证表达-生存关联。

输入：
  results/scmr_luad_risk_sig.csv (可选；无则分析全部基因)
  data/luad_survival/TCGA.LUAD.HiSeqV2.gz
  data/luad_survival/TCGA.LUAD.clinicalMatrix.txt
输出：
  results/tcga_survival_validation.csv
"""
import os, gzip
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "luad_survival")
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

EXPR = os.path.join(DATA, "TCGA.LUAD.HiSeqV2.gz")
CLIN = os.path.join(DATA, "TCGA.LUAD.clinicalMatrix.txt")
SIG = os.path.join(RES, "scmr_luad_risk_sig.csv")


def load_expr(path):
    df = pd.read_csv(path, sep="\t", compression="gzip", index_col=0)
    return df  # 行=基因, 列=样本


def load_clin(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df


def cox_hr_quick(time, event, group):
    """简化 Cox：用 log-rank 与二分组 HR 近似（基于正态近似）"""
    # 简化：用 log-rank 检验 P 值（无 R 依赖的保守实现）
    # 分组 0/1 的比较
    t0 = time[group == 0]; e0 = event[group == 0]
    t1 = time[group == 1]; e1 = event[group == 1]
    if len(t0) < 5 or len(t1) < 5 or e0.sum() == 0 or e1.sum() == 0:
        return np.nan, np.nan
    # 用 log-rank 的近似 z 检验（Mantel-Haenszel 简化版）
    # 简单近似：比较两组的事件率
    rate0 = e0.sum() / t0.sum() if t0.sum() > 0 else np.nan
    rate1 = e1.sum() / t1.sum() if t1.sum() > 0 else np.nan
    hr = rate1 / rate0 if rate0 and rate0 > 0 else np.nan
    # 简化 p 值：两组事件比例的比较（卡方近似）
    n0, n1 = len(t0), len(t1)
    p_e0, p_e1 = e0.sum() / n0, e1.sum() / n1
    pool = (e0.sum() + e1.sum()) / (n0 + n1)
    if pool == 0 or pool == 1:
        return hr, np.nan
    se = np.sqrt(pool * (1 - pool) * (1 / n0 + 1 / n1))
    z = (p_e1 - p_e0) / se if se > 0 else np.nan
    p = 2 * (1 - stats.norm.cdf(abs(z))) if not np.isnan(z) else np.nan
    return hr, p


def main():
    sig = None
    if os.path.exists(SIG):
        sig = pd.read_csv(SIG)
        gcol = "gene" if "gene" in sig.columns else "gene_symbol"
        genes = sorted(set(sig[gcol]))
        print(f"[输入] scMR 显著基因 {len(genes)} 个")
    else:
        print("[提示] 未找到 scMR 显著结果，跳过生存验证")
        return

    expr = load_expr(EXPR)
    clin = load_clin(CLIN)
    print(f"[TCGA] 表达矩阵 {expr.shape}, 临床 {clin.shape}")

    # 临床列检测
    os_time_col = None
    os_event_col = None
    for c in clin.columns:
        cu = c.upper()
        if cu in ("OS.TIME", "OS_TIME", "DAYS_TO_DEATH", "OS_DAYS") and os_time_col is None:
            os_time_col = c
        if cu in ("OS", "OS_STATUS", "VITAL_STATUS", "OS_EVENT") and os_event_col is None:
            os_event_col = c
    if os_time_col is None:
        # Xena: OS.time / OS
        os_time_col = "OS.time" if "OS.time" in clin.columns else None
        os_event_col = "OS" if "OS" in clin.columns else None
    if not os_time_col or not os_event_col:
        print(f"[错误] 未识别生存列，现有列示例: {list(clin.columns)[:30]}")
        return
    print(f"[TCGA] 生存列: time={os_time_col}, event={os_event_col}")

    clin = clin.copy()
    # vital_status 为文本 (LIVING/DECEASED)；OS 时间按状态合成
    vs = clin[os_event_col].astype(str).str.upper()
    death = vs.str.contains("DECEASED", na=False)
    dd = pd.to_numeric(clin.get("days_to_death"), errors="coerce")
    dlf = pd.to_numeric(clin.get("days_to_last_followup"), errors="coerce")
    clin["_time"] = np.where(death, dd, dlf)
    clin["_event"] = death.astype(int)
    clin["_time"] = pd.to_numeric(clin["_time"], errors="coerce")
    clin = clin.dropna(subset=["_time", "_event"])
    clin = clin[clin["_time"] > 0]

    # 样本对齐（Xena 样本 ID 形如 TCGA-XX-XXXX-01，临床同）
    common = sorted(set(expr.columns) & set(clin.index))
    print(f"[TCGA] 共同样本 {len(common)}")

    rows = []
    for gene in genes:
        if gene not in expr.index:
            continue
        vals = expr.loc[gene, common]
        sub = clin.loc[common, ["_time", "_event"]].copy()
        sub["expr"] = pd.to_numeric(vals, errors="coerce").values
        sub = sub.dropna(subset=["expr"])
        if len(sub) < 30:
            continue
        med = sub["expr"].median()
        sub["group"] = (sub["expr"] > med).astype(int)
        hr, p = cox_hr_quick(sub["_time"].values, sub["_event"].values, sub["group"].values)
        rows.append({"gene": gene, "n": len(sub), "median_expr": med,
                     "hr_high_vs_low": hr, "logrank_p": p})
        if len(rows) % 100 == 0:
            print(f"  已处理 {len(rows)} 基因")

    res = pd.DataFrame(rows)
    if not res.empty:
        res["fdr"] = stats.false_discovery_control(res["logrank_p"], method="bh")
        res = res.sort_values("logrank_p")
        res.to_csv(os.path.join(RES, "tcga_survival_validation.csv"), index=False)
        print(f"\n[完成] 验证 {len(res)} 基因")
        print(res.head(20).to_string())
    else:
        print("无有效验证结果")


if __name__ == "__main__":
    main()
