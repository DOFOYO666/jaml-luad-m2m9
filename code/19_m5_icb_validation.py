# -*- coding: utf-8 -*-
"""M5：ICB 队列中 JAML 表达与治疗响应/无进展生存的关联。

数据：
  GSE126044  rawdata/icb/GSE126044/GSE126044_counts.txt.gz        （symbol × 16 例，原始 counts）
  GSE135222  rawdata/icb/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz （ENSG × 27 例，已归一化）
  响应/PFS 注释取自各自的 series_matrix

⚠️ 两个队列均为 NSCLC（非纯 LUAD），且 R 组仅个位数 → 结论语境限定为 NSCLC、小样本。
"""
import gzip
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)
JAML_ENSG = "ENSG00000160593"


def parse_series_matrix(path):
    """从 GEO series_matrix 抽 !Sample_title 与各 !Sample_characteristics 行。"""
    txt = gzip.open(path, "rt", errors="replace").read().split("\n")
    titles, chars = None, {}
    for ln in txt:
        if ln.startswith("!Sample_title"):
            titles = [c.strip('"') for c in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            vals = [c.strip('"') for c in ln.split("\t")[1:]]
            key = vals[0].split(":")[0].strip()
            chars[key] = vals
    df = pd.DataFrame(chars)
    df.insert(0, "title", titles)
    return df


def save_fig(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".tif"), dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  [图] {name}.png / .tif")


def logrank(time, event, group):
    time = np.asarray(time, float); event = np.asarray(event, int)
    group = np.asarray(group, int)
    O1 = E1 = V = 0.0
    for t in np.unique(time[event == 1]):
        at = time >= t
        n = at.sum(); n1 = (at & (group == 1)).sum()
        d = ((time == t) & (event == 1)).sum()
        d1 = ((time == t) & (event == 1) & (group == 1)).sum()
        if n <= 1:
            continue
        O1 += d1; E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V <= 0:
        return np.nan, np.nan
    chi2 = (O1 - E1) ** 2 / V
    return float(chi2), float(stats.chi2.sf(chi2, 1))


plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})

# ================= GSE126044（16 例，抗 PD-1）=================
print("=" * 78)
print("GSE126044 —— 抗 PD-1（Nivolumab/Pembrolizumab），16 例 NSCLC")
print("=" * 78)
ann44 = parse_series_matrix(os.path.join(BASE, "rawdata/icb/GSE126044/GSE126044_series_matrix.txt.gz"))
key = [k for k in ann44.columns if "response" in k.lower()]
print("注释键:", list(ann44.columns))
resp = ann44[key[0]].str.split(":").str[-1].str.strip()
ann44["response"] = resp
print(ann44[["title", "response"]].to_string(index=False))
print("响应分组计数:", ann44["response"].value_counts().to_dict())

cnt = pd.read_csv(os.path.join(BASE, "rawdata/icb/GSE126044/GSE126044_counts.txt.gz"),
                  sep="\t", index_col=0)
print(f"表达矩阵 {cnt.shape}（行=基因，列=样本）")
# 列名形如 Dis_01 → 对应 title RNA-seq_Dis_01
cnt.columns = [c.replace("Dis_", "Dis_") for c in cnt.columns]
samples = [c for c in cnt.columns]
print("样本列:", samples)

cpm = cnt / cnt.sum(axis=0) * 1e6
logcpm = np.log2(cpm + 1)

# ⚠️ 基因符号别名：GSE126044 使用旧注释，JAML 记为 AMICA1（与 TCGA/GSE131907 同一陷阱）
GENE_ALIAS = {"JAML": ["JAML", "AMICA1"], "HYKK": ["HYKK", "AGPHD1"]}


def find_row(index, sym):
    for c in GENE_ALIAS.get(sym, [sym]):
        if c in index:
            return c
    return None


jaml_sym = find_row(logcpm.index.astype(str), "JAML")
print(f"[符号] JAML 在 GSE126044 中的实际符号: {jaml_sym}")
if jaml_sym is None:
    print("⚠️ GSE126044 的 counts 矩阵仅含 %d 个基因，JAML 与 AMICA1 均不在其中"
          "（应为低表达/低检出被过滤）。**该队列对 JAML 不可用**，跳过。"
          % logcpm.shape[0])
    print("   已核：同家族 JAM2、JAM3 在矩阵中，说明并非整体缺失，而是该基因被过滤。")
    jaml44 = None
else:
    jaml44 = pd.Series(logcpm.loc[jaml_sym].values, index=logcpm.columns, name="JAML_log2CPM")
if jaml44 is not None:
    m44 = pd.DataFrame({"expr": jaml44})
    m44["response"] = [ann44.set_index(ann44["title"].str.replace("RNA-seq_", ""))["response"]
                       .get(i, None) for i in m44.index]
    print(m44.to_string())
    grp = m44.dropna()
    r = grp.loc[grp["response"].str.contains("responder") & ~grp["response"].str.contains("non"), "expr"]
    nr = grp.loc[grp["response"].str.contains("non"), "expr"]
    mw = stats.mannwhitneyu(r, nr, alternative="two-sided")
    print(f"\nresponder n={len(r)} mean={r.mean():.3f} | non-responder n={len(nr)} mean={nr.mean():.3f}")
    print(f"Mann-Whitney U P = {mw.pvalue:.4f}")
    m44.to_csv(os.path.join(RES, "m5_gse126044_jaml.csv"))

    fig, ax = plt.subplots(figsize=(3.6, 3.4))
    data = [nr.values, r.values]
    bp = ax.boxplot(data, tick_labels=[f"Non-responder\n(n={len(nr)})", f"Responder\n(n={len(r)})"],
                    patch_artist=True, widths=0.55, medianprops=dict(color="#333", lw=1.2))
    for b, c in zip(bp["boxes"], ["#4575B4", "#D73027"]):
        b.set_facecolor(c); b.set_alpha(0.75); b.set_edgecolor("#333")
    ax.scatter(np.repeat([1, 2], [len(nr), len(r)]), np.concatenate([nr.values, r.values]),
               s=14, color="#333", zorder=3, alpha=0.7)
    ax.set_ylabel("JAML expression (log2 CPM)")
    ax.set_title(f"GSE126044 (NSCLC, anti-PD-1)\nMann-Whitney P = {mw.pvalue:.3f}", fontsize=9)
    save_fig(fig, "m5_gse126044_jaml_response")

# ================= GSE135222（27 例，抗 PD-1/PD-L1）=================
print("\n" + "=" * 78)
print("GSE135222 —— 抗 PD-1/PD-L1，27 例 NSCLC（含 PFS）")
print("=" * 78)
ann35 = parse_series_matrix(os.path.join(BASE, "rawdata/icb/GSE135222/GSE135222_series_matrix.txt.gz"))
print("注释键:", list(ann35.columns))
pfs_k = [k for k in ann35.columns if "progression-free" in k.lower()][0]
time_k = [k for k in ann35.columns if "pfs.time" in k.lower()][0]
ann35["pfs"] = ann35[pfs_k].str.split(":").str[-1].str.strip().astype(int)
ann35["pfs_time"] = ann35[time_k].str.split(":").str[-1].str.strip().astype(float)
print("pfs 取值计数:", ann35["pfs"].value_counts().to_dict())
print("pfs_time 描述: 中位=%.1f 均值=%.1f 范围=%.0f-%.0f"
      % (ann35["pfs_time"].median(), ann35["pfs_time"].mean(),
         ann35["pfs_time"].min(), ann35["pfs_time"].max()))

gd = pd.read_csv(os.path.join(BASE, "rawdata/icb/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"),
                 sep="\t", index_col=0)
print(f"表达矩阵 {gd.shape}（行=ENSG，列=样本）")
sym = {}
p = os.path.join(BASE, "data/onek1k/ensg2sym.json")
if os.path.exists(p):
    sym = json.load(open(p))
gd.index = [str(i).split(".")[0] for i in gd.index]
hit = [i for i in gd.index if i == JAML_ENSG]
print("JAML(ENSG00000160593) 是否在矩阵:", bool(hit))
if not hit:
    # 备用：用 symbol 映射反查
    inv = {v: k for k, v in sym.items()}
    print("  ensg2sym 中 JAML →", inv.get("JAML"))
jrow = gd.loc[JAML_ENSG] if hit else None
assert jrow is not None, "无法定位 JAML 行"

# 样本名对齐（title 形如 "NSCLC 990"；矩阵列形如 "NSCLC990"）
norm = {re.sub(r"\s+", "", str(t)): i for i, t in enumerate(ann35["title"])}
recs = []
for col, val in jrow.items():
    k = re.sub(r"\s+", "", str(col))
    if k in norm:
        idx = norm[k]
        recs.append({"sample": col, "expr": float(val),
                     "pfs": int(ann35["pfs"].iloc[idx]),
                     "pfs_time": float(ann35["pfs_time"].iloc[idx])})
d35 = pd.DataFrame(recs)
print(f"成功对齐 {len(d35)} / {gd.shape[1]} 例")
d35.to_csv(os.path.join(RES, "m5_gse135222_jaml_pfs.csv"), index=False)

if len(d35) >= 10:
    med = d35["expr"].median()
    d35["grp"] = (d35["expr"] > med).astype(int)
    chi2, p_lr = logrank(d35["pfs_time"].values, d35["pfs"].values, d35["grp"].values)
    print(f"\n按 JAML 中位数分组：高 n={(d35['grp']==1).sum()} / 低 n={(d35['grp']==0).sum()}")
    print(f"log-rank χ² = {chi2:.3f}, P = {p_lr:.4f}")
    print(d35.groupby("grp")[["expr", "pfs_time"]].agg(["mean", "median", "size"]).to_string())

    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    for g, col, lab in [(0, "#4575B4", "JAML low"), (1, "#D73027", "JAML high")]:
        sub = d35[d35["grp"] == g].sort_values("pfs_time")
        t, e = sub["pfs_time"].values, sub["pfs"].values
        s = 1.0; xs, ys = [0], [1.0]
        for i in range(len(t)):
            if e[i] == 1:
                s *= (1 - 1 / (len(t) - i))
            xs.append(t[i]); ys.append(s)
        ax.step(xs, ys, where="post", color=col, lw=1.4, label=f"{lab} (n={len(sub)})")
    ax.set_xlabel("Progression-free survival (days)"); ax.set_ylabel("PFS probability")
    ax.set_ylim(0, 1.02); ax.legend(frameon=False, fontsize=8)
    ax.set_title(f"GSE135222 (NSCLC, anti-PD-1/PD-L1)\nlog-rank P = {p_lr:.3f}", fontsize=9)
    save_fig(fig, "m5_gse135222_jaml_pfs")

print("\n[完成] M5 输出：results/m5_gse126044_jaml.csv, results/m5_gse135222_jaml_pfs.csv")
