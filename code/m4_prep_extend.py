# -*- coding: utf-8 -*-
"""
M4 输入扩展（第二轮补齐）：把 base GRN 中有 motif 证据、但不在当前 3013 基因里的
TF 从原始 29,634 基因 UMI 矩阵中取出，并入输入，使"因不在矩阵而跳过"的 TF 归零。

设计要点：
  * 只解析目标基因行（其余行只找第一个 tab），避免重复 m4_prep_bounded.py 的
    全量解析开销；单趟扫描。
  * 基因组与细胞序严格与现有 ext h5ad 对齐（按 barcode 匹配，全部命中才继续）。
  * 不做任何静默回退：目标基因若在矩阵中不存在，记录并显式报告。
"""
import os
import sys
import gzip
import json
import time

import numpy as np
import anndata as ad
import scipy.sparse as sp

D_ROOT = r"D:\workbuddy工作空间\JAML深度研究"
RES = os.path.join(D_ROOT, "results")
MAT = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化\rawdata\scRNA\GSE131907\GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
BASE = os.path.join(RES, "CD4T_celloracle_ext.h5ad")
OUT = os.path.join(RES, "CD4T_celloracle_ext2.h5ad")
QC = os.path.join(RES, "m4_prep_extend_qc.json")

# 需补的 61 个 TF（在 base GRN 的 TF 列中、但不在当前 3013 基因内）
TARGETS = ['MYB', 'RORC', 'EOMES', 'IRF6', 'IRF5', 'ZBTB7B', 'GFI1', 'BCL6', 'IRF8',
           'SPI1', 'CEBPA', 'ATF3', 'JDP2', 'MAFG', 'KLF4', 'KLF10', 'SP2', 'SP3', 'SP4',
           'EGR2', 'EGR3', 'RORB', 'NR1H3', 'STAT2', 'STAT5A', 'STAT6', 'SMAD2', 'SMAD3',
           'SMAD4', 'RUNX1', 'RUNX2', 'TCF3', 'TCF4', 'TCF12', 'TCF7L2', 'FOXO3', 'FOXO4',
           'IKZF4', 'ZEB1', 'SNAI1', 'SNAI2', 'TWIST1', 'ID1', 'ID4', 'HES1', 'HEY1',
           'NR2F6', 'NR2F1', 'NR2F2', 'ETS2', 'ELF4', 'FLI1', 'ERG', 'ETV6', 'ETV5',
           'GABPA', 'BHLHE41', 'ATF6', 'CREB1', 'ATF2', 'NFIL3']

# 额外强制的对照基因（用于确认扫描逻辑正确：它们应在矩阵中且表达模式已知）
SPIKE = ["JAML", "CXADR", "CD3E", "EPCAM"]


def avail_gb():
    import ctypes

    class M(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

    m = M(); m.dwLength = ctypes.sizeof(M)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ullAvailPhys / 1073741824


def main():
    t0 = time.time()
    a = avail_gb()
    print(f"[内存] 起始可用 {a:.2f} GB")
    if a < 0.40:
        sys.exit(f"可用内存仅 {a:.2f} GB < 0.40 GB，拒绝启动以免 OOM")

    print(f"[1] 载入基准 {BASE}")
    base = ad.read_h5ad(BASE)
    print(f"    shape {base.shape}；亚型 {sorted(set(map(str, base.obs['cell_subtype'])))}")
    keys = list(map(str, base.obs_names))
    have = set(map(str, base.var_names))
    for g in SPIKE:
        if g not in have:
            print(f"    [警告] 对照基因 {g} 不在基准中（仅用于校验，不影响流程）")

    # ---------------- 读表头，建立 barcode -> 列号 ----------------
    print("[2] 读取原始矩阵表头")
    with gzip.open(MAT, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    cells = header[1:]
    print(f"    列数 {len(cells)}；表头示例 {cells[:2]}")
    pos = {c: i for i, c in enumerate(cells)}
    missing_bc = [k for k in keys if k not in pos]
    if missing_bc:
        sys.exit(f"有 {len(missing_bc)} 个 barcode 在矩阵表头中找不到，例如 {missing_bc[:3]} —— 停止")
    col_idx = np.array([pos[k] for k in keys], dtype=np.int64)
    print(f"    barcode 全部命中 {len(keys)}/{len(keys)}；列号范围 {col_idx.min()}–{col_idx.max()}")

    # ---------------- 单趟扫描 ----------------
    want = set(TARGETS)
    print(f"[3] 单趟扫描 {MAT}")
    print(f"    目标基因 {len(want)} 个 + 对照 {len(SPIKE)} 个")
    found = {}
    n_gene = 0
    t_scan = time.time()
    with gzip.open(MAT, "rt", encoding="utf-8", errors="replace") as f:
        f.readline()                                     # 跳过表头
        for line in f:
            tab = line.find("\t")
            if tab <= 0:
                continue
            name = line[:tab]
            n_gene += 1
            if name in want or name in SPIKE:
                v = np.fromstring(line[tab + 1:], dtype=np.float32, sep="\t")
                if v.size != len(cells):
                    sys.exit(f"基因 {name} 解析出 {v.size} 个值，期望 {len(cells)} —— 格式异常，停止")
                found[name] = v[col_idx].copy()
                print(f"      [命中 {len(found):3d}] {name:10s} 阳性率 "
                      f"{100*(v[col_idx] > 0).mean():5.1f}%  均值 {v[col_idx].mean():7.3f}"
                      f"  ({time.time()-t_scan:.0f}s)", flush=True)
            if n_gene % 5000 == 0:
                print(f"      ...已扫描 {n_gene} 基因 ({time.time()-t_scan:.0f}s)", flush=True)
    print(f"[3] 扫描完成：{n_gene} 基因，用时 {time.time()-t_scan:.0f}s；命中 {len(found)}")
    print(f"[内存] 扫描后可用 {avail_gb():.2f} GB")

    hit = [g for g in TARGETS if g in found]
    miss = [g for g in TARGETS if g not in found]
    print(f"\n[3b] 目标 TF 命中 {len(hit)}/{len(TARGETS)}；缺失 {len(miss)} 个: {miss}")

    if not hit:
        sys.exit("没有任何目标 TF 被取到 —— 停止")

    # ---------------- 合并 ----------------
    print("[4] 合并到输入")
    Xb = base.X if sp.issparse(base.X) else sp.csr_matrix(base.X)
    Xn = np.stack([found[g] for g in hit], axis=0).astype(np.float32)   # 基因 × 细胞
    X = sp.hstack([Xb, sp.csr_matrix(Xn.T)], format="csr").astype(np.float32)

    new = ad.AnnData(X)
    new.obs = base.obs.copy()
    new.obs_names = base.obs_names.copy()
    new.var_names = list(map(str, base.var_names)) + hit
    new.layers["raw_count"] = X.copy()
    for k in base.obsm.keys():
        new.obsm[k] = np.asarray(base.obsm[k])
    new.write_h5ad(OUT, compression="gzip")
    print(f"[输出] {OUT}  ({new.shape[0]} 细胞 × {new.shape[1]} 基因)")

    # ---------------- 校验 ----------------
    print("[5] 校验")
    RC = new.layers["raw_count"]
    vn = list(map(str, new.var_names))
    for g in ["JAML"] + hit[:5]:
        col = RC[:, vn.index(g)]
        v = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
        print(f"    {g:10s} 阳性率 {100*(v > 0).mean():5.1f}%  均值 {v.mean():7.3f}  最大 {v.max():8.1f}")
    # 整体一致性：JAML 应与基准完全一致
    b_col = base.layers["raw_count"][:, list(map(str, base.var_names)).index("JAML")]
    b_v = np.asarray(b_col.todense()).ravel() if sp.issparse(b_col) else np.asarray(b_col).ravel()
    n_col = RC[:, vn.index("JAML")]
    n_v = np.asarray(n_col.todense()).ravel() if sp.issparse(n_col) else np.asarray(n_col).ravel()
    same = bool(np.array_equal(b_v, n_v))
    print(f"    JAML 与基准完全一致: {same}")
    if not same:
        sys.exit("JAML 列与基准不一致 —— 合并出错，停止")

    qc = {
        "base": os.path.basename(BASE),
        "output": os.path.basename(OUT),
        "shape": list(new.shape),
        "n_genes_added": len(hit),
        "genes_added": hit,
        "genes_missing_in_matrix": miss,
        "n_genes_scanned": int(n_gene),
        "jaml_identical_to_base": same,
        "elapsed_sec": round(time.time() - t0, 1),
        "avail_gb_end": round(avail_gb(), 2),
    }
    with open(QC, "w", encoding="utf-8") as f:
        json.dump(qc, f, ensure_ascii=False, indent=2)
    print(f"[输出] {QC}")
    print(f"\n总耗时 {time.time()-t0:.0f}s\nDONE")


if __name__ == "__main__":
    main()
