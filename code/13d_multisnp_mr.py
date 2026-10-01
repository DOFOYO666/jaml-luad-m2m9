# -*- coding: utf-8 -*-
"""13c_ld_prune_mr.py - M1 正式版：LD 剪枝(r²<0.3) + 等位基因方向校正的多 SNP MR
改进点（对比 13b 初步版）：
  1. 等位基因方向校正：GWAS effect allele 与 OneK1K alt 对齐（含互补链翻转），消除 β 抵消
  2. LD 剪枝：1000G Phase3 v5b（本地 bgzip+tbi，纯 Python tabix）计算 r²，贪心剪枝 r²<0.3
  3. 剪枝后做 IVW / MR-Egger / 加权中位数
LD 参考：D:/WorkBuddy/单细胞数据孟德尔随机化/data/ld_ref/
"""
import gzip
import io
import os
import pickle
import re
import struct
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
ALL_DIR = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/onek1k/all"
LD_DIR = r"D:/WorkBuddy/单细胞数据孟德尔随机化/data/ld_ref38"
GWAS = os.path.join(BASE, "data", "luad_gwas", "GCST004744.h.tsv.gz")
SIG = os.path.join(RES, "scmr_luad_risk_sig.csv")
CACHE = os.path.join(BASE, ".workbuddy", "tmp", "qtd_cache")
WINDOW = 250_000
R2_THRESH = 0.3
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# ---------------- 本地 tabix（BGZF + tbi 区域读取） ----------------

def _file_range(path, start, end):
    with open(path, "rb") as f:
        f.seek(start)
        return f.read(end - start)

def parse_tbi(tbi_path):
    with open(tbi_path, "rb") as f:
        raw = f.read()
    dec = gzip.GzipFile(fileobj=io.BytesIO(raw)).read()
    off = 0
    def u32():
        nonlocal off
        v = struct.unpack("<I", dec[off:off+4])[0]; off += 4; return v
    def i32():
        nonlocal off
        v = struct.unpack("<i", dec[off:off+4])[0]; off += 4; return v
    def u64():
        nonlocal off
        v = struct.unpack("<Q", dec[off:off+8])[0]; off += 8; return v
    if dec[off:off+4] != b"TBI\x01":
        raise ValueError("bad tbi")
    off += 4
    n_ref = i32(); i32(); i32(); i32(); i32(); i32(); i32()  # n_ref + fmt/col_seq/col_beg/col_end/meta/skip
    l_nm = i32()
    names = dec[off:off+l_nm].decode("utf-8", "replace").split("\x00"); off += l_nm
    refs = {}
    for _ in range(n_ref):
        name = names[_] if _ < len(names) else str(_+1)
        bins = {}
        n_bin = i32()
        for _b in range(n_bin):
            bin_id = u32(); n_chunk = i32()
            chunks = [(u64(), u64()) for _c in range(n_chunk)]
            bins[bin_id] = chunks
        n_intv = i32()
        intervals = [u64() for _i in range(n_intv)]
        refs[name] = {"bins": bins, "intervals": intervals}
    return refs

def reg2bins(beg, end):
    lst = []
    if beg >= end: return lst
    end -= 1
    if beg >> 14 == end >> 14:
        lst.append(((1 << 15) - 1) // 7 + (beg >> 14)); beg >>= 14; end >>= 14
    if beg >> 17 == end >> 17:
        lst.append(((1 << 12) - 1) // 7 + (beg >> 17)); beg >>= 17; end >>= 17
    if beg >> 20 == end >> 20:
        lst.append(((1 << 9) - 1) // 7 + (beg >> 20)); beg >>= 20; end >>= 20
    if beg >> 23 == end >> 23:
        lst.append(((1 << 6) - 1) // 7 + (beg >> 23)); beg >>= 23; end >>= 23
    if beg >> 26 == end >> 26:
        lst.append(((1 << 3) - 1) // 7 + (beg >> 26))
    lst.append(0)
    return lst

def fetch_region(vcf_path, tbi_path, chrom, beg, end):
    """用 intervals（16kb 桶线性偏移）读取区域行——不依赖 bin 语义，稳妥正确"""
    refs = parse_tbi(tbi_path)
    ref = refs.get(chrom) or refs.get("chr" + chrom)
    if ref is None: return []
    iv = ref["intervals"]
    i1 = beg >> 14
    i2 = (end - 1) >> 14
    if i1 >= len(iv): return []
    # 起点：i1 起第一个非零 interval 偏移
    o1 = 0
    for j in range(i1, len(iv)):
        if iv[j] > 0:
            o1 = iv[j]; break
    if o1 == 0: return []
    pos = o1 >> 16
    skip = o1 & 0xFFFF
    out = bytearray()
    stop = False
    while not stop:
        raw = _file_range(vcf_path, pos, pos + 4 * 1024 * 1024)
        if not raw: break
        p = 0
        got = False
        while p + 18 <= len(raw):
            if raw[p:p+2] != b"\x1f\x8b": break
            bsize = int.from_bytes(raw[p+16:p+18], "little") + 1
            if p + bsize > len(raw): break
            try:
                dec = gzip.GzipFile(fileobj=io.BytesIO(raw[p:p+bsize])).read()
            except Exception:
                break
            if skip > 0:
                if skip >= len(dec): skip -= len(dec)
                else:
                    out += dec[skip:]; skip = 0
            else:
                out += dec
                # 提前停止：本块首行位置已越过区域终点
                nl = dec.find(b"\n")
                if nl > 0:
                    fl = dec[:nl]
                    if not fl.startswith(b"#"):
                        ff = fl.split(b"\t")
                        if len(ff) > 1:
                            try:
                                if int(ff[1]) > end:
                                    p += bsize
                                    stop = True
                                    break
                            except ValueError:
                                pass
            p += bsize; got = True
        if not got: break
        pos += p
    rows = []
    for line in out.split(b"\n"):
        if not line or line.startswith(b"#"): continue
        f = line.decode("utf-8", "replace").split("\t")
        if len(f) < 5: continue
        try: pp = int(f[1])
        except ValueError: continue
        if beg <= pp < end:
            rows.append(f)
    return rows

# ---------------- 数据读取 ----------------

def load_gwas_rsid():
    keep = []
    for ch in pd.read_csv(GWAS, sep="\t", compression="gzip",
                          usecols=["hm_rsid", "hm_beta", "standard_error", "p_value",
                                   "hm_effect_allele", "hm_other_allele"],
                          dtype={"hm_rsid": str}, chunksize=3_000_000):
        ch = ch.rename(columns={"hm_rsid": "rsid", "hm_beta": "beta",
                                "standard_error": "se", "p_value": "p",
                                "hm_effect_allele": "ea", "hm_other_allele": "oa"})
        for c in ["beta", "se", "p"]:
            ch[c] = pd.to_numeric(ch[c], errors="coerce")
        ch = ch.dropna(subset=["beta", "se"]).drop_duplicates("rsid")
        keep.append(ch[["rsid", "beta", "se", "p", "ea", "oa"]])
    g2 = pd.concat(keep, ignore_index=True).drop_duplicates("rsid", keep="first")
    return dict(zip(g2["rsid"], zip(g2["beta"], g2["se"], g2["p"], g2["ea"], g2["oa"])))


def allele_sign(ref, alt, ea, oa):
    """GWAS effect allele(ea) 与 eQTL alt 对齐：返回 β 符号 +1/-1/None"""
    if ea == alt: return 1
    if ea == ref: return -1
    if ea in COMP and COMP[ea] == alt: return 1
    if ea in COMP and COMP[ea] == ref: return -1
    return None


def ivw(b_mr, se_mr):
    """标准 Wald-ratio IVW：每 SNP 的比值 β_GW/β_eQTL 加权平均"""
    w = 1 / se_mr ** 2
    return (b_mr * w).sum() / w.sum(), np.sqrt(1 / w.sum())


def egger(x, y, se):
    w = 1 / se ** 2
    X = np.column_stack([np.ones(len(x)), x])
    W = np.diag(w)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
    resid = y - X @ beta
    sigma2 = (resid ** 2 * w).sum() / (len(x) - 2)
    cov = sigma2 * np.linalg.inv(X.T @ W @ X)
    return beta[1], np.sqrt(cov[1, 1]), beta[0], np.sqrt(cov[0, 0])


def weighted_median(b, w):
    w = w / w.sum()
    order = np.argsort(b)
    cw = np.cumsum(w[order])
    return b[order][np.searchsorted(cw, 0.5)]


def do_mr(cl, gene, qtd, n_region, n_match, tag):
    """对剪枝后 SNP 集合做 IVW（Wald-ratio）/ MR-Egger / 加权中位数"""
    b_ivw, se_ivw = ivw(cl["b_mr"].values, cl["se_mr"].values)
    p_ivw = 2 * (1 - stats.norm.cdf(abs(b_ivw) / se_ivw))
    try:
        slope, se_slope, intercept, se_intercept = egger(cl["beta"].values, cl["beta_gw"].values, cl["se_gw"].values)
        p_egger = 2 * (1 - stats.norm.cdf(abs(slope) / se_slope))
        p_int = 2 * (1 - stats.norm.cdf(abs(intercept) / se_intercept))
    except Exception:
        slope = se_slope = intercept = se_intercept = p_egger = p_int = np.nan
    wmed = weighted_median(cl["b_mr"].values, 1 / cl["se_mr"].values ** 2)
    return {"gene": gene, "qtd": qtd, "tag": tag, "n_snp_region": n_region, "n_matched": n_match,
            "n_snp": len(cl), "r2_thresh": R2_THRESH,
            "beta_ivw": b_ivw, "p_ivw": p_ivw, "or_ivw": np.exp(b_ivw),
            "beta_egger": slope, "p_egger": p_egger,
            "egger_intercept": intercept, "p_intercept": p_int,
            "wmedian_beta": wmed, "or_wmedian": np.exp(wmed)}


def ld_prune(positions, pvalues, chrom, gstart, gend):
    """1000G GRCh38 vcf 按【位置】匹配（vcf rsid 列为 .，用坐标匹配）算 r² → 贪心剪枝 r²<0.3。
    返回剪枝后保留的 pos 列表。"""
    vcf = os.path.join(LD_DIR, f"ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
    tbi = vcf + ".tbi"
    if not os.path.exists(vcf):
        print(f"    [LD] chr{chrom} vcf 缺失", flush=True)
        return None
    rows = fetch_region(vcf, tbi, chrom, gstart, gend)
    pos2row = {}
    for f in rows:
        if len(f) < 10: continue
        try: pp = int(f[1])
        except ValueError: continue
        if pp not in pos2row:
            pos2row[pp] = f
    # 剂量矩阵（ALT 剂量 0/1/2）——只保留 vcf 可匹配、双等位、低缺失 SNP
    keep_pos, keep_pvals, X = [], [], []
    for pp, pv in zip(positions, pvalues):
        f = pos2row.get(pp)
        if f is None: continue
        if "," in f[4]: continue  # 多等位跳过
        dos = []
        miss = 0
        for g in f[9:]:
            if g.startswith(".") or len(g) < 3:
                dos.append(np.nan); miss += 1
                continue
            try:
                dos.append(int(g[0]) + int(g[2]))
            except (ValueError, IndexError):
                dos.append(np.nan); miss += 1
        if miss / len(f[9:]) > 0.05: continue
        X.append(dos); keep_pos.append(pp); keep_pvals.append(pv)
    if len(keep_pos) < 2:
        return keep_pos
    X = np.array(X, dtype=float)
    # 缺失置列均值
    col_mean = np.nanmean(X, axis=1)
    for i in range(X.shape[0]):
        nan_mask = np.isnan(X[i])
        if nan_mask.any():
            X[i][nan_mask] = col_mean[i]
    # r² 矩阵
    Z = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True)
    R = Z @ Z.T / X.shape[1]
    r2 = R ** 2
    # 贪心剪枝（按 eQTL p 升序）——索引与 keep_pos 对齐
    n = len(keep_pos)
    order = np.argsort(keep_pvals)
    kept = []
    excluded = set()
    for i in order:
        if i in excluded: continue
        kept.append(i)
        for j in range(n):
            if j != i and j not in excluded and r2[i, j] > R2_THRESH:
                excluded.add(j)
    return [keep_pos[i] for i in kept]


def main():
    sig = pd.read_csv(SIG)
    top = sig.sort_values("F", ascending=False).drop_duplicates("gene_symbol")
    gq = dict(zip(top["gene_symbol"], top["qtd"]))
    print(f"[输入] {len(gq)} 基因", flush=True)
    gwas = load_gwas_rsid()
    print(f"[GWAS] rsid 条目 {len(gwas)}", flush=True)

    coords = {}
    with gzip.open(os.path.join(BASE, "data", "onek1k", "gencode.v38.annotation.gtf.gz"), "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            m = re.search(r'gene_name "([^"]+)"', line)
            gid = re.search(r'gene_id "([^"]+)"', line)
            if not (m and gid and m.group(1) in gq): continue
            p = line.strip().split("\t")
            if len(p) < 9 or p[2] != "gene": continue
            coords[m.group(1)] = (p[0].replace("chr", ""), int(p[3]), int(p[4]), gid.group(1).split(".")[0])
    print(f"[坐标] {len(coords)}", flush=True)

    by_qtd = defaultdict(list)
    for g, q in gq.items():
        if g in coords: by_qtd[q].append(g)

    results = []
    for qtd, genes in by_qtd.items():
        fp = os.path.join(ALL_DIR, f"{qtd}.all.tsv.gz")
        if not os.path.exists(fp):
            print(f"  QTD {qtd} 缺失", flush=True); continue
        targets = {coords[g][3]: g for g in genes}
        region = defaultdict(list)
        print(f"[读取] {qtd} ({len(genes)} 基因)...", flush=True)
        with gzip.open(fp, "rt", errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(hdr)}
            gi, ci, pi = idx["molecular_trait_id"], idx["chromosome"], idx["position"]
            ri, ai, vi = idx["ref"], idx["alt"], idx["variant"]
            bi, si, pvi, rsi = idx["beta"], idx["se"], idx["pvalue"], idx["rsid"]
            for line in f:
                c = line.split("\t")
                ensg = c[gi]
                if ensg not in targets: continue
                try:
                    ch, pos = c[ci], int(c[pi])
                except ValueError:
                    continue
                g = targets[ensg]
                chrom, gs, ge, _ = coords[g]
                if ch == chrom and gs - WINDOW <= pos <= ge + WINDOW:
                    region[g].append((pos, c[ri].strip(), c[ai].strip(),
                                      float(c[bi]), float(c[si]), float(c[pvi]), c[rsi].strip()))
        for gene in genes:
            # QTD 区域缓存（重跑免全文件扫描）
            cache_fp = os.path.join(CACHE, f"{qtd}_{gene}.pkl")
            if os.path.exists(cache_fp):
                with open(cache_fp, "rb") as f:
                    df = pickle.load(f)
            else:
                if gene not in region:
                    print(f"  {gene}: 无区域 SNP", flush=True); continue
                df = pd.DataFrame(region[gene], columns=["pos", "ref", "alt", "beta", "se", "pvalue", "rsid"])
                os.makedirs(CACHE, exist_ok=True)
                with open(cache_fp, "wb") as f:
                    pickle.dump(df, f)
            df = df[df["pvalue"] < 0.05].copy()
            df["rsid2"] = df["rsid"].where(df["rsid"] != ".", np.nan)
            m = df.dropna(subset=["rsid2"]).copy()
            m["gwas"] = m["rsid2"].map(gwas)
            m = m.dropna(subset=["gwas"]).copy()
            if len(m) < 2:
                print(f"  {gene}: rsid 匹配 <2", flush=True); continue
            # 方向校正
            signs = []
            for _, r in m.iterrows():
                ea, oa = r["gwas"][3], r["gwas"][4]
                signs.append(allele_sign(r["ref"], r["alt"], str(ea), str(oa)))
            m["sign"] = signs
            m = m[m["sign"].notna()].copy()
            if len(m) < 2:
                print(f"  {gene}: 方向可对齐 <2", flush=True); continue
            m["beta_gw"] = [x[0] * s for x, s in zip(m["gwas"], m["sign"])]
            m["se_gw"] = [x[1] for x in m["gwas"]]
            m = m.dropna(subset=["beta_gw", "se_gw"])
            chrom, gs, ge, _ = coords[gene]
            # LD 剪枝（按位置匹配 GRCh38 1000G）
            kept = ld_prune(list(m["pos"]), list(m["pvalue"]), chrom, gs - WINDOW, ge + WINDOW)
            if kept is None:
                print(f"  {gene}: LD 参考不可用", flush=True); continue
            cl = m[m["pos"].isin(kept)].copy()
            if len(cl) < 2:
                print(f"  {gene}: 剪枝后 <2 SNP", flush=True); continue
            # Wald ratio 与 F 统计量（官方 SE）
            cl["b_mr"] = cl["beta_gw"] / cl["beta"]
            cl["se_mr"] = cl["se_gw"] / cl["beta"].abs()
            cl["F"] = (cl["beta"] / cl["se"]) ** 2
            cl = cl.dropna(subset=["b_mr", "se_mr"])
            # 版本1：剪枝后全部显著 SNP；版本2：F>10 严格
            out_row = do_mr(cl, gene, qtd, len(df), len(m), "all")
            cl_f10 = cl[cl["F"] > 10].copy()
            if len(cl_f10) >= 2:
                out_row2 = do_mr(cl_f10, gene, qtd, len(df), len(m), "f10")
            else:
                out_row2 = {"gene": gene, "qtd": qtd, "tag": "f10", "n_snp": 0,
                            "beta_ivw": np.nan, "p_ivw": np.nan, "or_ivw": np.nan,
                            "beta_egger": np.nan, "p_egger": np.nan,
                            "egger_intercept": np.nan, "p_intercept": np.nan,
                            "wmedian_beta": np.nan, "or_wmedian": np.nan}
            results.append(out_row); results.append(out_row2)
            # 注意：out_row["or_ivw"] 已经是对数尺度 β 取 exp 后的值（见 do_mr 的
            # "or_ivw": np.exp(b_ivw)），此处**直接引用**，不得再套一层 np.exp()——
            # 早期版本误写成 np.exp(out_row['or_ivw'])，导致日志中的 OR 被双重指数化
            # 而虚高（如 JAML 真实 OR=1.11 被打印成 3.05）。CSV 输出一直是正确的。
            _or_all = out_row["or_ivw"]
            _or_f10 = out_row2["or_ivw"] if out_row2["n_snp"] > 0 else float("nan")
            print(f"  {gene}: ALL{len(cl)}SNP IVW OR={_or_all:.3f} P={out_row['p_ivw']:.2g} "
                  f"| F>10:{len(cl_f10)}SNP OR={_or_f10:.3f} "
                  f"P={out_row2['p_ivw']:.2g} 截距P={out_row['p_intercept']:.2g}", flush=True)
    out = pd.DataFrame(results)
    out.to_csv(os.path.join(RES, "multisnp_mr_results_ld.csv"), index=False)
    print(f"\n[完成] {len(out)} 基因; results/multisnp_mr_results_ld.csv", flush=True)


if __name__ == "__main__":
    main()
