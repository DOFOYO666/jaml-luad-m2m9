# -*- coding: utf-8 -*-
"""
01_inspect_data.py - 检查已下载数据的格式（列名、行数、示例）
在运行主分析前执行，确认各文件的实际列结构。
"""
import gzip, os, sys, io

DATA = os.path.join(os.path.dirname(__file__), "..", "data")

def peek(path, n=3, is_gzip=None):
    if not os.path.exists(path):
        print(f"[MISSING] {path}")
        return
    size_mb = os.path.getsize(path) / 1024 / 1024
    print(f"[FILE] {os.path.basename(path)}  {size_mb:.1f} MB")
    try:
        opener = gzip.open if (is_gzip or path.endswith(".gz")) else open
        with opener(path, "rt", errors="replace") as f:
            for i in range(n):
                line = f.readline()
                if not line:
                    break
                print(f"  line{i}: {line.rstrip()[:400]}")
    except Exception as e:
        print(f"  ERROR reading: {e}")

def count_rows(path, is_gzip=True):
    opener = gzip.open if is_gzip else open
    n = 0
    with opener(path, "rt", errors="replace") as f:
        for _ in f:
            n += 1
    return n

def main():
    base = os.path.join(DATA, "onek1k")
    zip_path = os.path.join(base, "OneK1K_TensorQTL_top_eQTL_summary.zip")
    if os.path.exists(zip_path):
        import zipfile
        try:
            z = zipfile.ZipFile(zip_path)
            print("=== OneK1K zip contents ===")
            for i in z.infolist():
                print(f"  {i.filename}  {i.file_size/1024:.0f} KB")
        except Exception as e:
            print("zip error:", e)
    else:
        # 如果已解压，查看解压目录
        for fn in sorted(os.listdir(base)):
            fp = os.path.join(base, fn)
            if os.path.isfile(fp) and not fn.endswith(".zip"):
                print(f"=== {fn} ===")
                peek(fp)

    gwas = os.path.join(DATA, "luad_gwas", "GCST004744.h.tsv.gz")
    print("\n=== LUAD GWAS (harmonised) head ===")
    peek(gwas)

    tcga_expr = os.path.join(DATA, "luad_survival", "TCGA.LUAD.HiSeqV2.gz")
    tcga_clin = os.path.join(DATA, "luad_survival", "TCGA.LUAD.clinicalMatrix.txt")
    print("\n=== TCGA-LUAD expression head ===")
    peek(tcga_expr)
    print("\n=== TCGA-LUAD clinical head ===")
    peek(tcga_clin, n=2)

    eqtlgen = os.path.join(DATA, "eqtlgen", "cis-eQTLs_full_20180905.txt.gz")
    print("\n=== eQTLGen head ===")
    peek(eqtlgen)

if __name__ == "__main__":
    main()
