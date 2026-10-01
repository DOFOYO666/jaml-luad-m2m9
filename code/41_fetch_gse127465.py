# -*- coding: utf-8 -*-
"""稳健续传下载 GSE127465 表达矩阵（含 Range 语义校验）。

上一版失败的原因（已定位）：
  用 `Range: bytes=N-` 续传时**没有检查响应码**。当服务器忽略 Range、返回
  完整文件（200）而我方按分片追加（ab），就会把整份数据再写一遍，
  造成文件超长且中间内容错位（实测末尾多出 4,678,900 B，gzip 报
  `invalid literal/length code`）。
本版强制：
  1) 只有 **206 Partial Content** 才按分片追加；返回 200 则从 0 重写。
  2) 单轮请求的超时/中断只影响该轮，不影响已有字节。
  3) 收尾做 **完整 gzip 解压校验**，并核对 MM 头里的 NNZ 与数据行数是否一致。

用法：python 41_fetch_gse127465.py [--fresh]
      --fresh 表示忽略本地文件、从头下载
"""
import gzip
import os
import sys
import time
import urllib.request

URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/"
       "GSE127465_human_counts_normalized_54773x41861.mtx.gz")
DST = os.path.join("rawdata", "scRNA", "GSE127465",
                   "GSE127465_human_counts_normalized_54773x41861.mtx.gz")
FRESH = "--fresh" in sys.argv
MAX_ROUNDS = 500
STALL_LIMIT = 30

os.makedirs(os.path.dirname(DST), exist_ok=True)
if FRESH and os.path.isfile(DST):
    bak = DST + ".bad"
    i = 0
    while os.path.exists(bak):
        i += 1
        bak = DST + f".bad{i}"
    os.replace(DST, bak)
    print(f"[重置] 旧文件移至 {os.path.basename(bak)}", flush=True)


def size():
    return os.path.getsize(DST) if os.path.isfile(DST) else 0


def head_total():
    try:
        req = urllib.request.Request(URL, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            return int(r.headers.get("Content-Length") or 0)
    except Exception as e:
        print(f"[HEAD 失败] {type(e).__name__}: {e}", flush=True)
        return 0


TOTAL = head_total()
print(f"[目标] {URL}", flush=True)
print(f"[大小] 服务器 {TOTAL} B ({TOTAL/1048576:.2f} MiB)", flush=True)
print(f"[起始] 本地 {size()/1048576:.2f} MiB", flush=True)

stall = 0
prev = size()
n_rewrite = 0

for rnd in range(1, MAX_ROUNDS + 1):
    have = size()
    if TOTAL and have >= TOTAL:
        break

    req = urllib.request.Request(URL)
    if have:
        req.add_header("Range", f"bytes={have}-")
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            code = resp.status
            if have and code != 206:
                # 服务器忽略了 Range，返回完整文件 → 必须从 0 重写，否则内容错位
                print(f"[轮 {rnd}] 服务器返回 {code}（非 206）→ 从头重写", flush=True)
                mode, have = "wb", 0
                n_rewrite += 1
            else:
                mode = "ab" if have else "wb"
            with open(DST, mode) as f:
                while True:
                    chunk = resp.read(1 << 16)
                    if not chunk:
                        break
                    f.write(chunk)
    except Exception as e:
        print(f"[轮 {rnd}] 中断: {type(e).__name__}: {str(e)[:70]}", flush=True)
        time.sleep(2)

    now = size()
    got = now - prev
    print(f"[轮 {rnd}] {now/1048576:.2f} MiB (+{got/1024:.0f} KB)", flush=True)
    if got <= 0:
        stall += 1
        if stall >= STALL_LIMIT:
            print(f"[放弃] 连续 {stall} 轮无增长", flush=True)
            break
    else:
        stall = 0
    prev = now

final = size()
print(f"\n[结束] {final/1048576:.2f} MiB（重写次数 {n_rewrite}）", flush=True)

# ---------------- 完整性校验 ----------------
ok = False
try:
    n = 0
    with gzip.open(DST, "rt", errors="strict") as f:
        l1 = f.readline().rstrip()
        f.readline()
        l3 = f.readline().rstrip()
        nnz = int(l3.split()[2])
        for _ in f:
            n += 1
    print(f"[校验] MM 头声明 NNZ = {nnz:,}；实际数据行 = {n:,}", flush=True)
    ok = (n == nnz)
    print(f"[校验] gzip 可完整解压，行数{'一致 ✓' if ok else '不一致 ✗'}", flush=True)
except Exception as e:
    print(f"[校验] 失败: {type(e).__name__}: {str(e)[:100]}", flush=True)

print("VERDICT:", "OK" if ok else "FAIL", flush=True)
