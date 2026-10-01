# -*- coding: utf-8 -*-
"""纯 Python tabix(.tbi) 区域查询：通过 HTTP range 从远程 BGZF 压缩的 TSV 提取区域行。
无需 pysam，适用于 eQTL Catalogue / GTEx sumstats 等标准 tabix 文件。
用法:
    from tools_tabix import TabixQuery
    tq = TabixQuery(base_url, tbi_path)   # tbi_path 为本地已下载的 .tbi
    rows = tq.fetch('chr1', 1000, 200000, chrom_col=13, pos_col=14, header=1)
"""
import gzip
import io
import struct
import urllib.request

__all__ = ["TabixQuery", "fetch_region_virtual", "parse_tbi"]


def _http_range(url: str, start: int, end: int, timeout=180) -> bytes:
    req = urllib.request.Request(
        url, headers={"User-Agent": "Mozilla/5.0", "Range": "bytes=%d-%d" % (start, end - 1)})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _bgzf_decompress(raw: bytes) -> bytes:
    return gzip.GzipFile(fileobj=io.BytesIO(raw)).read()


def parse_tbi(tbi_path: str) -> dict:
    with open(tbi_path, "rb") as f:
        raw = f.read()
    dec = _bgzf_decompress(raw)
    off = 0

    def u32():
        nonlocal off
        v = struct.unpack("<I", dec[off:off + 4])[0]; off += 4
        return v

    def i32():
        nonlocal off
        v = struct.unpack("<i", dec[off:off + 4])[0]; off += 4
        return v

    def u64():
        nonlocal off
        v = struct.unpack("<Q", dec[off:off + 8])[0]; off += 8
        return v

    magic = dec[off:off + 4]; off += 4
    if magic != b"TBI\x01":
        raise ValueError("bad tabix magic: %r" % magic)
    n_ref = i32()
    _fmt = i32()
    col_seq = i32()
    col_beg = i32()
    _col_end = i32()
    _meta = i32()
    _skip = i32()
    l_nm = i32()
    names = dec[off:off + l_nm].decode("utf-8", "replace").split("\x00"); off += l_nm
    refs = {}
    for _ in range(n_ref):
        name = names[_] if _ < len(names) else str(_ + 1)
        bins = {}
        n_bin = i32()
        for _b in range(n_bin):
            bin_id = u32()
            n_chunk = i32()
            chunks = [(u64(), u64()) for _c in range(n_chunk)]
            bins[bin_id] = chunks
        n_intv = i32()
        intervals = [u64() for _i in range(n_intv)]
        refs[name] = {"bins": bins, "intervals": intervals}
    return {"col_seq": col_seq, "col_beg": col_beg, "refs": refs}


def reg2bins(beg: int, end: int):
    lst = []
    if beg >= end:
        return lst
    end -= 1
    if beg >> 14 == end >> 14:
        lst.append(((1 << 15) - 1) // 7 + (beg >> 14))
        beg >>= 14; end >>= 14
    if beg >> 17 == end >> 17:
        lst.append(((1 << 12) - 1) // 7 + (beg >> 17))
        beg >>= 17; end >>= 17
    if beg >> 20 == end >> 20:
        lst.append(((1 << 9) - 1) // 7 + (beg >> 20))
        beg >>= 20; end >>= 20
    if beg >> 23 == end >> 23:
        lst.append(((1 << 6) - 1) // 7 + (beg >> 23))
        beg >>= 23; end >>= 23
    if beg >> 26 == end >> 26:
        lst.append(((1 << 3) - 1) // 7 + (beg >> 26))
        beg >>= 26; end >>= 26
    lst.append(0)
    return lst


def merge_chunks(chunks):
    if not chunks:
        return []
    chunks = sorted(chunks)
    out = [list(chunks[0])]
    for b, e in chunks[1:]:
        if b <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([b, e])
    return [tuple(c) for c in out]


def fetch_region_virtual(url: str, vb: int, ve: int, skip_bytes: int = 0,
                         max_block_bytes: int = 16 * 1024 * 1024, timeout=120) -> bytes:
    """按虚拟偏移区间 [vb, ve) 提取解压数据。
    - 虚拟偏移: 高48位=BGZF块起始文件偏移, 低16位=块内解压偏移
    - skip_bytes: 跳过起始块内前 N 个解压字节（来自 vb 低16位）
    实现: 从 vb>>16 的块边界开始，按 BGZF 块(BC 字段 BSIZE)逐步 range 请求并解压。
    """
    pos = vb >> 16
    end_file = (ve >> 16) + 1  # 多取一块，确保覆盖
    skip = skip_bytes
    out = bytearray()
    while pos < end_file:
        req_end = min(pos + max_block_bytes, end_file + max_block_bytes)
        try:
            raw = _http_range(url, pos, req_end, timeout=timeout)
        except Exception:
            break
        p = 0
        got_block = False
        while p + 18 <= len(raw):
            if raw[p:p + 2] != b"\x1f\x8b":
                break
            xlen = int.from_bytes(raw[p + 10:p + 12], "little")
            bsize = int.from_bytes(raw[p + 16:p + 18], "little") + 1
            if p + bsize > len(raw):
                break  # 块不完整，需要更大请求
            block = raw[p:p + bsize]
            try:
                dec = _bgzf_decompress(block)
            except Exception:
                break
            if skip > 0:
                if skip >= len(dec):
                    skip -= len(dec)
                else:
                    out += dec[skip:]
                    skip = 0
            else:
                out += dec
            p += bsize
            got_block = True
        if not got_block:
            # 请求内无完整块：扩大范围重试一次，否则终止
            req_end2 = min(end_file + 3 * max_block_bytes, pos + 4 * max_block_bytes)
            if req_end2 <= req_end:
                break
            try:
                raw = _http_range(url, pos, req_end2, timeout=timeout)
            except Exception:
                break
            p = 0
            while p + 18 <= len(raw):
                if raw[p:p + 2] != b"\x1f\x8b":
                    break
                xlen = int.from_bytes(raw[p + 10:p + 12], "little")
                bsize = int.from_bytes(raw[p + 16:p + 18], "little") + 1
                if p + bsize > len(raw):
                    break
                block = raw[p:p + bsize]
                try:
                    dec = _bgzf_decompress(block)
                except Exception:
                    break
                if skip > 0:
                    if skip >= len(dec):
                        skip -= len(dec)
                    else:
                        out += dec[skip:]
                        skip = 0
                else:
                    out += dec
                p += bsize
            if p == 0:
                break
        pos += p
    return bytes(out)


class TabixQuery:
    def __init__(self, data_url: str, tbi_path: str):
        self.url = data_url
        self.meta = parse_tbi(tbi_path)

    def _contig(self, chrom: str) -> dict:
        refs = self.meta["refs"]
        if chrom in refs:
            return refs[chrom]
        if chrom.startswith("chr"):
            return refs.get(chrom[3:])
        return refs.get("chr" + chrom)

    def fetch(self, chrom: str, beg: int, end: int, chrom_col: int, pos_col: int,
              header_rows: int = 1, timeout=180):
        ref = self._contig(chrom)
        if ref is None:
            return []
        bins = reg2bins(beg, end)
        chunks = []
        for b in bins:
            chunks.extend(ref["bins"].get(b, []))
        if not chunks:
            iv = ref["intervals"]
            if iv:
                i1 = min(beg >> 14, len(iv) - 1)
                i2 = min((end - 1) >> 14, len(iv) - 1)
                o1 = iv[i1]
                o2 = iv[i2]
                if o2 == 0:
                    j = i2
                    while j >= 0 and iv[j] == 0:
                        j -= 1
                    o2 = iv[j] if j >= 0 else o1
                if o1 and o2 and o2 >= o1:
                    chunks = [(o1, o2)]
        merged = merge_chunks(chunks)
        rows = []
        for s, e in merged:
            try:
                text = fetch_region_virtual(self.url, s, e, skip_bytes=s & 0xFFFF, timeout=timeout)
            except Exception:
                continue
            for line in text.split(b"\n"):
                if not line:
                    continue
                if line.startswith(b"#"):
                    continue
                fields = line.decode("utf-8", "replace").split("\t")
                if len(fields) < max(chrom_col, pos_col):
                    continue
                c = fields[chrom_col - 1]
                if c != chrom and c != chrom.replace("chr", ""):
                    continue
                try:
                    p = int(fields[pos_col - 1])
                except ValueError:
                    continue
                if beg <= p < end:
                    rows.append(fields)
        return rows
