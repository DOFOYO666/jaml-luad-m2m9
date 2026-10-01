# -*- coding: utf-8 -*-
"""把《JAML机制课题研究》的三份 Markdown 文档转换为 docx（便于申报填写与打印）。

复用 scripts/49_build_new_docx.py 中的转换器（三线表、字体约定一致）。

用法：python 52_build_ke_docx.py
"""
import importlib.util
import os
import sys

ROOT = r"C:\Users\86159\WorkBuddy\单细胞数据孟德尔随机化"
KE = os.path.join(ROOT, "JAML机制课题研究")


def load_converter():
    p = os.path.join(ROOT, "scripts", "49_build_new_docx.py")
    spec = importlib.util.spec_from_file_location("md2docx", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["md2docx"] = mod
    spec.loader.exec_module(mod)
    return mod


DOCS = [
    ("01_课题立项书_JAML机制研究.md", "01_课题立项书_JAML机制研究.docx", "cn"),
    ("02_预实验方案.md", "02_预实验方案.docx", "cn"),
    ("03_前期工作基础.md", "03_前期工作基础.docx", "cn"),
]

if __name__ == "__main__":
    m = load_converter()
    for src_name, out_name, lang in DOCS:
        cfg = dict(m.CFG[lang])
        src = os.path.join(KE, src_name)
        dst = os.path.join(KE, out_name)
        m.convert(src, dst, cfg)
        print(f"[docx] {out_name}  ({os.path.getsize(dst)/1024:.0f} KB)")
    print("DONE")
