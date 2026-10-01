# -*- coding: utf-8 -*-
"""验证本机自建的 velocyto 精简包（M4 CellOracle 的前置条件）。

背景
----
CellOracle 0.20.0 依赖 velocyto 的三个子模块：
  diffusion.Diffusion / estimation.colDeltaCor* / serialization.dump_hdf5,load_hdf5
velocyto 在 PyPI 上**没有任何平台的 wheel**（全部只有 sdist），而它的 sdist 只附带
Cython 0.29.2 生成的 C 代码，在 CPython 3.13 + numpy 2.x 下无法编译：
  ① `Py_SIZE(list) = ...` —— Py_SIZE 自 3.11 起不再是左值
  ② `_PyList_Extend` —— 私有 API 无声明
  ③ `PyArray_Descr.subarray` —— numpy 2.0 起该成员已 opaque
处理：用 Cython 3.3.0 从 speedboosted.pyx 重新生成 C，再用 Rtools44 的
x86_64-w64-mingw32 GCC 13.3.0 编译，链接 python313.lib。**未手改任何算法代码**。

本脚本做三件事：
  1. 包完整性检查
  2. **数值正确性**：speedboosted 的 colDeltaCor 与 numpy 参考实现逐元素比对
     （这才是关键——编译通过不代表算得对）
  3. 串行 / 并行（OpenMP）结果一致性

用法：python scripts/40_verify_velocyto_minimal.py
"""
import os
import sys
import sysconfig

import numpy as np

SP = sysconfig.get_paths()["purelib"]          # .../Lib/site-packages
PKG = os.path.join(SP, "velocyto")
EXPECTED = ["__init__.py", "diffusion.py", "estimation.py", "serialization.py",
            "speedboosted.cp313-win_amd64.pyd"]

fails = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  {detail}" if detail else ""))
    if not ok:
        fails.append(label)


print("=" * 74)
print("1. 包完整性")
print("=" * 74)
check("velocyto 目录存在", os.path.isdir(PKG), PKG)
if os.path.isdir(PKG):
    have = set(os.listdir(PKG))
    for f in EXPECTED:
        check(f"包含 {f}", f in have)

print()
print("=" * 74)
print("2. 导入（含 6 个 speedboosted 包装函数）")
print("=" * 74)
try:
    import velocyto
    print(f"  velocyto.__version__ = {velocyto.__version__}")
    check("velocyto 可导入", True)
except Exception as e:
    check("velocyto 可导入", False, f"{type(e).__name__}: {e}")
    sys.exit(1)

from velocyto.diffusion import Diffusion                                    # noqa: E402
from velocyto.serialization import dump_hdf5, load_hdf5                     # noqa: E402
from velocyto.estimation import (colDeltaCor, colDeltaCorLog10,             # noqa: E402
                                 colDeltaCorLog10partial, colDeltaCorpartial,
                                 colDeltaCorSqrt, colDeltaCorSqrtpartial)
from velocyto import speedboosted as sb                                     # noqa: E402

check("Diffusion 类", isinstance(Diffusion, type))
check("dump_hdf5/load_hdf5 可调用", callable(dump_hdf5) and callable(load_hdf5))
for fn in [colDeltaCor, colDeltaCorLog10, colDeltaCorLog10partial,
           colDeltaCorpartial, colDeltaCorSqrt, colDeltaCorSqrtpartial]:
    check(f"estimation.{fn.__name__}", callable(fn))

print()
print("=" * 74)
print("3. 数值正确性：colDeltaCor  vs  numpy 参考实现")
print("=" * 74)


def ref_colDeltaCor(e, d):
    """按 speedboosted.pyx 的算法逐字转写（列 c 与所有列的相关）。

    rm[c, i] = corr( e[:, i] - e[:, c], d[:, c] )
    """
    _, cols = e.shape
    rm = np.zeros((cols, cols))
    for c in range(cols):
        A = e - e[:, [c]]                       # A[:, i] = e[:, i] - e[:, c]
        Am = A - A.mean(axis=0, keepdims=True)
        b = d[:, c]
        bm = b - b.mean()
        num = Am.T @ bm
        den = np.linalg.norm(Am, axis=0) * np.linalg.norm(bm)
        rm[c] = num / den
    return rm


rng = np.random.default_rng(20260930)
rows, cols = 5, 60                              # 5 基因 × 60 细胞（按 pyx：rows=基因, cols=细胞）
e = rng.normal(size=(rows, cols)) * 2 + 1
d = rng.normal(size=(rows, cols))

ref = ref_colDeltaCor(e, d)
got = np.asarray(colDeltaCor(e, d, threads=1))
print(f"  形状: speedboosted {got.shape}  vs  reference {ref.shape}")
check("输出形状一致", got.shape == ref.shape)

if got.shape == ref.shape:
    gn, rn = np.isnan(got), np.isnan(ref)
    print(f"  NaN 出现：speedboosted {int(gn.sum())} 处，reference {int(rn.sum())} 处；"
          f"不一致 {int((gn ^ rn).sum())} 处")
    # i == c 时位移向量恒为 0，相关系数未定义 → 双方都应是 NaN。这是算法固有行为，
    # 关键是"出现位置一致"，而非"没有 NaN"。
    check("NaN 出现位置一致（i==c 退化为 0/0）", int((gn ^ rn).sum()) == 0,
          f"不一致 {int((gn ^ rn).sum())} 处")
    ok = ~(gn | rn)
    diff = np.abs(got[ok] - ref[ok])
    print(f"  非 NaN 元素 n = {int(ok.sum())}：最大绝对差 = {diff.max():.3e}，"
          f"平均 = {diff.mean():.3e}")
    check("逐元素一致 (< 1e-10)", diff.max() < 1e-10, f"max|Δ| = {diff.max():.2e}")
    vmin, vmax = float(np.nanmin(got)), float(np.nanmax(got))
    check("取值落在 [-1, 1]", vmin >= -1 - 1e-9 and vmax <= 1 + 1e-9,
          f"范围 [{vmin:.6f}, {vmax:.6f}]")

print()
print("=" * 74)
print("4. OpenMP 并行一致性（threads 1 vs 8）")
print("=" * 74)
a1 = np.asarray(colDeltaCor(e, d, threads=1))
a8 = np.asarray(colDeltaCor(e, d, threads=8))
okp = ~(np.isnan(a1) | np.isnan(a8))
check("并行/串行 NaN 模式一致", int((np.isnan(a1) ^ np.isnan(a8)).sum()) == 0)
d18 = float(np.abs(a1[okp] - a8[okp]).max())
print(f"  串行 vs 8 线程（n = {int(okp.sum())}）最大差 = {d18:.3e}")
check("并行结果与串行一致 (< 1e-12)", d18 < 1e-12, f"max|Δ| = {d18:.2e}")

print()
print("=" * 74)
print("5. colDeltaCorLog10 自洽性（与手算 log10 位移比对）")
print("=" * 74)
# psc 取 0.5（而非 1.0）：这样 tmp==0 时 pyx 的 else 分支给 -log10(0.5)=+0.301，
# 与 np.sign 语义（给 0）明显可分 —— 能真正检验分支是否被照抄。
psc = 0.5
e2 = np.abs(rng.normal(size=(4, 40))) + 0.5      # 正值，便于 log10
d2 = rng.normal(size=(4, 40))
got10 = np.asarray(colDeltaCorLog10(e2, d2, threads=1, psc=psc))


def ref_colDeltaCorLog10(e, d, psc):
    """按 speedboosted.pyx 的 x_colDeltaCorLog10 逐字转写。

    pyx 逻辑:  if tmp > 0: A = log10(tmp + psc)
               else:       A = -log10(-tmp + psc)
    tmp == 0 走 else 分支（得 -log10(psc)）—— 用 np.where 精确复现分支语义，
    而不是用 np.sign（后者在 0 处给 0，会掩盖差异）。
    """
    _, cols = e.shape
    rm = np.zeros((cols, cols))
    for c in range(cols):
        tmp = e - e[:, [c]]
        A = np.where(tmp > 0,
                     np.log10(np.abs(tmp) + psc),
                     -np.log10(np.abs(tmp) + psc))
        Am = A - A.mean(axis=0, keepdims=True)
        b = d[:, c]
        bm = b - b.mean()
        rm[c] = (Am.T @ bm) / (np.linalg.norm(Am, axis=0) * np.linalg.norm(bm))
    return rm


ref10 = ref_colDeltaCorLog10(e2, d2, psc)
gn, rn = np.isnan(got10), np.isnan(ref10)
diag = np.eye(got10.shape[0], dtype=bool)      # i == c：位移恒为 0，相关系数无定义
off = ~diag
n_bad_off = int(((gn ^ rn) & off).sum())
check("非对角线：NaN 位置一致", n_bad_off == 0, f"不一致 {n_bad_off} 处")
ok10 = off & ~(gn | rn)
d10 = float(np.abs(got10[ok10] - ref10[ok10]).max())
print(f"  psc = {psc}；非对角线有效元素 n = {int(ok10.sum())}；最大绝对差 = {d10:.3e}")
check("log10 版（非对角线）逐元素一致 (< 1e-10)", d10 < 1e-10, f"max|Δ| = {d10:.2e}")
n_diag_diff = int((gn ^ rn).sum()) - n_bad_off
print(f"  [披露] 对角线（i == c，位移向量恒为 0、相关系数无定义）：共 {int(diag.sum())} 个位置，"
      f"两实现不一致 {n_diag_diff} 个。\n"
      f"         参考实现给 NaN；speedboosted 在 -ffast-math 下把 1/sqrt(数值噪声) 放大为有限值。\n"
      f"         对角线无定义，任何实现都不得作为结果使用（CellOracle 亦不取对角线）。")
check("参考实现的对角线为 NaN（确认其为退化情形而非计算差异）",
      bool(np.isnan(ref10[diag]).all()),
      f"对角线 NaN {int(np.isnan(ref10[diag]).sum())}/{int(diag.sum())}")

print()
print("=" * 74)
print("6. CellOracle 可导入性")
print("=" * 74)
try:
    import celloracle as co
    check("import celloracle", True, f"version {co.__version__}")
    for attr in ["Oracle", "data", "trajectory"]:
        check(f"celloracle.{attr} 存在", hasattr(co, attr))
except Exception as ex:
    import traceback
    traceback.print_exc()
    check("import celloracle", False, f"{type(ex).__name__}: {str(ex)[:160]}")

print()
print("=" * 74)
if fails:
    print(f"结果：{len(fails)} 项未通过")
    for f in fails:
        print("   -", f)
    sys.exit(1)
print("结果：全部通过 —— 自建 velocyto 精简包可替代官方包（数值已验证）")
print("=" * 74)
