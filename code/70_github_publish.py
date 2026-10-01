# -*- coding: utf-8 -*-
"""发布助手：把 GitHub 建仓推送、以及拿到 DOI 后的"写回 + 重建"两个机械环节各变成一条命令。

用法（三条命令对应手册的 §3 / §6，其余步骤只能在网页上点）

  1) 先自检，看看缺什么（不改任何东西）：
       python scripts/70_github_publish.py --check

  1b) 校验归档没损坏（不用 sha256sum，本机没有这个命令）：
       python scripts/70_github_publish.py --verify-tar

  2) 网页建好空仓库之后，一条命令完成"配身份 + 加 remote + 改名 main + 推送"：
       python scripts/70_github_publish.py --url https://github.com/<你的用户名>/jaml-luad-m2m9
     （推送时若弹出 Git Credential Manager 窗口，用户名填你的 GitHub 用户名，密码处粘贴 PAT）

  3) 从 Zenodo 拿到 DOI 之后，一条命令完成"写回稿件 + 重建全部交付物 + 再推送一次"：
       python scripts/70_github_publish.py --url https://github.com/<你的用户名>/jaml-luad-m2m9 \
            --doi 10.5281/zenodo.1234567 --writeback

本脚本**不联网建仓库、不申请 DOI**（那需要你的账号），只把本机能做的那部分做完。
"""
import hashlib
import os
import subprocess
import sys

sys.stdout.reconfigure(encoding="utf-8")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.join(BASE, "JAML_M2M9_code_release")
PY = sys.executable

GIT_NAME = "Huayong Liu"
GIT_EMAIL = "xingxinghuoshu@163.com"

OK, WARN, BAD = "  [OK]  ", "  [注意] ", "  [缺]  "


def log(*a):
    print(*a, flush=True)


def run(args, cwd=None, capture=True, check=False):
    r = subprocess.run(args, cwd=cwd, capture_output=capture, text=True)
    if check and r.returncode != 0:
        raise SystemExit("命令失败：%s\n%s" % (" ".join(args), (r.stderr or r.stdout or "").strip()[:400]))
    return r


def git_ok():
    try:
        run(["git", "--version"], check=True)
        return True
    except Exception:
        return False


def current_branch():
    r = run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=REPO)
    return r.stdout.strip() if r.returncode == 0 else "?"


# --------------------------------------------------------------------- check --
def do_check():
    log("=" * 78)
    log("发布前自检（只看不改）")
    log("=" * 78)

    log("\n[1] 本机工具")
    if git_ok():
        log(OK + "git 已安装：%s" % run(["git", "--version"]).stdout.strip())
    else:
        log(BAD + "git 未安装 —— 先装 Git for Windows：https://git-scm.com/download/win")

    log("\n[2] git 全局身份（提交时会用）")
    name = run(["git", "config", "--global", "user.name"]).stdout.strip()
    mail = run(["git", "config", "--global", "user.email"]).stdout.strip()
    log((OK if name else WARN) + "user.name  = %s" % (name or "未配置（步骤 2 会自动配）"))
    log((OK if mail else WARN) + "user.email = %s" % (mail or "未配置（步骤 2 会自动配）"))

    log("\n[3] 本地发布树")
    if not os.path.isdir(REPO):
        log(BAD + "找不到 %s —— 先跑 python scripts/66_release_code.py" % REPO)
        return
    log(OK + "目录存在：%s" % REPO)
    if os.path.isdir(os.path.join(REPO, ".git")):
        log(OK + "已是 git 仓库，当前分支 %s" % current_branch())
    else:
        log(BAD + "还没有 .git —— 先跑 python scripts/66_release_code.py")
        return
    head = run(["git", "log", "--oneline", "-1"], cwd=REPO).stdout.strip()
    tracked = run(["git", "ls-files"], cwd=REPO).stdout.splitlines()
    log(OK + "最新提交：%s" % head)
    log(OK + "跟踪文件：%d 个" % len(tracked))
    bad = [l for l in tracked if l.startswith("_superseded/")]
    log((OK if not bad else BAD) + "回收区未上传：%s" % ("是（0 个）" if not bad else "%d 个会被推上去！" % len(bad)))
    ls = os.path.join(BASE, "JAML_M2M9_code_release.tar.gz.sha256")
    log((OK if os.path.exists(ls) else WARN) + "归档校验文件：%s"
        % (open(ls, encoding="utf-8").read().strip()[:24] + "…" if os.path.exists(ls) else "缺（跑 66 生成）"))

    log("\n[4] 远程仓库")
    rem = run(["git", "remote", "-v"], cwd=REPO).stdout.strip()
    if rem:
        log(OK + "已配置 remote：")
        for line in rem.splitlines():
            log("         " + line)
    else:
        log(WARN + "尚未配置 remote —— 网页建好空仓库后跑步骤 2 的命令即可")

    log("\n[5] 网页操作是否已完成的判断依据")
    log("      · 仓库页能看到 %d 个文件 → 推送已完成" % len(tracked))
    log("      · Zenodo → Settings → GitHub 里该仓库开关为 ON → 可以发 Release 了")
    log("      · Zenodo 出现 DOI 形如 10.5281/zenodo.1234567 → 可以跑 --writeback")

    log("\n自检结束。缺项按上面提示处理，其余按手册步骤走。")
    log("=" * 78)


# --------------------------------------------------------------- push to git --
def do_push(url):
    if not git_ok():
        raise SystemExit("本机没有 git，先安装 Git for Windows")

    log("=" * 78)
    log("步骤 2：配置身份 → 加 remote → 切到 main → 推送")
    log("=" * 78)

    r = run(["git", "config", "--global", "user.name"])
    if not r.stdout.strip():
        run(["git", "config", "--global", "user.name", GIT_NAME], check=True)
        log(OK + "已设置 user.name = %s" % GIT_NAME)
    else:
        log(OK + "user.name 已是 %s（未改动）" % r.stdout.strip())
    r = run(["git", "config", "--global", "user.email"])
    if not r.stdout.strip():
        run(["git", "config", "--global", "user.email", GIT_EMAIL], check=True)
        log(OK + "已设置 user.email = %s" % GIT_EMAIL)
    else:
        log(OK + "user.email 已是 %s（未改动）" % r.stdout.strip())

    has = run(["git", "remote"]).stdout.split()
    if has:
        run(["git", "remote", "set-url", "origin", url], cwd=REPO, check=True)
        log(OK + "origin 地址已更新为 %s" % url)
    else:
        run(["git", "remote", "add", "origin", url], cwd=REPO, check=True)
        log(OK + "已添加 origin = %s" % url)

    log("      → 切到 main 分支并推送（首次会要求认证：用户名=GitHub 用户名，密码=粘贴 PAT）")
    run(["git", "branch", "-M", "main"], cwd=REPO, check=True)
    r = subprocess.run(["git", "push", "-u", "origin", "main"], cwd=REPO)

    if r.returncode != 0:
        log("\n" + BAD + "推送失败。按下面三条自查：")
        log("      1. 报 'password authentication was removed' → 密码处要粘贴 PAT（不是账号密码）")
        log("      2. 报 'Permission denied (publickey)' → 你用的是 SSH 地址，但没配公钥；改回 https:// 地址")
        log("      3. 报 'remote origin already exists' / 'rejected' → 建仓时勾了 README，请删掉该仓库重建一个空仓库")
        raise SystemExit(1)

    local = run(["git", "rev-parse", "HEAD"], cwd=REPO).stdout.strip()
    remote = run(["git", "ls-remote", "--heads", "origin", "main"]).stdout.strip()
    ok = remote.startswith(local)
    log("\n" + (OK if ok else WARN) + "本地 HEAD = %s" % local[:12])
    log((OK if ok else WARN) + "远程 main = %s" % (remote.split()[0][:12] if remote else "(读不到)"))
    if ok:
        log(OK + "推送成功。下一步：打开仓库网页确认文件数，然后去 Zenodo（手册 §5 -> 5a）")
        log("\n      接着用这条命令把仓库地址写进发布元数据，并同步到远程：")
        log('      python scripts/70_github_publish.py --url %s --record-only' % url)
    log("=" * 78)


# ----------------------------------------------------------- record-only url --
def do_record_only(url):
    log("把仓库地址写入发布元数据（CITATION.cff / README），并推送该更新")
    run([PY, os.path.join(BASE, "scripts", "66_release_code.py"), "--url", url], check=True)
    run(["git", "add", "-A"], cwd=REPO, check=True)
    r = run(["git", "-c", "user.name=%s" % GIT_NAME, "-c", "user.email=%s" % GIT_EMAIL,
             "commit", "-q", "-m", "Record repository URL in metadata"], cwd=REPO)
    if r.returncode != 0:
        log(WARN + "没有需要提交的改动（可能已记录过）")
    else:
        log(OK + "已提交仓库地址")
    run(["git", "push"], cwd=REPO, check=True)
    log(OK + "已推送。下一步：Zenodo 开开关 → 发 Release v1.0.0 → 拿 DOI")


# ------------------------------------------------------------------ writeback --
def do_writeback(url, doi):
    log("=" * 78)
    log("步骤 6：把 DOI 写回稿件并重建全部交付物（顺序 67 → 49 → 64 → 65 → 66）")
    log("=" * 78)
    seq = [
        ("67 写回中英母稿的可用性声明", [PY, os.path.join(BASE, "scripts", "67_finalise_review2.py"),
                                 "--doi", doi, "--url", url]),
        ("49 重建中英母稿 docx", [PY, os.path.join(BASE, "scripts", "49_build_new_docx.py")]),
        ("64 重建 JTM 投稿包", [PY, os.path.join(BASE, "scripts", "64_build_jtm_pkg.py")]),
        ("65 重建 JTM 的 Word 稿", [PY, os.path.join(BASE, "scripts", "65_build_jtm_docx.py")]),
        ("66 重建代码发布树（并把 URL 与 DOI 记进元数据）",
         [PY, os.path.join(BASE, "scripts", "66_release_code.py"), "--url", url, "--doi", doi]),
    ]
    for label, cmd in seq:
        log("\n--- %s" % label)
        r = subprocess.run(cmd, cwd=BASE)
        if r.returncode != 0:
            raise SystemExit("这一步失败：%s —— 修好后可只重跑这一条" % label)

    log("\n--- 把更新后的发布树推到远程")
    run(["git", "push"], cwd=REPO, check=True)

    log("\n" + OK + "全部完成。请核对：")
    log("      · M2-M9深化研究稿/manuscript_EN.md 的可用性声明含 %s" % url)
    log("      · 同文件与 manuscript_CN.md 含 DOI %s" % doi)
    log("      · JTM投稿_M2M9/ 里 4 个文件（md / docx / README / BUILD_LOG）时间戳已更新")
    log("=" * 78)


# ------------------------------------------------------------------ verify tar --
def do_verify_tar():
    """校验归档本身（不依赖 sha256sum —— 本机 Git Bash 里没有这个命令）。"""
    tar = os.path.join(BASE, "JAML_M2M9_code_release.tar.gz")
    side = tar + ".sha256"
    log("=" * 78)
    log("校验归档 JAML_M2M9_code_release.tar.gz")
    log("=" * 78)
    for p in (tar, side):
        if not os.path.exists(p):
            raise SystemExit("找不到 %s —— 先跑：python scripts/66_release_code.py" % p)
    want = open(side, encoding="utf-8").read().split()[0]
    got = hashlib.sha256(open(tar, "rb").read()).hexdigest()
    size = os.path.getsize(tar)
    log("  文件大小 : %d 字节（%.2f MB）" % (size, size / 1e6))
    log("  应为哈希 : %s" % want)
    log("  实际哈希 : %s" % got)
    if want != got:
        raise SystemExit("  [失败] 与 sidecar 不一致 —— 归档被改动过，重新跑 scripts/66_release_code.py")
    log("  [OK]  一致，归档完好，可以上传 Zenodo")
    log("=" * 78)


# ------------------------------------------------------------------ final check --
def do_final_check():
    """发布后的核对（合并手册第 9、10 步）。

    只用 Python 实现：不依赖 `grep`/`head` —— 那两个命令在 PowerShell 里没有，
    在中文路径下还容易出编码问题（这也是手册初版会失败的原因之一）。
    """
    log("=" * 78)
    log("发布后核对：DOI 是否写回 + 交付物是否重建 + 许可是否一致")
    log("=" * 78)

    log("\n[A] 稿件里的仓库地址与 DOI")
    targets = [("M2-M9深化研究稿/manuscript_EN.md", "英文母稿"),
               ("M2-M9深化研究稿/manuscript_CN.md", "中文母稿"),
               ("JTM投稿_M2M9/manuscript_JTM.md", "JTM 稿")]
    doi_found = None
    for rel, name in targets:
        p = os.path.join(BASE, rel)
        if not os.path.exists(p):
            log(BAD + "%-8s 找不到 %s" % (name, rel))
            continue
        t = open(p, encoding="utf-8").read()
        i = t.find("10.5281/zenodo.")
        have_url = "github.com/" in t
        if i == -1:
            log(WARN + "%-8s 还没有 DOI（说明 --writeback 尚未跑：此时稿件里还是"
                       "'正存入公共仓库…'的旧话）" % name)
        else:
            doi_found = t[i:i + 60].split()[0].rstrip(").,;，。")
            log(OK + "%-8s DOI = %s   仓库地址 = %s"
                % (name, doi_found, "已写入" if have_url else "未见"))

    log("\n[B] 四个交付物是否刚重建过（时间戳应为刚才那一刻）")
    import datetime
    arts = [("M2-M9深化研究稿/01_Manuscript_EN.docx", "英文母稿 Word"),
            ("M2-M9深化研究稿/01_Manuscript_CN.docx", "中文母稿 Word"),
            ("JTM投稿_M2M9/01_Manuscript_JTM.docx", "JTM Word 稿"),
            ("JTM投稿_M2M9/manuscript_JTM.md", "JTM Markdown")]
    for rel, name in arts:
        p = os.path.join(BASE, rel)
        if os.path.exists(p):
            mt = datetime.datetime.fromtimestamp(os.path.getmtime(p))
            log(OK + "%-14s %s" % (name, mt.strftime("%Y-%m-%d %H:%M:%S")))
        else:
            log(BAD + "%-14s 不存在：%s" % (name, rel))

    log("\n[C] 许可三处是否一致")
    def first_line(p):
        try:
            return open(p, encoding="utf-8").read().splitlines()[0].strip()
        except Exception:
            return "(读不到)"

    def find_line(p, needle):
        try:
            for line in open(p, encoding="utf-8").read().splitlines():
                if needle in line:
                    return line.strip()
        except Exception:
            pass
        return "(未找到)"

    a = first_line(os.path.join(REPO, "LICENSE"))
    b = find_line(os.path.join(REPO, "CITATION.cff"), "license:")
    c = find_line(os.path.join(REPO, ".zenodo.json"), "license")
    d = find_line(os.path.join(REPO, "LICENSE"), "Copyright")
    log(OK + "LICENSE      : %s" % a)
    log(OK + "CITATION.cff : %s" % b)
    log(OK + ".zenodo.json : %s" % c)
    log(OK + "著作权人     : %s" % d)
    same = a.lower().startswith("mit") and "MIT" in b and "MIT" in c
    log("\n" + (OK if same else BAD) + "三处许可一致：%s" % ("是，MIT" if same else "否，请人工确认"))
    if doi_found:
        log("\n" + OK + "全部就绪：DOI %s 已写回，可以投稿" % doi_found)
    log("=" * 78)


def main():
    if "--final-check" in sys.argv or "--license-check" in sys.argv:
        do_final_check()
        return
    if "--verify-tar" in sys.argv:
        do_verify_tar()
        return
    if "--check" in sys.argv:
        do_check()
        return
    url = sys.argv[sys.argv.index("--url") + 1] if "--url" in sys.argv else None
    doi = sys.argv[sys.argv.index("--doi") + 1] if "--doi" in sys.argv else None
    if "--writeback" in sys.argv:
        if not (url and doi):
            raise SystemExit("--writeback 需要同时给 --url 和 --doi")
        do_writeback(url, doi)
    elif "--record-only" in sys.argv:
        if not url:
            raise SystemExit("--record-only 需要 --url")
        do_record_only(url)
    elif url:
        do_push(url)
    else:
        log(__doc__)


if __name__ == "__main__":
    main()
