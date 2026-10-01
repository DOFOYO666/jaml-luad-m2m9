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

  2b) 若推送报"连接被重置"（大陆网络常见）——先看网络诊断，再给 git 配代理：
       python scripts/70_github_publish.py --net
       python scripts/70_github_publish.py --url https://github.com/<你的用户名>/jaml-luad-m2m9 --proxy http://127.0.0.1:29290
       python scripts/70_github_publish.py --proxy none          # 用完清掉

  3) 从 Zenodo 拿到 DOI 之后，一条命令完成"写回稿件 + 重建全部交付物 + 再推送一次"：
       python scripts/70_github_publish.py --url https://github.com/<你的用户名>/jaml-luad-m2m9 \
            --doi 10.5281/zenodo.1234567 --writeback

本脚本**不联网建仓库、不申请 DOI**（那需要你的账号），只把本机能做的那部分做完。
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time


def _init_console():
    """Windows 控制台不要强行把 stdout 改成 UTF-8。

    PEP 528 已让中文在控制台正确显示；强行改编码反而会让 5.1 的控制台把 UTF-8 字节
    按 GBK 解码，显示成"鍙戝竷鍓嶈嚜妗€"这类乱码（实测）。
    非 Windows（本项目的沙箱 bash）才设 UTF-8。
    """
    try:
        if os.name == "nt":
            sys.stdout.reconfigure(errors="replace")
        else:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


_init_console()

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.join(BASE, "JAML_M2M9_code_release")
PY = sys.executable

GIT_CANDIDATES = [
    r"C:\Users\86159\.workbuddy\vendor\PortableGit\cmd\git.exe",
    r"C:\Users\86159\.workbuddy\vendor\PortableGit\mingw64\bin\git.exe",
    r"C:\Program Files\Git\cmd\git.exe",
    r"C:\Program Files (x86)\Git\cmd\git.exe",
    os.path.join(os.environ.get("LOCALAPPDATA", ""), "Programs", "Git", "cmd", "git.exe"),
]


def find_git():
    """返回 (git 可执行文件, 来源说明)。

    本机**没有装 Git for Windows**，git 来自 WorkBuddy 自带的 PortableGit，它的目录不在
    PATH 里 —— 所以在 PowerShell 窗口里敲 `git` 会直接返回 9009（找不到命令）。
    这里解析出绝对路径，后面所有 git 调用都用它。
    """
    p = shutil.which("git")
    if p:
        return p, "PATH"
    for c in GIT_CANDIDATES:
        if c and os.path.isfile(c):
            return c, "绝对路径（不在 PATH 里）"
    return None, "未找到"


GIT, GIT_SOURCE = find_git()

GIT_NAME = "Huayong Liu"
GIT_EMAIL = "xingxinghuoshu@163.com"

OK, WARN, BAD = "  [OK]  ", "  [注意] ", "  [缺]  "


def log(*a):
    print(*a, flush=True)


def net_probe(quiet=False):
    """纯 Python 探测网络可达性与系统代理。

    不用 ping / curl / Test-NetConnection —— 那些在别的 shell 里没有；
    socket 直连在 PowerShell / Git Bash / cmd 下行为一致。
    """
    import socket
    log("\n[网络] 目标站点 443 端口连通性")
    results = {}
    for host in ("github.com", "zenodo.org"):
        try:
            t0 = time.time()
            with socket.create_connection((host, 443), timeout=8):
                dt = time.time() - t0
            results[host] = True
            log(OK + "%-12s 可连通（%.1f 秒）" % (host, dt))
        except Exception as e:
            results[host] = False
            log(BAD + "%-12s 连不上：%s" % (host, "%s: %s" % (type(e).__name__, str(e)[:60])))
    log("\n[网络] 系统代理设置")
    proxy = None
    if os.name == "nt":
        try:
            import winreg
            k = winreg.OpenKey(winreg.HKEY_CURRENT_USER,
                               r"Software\Microsoft\Windows\CurrentVersion\Internet Settings")
            enable = winreg.QueryValueEx(k, "ProxyEnable")[0]
            try:
                server = winreg.QueryValueEx(k, "ProxyServer")[0]
            except FileNotFoundError:
                server = "(未设置)"
            try:
                pac = winreg.QueryValueEx(k, "AutoConfigURL")[0]
            except FileNotFoundError:
                pac = "(未设置)"
            log(OK + "ProxyEnable = %s    ProxyServer = %s    AutoConfigURL = %s"
                % (enable, server, pac))
            if enable and server and server != "(未设置)":
                s = str(server)
                if not s.lower().startswith("http"):
                    s = "http://" + s
                proxy = s          # 归一化后再用，避免拼出 http://http://...
        except Exception as e:
            log(WARN + "读系统代理失败：%s" % e)
    gp = run([_g(), "config", "--global", "--get", "http.proxy"]).stdout.strip()
    log((OK if gp else WARN) + "git 自己的 http.proxy = %s" % (gp or "(未设置)"))
    if proxy and not gp:
        log("\n" + WARN + "关键线索：系统有代理（%s）但 git 没配 → git 是直连，所以被重置。" % proxy)
        log("      正解：给 git 配上同一个代理后重试（把端口换成你实际看到的）：")
        log('      <本脚本> --url <仓库地址> --proxy %s' % proxy)
        log("      用完后清掉： <本脚本> --proxy none")
    if not results.get("github.com") and not proxy:
        log("\n" + BAD + "github.com 直连不通且没有系统代理 → 属于网络限制，三条出路：")
        log("      ① 开加速/代理后重试（推荐，见上面那条命令）")
        log("      ② 换网络（手机热点常常能通）后再 push")
        log("      ③ 放弃 GitHub，改用 Zenodo 直传上传归档拿 DOI（手册附录 B，不需要 push）")
    return results


def _set_git_proxy(value):
    """给 git 配/清 http(s) 代理。value 为 'none'/'' 时清除。"""
    if value in (None, "", "none", "None"):
        for k in ("http.proxy", "https.proxy"):
            run([_g(), "config", "--global", "--unset", k])
        log(OK + "已清除 git 的 http/https 代理")
        return None
    for k in ("http.proxy", "https.proxy"):
        run([_g(), "config", "--global", k, value], check=True)
    log(OK + "已给 git 设置代理：%s（http 与 https 都设了）" % value)
    return value


def _g():
    """git 可执行文件；找不到时退回裸 "git"（由 run() 给出友好提示）。"""
    return GIT or "git"


def _norm_url(u):
    """清洗用户粘贴的仓库地址。

    实测用户会从占位符 `<你的用户名>` 里把尖括号一起带出来，也会带上 `.git` 后缀；
    这两者在 cmd 里 `<` 还是**输入重定向**、会直接报"文件名、目录名或卷标语法不正确"。
    """
    if not u:
        return u
    v = u.strip().strip("<>").strip()
    while v.endswith("/"):
        v = v[:-1]
    if v.endswith(".git"):
        v = v[:-4]
    return v


def _norm_doi(d):
    """清洗并粗校验 DOI。"""
    if not d:
        return d
    v = d.strip().strip("<>").strip().rstrip(".,;，。")
    if v.lower().startswith("doi:"):
        v = v[4:].strip()
    if v.lower().startswith("https://doi.org/"):
        v = v[len("https://doi.org/"):]
    if not re.match(r"^10\.\d{4,9}/\S+$", v):
        log(WARN + "'%s' 看起来不像合法 DOI（应为 10.xxxx/…）；仍会按你给的值写入，请自行确认" % v)
    return v


def _norm_url(u):
    """清洗用户粘贴的仓库地址。

    实测用户会从占位符 `<你的用户名>` 里把尖括号一起带出来，也会带上 `.git` 后缀；
    这两者在 cmd 里 `<` 还是**输入重定向**、会直接报"文件名、目录名或卷标语法不正确"。
    """
    if not u:
        return u
    v = u.strip().strip("<>").strip()
    while v.endswith("/"):
        v = v[:-1]
    if v.endswith(".git"):
        v = v[:-4]
    return v


def _norm_doi(d):
    """清洗并粗校验 DOI。"""
    if not d:
        return d
    v = d.strip().strip("<>").strip().rstrip(".,;，。")
    if v.lower().startswith("doi:"):
        v = v[4:].strip()
    if v.lower().startswith("https://doi.org/"):
        v = v[len("https://doi.org/"):]
    if not re.match(r"^10\.\d{4,9}/\S+$", v):
        log(WARN + "'%s' 看起来不像合法 DOI（应为 10.xxxx/…）；仍会按你给的值写入，请自行确认" % v)
    return v


class _Missing:
    """占位返回值：让调用方按 returncode != 0 处理，而不是抛出 FileNotFoundError 堆栈。"""

    returncode = 127
    stdout = ""
    stderr = ""


def run(args, cwd=None, capture=True, check=False):
    try:
        r = subprocess.run(args, cwd=cwd, capture_output=capture, text=True)
    except FileNotFoundError:
        r = _Missing()
        r.stderr = "找不到可执行文件：%s" % args[0]
        log(WARN + r.stderr + "（见自检 [1] 提示）")
    if check and r.returncode != 0:
        raise SystemExit("命令失败：%s\n%s"
                         % (" ".join(str(a) for a in args),
                            (r.stderr or r.stdout or "").strip()[:400]))
    return r


def git_ok():
    try:
        run([_g(), "--version"], check=True)
        return True
    except Exception:
        return False


def current_branch():
    r = run([_g(), "rev-parse", "--abbrev-ref", "HEAD"], cwd=REPO)
    return r.stdout.strip() if r.returncode == 0 else "?"


# --------------------------------------------------------------------- check --
def do_check():
    log("=" * 78)
    log("发布前自检（只看不改）")
    log("=" * 78)

    log("\n[1] 本机工具")
    log(OK + "python : %s" % sys.executable)
    if GIT:
        ver = run([_g(), "--version"]).stdout.strip()
        log(OK + "git    : %s  [%s]  %s" % (GIT, GIT_SOURCE, ver))
    else:
        log(BAD + "git    : 未找到 —— 装一个 Git for Windows（https://git-scm.com/download/win），"
                 "或在本机重新生成发布树后重试")

    log("\n[2] git 全局身份（提交时会用）")
    name = run([_g(), "config", "--global", "user.name"]).stdout.strip()
    mail = run([_g(), "config", "--global", "user.email"]).stdout.strip()
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
    head = run([_g(), "log", "--oneline", "-1"], cwd=REPO).stdout.strip()
    tracked = run([_g(), "ls-files"], cwd=REPO).stdout.splitlines()
    log(OK + "最新提交：%s" % head)
    log(OK + "跟踪文件：%d 个" % len(tracked))
    bad = [l for l in tracked if l.startswith("_superseded/")]
    log((OK if not bad else BAD) + "回收区未上传：%s" % ("是（0 个）" if not bad else "%d 个会被推上去！" % len(bad)))
    ls = os.path.join(BASE, "JAML_M2M9_code_release.tar.gz.sha256")
    log((OK if os.path.exists(ls) else WARN) + "归档校验文件：%s"
        % (open(ls, encoding="utf-8").read().strip()[:24] + "…" if os.path.exists(ls) else "缺（跑 66 生成）"))

    log("\n[4] 远程仓库")
    rem = run([_g(), "remote", "-v"], cwd=REPO).stdout.strip()
    if rem:
        log(OK + "已配置 remote：")
        for line in rem.splitlines():
            log("         " + line)
    else:
        log(WARN + "尚未配置 remote —— 网页建好空仓库后跑步骤 2 的命令即可")

    net_probe()

    log("\n[5] 网页操作是否已完成的判断依据")
    log("      · 仓库页能看到 %d 个文件 → 推送已完成" % len(tracked))
    log("      · Zenodo → Settings → GitHub 里该仓库开关为 ON → 可以发 Release 了")
    log("      · Zenodo 出现 DOI 形如 10.5281/zenodo.1234567 → 可以跑 --writeback")

    log("\n自检结束。缺项按上面提示处理，其余按手册步骤走。")
    log("=" * 78)


# --------------------------------------------------------------- push to git --
def do_push(url, proxy=None):
    if not git_ok():
        raise SystemExit("本机没有 git，先安装 Git for Windows")

    log("=" * 78)
    log("步骤 2：配置身份 → 加 remote → 切到 main → 推送")
    log("=" * 78)
    if proxy:
        _set_git_proxy(proxy)

    r = run([_g(), "config", "--global", "user.name"])
    if not r.stdout.strip():
        run([_g(), "config", "--global", "user.name", GIT_NAME], check=True)
        log(OK + "已设置 user.name = %s" % GIT_NAME)
    else:
        log(OK + "user.name 已是 %s（未改动）" % r.stdout.strip())
    r = run([_g(), "config", "--global", "user.email"])
    if not r.stdout.strip():
        run([_g(), "config", "--global", "user.email", GIT_EMAIL], check=True)
        log(OK + "已设置 user.email = %s" % GIT_EMAIL)
    else:
        log(OK + "user.email 已是 %s（未改动）" % r.stdout.strip())

    # 注意：必须带 cwd=REPO —— 否则在非仓库目录里 `git remote` 会失败并返回空，
    # 脚本就会误判"没有 remote"而去 add，撞上 "remote origin already exists."（实测踩过）。
    has = run([_g(), "remote"], cwd=REPO).stdout.split()
    if has:
        run([_g(), "remote", "set-url", "origin", url], cwd=REPO, check=True)
        log(OK + "origin 地址已更新为 %s" % url)
    else:
        r = run([_g(), "remote", "add", "origin", url], cwd=REPO)
        if r.returncode != 0 and "already exists" in (r.stderr or r.stdout or ""):
            # 双保险：万一还是判定错了，也不要中止，直接改地址
            run([_g(), "remote", "set-url", "origin", url], cwd=REPO, check=True)
            log(OK + "origin 已存在，改为更新地址 = %s" % url)
        else:
            log(OK + "已添加 origin = %s" % url)

    log("      → 切到 main 分支并推送（首次会要求认证：用户名=GitHub 用户名，密码=粘贴 PAT）")
    run([_g(), "branch", "-M", "main"], cwd=REPO, check=True)
    r = subprocess.run([_g(), "push", "-u", "origin", "main"], cwd=REPO)

    if r.returncode != 0:
        # 先把网络层查清楚，再给建议 —— 连接被重置 vs 认证失败 vs 仓库冲突，处理方式完全不同
        net_probe()
        log("\n" + BAD + "推送失败。按报错原文对号入座：")
        log("      含 'Connection was reset' / 'unable to access' / 'Failed to connect' / 超时"
            " → 网络问题（不是账号问题）：看上面 [网络] 两节，多半要给 git 配代理，"
            "或换网络，或改用 Zenodo 直传（手册附录 B）")
        log("      含 'password authentication was removed' → 密码处要粘贴 PAT（不是账号密码）")
        log("      含 'Permission denied (publickey)' → 你用了 SSH 地址但没配公钥；改回 https:// 地址")
        log("      含 'rejected' / 'fetch first' → 建仓时勾了 README；删掉该仓库，重建一个空仓库")
        log("      含 'Repository not found' → 仓库名或用户名写错，或仓库是 Private 而令牌没勾 repo 权限")
        raise SystemExit(1)

    local = run([_g(), "rev-parse", "HEAD"], cwd=REPO).stdout.strip()
    remote = run([_g(), "ls-remote", "--heads", "origin", "main"]).stdout.strip()
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
    run([_g(), "add", "-A"], cwd=REPO, check=True)
    r = run([_g(), "-c", "user.name=%s" % GIT_NAME, "-c", "user.email=%s" % GIT_EMAIL,
             "commit", "-q", "-m", "Record repository URL in metadata"], cwd=REPO)
    if r.returncode != 0:
        log(WARN + "没有需要提交的改动（可能已记录过）")
    else:
        log(OK + "已提交仓库地址")
    run([_g(), "push"], cwd=REPO, check=True)
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
    run([_g(), "push"], cwd=REPO, check=True)

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
    if "--net" in sys.argv:
        net_probe()
        return
    if "--proxy" in sys.argv:
        _set_git_proxy(sys.argv[sys.argv.index("--proxy") + 1])
        if not url:
            log("\n下一步：重跑带 --url 的推送命令")
            return
    proxy = sys.argv[sys.argv.index("--proxy") + 1] if "--proxy" in sys.argv else None

    # 统一清洗用户输入（去掉 <...>、结尾 /、.git、doi: 前缀等），并告知被改了什么
    if url:
        u2 = _norm_url(url)
        if u2 != url:
            log("注意：仓库地址已规范化为 %s" % u2)
            url = u2
    if doi:
        d2 = _norm_doi(doi)
        if d2 != doi:
            log("注意：DOI 已规范化为 %s" % d2)
            doi = d2

    # 统一清洗用户输入（去掉 <...>、结尾 /、.git、doi: 前缀等），并告知被改了什么
    if url:
        u2 = _norm_url(url)
        if u2 != url:
            log("注意：仓库地址已规范化为 %s" % u2)
            url = u2
    if doi:
        d2 = _norm_doi(doi)
        if d2 != doi:
            log("注意：DOI 已规范化为 %s" % d2)
            doi = d2
    if "--writeback" in sys.argv:
        if not (url and doi):
            raise SystemExit("--writeback 需要同时给 --url 和 --doi")
        do_writeback(url, doi)
    elif "--record-only" in sys.argv:
        if not url:
            raise SystemExit("--record-only 需要 --url")
        do_record_only(url)
    elif url:
        do_push(url, proxy)
    else:
        log(__doc__)


if __name__ == "__main__":
    main()
