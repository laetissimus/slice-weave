"""Network audit of SliceWeave.

Static (default): lists every place in the program's own code (sliceweave/ and the executable
launcher) that could open a network connection, start another program, or load a resource
from a URL, and fails if any is not on the reviewed list below or if an external
(non-loopback) address appears in the code.

    python tools/audit_network.py

Trace: summarizes a SLICEWEAVE_NETTRACE file written while the program ran, and fails if any
network event targeted anything other than this computer.

    python tools/audit_network.py --trace trace.jsonl
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROGRAM = sorted((ROOT / "sliceweave").rglob("*.py")) + [ROOT / "packaging" / "launcher.py"]
WEB = sorted((ROOT / "sliceweave" / "web").glob("*"))

PY_PATTERNS = {
    "socket": r"\bimport socket\b|\bsocket\.",
    "urllib.request": r"urllib\.request|urlopen\(",
    "http.client": r"http\.client|HTTPConnection",
    "requests/httpx/aiohttp": r"\bimport (requests|httpx|aiohttp)\b|\b(requests|httpx|aiohttp)\.(get|post|put|Session|Client)",
    "ftp/smtp/telnet/xmlrpc": r"ftplib|smtplib|telnetlib|xmlrpc",
    "websocket": r"websocket",
    "asyncio network": r"open_connection|create_connection|start_server",
    "server bind": r"uvicorn\.run|host=",
    "browser launch": r"webbrowser\.",
    "external program": r"subprocess\.|os\.system|os\.popen|Popen\(",
}
JS_PATTERNS = {
    "fetch": r"\bfetch\(",
    "XMLHttpRequest": r"XMLHttpRequest",
    "WebSocket/EventSource": r"WebSocket|EventSource",
    "sendBeacon": r"sendBeacon",
    "importScripts": r"importScripts",
    "external resource tag": r"<(script|link|img|iframe)[^>]+(src|href)=[\"']?(https?:)?//",
    "css import/url": r"@import|url\(\s*[\"']?(https?:)?//",
}
URL = re.compile(r"(?:https?|wss?|ftp)://([A-Za-z0-9.\-\[\]:]+)")
LOCAL = re.compile(r"^(127\.0\.0\.1|localhost|\[::1\])(:\d+)?$")

# Reviewed uses: (file name, category) -> why it is safe
REVIEWED = {
    ("netguard.py", "socket"): "the guard itself: inspects socket events, opens nothing",
    ("netguard.py", "urllib.request"): "the guard itself: blocks urllib requests to non-local hosts",
    ("netguard.py", "http.client"): "the guard itself: blocks http.client connections to non-local hosts",
    ("netguard.py", "browser launch"): "the guard itself: allows only local URLs to be opened",
    ("netguard.py", "external program"): "the guard itself: records external programs in the trace",
    ("server.py", "server bind"): "uvicorn bound to 127.0.0.1 only; no option to change it",
    ("server.py", "browser launch"): "opens the local login link http://127.0.0.1:<port>/login in the user's browser",
    ("security.py", "external program"): "local OS checks: fdesetup / PowerShell BitLocker property / manage-bde / df / lsblk",
    ("cli.py", "server bind"): "false positive: argparse help text",
    ("app.js", "fetch"): "same-origin only: every call goes through api() with a path starting /api/",
}


def scan_static() -> int:
    findings, problems = [], []
    for f in PROGRAM + WEB:
        if not f.is_file() or f.suffix not in (".py", ".js", ".html", ".css"):
            continue
        text = f.read_text(encoding="utf-8", errors="replace")
        pats = PY_PATTERNS if f.suffix == ".py" else JS_PATTERNS
        for cat, pat in pats.items():
            hits = [n + 1 for n, line in enumerate(text.splitlines()) if re.search(pat, line)]
            if hits:
                why = REVIEWED.get((f.name, cat))
                findings.append((f.relative_to(ROOT), cat, hits, why))
                if why is None:
                    problems.append(f"{f.relative_to(ROOT)}: {cat} at lines {hits} is not reviewed")
        for n, line in enumerate(text.splitlines(), 1):
            for m in URL.finditer(line):
                host = m.group(1).rstrip(":")
                if not LOCAL.match(host):
                    problems.append(f"{f.relative_to(ROOT)}:{n}: external address {m.group(0)}")
    if any(f.name == "app.js" for f in WEB):
        js = (ROOT / "sliceweave" / "web" / "app.js").read_text(encoding="utf-8")
        calls = re.findall(r"""api\(\s*"(?:GET|POST|PUT|DELETE)"\s*,\s*([`'"])(.*?)\1""", js)
        bad = [p for _, p in calls if not p.startswith("/api/")]
        if bad:
            problems.append(f"app.js: api() calls to non-API paths {bad}")
        fetches = re.findall(r"\bfetch\(([^,)]*)", js)
        if any(a.strip() != "path" for a in fetches):
            problems.append(f"app.js: fetch() outside the api() helper: {fetches}")
        csp = (ROOT / "sliceweave" / "server.py").read_text(encoding="utf-8")
        if "default-src 'self'" not in csp or "connect-src 'self'" not in csp:
            problems.append("server.py: Content-Security-Policy no longer restricts the page to this server")
    print("Network-relevant code in SliceWeave (program files only):")
    for path, cat, hits, why in findings:
        print(f"  {str(path):34s} {cat:24s} lines {hits}\n      -> {why or 'NOT REVIEWED'}")
    print(f"\n{len(PROGRAM)} Python files and {len(WEB)} web files scanned.")
    if problems:
        print("PROBLEMS:")
        for p in problems:
            print("  " + p)
        return 1
    print("Result: no code that sends data off this computer; the only external programs are the local tools listed.")
    return 0


def summarize_trace(path: str) -> int:
    rows = [json.loads(line) for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    by_event = Counter(r["event"] for r in rows)
    targets = Counter()
    nonlocal_rows = []
    for r in rows:
        t = r["target"]
        host = t[0] if isinstance(t, list) and t else t
        if r["event"].startswith("subprocess"):
            targets[f"program: {t}"] += 1
            continue
        targets[f"{r['event']}: {host}"] += 1
        if not r["allowed"]:
            nonlocal_rows.append(r)
    print(f"Trace {path}: {len(rows)} network-related events")
    for k, v in sorted(by_event.items()):
        print(f"  {k:62s} {v}")
    print("Targets:")
    for k, v in sorted(targets.items()):
        print(f"  {k:62s} {v}")
    if nonlocal_rows:
        print(f"BLOCKED non-local attempts: {len(nonlocal_rows)}")
        for r in nonlocal_rows[:20]:
            print("  ", r)
        return 1
    print("Result: every network event targeted this computer (loopback); no attempt to reach anything else.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--trace":
        sys.exit(summarize_trace(sys.argv[2]))
    sys.exit(scan_static())
