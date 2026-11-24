"""In-process network guard: the workbench may talk to this computer (loopback) and nothing else.

Installed at the start of every `sliceweave` command and of the Windows executable. It uses
Python audit hooks (PEP 578), which fire for every socket connection, DNS lookup, datagram
send, socket bind, urllib/http.client request and browser launch made by Python code in this
process, including code inside third-party libraries. Anything that is not loopback is
refused with PermissionError before any packet is sent.

Verification: set SLICEWEAVE_NETTRACE=<file> and every network-related event (allowed or
blocked) plus every external program started is written to that file as JSON lines.
`python tools/audit_network.py --trace <file>` summarizes it.

Limits (see SECURITY.md): external programs (PowerShell for the BitLocker check) run in
their own processes and are not covered; native code that opens sockets without going through
Python would not be seen. None of the bundled native libraries (numpy, scipy, pydicom, GDCM)
has networking code paths that SliceWeave uses.
"""
from __future__ import annotations

import ipaddress
import json
import os
import sys
import threading
import time

_installed = False
_lock = threading.Lock()
_trace_path = None
blocked: list[dict] = []
counts = {"allowed": 0, "blocked": 0}
LOCAL_NAMES = {"localhost", "localhost.localdomain", "ip6-localhost", "ip6-loopback"}


class NetworkBlocked(PermissionError):
    pass


def _is_local_host(host) -> bool:
    if host is None:
        return True  # getaddrinfo(None, port) = local wildcard, used for binding
    if isinstance(host, bytes):
        host = host.decode("ascii", "replace")
    host = str(host).strip().strip("[]").lower()
    if host in ("",) or host in LOCAL_NAMES:
        return True
    if "%" in host:
        host = host.split("%", 1)[0]
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


def _is_local_address(addr) -> bool:
    if isinstance(addr, (str, bytes)):  # AF_UNIX path: local by definition
        return True
    if isinstance(addr, tuple) and addr:
        return _is_local_host(addr[0])
    return False


def _is_local_url(url) -> bool:
    from urllib.parse import urlsplit
    try:
        u = urlsplit(str(url))
    except ValueError:
        return False
    if u.scheme in ("file", ""):
        return True
    return _is_local_host(u.hostname)


def _record(event: str, target, allowed: bool):
    entry = {"t": round(time.time(), 3), "event": event, "target": _safe(target), "allowed": allowed}
    with _lock:
        counts["allowed" if allowed else "blocked"] += 1
        if not allowed:
            blocked.append(entry)
        if _trace_path:
            try:
                with open(_trace_path, "a", encoding="utf-8") as f:
                    f.write(json.dumps(entry) + "\n")
            except OSError:
                pass


def _safe(x):
    try:
        if isinstance(x, (str, int, float)) or x is None:
            return x
        if isinstance(x, bytes):
            return x.decode("ascii", "replace")
        if isinstance(x, (tuple, list)):
            return [_safe(v) for v in x][:6]
        return str(x)[:200]
    except Exception:  # noqa: BLE001
        return "?"


def _deny(event, target):
    _record(event, target, False)
    msg = (f"SliceWeave network guard: blocked {event} to {_safe(target)!r}. "
           "This program works offline and never sends data off this computer.")
    try:
        sys.stderr.write(msg + "\n")
    except Exception:  # noqa: BLE001
        pass
    raise NetworkBlocked(msg)


def _hook(event: str, args):
    # fast path: almost every audited event is irrelevant
    c = event[0]
    if c not in "suhw":
        return
    if event in ("socket.connect", "socket.sendto", "socket.sendmsg"):
        addr = args[1] if len(args) > 1 else None
        if event == "socket.sendmsg" and addr is None:
            return  # sendmsg on a connected socket: the connect was already checked
        if _is_local_address(addr):
            _record(event, addr, True)
        else:
            _deny(event, addr)
    elif event == "socket.bind":
        addr = args[1] if len(args) > 1 else None
        if isinstance(addr, tuple) and addr and str(addr[0]) in ("", "0.0.0.0", "::"):
            _deny(event, addr)  # listening on all interfaces would expose the server
        elif _is_local_address(addr):
            _record(event, addr, True)
        else:
            _deny(event, addr)
    elif event in ("socket.getaddrinfo", "socket.gethostbyname", "socket.gethostbyname_ex",
                   "socket.gethostbyaddr"):
        host = args[0] if args else None
        if _is_local_host(host):
            _record(event, host, True)
        else:
            _deny(event, host)
    elif event == "socket.getnameinfo":
        sa = args[0] if args else None
        if _is_local_address(sa):
            _record(event, sa, True)
        else:
            _deny(event, sa)
    elif event == "urllib.Request":
        url = args[0] if args else None
        if _is_local_url(url):
            _record(event, url, True)
        else:
            _deny(event, url)
    elif event == "http.client.connect":
        host = args[1] if len(args) > 1 else None
        if _is_local_host(host):
            _record(event, host, True)
        else:
            _deny(event, host)
    elif event == "webbrowser.open":
        url = args[0] if args else None
        if _is_local_url(url):
            _record(event, url, True)
        else:
            _deny(event, url)
    elif event == "subprocess.Popen" and _trace_path:
        exe = args[0] if args else None
        argv = args[1] if len(args) > 1 else None
        name = os.path.basename(str(exe or (argv[0] if argv else "")))
        _record("subprocess.Popen (separate process; not covered by this guard)", name, True)


def install(trace_path: str | None = None) -> bool:
    """Install once per process. Returns True if the guard is active."""
    global _installed, _trace_path
    if _installed:
        return True
    _trace_path = trace_path or os.environ.get("SLICEWEAVE_NETTRACE") or None
    sys.addaudithook(_hook)
    _installed = True
    return True


def status() -> dict:
    return {"active": _installed, "allowed_local_events": counts["allowed"], "blocked": counts["blocked"],
            "trace": bool(_trace_path)}
