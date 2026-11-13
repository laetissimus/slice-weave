"""Offline guarantees: the network guard and a static scan of the program code."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

GUARD_CHECK = r"""
import json, socket, sys, urllib.request
sys.path.insert(0, %r)
from sliceweave import netguard
netguard.install()
r = {}
def attempt(name, fn):
    try:
        fn(); r[name] = "allowed"
    except PermissionError:
        r[name] = "blocked"
    except OSError as e:
        r[name] = "allowed (then %%s)" %% type(e).__name__
attempt("connect_public_ip", lambda: socket.create_connection(("93.184.215.14", 80), 2))
attempt("dns_public_name", lambda: socket.getaddrinfo("example.com", 443))
attempt("udp_public", lambda: socket.socket(socket.AF_INET, socket.SOCK_DGRAM).sendto(b"x", ("8.8.8.8", 53)))
attempt("urllib_public", lambda: urllib.request.urlopen("https://example.com", timeout=2))
attempt("bind_all_interfaces", lambda: socket.socket().bind(("0.0.0.0", 0)))
srv = socket.socket(); srv.bind(("127.0.0.1", 0)); srv.listen()
attempt("loopback_connect", lambda: socket.create_connection(srv.getsockname(), 2).close())
attempt("dns_localhost", lambda: socket.getaddrinfo("localhost", 80))
print(json.dumps(r))
"""


def test_guard_blocks_everything_but_loopback(tmp_path):
    import os
    env = {k: v for k, v in os.environ.items() if k != "SLICEWEAVE_NETTRACE"}  # keep deliberate probes out of program traces
    out = subprocess.run([sys.executable, "-c", GUARD_CHECK % str(ROOT)], capture_output=True, text=True, timeout=60, env=env)
    r = json.loads(out.stdout.strip().splitlines()[-1])
    assert r["connect_public_ip"] == "blocked"
    assert r["dns_public_name"] == "blocked"
    assert r["udp_public"] == "blocked"
    assert r["urllib_public"] == "blocked"
    assert r["bind_all_interfaces"] == "blocked"
    assert r["loopback_connect"] == "allowed"
    assert r["dns_localhost"].startswith("allowed")


def test_static_scan_finds_no_outbound_code():
    out = subprocess.run([sys.executable, str(ROOT / "tools" / "audit_network.py")], capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stdout + out.stderr
