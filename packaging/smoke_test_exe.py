"""Test a built executable the way a user would run it (used by the Windows build workflow).

  python packaging/smoke_test_exe.py dist/SliceWeave.exe [--full]

Checks: --version, doctor, the server (login link, session, security headers, case creation),
and with --full the end-to-end self-test run by the executable itself."""
import json
import os
import queue
import re
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import requests


def main():
    exe = str(Path(sys.argv[1]).resolve())
    full = "--full" in sys.argv
    data = tempfile.mkdtemp(prefix="sliceweave-smoke-")
    env = dict(os.environ, SLICEWEAVE_DATA=data, PYTHONUNBUFFERED="1")
    results = {}

    def run(args, timeout):
        t = time.time()
        r = subprocess.run([exe] + args, capture_output=True, text=True, timeout=timeout, env=env)
        print(f"$ {Path(exe).name} {' '.join(args)}  (exit {r.returncode}, {time.time() - t:.0f} s)")
        print(r.stdout[-4000:])
        if r.stderr.strip():
            print(r.stderr[-4000:])
        return r

    r = run(["--version"], 300)
    results["version"] = r.returncode == 0 and "SliceWeave" in r.stdout
    r = run(["doctor"], 300)
    results["doctor"] = r.returncode == 0 and "numpy" in r.stdout and "MISSING" not in r.stdout
    results["network_guard_active"] = "Network guard: ACTIVE" in r.stdout
    # DICOM decoding inside the executable: JPEG-LS sample must match the uncompressed one
    dd = Path(__file__).resolve().parents[1] / "tests" / "data" / "dicom"
    r = run(["dicom-info", str(dd / "MR_small.dcm"), str(dd / "MR_small_jpeg_ls_lossless.dcm")], 300)
    sums = re.findall(r"pixels sha256 ([0-9a-f]+)", r.stdout)
    results["jpeg_dicom_decoding"] = len(sums) == 2 and sums[0] == sums[1]
    # server
    port = 8791
    p = subprocess.Popen([exe, "serve", "--no-browser", "--port", str(port)], stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, env=env)
    q: queue.Queue = queue.Queue()
    threading.Thread(target=lambda: [q.put(line) for line in p.stdout], daemon=True).start()
    url, deadline, out = None, time.time() + 300, []
    while time.time() < deadline and url is None:
        try:
            line = q.get(timeout=1)
        except queue.Empty:
            continue
        out.append(line)
        m = re.search(r"(http://127\.0\.0\.1:\d+/login\?token=\S+)", line)
        if m:
            url = m.group(1)
    print("".join(out).replace(url or "§", "<login link>"))
    ok = False
    if url:
        s = requests.Session()
        s.trust_env = False
        for _ in range(60):
            try:
                if s.get(f"http://127.0.0.1:{port}/api/cases", timeout=2).status_code == 401:
                    break
            except requests.ConnectionError:
                time.sleep(1)
        r1 = s.get(url, allow_redirects=False, timeout=10)
        sess = s.get(f"http://127.0.0.1:{port}/api/session", timeout=10)
        page = s.get(f"http://127.0.0.1:{port}/", timeout=10)
        js = s.get(f"http://127.0.0.1:{port}/static/app.js", timeout=10)
        csrf = sess.json().get("csrf", "") if sess.ok else ""
        made = s.post(f"http://127.0.0.1:{port}/api/cases", json={"case_id": "SMOKE-01"}, headers={"X-SW-CSRF": csrf}, timeout=10)
        nocsrf = s.post(f"http://127.0.0.1:{port}/api/cases", json={"case_id": "SMOKE-02"}, timeout=10)
        ok = (r1.status_code == 303 and sess.ok and "default-src 'self'" in page.headers.get("content-security-policy", "")
              and js.ok and made.ok and nocsrf.status_code == 403)
        results["server_encryption_check"] = sess.json().get("encryption") if sess.ok else None
        results["server_network_guard_active"] = bool(sess.ok and sess.json().get("network_guard", {}).get("active"))
    results["server"] = ok
    if os.name == "nt":
        subprocess.run(["taskkill", "/F", "/T", "/PID", str(p.pid)], capture_output=True)
    else:
        p.terminate()
    if full:
        r = run(["selftest"], 3600)
        results["selftest"] = r.returncode == 0 and "self-test PASSED" in r.stdout
    print(json.dumps(results, indent=1))
    rep = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--report=")), "smoke-report.json")
    Path(rep).write_text(json.dumps(results, indent=1), encoding="utf-8")
    sys.exit(0 if all(v for k, v in results.items() if k != "server_encryption_check") else 1)


if __name__ == "__main__":
    main()
