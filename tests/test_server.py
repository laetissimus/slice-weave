"""Security controls of the local server, and an end-to-end case through the API."""
import io
import socket
import threading
import time
import zipfile

import pytest
import requests

from sliceweave import config
from sliceweave.server import create_app


@pytest.fixture(scope="module")
def srv(tmp_path_factory):
    import uvicorn
    s = config.Settings(data_dir=str(tmp_path_factory.mktemp("d"))).ensure()
    sock = socket.socket()
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()
    app = create_app(s, "tok-xyz", port)
    server = uvicorn.Server(uvicorn.Config(app, host="127.0.0.1", port=port, log_level="error"))
    t = threading.Thread(target=server.run, daemon=True)
    t.start()
    for _ in range(50):
        try:
            requests.get(f"http://127.0.0.1:{port}/", timeout=0.5)
            break
        except requests.ConnectionError:
            time.sleep(0.1)
    yield f"http://127.0.0.1:{port}"
    server.should_exit = True


def _session():
    s = requests.Session()
    s.trust_env = False
    return s


def _login(srv):
    s = _session()
    s.get(srv + "/login?token=tok-xyz", allow_redirects=False)
    csrf = s.get(srv + "/api/session").json()["csrf"]
    s.headers["X-SW-CSRF"] = csrf
    return s


def test_requires_login(srv):
    s = _session()
    assert s.get(srv + "/api/cases").status_code == 401
    assert s.get(srv + "/login?token=wrong", allow_redirects=False).status_code == 401


def test_host_header_checked(srv):
    assert _session().get(srv + "/", headers={"Host": "evil.example"}).status_code == 403


def test_login_csrf_origin(srv):
    s = _session()
    r = s.get(srv + "/login?token=tok-xyz", allow_redirects=False)
    assert r.status_code == 303 and "HttpOnly" in r.headers["set-cookie"] and "SameSite=strict" in r.headers["set-cookie"]
    csrf = s.get(srv + "/api/session").json()["csrf"]
    assert s.post(srv + "/api/cases", json={"case_id": "FET-T-9"}).status_code == 403
    assert s.post(srv + "/api/cases", json={"case_id": "FET-T-9"},
                  headers={"X-SW-CSRF": csrf, "Origin": "https://evil.example"}).status_code == 403
    assert s.post(srv + "/api/cases", json={"case_id": "MRN12345678"}, headers={"X-SW-CSRF": csrf}).status_code == 400
    assert s.post(srv + "/api/cases", json={"case_id": "FET-T-9"}, headers={"X-SW-CSRF": csrf}).status_code == 200
    r = s.get(srv + "/")
    assert "default-src 'self'" in r.headers["content-security-policy"]


def test_upload_import_plan_layers(srv, tmp_path):
    """NIfTI stacks uploaded through the API become a planned case with viewer layers."""
    import numpy as np

    from sliceweave.nifti import write_nifti
    s = _login(srv)
    assert s.post(srv + "/api/cases", json={"case_id": "API-T-1"}).ok
    rng = np.random.default_rng(1)
    files = []
    for n, (perm, name) in enumerate([((0, 1, 2), "ax"), ((0, 2, 1), "cor"), ((2, 0, 1), "sag")]):
        d = np.zeros((40, 40, 10), np.float32)
        d[8:32, 8:32, 2:8] = 100 + rng.random((24, 24, 6)) * 10
        A = np.eye(4)
        sp = np.array([1.0, 1.0, 4.0])
        A[:3, :3] = np.eye(3)[:, perm] * sp
        A[:3, 3] = -(A[:3, :3] @ (np.array(d.shape) - 1) / 2)
        p = tmp_path / f"{name}_haste.nii.gz"
        write_nifti(p, d, A, "test")
        files.append(("files", (p.name, p.read_bytes(), "application/octet-stream")))
    assert s.post(srv + "/api/cases/API-T-1/upload", files=files).ok
    assert s.post(srv + "/api/cases/API-T-1/import").ok
    for _ in range(120):
        j = s.get(srv + "/api/cases/API-T-1/job").json()
        if j["status"] in ("done", "error"):
            break
        time.sleep(0.5)
    assert j["status"] == "done", j
    case = s.get(srv + "/api/cases/API-T-1").json()
    assert len(case["stacks"]) == 3 and case["plan"]["mode"] in ("adult", "fetal")
    idx = s.get(srv + "/api/cases/API-T-1/grid").json()
    lid = idx["layers"][0]["id"]
    r = s.get(srv + f"/api/cases/API-T-1/layer?lid={lid}")
    assert r.ok and len(r.content) == int(np.prod(idx["grid"]["shape"]))


def test_no_result_no_download(srv):
    s = _login(srv)
    s.post(srv + "/api/cases", json={"case_id": "API-T-2"})
    assert s.get(srv + "/api/cases/API-T-2/download/nifti").status_code == 400
