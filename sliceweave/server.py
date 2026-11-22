"""Local web server (Starlette + uvicorn) for the SliceWeave viewer.

Security model (see SECURITY.md), shared with the GBM Infiltration Workbench:
* binds to the loopback interface only (127.0.0.1); there is no option to listen elsewhere;
* a random one-time login token is printed at start-up; opening the login link sets an
  HttpOnly, SameSite=Strict session cookie; every other request without it is refused;
* state-changing requests also need a CSRF header and, when the browser sends one, a matching
  Origin; the Host header must be 127.0.0.1/localhost on our port (DNS rebinding);
* strict Content-Security-Policy: the page loads nothing from the internet;
* uploads are staged in the data folder's temp area and deleted after import;
* the network guard (netguard.py) refuses every non-loopback connection of this process.
"""
from __future__ import annotations

import hmac
import re
import secrets
import shutil
import tempfile
import threading
import time
import traceback
from pathlib import Path

from starlette.applications import Starlette
from starlette.concurrency import run_in_threadpool
from starlette.requests import Request
from starlette.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, Response
from starlette.routing import Mount, Route
from starlette.staticfiles import StaticFiles

from . import export, hooks, ingest, layers, netguard, security, workflow
from .config import NAME, RESEARCH_BANNER, VERSION, Settings
from .deid import DeidError, validate_case_id
from .geometry import Grid
from .nifti import read_nifti
from .recon import STEPS
from .store import CaseStore, StoreError, clean_tmp

WEB = Path(__file__).parent / "web"
CSP = ("default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data: blob:; "
       "connect-src 'self'; font-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; "
       "form-action 'self'")
UNSAFE = {"POST", "PUT", "DELETE", "PATCH"}
CHUNK = 1 << 20


class Jobs:
    """One background worker: heavy work runs one job at a time."""

    def __init__(self):
        self.state: dict[str, dict] = {}
        self.lock = threading.Lock()
        self.queue: list = []
        self.cv = threading.Condition(self.lock)
        threading.Thread(target=self._loop, daemon=True).start()

    def submit(self, key: str, kind: str, fn, *args, **kw):
        with self.lock:
            if self.state.get(key, {}).get("status") in ("queued", "running"):
                raise RuntimeError("a job is already queued or running for this case")
            st = {"status": "queued", "kind": kind, "started": None, "error": None, "cancel": False,
                  "steps": {k: {"status": "pending", "msg": ""} for k, _ in STEPS} if kind == "run" else {}}
            self.state[key] = st
            self.queue.append((key, fn, args, kw))
            self.cv.notify()

    def progress_fn(self, key):
        def progress(step, frac, msg=""):
            st = self.state.get(key)
            if not st:
                return
            steps = st["steps"]
            seen = False
            for k, _ in STEPS:  # earlier steps are done once a later one starts
                if k == step:
                    seen = True
                    steps[k]["status"] = "done" if frac is not None and frac >= 1.0 else "running"
                    steps[k]["msg"] = msg or steps[k]["msg"]
                elif not seen and steps[k]["status"] in ("pending", "running"):
                    steps[k]["status"] = "done"
        return progress

    def cancel(self, key):
        st = self.state.get(key)
        if st and st["status"] in ("queued", "running"):
            st["cancel"] = True

    def _loop(self):
        while True:
            with self.cv:
                while not self.queue:
                    self.cv.wait()
                key, fn, args, kw = self.queue.pop(0)
            st = self.state[key]
            st.update(status="running", started=time.time())
            try:
                st["result"] = fn(*args, **kw)
                st["status"] = "done"
                for v in st["steps"].values():
                    v["status"] = "done"
            except Exception as e:  # noqa: BLE001 - surfaced to the UI
                st.update(status="cancelled" if str(e) == "cancelled" else "error", error=str(e),
                          trace=traceback.format_exc()[-2000:])
                for v in st["steps"].values():
                    if v["status"] == "running":
                        v["status"] = "error"
            st["finished"] = time.time()


def create_app(settings: Settings, token: str, port: int) -> Starlette:
    store = CaseStore(settings)
    clean_tmp(settings)
    spool = settings.tmp_dir / "spool"
    spool.mkdir(parents=True, exist_ok=True)
    tempfile.tempdir = str(spool)  # multipart spooling stays inside the data folder
    jobs = Jobs()
    session_secret = secrets.token_urlsafe(32)
    csrf_secret = secrets.token_urlsafe(32)
    allowed_hosts = {f"127.0.0.1:{port}", f"localhost:{port}"}
    allowed_origins = {f"http://127.0.0.1:{port}", f"http://localhost:{port}"}
    staging_root = settings.tmp_dir / "staging"

    class Guard:
        def __init__(self, app):
            self.app = app

        async def __call__(self, scope, receive, send):
            if scope["type"] != "http":
                return await self.app(scope, receive, send)
            req = Request(scope)
            if req.headers.get("host", "") not in allowed_hosts:
                return await JSONResponse({"error": "bad host"}, 403)(scope, receive, send)
            path = scope["path"]
            if path != "/login":
                if not hmac.compare_digest(req.cookies.get("sw_session", ""), session_secret):
                    if path.startswith("/api/"):
                        return await JSONResponse({"error": "not logged in"}, 401)(scope, receive, send)
                    return await HTMLResponse(LOGIN_HELP, 401)(scope, receive, send)
                if scope["method"] in UNSAFE:
                    if not hmac.compare_digest(req.headers.get("x-sw-csrf", ""), csrf_secret):
                        return await JSONResponse({"error": "missing CSRF token"}, 403)(scope, receive, send)
                    origin = req.headers.get("origin")
                    if origin and origin not in allowed_origins:
                        return await JSONResponse({"error": "bad origin"}, 403)(scope, receive, send)
                    sfs = req.headers.get("sec-fetch-site")
                    if sfs and sfs not in ("same-origin", "none"):
                        return await JSONResponse({"error": "cross-site request"}, 403)(scope, receive, send)

            async def send_wrapped(msg):
                if msg["type"] == "http.response.start":
                    h = list(msg.get("headers", []))
                    h += [(b"content-security-policy", CSP.encode()), (b"x-content-type-options", b"nosniff"),
                          (b"referrer-policy", b"no-referrer"), (b"x-frame-options", b"DENY"),
                          (b"cross-origin-opener-policy", b"same-origin"),
                          (b"cross-origin-resource-policy", b"same-origin")]
                    if path.startswith("/api/"):
                        h.append((b"cache-control", b"no-store"))
                    msg["headers"] = h
                await send(msg)
            return await self.app(scope, receive, send_wrapped)

    def ok(obj=None, code=200):
        return JSONResponse(obj if obj is not None else {"ok": True}, code)

    def err(msg, code=400):
        return JSONResponse({"error": msg}, code)

    def cid_of(req) -> str:
        cid = req.path_params["cid"]
        store.dir(cid)
        if not store.exists(cid):
            raise StoreError("no such case")
        return cid

    async def body_json(req):
        try:
            return await req.json()
        except Exception:  # noqa: BLE001
            return {}

    def job_view(cid):
        st = jobs.state.get(cid)
        v = ({"status": st["status"], "kind": st["kind"], "error": st.get("error"), "steps": st["steps"],
              "elapsed": round(time.time() - st["started"], 1) if st.get("started") else 0}
             if st else {"status": "idle"})
        v["log"] = store.read_log(cid, 150) if store.exists(cid) else []
        return v

    # ------------------------------------------------------------ routes
    async def login(req: Request):
        if not hmac.compare_digest(req.query_params.get("token", ""), token):
            return HTMLResponse(LOGIN_HELP, 401)
        r = RedirectResponse("/", 303)
        r.set_cookie("sw_session", session_secret, httponly=True, samesite="strict", path="/")
        return r

    async def index(req):
        return FileResponse(WEB / "index.html", headers={"cache-control": "no-store"})

    async def session(req):
        over = [c["case_id"] for c in store.list() if c["age_days"] is not None and c["age_days"] > settings.retention_days]
        return ok({"csrf": csrf_secret, "version": VERSION, "name": NAME, "banner": RESEARCH_BANNER,
                   "data_dir": str(settings.root), "retention_days": settings.retention_days, "over_retention": over,
                   "encryption": security.disk_encryption_status(settings.root), "network_guard": netguard.status(),
                   "hooks": hooks.status(), "steps": STEPS})

    async def cases(req):
        if req.method == "GET":
            return ok({"cases": store.list()})
        b = await body_json(req)
        cid = (b.get("case_id") or "").strip()
        problems = validate_case_id(cid, None, bool(b.get("allow_digits")))
        if problems:
            return err(" ".join(problems))
        if store.exists(cid):
            return err("case already exists")
        store.create(cid)
        return ok({"case_id": cid})

    async def case(req):
        cid = cid_of(req)
        if req.method == "DELETE":
            store.delete(cid)
            shutil.rmtree(staging_root / cid, ignore_errors=True)
            return ok()
        c = store.load(cid)
        c["staged"] = _staged(cid)
        c["job"] = job_view(cid)
        return ok(c)

    def _staged(cid):
        d = staging_root / cid
        files = [f for f in d.rglob("*") if f.is_file()] if d.exists() else []
        return {"files": len(files), "bytes": sum(f.stat().st_size for f in files)}

    async def upload(req: Request):
        cid = cid_of(req)
        dest = staging_root / cid
        dest.mkdir(parents=True, exist_ok=True)
        limit = settings.max_upload_gb * 1024 ** 3
        used = sum(f.stat().st_size for f in dest.rglob("*") if f.is_file())
        n = 0
        async with req.form(max_files=100000, max_fields=1000) as form:
            for item in form.getlist("files"):
                if not hasattr(item, "read"):
                    continue
                name = Path((item.filename or "").replace("\\", "/")).name
                low = name.lower()
                if low.endswith((".nii", ".nii.gz", ".json")):
                    # keep a sanitised base name so that a JSON sidecar pairs with its NIfTI and
                    # "mask" files are recognised; the staging folder is deleted after import
                    safe = re.sub(r"\d{5,}", "n", re.sub(r"[^A-Za-z0-9._-]", "_", name))[:80]
                    target = dest / safe
                    k = 1
                    while target.exists():
                        target = dest / f"u{k}" / safe
                        k += 1
                    target.parent.mkdir(parents=True, exist_ok=True)
                elif low.endswith(".zip"):
                    target = dest / f"z{n:05d}_{secrets.token_hex(3)}.zip"
                else:
                    target = dest / f"f{secrets.token_hex(6)}.dcm"
                with open(target, "wb") as out:
                    while True:
                        chunk = await item.read(CHUNK)
                        if not chunk:
                            break
                        used += len(chunk)
                        if used > limit:
                            out.close()
                            shutil.rmtree(dest, ignore_errors=True)
                            return err(f"upload exceeds {settings.max_upload_gb} GB", 413)
                        out.write(chunk)
                n += 1
        return ok({"received": n, "staged": _staged(cid)})

    async def clear_staging(req):
        shutil.rmtree(staging_root / cid_of(req), ignore_errors=True)
        return ok()

    async def do_import(req):
        cid = cid_of(req)
        src = staging_root / cid
        if not src.exists() or not any(src.rglob("*")):
            return err("nothing staged: add DICOM folders, zips or NIfTI files first")

        def work():
            tmp = settings.tmp_dir / f"unzip_{cid}"
            try:
                out = ingest.ingest(store, cid, [src], tmp, log=lambda m: store.log(cid, m))
                if store.load(cid)["stacks"]:
                    workflow.make_plan(store, cid)
                    _clear_viewer(cid)
                return out
            finally:
                shutil.rmtree(src, ignore_errors=True)
                shutil.rmtree(tmp, ignore_errors=True)
        try:
            jobs.submit(cid, "import", work)
        except RuntimeError as e:
            return err(str(e), 409)
        return ok({"queued": True})

    def _clear_viewer(cid, prefix=""):
        vd = store.dir(cid) / "work" / "viewer"
        if vd.exists():
            for f in vd.glob(prefix + "*"):
                f.unlink(missing_ok=True)

    async def stacks(req):
        cid = cid_of(req)
        b = await body_json(req)
        c = await run_in_threadpool(workflow.update_stacks, store, cid, b.get("changes") or [])
        return ok({"stacks": c["stacks"], "plan": c.get("plan")})

    async def plan(req):
        cid = cid_of(req)
        b = await body_json(req)
        if req.method == "POST":
            p = await run_in_threadpool(workflow.make_plan, store, cid, b.get("mode"))
        else:
            p = await run_in_threadpool(workflow.update_plan, store, cid, b)
        g = await run_in_threadpool(workflow.planned_grid, store, cid)
        return ok({"plan": p, "grid": g.to_dict(), "voxels_m": round(g.size / 1e6, 2)})

    def _grid(cid, which) -> Grid:
        c = store.load(cid)
        if which == "result":
            if not c.get("result"):
                raise StoreError("no reconstruction yet")
            return Grid.from_dict(c["result"]["grid"])
        if not c.get("stacks"):
            raise StoreError("import stacks first")
        ref_id = (c.get("plan") or {}).get("reference") or c["stacks"][0]["id"]
        g = c.get("overview_grid") or {}
        if g.get("ref") != ref_id:
            _, entries, sts = workflow.load_stacks(store, cid, c)
            ref = sts[[e["id"] for e in entries].index(ref_id)]
            grid = layers.overview_grid(sts, ref)
            c = store.load(cid)
            c["overview_grid"] = {"ref": ref_id, "grid": grid.to_dict()}
            store.save(cid, c)
            return grid
        return Grid.from_dict(g["grid"])

    async def grid_index(req):
        cid = cid_of(req)
        which = req.query_params.get("which", "overview")
        g = await run_in_threadpool(_grid, cid, which)
        c = store.load(cid)
        lay = []
        if which == "result":
            lay += [{"id": "srr", "title": "Super-resolution", "group": "Result"},
                    {"id": "average", "title": "Average (no SRR)", "group": "Result"}]
        for e in c.get("stacks", []):
            if which == "result" and not e.get("include"):
                continue
            lay.append({"id": e["id"], "title": f"{e['label']}" + ("" if e.get("include") else " (not used)"),
                        "group": "Input stacks (as acquired)", "plane": e["plane"]})
        import numpy as np
        return ok({"which": which, "grid": g.to_dict(), "inv": np.round(g.inv, 8).tolist(), "layers": lay,
                   "key": g.key()})

    def _layer_bytes(cid, which, lid):
        g = _grid(cid, which)
        vd = store.dir(cid) / "work" / "viewer"
        key = f"{'res' if which == 'result' else 'ov'}_{g.key()}_{lid}"
        hit = layers.read_layer(vd, key)
        if hit:
            return hit
        c = store.load(cid)
        if lid in ("srr", "average"):
            v = read_nifti(store.dir(cid) / "output" / f"{lid}.nii.gz")
            vol = v.data if v.shape == tuple(g.shape) else layers.resample(v.data, Grid(v.shape, v.affine), g)
        else:
            e = next((x for x in c.get("stacks", []) if x["id"] == lid), None)
            if e is None:
                raise StoreError("no such layer")
            vol = layers.stack_on_grid(store.load_stack(cid, e), g)
        layers.write_layer(vd, key, vol)
        return layers.read_layer(vd, key)

    async def layer(req):
        cid = cid_of(req)
        which = req.query_params.get("which", "overview")
        lid = req.query_params.get("lid", "")
        if not re.fullmatch(r"[a-z0-9_]{1,16}", lid):
            return err("bad layer id")
        data, meta = await run_in_threadpool(_layer_bytes, cid, which, lid)
        return Response(data, media_type="application/octet-stream",
                        headers={"x-sw-offset": str(meta["offset"]), "x-sw-scale": str(meta["scale"]),
                                 "x-sw-window": f"{meta['window'][0]},{meta['window'][1]}",
                                 "x-sw-window-brain": ",".join(str(v) for v in meta.get("window_brain", meta["window"]))})

    async def run_case(req):
        cid = cid_of(req)
        c = store.load(cid)
        if not c.get("stacks"):
            return err("import stacks first")
        try:
            jobs.submit(cid, "run", workflow.run, store, cid, progress=jobs.progress_fn(cid),
                        cancel=lambda: jobs.state.get(cid, {}).get("cancel", False))
        except RuntimeError as e:
            return err(str(e), 409)
        return ok({"queued": True})

    async def cancel(req):
        jobs.cancel(cid_of(req))
        return ok()

    async def job(req):
        return ok(job_view(cid_of(req)))

    async def report(req):
        return ok(await run_in_threadpool(export.report_json, store, cid_of(req)))

    async def reorient(req):
        cid = cid_of(req)
        b = await body_json(req)
        r = await run_in_threadpool(workflow.reorient, store, cid, b.get("rx", 0), b.get("ry", 0), b.get("rz", 0))
        return ok({"result": r})

    async def download(req):
        cid = cid_of(req)
        kind = req.path_params["kind"]
        c = store.load(cid)
        d = store.dir(cid) / "output"
        if kind in ("nifti", "average"):
            p = export.output_paths(store, cid)["srr" if kind == "nifti" else "average"]
            name = export.nifti_name(c) if kind == "nifti" else export.nifti_name(c).replace("_sliceweave_", "_average_")
            return FileResponse(p, filename=name, media_type="application/gzip")
        if kind == "report":
            return FileResponse(d / "report.json", filename=f"{cid}_report.json", media_type="application/json")
        if kind == "csv":
            return FileResponse(d / "slices.csv", filename=f"{cid}_slices.csv", media_type="text/csv")
        if kind == "bundle":
            data = await run_in_threadpool(export.bundle_zip, store, cid)
            return _zip(data, f"{cid}_sliceweave_bundle.zip")
        if kind == "dicom":
            which = req.query_params.get("which", "srr")
            data = await run_in_threadpool(export.dicom_zip, store, cid, which if which in ("srr", "average") else "srr")
            return _zip(data, f"{cid}_sliceweave_dicom.zip")
        return err("unknown export")

    async def dicom_attached(req):
        """DICOM export with the original study's patient/study attributes, read from one DICOM
        file chosen now. The file and the attributes are not stored; the zip is streamed back."""
        cid = cid_of(req)
        raw = None
        async with req.form(max_files=2, max_fields=10) as form:
            f = form.get("reference")
            if f is not None and hasattr(f, "read"):
                raw = await f.read()
            which = form.get("which") or "srr"
        if not raw:
            return err("choose one DICOM file of the original study")
        try:
            attrs = export.reference_attributes(raw)
        except Exception as e:  # noqa: BLE001
            return err(f"not a usable DICOM file: {e}")
        data = await run_in_threadpool(export.dicom_zip, store, cid, which if which in ("srr", "average") else "srr", attrs)
        store.log(cid, "DICOM export attached to the original study (identifiers not stored)")
        return _zip(data, f"{cid}_sliceweave_dicom_study.zip")

    def _zip(data: bytes, name: str):
        return Response(data, media_type="application/zip",
                        headers={"content-disposition": f'attachment; filename="{name}"'})

    async def demo(req):
        b = await body_json(req)
        kind = b.get("kind") if b.get("kind") in ("fetal", "flair") else "adult"
        cid = f"DEMO-{kind.upper()}-{secrets.token_hex(2).upper()}"
        store.create(cid, simulated=True)

        def work():
            from .selftest import make_demo
            make_demo(store, cid, kind, log=lambda m: store.log(cid, m))
            workflow.make_plan(store, cid, kind)
            return {"case_id": cid}
        jobs.submit(cid, "import", work)
        return ok({"case_id": cid})

    async def purge(req):
        gone = store.purge_older_than(settings.retention_days)
        return ok({"deleted": gone})

    routes = [
        Route("/login", login), Route("/", index), Route("/api/session", session),
        Route("/api/cases", cases, methods=["GET", "POST"]),
        Route("/api/demo", demo, methods=["POST"]),
        Route("/api/cases/{cid}", case, methods=["GET", "DELETE"]),
        Route("/api/cases/{cid}/upload", upload, methods=["POST"]),
        Route("/api/cases/{cid}/staging", clear_staging, methods=["DELETE"]),
        Route("/api/cases/{cid}/import", do_import, methods=["POST"]),
        Route("/api/cases/{cid}/stacks", stacks, methods=["PUT"]),
        Route("/api/cases/{cid}/plan", plan, methods=["POST", "PUT"]),
        Route("/api/cases/{cid}/grid", grid_index),
        Route("/api/cases/{cid}/layer", layer),
        Route("/api/cases/{cid}/run", run_case, methods=["POST"]),
        Route("/api/cases/{cid}/cancel", cancel, methods=["POST"]),
        Route("/api/cases/{cid}/job", job),
        Route("/api/cases/{cid}/report", report),
        Route("/api/cases/{cid}/reorient", reorient, methods=["POST"]),
        Route("/api/cases/{cid}/download/{kind}", download),
        Route("/api/cases/{cid}/dicom-attached", dicom_attached, methods=["POST"]),
        Route("/api/retention/purge", purge, methods=["POST"]),
        Mount("/static", StaticFiles(directory=str(WEB)), name="static"),
    ]

    async def on_error(req, exc):
        return JSONResponse({"error": str(exc)}, 400)

    app = Starlette(routes=routes, exception_handlers={StoreError: on_error, ValueError: on_error,
                                                       DeidError: on_error})
    app.add_middleware(Guard)
    return app


LOGIN_HELP = """<!doctype html><meta charset="utf-8"><title>SliceWeave</title>
<body style="font-family:system-ui;margin:3em;max-width:40em;background:#0f1115;color:#e7e9ee">
<h2>SliceWeave</h2>
<p>Not logged in. Open the login link printed in the window that started SliceWeave
(it looks like <code>http://127.0.0.1:8766/login?token=...</code>).</p></body>"""


def serve(settings: Settings, port: int | None = None, open_browser: bool = True):
    import webbrowser

    import uvicorn
    netguard.install()
    port = port or settings.port
    token = secrets.token_urlsafe(24)
    app = create_app(settings, token, port)
    url = f"http://127.0.0.1:{port}/login?token={token}"
    print(f"\n  {NAME} {VERSION} - 3-plane HASTE to isotropic volume", flush=True)
    print("  " + RESEARCH_BANNER)
    print(f"  Data folder: {settings.root}")
    enc = security.disk_encryption_status(settings.root)
    print(f"  Disk encryption: {enc['status']} ({enc['detail']})")
    print("  Network: offline. Only this computer (127.0.0.1) can be contacted; anything else is blocked.")
    print(f"\n  Open this link (it contains your one-time login token):\n  {url}\n")
    print("  Keep this window open while you work. Close it (or press Ctrl+C) to stop.\n", flush=True)
    if open_browser:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    try:
        uvicorn.run(app, host="127.0.0.1", port=port, log_level="warning", access_log=False,
                    server_header=False, date_header=False)
    except KeyboardInterrupt:
        pass
    print(f"  {NAME} stopped.", flush=True)
