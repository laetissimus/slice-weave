"""Case folders: create, load, save, list, delete; stacks saved without any header field."""
from __future__ import annotations

import datetime as _dt
import json
import os
import shutil
import threading
from pathlib import Path

import numpy as np

from .config import VERSION, Settings
from .deid import CASE_ID_RE
from .stacks import Stack

_lock = threading.RLock()


class StoreError(ValueError):
    pass


def now_iso() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


class CaseStore:
    def __init__(self, settings: Settings):
        self.s = settings.ensure()

    # ------------------------------------------------------------ paths
    def dir(self, case_id: str) -> Path:
        if not CASE_ID_RE.match(case_id or ""):
            raise StoreError("invalid case ID")
        p = (self.s.cases_dir / case_id).resolve()
        if p.parent != self.s.cases_dir.resolve():
            raise StoreError("invalid case path")
        return p

    def exists(self, case_id: str) -> bool:
        try:
            return (self.dir(case_id) / "case.json").exists()
        except StoreError:
            return False

    # ------------------------------------------------------------- CRUD
    def create(self, case_id: str, simulated: bool = False) -> dict:
        d = self.dir(case_id)
        with _lock:
            if (d / "case.json").exists():
                raise StoreError(f"case {case_id} already exists")
            for sub in ("stacks", "work", "output", "export"):
                (d / sub).mkdir(parents=True, exist_ok=True)
            try:
                os.chmod(d, 0o700)
            except OSError:
                pass
            case = {"case_id": case_id, "created": now_iso(), "updated": now_iso(), "version": VERSION,
                    "simulated": simulated, "status": "created", "stacks": [], "plan": None, "result": None}
            self.save(case_id, case)
        return case

    def load(self, case_id: str) -> dict:
        p = self.dir(case_id) / "case.json"
        if not p.exists():
            raise StoreError(f"no case {case_id}")
        return json.loads(p.read_text(encoding="utf-8"))

    def save(self, case_id: str, case: dict):
        case["updated"] = now_iso()
        p = self.dir(case_id) / "case.json"
        tmp = p.with_suffix(".json.part")
        with _lock:
            tmp.write_text(json.dumps(case, indent=1, default=_json_default), encoding="utf-8")
            tmp.replace(p)

    def update(self, case_id: str, **fields) -> dict:
        with _lock:
            c = self.load(case_id)
            c.update(fields)
            self.save(case_id, c)
            return c

    def list(self) -> list[dict]:
        out = []
        for d in sorted(self.s.cases_dir.glob("*/case.json")):
            try:
                c = json.loads(d.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            out.append({"case_id": c.get("case_id"), "status": c.get("status"), "simulated": c.get("simulated", False),
                        "created": c.get("created", ""), "age_days": _age_days(c.get("created", "")),
                        "n_stacks": len(c.get("stacks", [])), "mode": (c.get("plan") or {}).get("mode")})
        return out

    def delete(self, case_id: str):
        d = self.dir(case_id)
        with _lock:
            if d.exists():
                shutil.rmtree(d)

    def purge_older_than(self, days: int) -> list[str]:
        gone = []
        for c in self.list():
            if c["age_days"] is not None and c["age_days"] > days:
                self.delete(c["case_id"])
                gone.append(c["case_id"])
        return gone

    def log(self, case_id: str, msg: str):
        p = self.dir(case_id) / "log.txt"
        with _lock, open(p, "a", encoding="utf-8") as f:
            f.write(f"{_dt.datetime.now().strftime('%H:%M:%S')} {msg}\n")

    def read_log(self, case_id: str, tail: int = 400) -> list[str]:
        p = self.dir(case_id) / "log.txt"
        if not p.exists():
            return []
        return p.read_text(encoding="utf-8").splitlines()[-tail:]

    # ----------------------------------------------------------- stacks
    def save_stack(self, case_id: str, sid: str, st: Stack):
        p = self.dir(case_id) / "stacks" / f"{sid}.npz"
        arrays = {"data": st.data.astype(np.float32), "slice_aff": st.slice_aff, "thickness": np.float64(st.thickness)}
        if st.mask is not None:
            arrays["mask"] = st.mask.astype(np.uint8)
        if st.meta.get("acq_order") is not None:
            arrays["acq_order"] = np.asarray(st.meta["acq_order"], np.int32)
        tmp = p.with_name(p.name + ".part.npz")
        np.savez_compressed(tmp, **arrays)
        tmp.replace(p)

    def load_stack(self, case_id: str, entry: dict) -> Stack:
        p = self.dir(case_id) / "stacks" / f"{entry['id']}.npz"
        z = np.load(p)
        meta = {"id": entry["id"]}
        if "acq_order" in z:
            meta["acq_order"] = [int(v) for v in z["acq_order"]]
        mask = z["mask"].astype(bool) if "mask" in z else None
        th = float(entry.get("thickness_override") or z["thickness"])
        return Stack(entry.get("label") or entry["id"], z["data"], z["slice_aff"], th, mask, entry.get("plane", ""), meta)


def _age_days(iso: str):
    try:
        t = _dt.datetime.fromisoformat(iso)
        return (_dt.datetime.now(t.tzinfo) - t).days
    except (ValueError, TypeError):
        return None


def _json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o).__name__)


def clean_tmp(settings: Settings):
    """Remove anything left in the temp folder (e.g. after a crash)."""
    t = settings.tmp_dir
    if t.exists():
        for p in t.iterdir():
            try:
                shutil.rmtree(p) if p.is_dir() else p.unlink()
            except OSError:
                pass
