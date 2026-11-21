"""Case workflow: plan (mode, reference stack, region, settings) -> reconstruction -> outputs.

Used by the local web app and by the command line.
"""
from __future__ import annotations

import csv
import io
import json
import time
from pathlib import Path

import numpy as np

from . import layers
from .config import VERSION
from .geometry import Grid, euler_deg_to_matrix, make_grid, nearest_ras_axes
from .ingest import suggest_reference
from .nifti import read_nifti, write_nifti
from .recon import FLAIR_SLICE_MSG, MODE_LABELS, MODES, ROI, Params, Reconstruction, default_roi, grid_for
from .store import CaseStore, StoreError, now_iso

EDITABLE = ("resolution", "motion", "cycles", "alpha", "regularizer", "outlier_rejection", "bias_correction",
            "ncc_threshold", "final_iterations", "psf_through", "mask_output")


# ------------------------------------------------------------------ loading
def load_stacks(store: CaseStore, cid: str, case: dict | None = None):
    case = case or store.load(cid)
    entries = case.get("stacks", [])
    stacks = [store.load_stack(cid, e) for e in entries]
    return case, entries, stacks


def included(entries, stacks):
    pairs = [(e, s) for e, s in zip(entries, stacks) if e.get("include")]
    return [p[0] for p in pairs], [p[1] for p in pairs]


# --------------------------------------------------------------------- plan
def make_plan(store: CaseStore, cid: str, mode: str | None = None) -> dict:
    case, entries, stacks = load_stacks(store, cid)
    if not entries:
        raise StoreError("import stacks first")
    prev = (case.get("plan") or {}).get("mode")
    mode = mode or prev or case.get("mode_hint") or "adult"
    if mode not in MODES:
        raise StoreError(f"mode must be one of {', '.join(MODES)}")
    if prev and mode != prev:
        _select_for_mode(entries, mode)
    inc_e, inc_s = included(entries, stacks)
    if not inc_e:
        raise StoreError("include at least one stack")
    ref_id = suggest_reference(inc_s, inc_e, mode)
    ref = [e["id"] for e in inc_e].index(ref_id)
    roi = default_roi(inc_s, ref, mode)
    p = Params.preset(mode)
    plan = {"mode": mode, "reference": ref_id, "params": {k: getattr(p, k) for k in EDITABLE},
            "roi": roi.to_dict(), "reference_suggested": ref_id, "roi_source": "automatic" if mode != "fetal"
            else ("brain mask" if inc_s[ref].mask is not None else "default - place it on the brain")}
    plan["warnings"] = plan_warnings(entries, mode)
    case["stacks"] = entries
    case["plan"] = plan
    case["status"] = "planned"
    store.save(cid, case)
    return plan


def update_plan(store: CaseStore, cid: str, changes: dict) -> dict:
    case = store.load(cid)
    plan = case.get("plan") or make_plan(store, cid)
    case = store.load(cid)
    if "mode" in changes and changes["mode"] != plan.get("mode"):
        plan = make_plan(store, cid, changes["mode"])
        case = store.load(cid)
    if "reference" in changes:
        ids = [e["id"] for e in case["stacks"] if e.get("include")]
        if changes["reference"] not in ids:
            raise StoreError("the reference must be an included stack")
        if changes["reference"] != plan["reference"]:
            plan["reference"] = changes["reference"]
            # the region and the viewer grids follow the reference stack's axes
            st = store.load_stack(cid, next(e for e in case["stacks"] if e["id"] == plan["reference"]))
            roi = ROI.from_dict(plan["roi"])
            roi.R = nearest_ras_axes(st.axes)
            plan["roi"] = roi.to_dict()
    for k, v in (changes.get("params") or {}).items():
        if k not in EDITABLE:
            raise StoreError(f"unknown setting {k}")
        v = _coerce(k, v)
        if k == "motion" and v == "slice" and plan.get("mode") == "flair":
            raise StoreError(FLAIR_SLICE_MSG)
        plan["params"][k] = v
    if "roi" in changes and changes["roi"]:
        r = changes["roi"]
        roi = ROI.from_dict(plan["roi"])
        if "center" in r:
            roi.center = np.array(r["center"], float)
        if "half" in r:
            roi.half = np.maximum(np.array(r["half"], float), 5.0)
        if "shape" in r:
            roi.shape = "ellipsoid" if r["shape"] == "ellipsoid" else "box"
        plan["roi"] = roi.to_dict()
        plan["roi_source"] = "user"
    plan["warnings"] = plan_warnings(case["stacks"], plan.get("mode", "adult"))
    case["plan"] = plan
    store.save(cid, case)
    return plan


def _select_for_mode(entries: list[dict], mode: str):
    """Switching to FLAIR selects the FLAIR stacks; switching back selects the single-shot ones
    (a reconstruction combines one kind of series)."""
    fl = [e for e in entries if e.get("flair_like")]
    ss = [e for e in entries if e.get("haste_like")]
    if mode == "flair" and len(fl) >= 1:
        for e in entries:
            e["include"] = bool(e.get("flair_like"))
    elif mode != "flair" and fl and ss:
        for e in entries:
            e["include"] = bool(e.get("haste_like"))


def plan_warnings(entries: list[dict], mode: str) -> list[str]:
    """Things that make the stacks unsuitable to combine, shown before the reconstruction runs."""
    inc = [e for e in entries if e.get("include")]
    out = []
    kinds = {e.get("sequence", "other") for e in inc}
    if "flair" in kinds and len(kinds) > 1:
        out.append("FLAIR and non-FLAIR series are both selected: one reconstruction must combine one contrast.")
    if mode == "flair" and inc and "flair" not in kinds:
        out.append("FLAIR mode, but no selected series looks like FLAIR (name or long inversion time).")
    if mode != "flair" and "flair" in kinds and len(kinds) == 1:
        out.append("The selected series look like FLAIR: choose the 2-D FLAIR mode (motion between series only).")
    for key, name, tol in (("te_ms", "echo time", 0.15), ("tr_ms", "repetition time", 0.25), ("ti_ms", "inversion time", 0.1)):
        vals = [float(e[key]) for e in inc if isinstance(e.get(key), (int, float)) and e[key] > 0]
        if len(vals) >= 2 and (max(vals) - min(vals)) > tol * max(vals):
            out.append(f"Selected series differ in {name} ({min(vals):.0f}-{max(vals):.0f} ms): contrast may differ "
                       "between planes.")
    if len(inc) < 3:
        out.append("Fewer than three stacks: through-plane detail cannot be recovered in every direction.")
    else:
        planes = {e.get("plane") for e in inc}
        if len(planes) < 3:
            out.append(f"Only {', '.join(sorted(p for p in planes if p))} planes selected: three orthogonal planes "
                       "are needed for isotropic detail.")
    return out


def update_stacks(store: CaseStore, cid: str, changes: list[dict]) -> dict:
    case = store.load(cid)
    by = {e["id"]: e for e in case["stacks"]}
    for ch in changes:
        e = by.get(ch.get("id"))
        if not e:
            continue
        if "include" in ch:
            e["include"] = bool(ch["include"])
        if "thickness_override" in ch:
            v = ch["thickness_override"]
            e["thickness_override"] = float(v) if v not in (None, "", 0) else None
    if not any(e.get("include") for e in case["stacks"]):
        raise StoreError("include at least one stack")
    plan = case.get("plan")
    if plan and not by.get(plan.get("reference"), {}).get("include"):
        plan["reference"] = next(e["id"] for e in case["stacks"] if e.get("include"))
    store.save(cid, case)
    return case


def _coerce(k, v):
    if k in ("cycles", "final_iterations"):
        return int(v)
    if k in ("outlier_rejection", "bias_correction", "mask_output"):
        return bool(v)
    if k == "alpha":
        return None if v in (None, "", "auto") else float(v)
    if k in ("motion", "regularizer"):
        allowed = {"motion": ("none", "stack", "slice"), "regularizer": ("huber", "tikhonov")}[k]
        if v not in allowed:
            raise StoreError(f"{k} must be one of {', '.join(allowed)}")
        return v
    return float(v)


def plan_objects(case: dict, entries, stacks):
    plan = case["plan"]
    inc = [bool(e.get("include")) for e in entries]
    ref = [e["id"] for e in entries].index(plan["reference"])
    p = Params.preset(plan["mode"], reference=ref, **{k: v for k, v in plan["params"].items() if k != "alpha"})
    p.alpha = plan["params"].get("alpha")
    p.reference = ref
    roi = ROI.from_dict(plan["roi"])
    return p, roi, inc


def planned_grid(store: CaseStore, cid: str) -> Grid:
    case, entries, stacks = load_stacks(store, cid)
    p, roi, inc = plan_objects(case, entries, stacks)
    return grid_for(roi, p, [s for s, i in zip(stacks, inc) if i])


# ---------------------------------------------------------------------- run
def run(store: CaseStore, cid: str, progress=None, cancel=None, log=None) -> dict:
    case, entries, stacks = load_stacks(store, cid)
    if not case.get("plan"):
        make_plan(store, cid)
        case, entries, stacks = load_stacks(store, cid)
    p, roi, inc = plan_objects(case, entries, stacks)

    def _log(m):
        store.log(cid, m)
        if log:
            log(m)
    _log(f"reconstruction started: mode {p.mode}, {sum(inc)} stacks, reference {entries[p.reference]['label']}, "
         f"{p.resolution} mm")
    t0 = time.time()
    rec = Reconstruction(stacks, roi, p, include=inc, log=_log, progress=progress, cancel=cancel)
    res = rec.run()
    d = store.dir(cid)
    out = d / "output"
    out.mkdir(exist_ok=True)
    desc = f"SliceWeave {VERSION} SRR {MODE_LABELS.get(p.mode, p.mode)} RESEARCH"
    write_nifti(out / "srr_native.nii.gz", res.volume, res.grid.affine, desc)
    write_nifti(out / "average_native.nii.gz", res.average, res.grid.affine, desc + " average")
    write_nifti(out / "srr.nii.gz", res.volume, res.grid.affine, desc)
    write_nifti(out / "average.nii.gz", res.average, res.grid.affine, desc + " average")
    report = res.report
    report.update({"version": VERSION, "case_id": cid, "finished": now_iso(), "simulated": case.get("simulated", False),
                   "stack_ids": [e["id"] for e, i in zip(entries, inc) if i],
                   "stack_labels": [e["label"] for e, i in zip(entries, inc) if i],
                   "wall_time_s": round(time.time() - t0, 1)})
    for s, lab in zip(report["stacks"], report["stack_labels"]):
        s["name"] = lab
    (out / "report.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    (out / "slices.csv").write_text(slices_csv(report), encoding="utf-8")
    np.savez_compressed(out / "transforms.npz", **{f"stack{i}": T for i, T in enumerate(res.transforms)})
    case = store.load(cid)
    case["result"] = {"grid": res.grid.to_dict(), "native_grid": res.grid.to_dict(), "finished": report["finished"],
                      "reorient": None, "summary": summary(report)}
    case["status"] = "reconstructed"
    store.save(cid, case)
    _clear_cache(store, cid, "res_")
    _log(f"reconstruction finished in {report['wall_time_s']} s")
    return case["result"]


def summary(report: dict) -> dict:
    st = report["stacks"]
    return {"stacks": len(st), "slices": sum(s["slices"] for s in st), "slices_used": sum(s["slices_used"] for s in st),
            "voxel_mm": report["grid"]["spacing"], "shape": report["grid"]["shape"],
            "alpha": report["alpha_used"], "seconds": report.get("wall_time_s")}


def slices_csv(report: dict) -> str:
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["stack", "slice", "included", "ncc", "intensity_scale", "rotation_deg_vs_stack_mean",
                "displacement_mm_vs_stack_mean", "pixels"])
    for s in report["stacks"]:
        for r in s["slices_detail"]:
            w.writerow([s["name"], r["slice"], int(r["included"]), r["ncc"], r["scale"], r["rotation_deg"],
                        r["displacement_mm"], r["pixels"]])
    return buf.getvalue()


def _clear_cache(store: CaseStore, cid: str, prefix: str):
    vd = store.dir(cid) / "work" / "viewer"
    if vd.exists():
        for f in vd.glob(prefix + "*"):
            f.unlink(missing_ok=True)


# ----------------------------------------------------------------- reorient
def reorient(store: CaseStore, cid: str, rx: float, ry: float, rz: float) -> dict:
    """Rotate the result (about the centre of the volume, in its own axes) and resample it with
    cubic splines, e.g. to show a fetal brain in standard planes. (0, 0, 0) restores the original."""
    case = store.load(cid)
    if not case.get("result"):
        raise StoreError("reconstruct first")
    d = store.dir(cid) / "output"
    native = Grid.from_dict(case["result"]["native_grid"])
    angles = [float(rx), float(ry), float(rz)]
    if not any(abs(a) > 1e-6 for a in angles):
        for name in ("srr", "average"):
            v = read_nifti(d / f"{name}_native.nii.gz")
            write_nifti(d / f"{name}.nii.gz", v.data, v.affine, f"SliceWeave {VERSION} SRR RESEARCH")
        case["result"]["grid"] = native.to_dict()
        case["result"]["reorient"] = None
    else:
        Rot = euler_deg_to_matrix(*angles)
        R_new = native.R @ Rot
        h = native.spacing
        half_native = (np.array(native.shape) - 1) / 2 * h
        # extents of the rotated box so nothing is cut
        corners = np.array([[sx, sy, sz] for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]) * half_native
        world = corners @ native.R.T
        half = np.abs(world @ R_new).max(axis=0)
        g = make_grid(R_new, native.center, half, h, fast=False)
        for name in ("srr", "average"):
            v = read_nifti(d / f"{name}_native.nii.gz")
            out = layers.resample(v.data, native, g, order=3)
            write_nifti(d / f"{name}.nii.gz", out, g.affine, f"SliceWeave {VERSION} SRR RESEARCH reoriented")
        case["result"]["grid"] = g.to_dict()
        case["result"]["reorient"] = {"rx": angles[0], "ry": angles[1], "rz": angles[2]}
    store.save(cid, case)
    _clear_cache(store, cid, "res_")
    store.log(cid, f"reorientation: {case['result']['reorient'] or 'none'}")
    return case["result"]
