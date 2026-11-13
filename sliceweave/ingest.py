"""Import: DICOM folders / zips / NIfTI files -> de-identified stacks in a case folder.

Policy (see SECURITY.md): original files are read from a staging folder that is deleted
afterwards; only pixel data, slice geometry, slice thickness, the order in which slices were
acquired and a whitelist of technical parameters are kept. Series descriptions are scrubbed
of the patient's own name/ID tokens, dates and long digit runs.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from . import deid, series
from .nifti import read_nifti
from .stacks import Stack
from .store import CaseStore

FETAL_HINTS = re.compile(r"fet|foet|fetus|fetal|pregn|uter|obstet|placent|gestat", re.I)


def ingest(store: CaseStore, case_id: str, paths: list, tmpdir: Path, log=None) -> dict:
    log = log or (lambda m: None)
    found = series.scan(paths, tmpdir=tmpdir, log=log)
    problems = deid.validate_case_id(case_id, found, allow_digits=False)
    problems = [p for p in problems if "found in the images" in p]
    if problems:
        raise deid.DeidError(" ".join(problems) + " Delete this case and use a coded ID.")
    phi = set()
    for s in found:
        phi |= deid._phi_tokens(s.phi)
    case = store.load(case_id)
    entries = list(case.get("stacks", []))
    n0 = len(entries)
    rejected, masks, fetal_hint = [], [], False
    order = sorted(found, key=lambda s: (s.kind != "dicom", s.series_number or 0, s.uid))
    for s in order:
        desc = deid.scrub(s.text() if s.kind == "dicom" else _safe_name(s.description), phi)
        if FETAL_HINTS.search(s.text()):
            fetal_hint = True
        if s.reject:
            rejected.append({"description": desc, "reason": s.reject})
            continue
        if s.kind == "nifti" and series.is_mask_name(s.description):
            masks.append(s)
            continue
        if (s.acq_type or "").upper() == "3D":
            rejected.append({"description": desc, "reason": "3-D acquisition (not a stack of 2-D slices)"})
            continue
        try:
            st = series.build_stack(s, desc)
        except Exception as e:  # noqa: BLE001 - reported per series
            rejected.append({"description": desc, "reason": f"could not build the stack: {e}"})
            continue
        if st.nz < 5:
            rejected.append({"description": desc, "reason": f"only {st.nz} slices"})
            continue
        sid = f"s{len(entries) + 1:02d}"
        dx, dy = st.pixel
        side = st.meta.pop("sidecar_text", "")
        haste = series.looks_like_haste(s) or any(h in side.lower() for h in series.HASTE_HINTS)
        if FETAL_HINTS.search(side):
            fetal_hint = True
        entries.append({
            "id": sid, "label": f"{st.plane} {s.series_number or len(entries) + 1}", "description": desc,
            "plane": st.plane, "shape": [int(v) for v in st.shape], "pixel_mm": [round(dx, 3), round(dy, 3)],
            "thickness_mm": round(float(st.thickness), 3), "spacing_mm": round(st.spacing, 3),
            "thickness_source": "DICOM SliceThickness" if s.kind == "dicom" and s.thickness else
            (st.meta.get("thickness_from") or "slice spacing (no thickness in the file)"),
            "source": s.kind, "haste_like": bool(haste), "include": True, "has_mask": False,
            "acq_order_known": st.meta.get("acq_order") is not None,
            "series_number": s.series_number, "te_ms": s.te, "tr_ms": s.tr, "field_t": s.field_t,
            "warnings": list(s.warnings) + ([st.meta["warning"]] if st.meta.get("warning") else []),
        })
        st.meta["sid"] = sid
        store.save_stack(case_id, sid, st)
        log(f"imported {entries[-1]['label']}: {st.shape[0]}x{st.shape[1]}x{st.nz}, pixel {dx:.2f} mm, "
            f"thickness {st.thickness:.2f} mm" + (" (HASTE-like)" if haste else ""))
    new = entries[n0:]
    # include only HASTE-like 2-D series when some are recognised
    if any(e["haste_like"] for e in new):
        for e in new:
            e["include"] = e["haste_like"]
    # brain masks (NIfTI) made by another tool: matched to a stack with the same grid
    for m in masks:
        _attach_mask(store, case_id, entries, Path(m.nifti_path), log)
    case = store.load(case_id)
    case["stacks"] = entries
    case["status"] = "imported" if entries else "created"
    if fetal_hint:
        case["mode_hint"] = "fetal"
    store.save(case_id, case)
    summary = {"imported": len(new), "rejected": rejected, "masks": len(masks), "fetal_hint": fetal_hint,
               "series_found": len(found)}
    store.log(case_id, f"import: {len(new)} stack(s), {len(rejected)} series not used, {len(masks)} mask file(s)")
    return summary


def _safe_name(name: str) -> str:
    n = re.sub(r"\d{5,}", "[n]", name or "")
    return n[:48]


def _attach_mask(store: CaseStore, case_id: str, entries: list, path: Path, log) -> bool:
    v = read_nifti(path)
    probe = Stack.from_volume("mask", v.data, v.affine)
    for e in entries:
        st = store.load_stack(case_id, e)
        if st.shape != probe.shape:
            continue
        if np.allclose(st.slice_aff[:, :3, :], probe.slice_aff[:, :3, :], atol=0.05):
            st.mask = probe.data > 0.5
            store.save_stack(case_id, e["id"], st)
            e["has_mask"] = True
            log(f"brain mask attached to {e['label']}")
            return True
    log(f"mask file {_safe_name(path.name)}: no stack with the same grid; ignored")
    return False


def motion_score(st: Stack, fg: np.ndarray | None = None) -> float:
    """Median correlation of neighbouring slices in the central part of the field of view:
    lower when the subject moved between slices (interleaved acquisition)."""
    nx, ny, nz = st.shape
    sl = (slice(nx // 4, 3 * nx // 4), slice(ny // 4, 3 * ny // 4))
    vals = []
    for s in range(nz - 1):
        a, b = st.data[sl[0], sl[1], s].ravel(), st.data[sl[0], sl[1], s + 1].ravel()
        if fg is not None:
            m = (fg[sl[0], sl[1], s] | fg[sl[0], sl[1], s + 1]).ravel()
            if m.sum() < 50:
                continue
            a, b = a[m], b[m]
        if a.std() > 0 and b.std() > 0:
            vals.append(float(np.corrcoef(a, b)[0, 1]))
    return float(np.median(vals)) if vals else 0.0


def suggest_reference(stacks: list[Stack], entries: list[dict], mode: str) -> str | None:
    cand = [(e, st) for e, st in zip(entries, stacks) if e.get("include")]
    if not cand:
        return None
    scores = {e["id"]: motion_score(st, st.mask) for e, st in cand}
    best = max(scores, key=scores.get)
    if mode != "fetal":
        # adults: the axial stack (the usual reading plane) unless none is included; the neighbour
        # correlation is not comparable between planes of a still head
        ax = [e["id"] for e, _ in cand if e["plane"] == "AX"]
        if ax:
            best = max(ax, key=scores.get)
    for e, _ in cand:
        e["motion_score"] = round(scores[e["id"]], 3)
    return best
