"""Exports: NIfTI, DICOM (de-identified, or attached to the original study) and a research bundle."""
from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path

import numpy as np

from . import dicom
from .config import RESEARCH_BANNER, VERSION
from .dicom_write import new_uid, write_series
from .nifti import read_nifti
from .store import CaseStore, StoreError

# attributes copied from a user-supplied DICOM file of the original study ("attach to study")
STUDY_KEYS = ["PatientName", "PatientID", "PatientBirthDate", "PatientSex", "PatientAge", "OtherPatientIDs",
              "StudyInstanceUID", "StudyDate", "StudyTime", "AccessionNumber", "StudyID", "StudyDescription",
              "ReferringPhysicianName", "FrameOfReferenceUID", "InstitutionName", "Manufacturer"]
EXTRA_TAGS = {"SoftwareVersions": (0x00181020, "LO"), "DerivationDescription": (0x00082111, "ST"),
              "SeriesDate": (0x00080021, "DA"), "PatientPosition": (0x00185100, "CS")}


def output_paths(store: CaseStore, cid: str) -> dict:
    d = store.dir(cid) / "output"
    p = {k: d / f"{k}.nii.gz" for k in ("srr", "average")}
    if not p["srr"].exists():
        raise StoreError("no reconstruction yet")
    return p


def nifti_name(case: dict) -> str:
    m = (case.get("plan") or {}).get("mode", "srr")
    vox = (case.get("result") or {}).get("grid", {}).get("spacing", 0)
    return f"{case['case_id']}_sliceweave_{m}_{vox:.2f}".replace(".", "p") + "mm.nii.gz"


def reference_attributes(raw: bytes) -> dict:
    """Patient/study attributes of a DICOM file of the original study (read in memory only)."""
    ds = dicom.read_bytes(raw, stop_before_pixels=True)
    out = {}
    for k in STUDY_KEYS:
        v = getattr(ds, k)
        if v not in (None, "", b""):
            out[k] = v
    if "StudyInstanceUID" not in out:
        raise StoreError("the file has no StudyInstanceUID: choose an image file of the original study")
    return out


def dicom_zip(store: CaseStore, cid: str, which: str = "srr", attach: dict | None = None) -> bytes:
    case = store.load(cid)
    p = output_paths(store, cid)[which]
    v = read_nifti(p)
    data = np.clip(v.data, 0, None)
    hi = float(np.percentile(data[data > 0], 99.9)) if (data > 0).any() else 1.0
    slope = max(hi / 4000.0, 1e-6)  # 12-bit-like range, unsigned
    mode = (case.get("plan") or {}).get("mode", "")
    vox = float(np.linalg.norm(v.affine[:3, 0]))
    label = "SRR" if which == "srr" else "AVG"
    hdr = {
        "SeriesDescription": f"SliceWeave {label} {mode} {vox:.2f}mm RESEARCH",
        "SeriesNumber": 9100 if which == "srr" else 9101,
        "ImageType": ["DERIVED", "SECONDARY", "SRR" if which == "srr" else "AVERAGE"],
        "Manufacturer": "SliceWeave (research software)",
        "ImageComments": f"{RESEARCH_BANNER}. Super-resolution reconstruction of 2-D stacks, SliceWeave {VERSION}.",
        "BodyPartExamined": "BRAIN", "BurnedInAnnotation": "NO",
        "WindowCenter": hi / 2, "WindowWidth": hi,
        EXTRA_TAGS["SoftwareVersions"]: f"SliceWeave {VERSION}",
        EXTRA_TAGS["DerivationDescription"]: "Motion-corrected super-resolution reconstruction from 2-D stacks",
    }
    if attach:
        hdr.update({k: v for k, v in attach.items() if k in STUDY_KEYS})
        hdr["SeriesInstanceUID"] = new_uid()
        if "FrameOfReferenceUID" not in attach:
            hdr["FrameOfReferenceUID"] = new_uid()
    else:
        uid = case.get("export_study_uid")
        if not uid:
            uid = new_uid()
            case["export_study_uid"] = uid
            store.save(cid, case)
        hdr.update({"PatientName": "", "PatientID": cid, "StudyInstanceUID": uid,
                    "StudyDescription": "SliceWeave research export", "StudyID": "SW"})
    buf = io.BytesIO()
    tmp = store.s.tmp_dir / f"dcm_{new_uid()[-12:]}"
    try:
        paths = write_series(data, v.affine, tmp, hdr, slope=slope, intercept=0.0)
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for f in paths:
                z.write(f, f"{label}/{f.name}")
            z.writestr("README.txt", _readme(case, attach is not None))
    finally:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)
    return buf.getvalue()


def bundle_zip(store: CaseStore, cid: str) -> bytes:
    case = store.load(cid)
    d = store.dir(cid) / "output"
    output_paths(store, cid)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name in ("srr.nii.gz", "average.nii.gz", "report.json", "slices.csv"):
            f = d / name
            if f.exists():
                z.write(f, f"{cid}/{name}")
        z.writestr(f"{cid}/README.txt", _readme(case, False))
    return buf.getvalue()


def _readme(case: dict, identified: bool) -> str:
    r = case.get("result") or {}
    plan = case.get("plan") or {}
    lines = [f"SliceWeave {VERSION} export - {RESEARCH_BANNER}", "",
             f"Case: {case['case_id']}" + ("  (SIMULATED DATA)" if case.get("simulated") else ""),
             f"Mode: {plan.get('mode')}   voxel: {r.get('grid', {}).get('spacing')} mm   grid: {r.get('grid', {}).get('shape')}",
             f"Reoriented: {r.get('reorient') or 'no'}", "",
             "srr: super-resolution reconstruction (motion-corrected, PSF deconvolution, edge-preserving regularisation).",
             "average: PSF-weighted mean of the same slices without deconvolution (baseline for comparison).",
             "report.json / slices.csv: settings, stack and slice registration, excluded slices, timing.", ""]
    if identified:
        lines += ["This DICOM series carries the patient and study attributes of the original study file chosen at",
                  "export time, so PACS files it with that study. Handle it as identifiable data."]
    else:
        lines += ["De-identified: PatientID is the coded case ID; no patient attributes are included."]
    lines += ["", "The volume is computed from the acquired slices and can contain reconstruction artefacts;",
              "always compare with the original 2-D series."]
    return "\n".join(lines) + "\n"


def report_json(store: CaseStore, cid: str) -> dict:
    p = store.dir(cid) / "output" / "report.json"
    if not p.exists():
        raise StoreError("no reconstruction yet")
    return json.loads(p.read_text(encoding="utf-8"))
