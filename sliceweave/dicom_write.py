"""Minimal DICOM writer (explicit VR little endian).

Used for (1) the research-labeled derived-series export and (2) building synthetic test
series. A de-identified export carries the coded case ID as PatientID and an empty
PatientName; an export "attached to the original study" copies the patient and study
attributes of a DICOM file the user supplies at export time (nothing is stored).
"""
from __future__ import annotations

import struct
import uuid
from pathlib import Path

import numpy as np

from .dicom import TAGS

MR_IMAGE_STORAGE = "1.2.840.10008.5.1.4.1.1.4"
EXPLICIT_LE = "1.2.840.10008.1.2.1"
IMPL_UID = "2.25.98437751235518392275810385470152239617"  # UUID-derived implementation UID
IMPL_NAME = "SLICEWEAVE_1"


def new_uid() -> str:
    return "2.25." + str(uuid.uuid4().int)


def _enc_value(vr: str, v) -> bytes:
    if vr in ("US", "SS", "UL", "SL", "FL", "FD"):
        fmt = {"US": "H", "SS": "h", "UL": "I", "SL": "i", "FL": "f", "FD": "d"}[vr]
        vals = v if isinstance(v, (list, tuple)) else [v]
        return struct.pack("<" + fmt * len(vals), *vals)
    if vr in ("OB", "OW", "UN"):
        b = bytes(v)
        return b + (b"\x00" if len(b) % 2 else b"")
    if isinstance(v, (list, tuple)):
        parts = []
        for x in v:
            parts.append(_num(x, vr))
        s = "\\".join(parts)
    else:
        s = _num(v, vr)
    b = s.encode("latin-1", "replace")
    if len(b) % 2:
        b += b"\x00" if vr == "UI" else b" "
    return b


def _num(x, vr):
    if vr == "DS" and isinstance(x, (int, float, np.floating, np.integer)):
        s = f"{float(x):.10g}"
        if len(s) > 16:
            s = f"{float(x):.8e}" if abs(float(x)) > 1e8 else f"{float(x):.10f}"[:16]
        return s
    if vr == "IS" and isinstance(x, (int, float, np.integer, np.floating)):
        return str(int(x))
    return str(x)


def _element(tag: int, vr: str, value) -> bytes:
    data = _enc_value(vr, value)
    g, e = tag >> 16, tag & 0xFFFF
    if vr in ("OB", "OD", "OF", "OL", "OW", "SQ", "UC", "UN", "UR", "UT"):
        return struct.pack("<HH2sHI", g, e, vr.encode(), 0, len(data)) + data
    return struct.pack("<HH2sH", g, e, vr.encode(), len(data)) + data


def write_file(path: Path, elements: dict, sop_class: str, sop_uid: str):
    """``elements``: {keyword or tag: value}; VRs come from :data:`gbm_wb.dicom.TAGS`."""
    items = []
    for k, v in elements.items():
        if v is None:
            continue
        if isinstance(k, tuple):
            tag, vr = k
        else:
            tag, vr = TAGS[k]
        if k == "PixelData":
            vr = "OW"
        items.append((tag, vr, v))
    items.append((0x00080016, "UI", sop_class))
    items.append((0x00080018, "UI", sop_uid))
    items.sort(key=lambda t: t[0])
    body = b"".join(_element(t, vr, v) for t, vr, v in items)
    meta_elems = [
        (0x00020001, "OB", b"\x00\x01"),
        (0x00020002, "UI", sop_class),
        (0x00020003, "UI", sop_uid),
        (0x00020010, "UI", EXPLICIT_LE),
        (0x00020012, "UI", IMPL_UID),
        (0x00020013, "SH", IMPL_NAME),
    ]
    meta = b"".join(_element(t, vr, v) for t, vr, v in meta_elems)
    meta = _element(0x00020000, "UL", len(meta)) + meta
    Path(path).write_bytes(b"\x00" * 128 + b"DICM" + meta + body)


def write_series(vol_data: np.ndarray, affine_ras: np.ndarray, out_dir: Path, header: dict,
                 slope: float = 1.0, intercept: float = 0.0, signed: bool = False, prefix: str = "IM") -> list[Path]:
    """Write a 3-D array (i, j, k) as a single-frame MR series, one file per k.

    ``header`` supplies study/series-level attributes (keywords). Pixel values
    are stored as (value - intercept) / slope in 16 bits."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lps = np.diag([-1.0, -1.0, 1.0, 1.0]) @ affine_ras
    sp = np.sqrt((lps[:3, :3] ** 2).sum(axis=0))
    row_dir = lps[:3, 0] / sp[0]
    col_dir = lps[:3, 1] / sp[1]
    study_uid = header.get("StudyInstanceUID") or new_uid()
    series_uid = header.get("SeriesInstanceUID") or new_uid()
    for_uid = header.get("FrameOfReferenceUID") or new_uid()
    stored = np.round((vol_data - intercept) / slope)
    if signed:
        stored = np.clip(stored, -32768, 32767).astype("<i2")
    else:
        stored = np.clip(stored, 0, 65535).astype("<u2")
    paths = []
    for k in range(vol_data.shape[2]):
        ipp = lps[:3, 3] + lps[:3, 2] * k
        plane = stored[:, :, k].T  # (rows=j, cols=i)
        el = dict(header)
        el.update({
            "StudyInstanceUID": study_uid, "SeriesInstanceUID": series_uid, "FrameOfReferenceUID": for_uid,
            "Modality": "MR", "InstanceNumber": k + 1,
            "ImagePositionPatient": [float(x) for x in ipp],
            "ImageOrientationPatient": [float(x) for x in np.concatenate([row_dir, col_dir])],
            "PixelSpacing": [float(sp[1]), float(sp[0])],
            "SliceThickness": float(header.get("SliceThickness") or np.linalg.norm(lps[:3, 2])),
            "SliceLocation": float(ipp[2]),
            "Rows": int(plane.shape[0]), "Columns": int(plane.shape[1]),
            "SamplesPerPixel": 1, "PhotometricInterpretation": "MONOCHROME2",
            "BitsAllocated": 16, "BitsStored": 16, "HighBit": 15, "PixelRepresentation": 1 if signed else 0,
            "RescaleSlope": slope, "RescaleIntercept": intercept,
            "PixelData": np.ascontiguousarray(plane).tobytes(),
        })
        el.setdefault("ImageType", ["ORIGINAL", "PRIMARY", "M", "ND"])
        el.setdefault("SliceLocation", float(ipp[2]))
        p = out_dir / f"{prefix}{k + 1:04d}.dcm"
        write_file(p, el, MR_IMAGE_STORAGE, new_uid())
        paths.append(p)
    return paths
