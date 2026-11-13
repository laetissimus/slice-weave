"""Find image series in a folder (DICOM and NIfTI) and turn each 2-D multi-slice series into a
:class:`sliceweave.stacks.Stack` that keeps every slice's own position (no resampling).

Identifiers read here (PatientName, PatientID, dates...) live only in memory
inside :class:`SeriesInfo.phi`; :mod:`gbm_wb.deid` uses them to scrub free text
and to check the coded case ID, then discards them.
"""
from __future__ import annotations

import os
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import dicom
from .dicom import DicomError
from .nifti import read_nifti
from .stacks import Stack, plane_label

SC_SOP_PREFIX = "1.2.840.10008.5.1.4.1.1.7"  # (multi-frame) secondary capture
MAX_ZIP_BYTES = 30 * 1024 ** 3
MAX_ZIP_FILES = 200_000


class IngestError(Exception):
    pass


@dataclass
class Frame:
    path: str
    frame: int  # index within a multi-frame file (0 for single-frame)
    ipp: np.ndarray
    iop: np.ndarray
    spacing: tuple  # (row spacing, column spacing)
    rows: int
    cols: int
    instance: int = 0
    thickness: float | None = None
    acq_time: float | None = None  # seconds since midnight; used only to order slices in time, never stored


def _tm_seconds(x) -> float | None:
    t = _s(x).strip()
    if len(t) < 6 or not t[:6].isdigit():
        return None
    try:
        return int(t[0:2]) * 3600 + int(t[2:4]) * 60 + float(t[4:])
    except ValueError:
        return None


@dataclass
class SeriesInfo:
    uid: str
    study_uid: str
    kind: str  # "dicom" | "nifti"
    description: str = ""
    protocol: str = ""
    sequence: str = ""
    scanning: str = ""
    variant: str = ""
    acq_type: str = ""
    te: float | None = None
    tr: float | None = None
    ti: float | None = None
    flip: float | None = None
    field_t: float | None = None
    manufacturer: str = ""
    model: str = ""
    contrast_tag: bool = False
    image_type: list = field(default_factory=list)
    bvalue: float | None = None
    echo: float | None = None
    series_number: int | None = None
    study_date: str = ""   # raw YYYYMMDD; converted to a day offset by deid, never stored
    series_time: str = ""  # raw HHMMSS; converted to an order rank, never stored
    frames: list = field(default_factory=list)
    nifti_path: str = ""
    n_slices: int = 0
    thickness: float | None = None
    in_plane: tuple | None = None
    orientation: str = ""
    reject: str = ""
    warnings: list = field(default_factory=list)
    phi: dict = field(default_factory=dict)
    sop_class: str = ""

    def text(self) -> str:
        return " ".join(x for x in (self.description, self.protocol, self.sequence) if x)


# ------------------------------------------------------------------ helpers
def _f(x):
    if isinstance(x, (list, tuple)):
        x = x[0] if x else None
    try:
        return float(x) if x not in (None, "") else None
    except (TypeError, ValueError):
        return None


def _s(x) -> str:
    if x is None:
        return ""
    if isinstance(x, (list, tuple)):
        return "\\".join(str(v) for v in x)
    if isinstance(x, bytes):
        return ""
    return str(x)


def _lst(x) -> list:
    if x is None or x == "":
        return []
    return [str(v).upper() for v in (x if isinstance(x, (list, tuple)) else [x])]


def _orientation_label(iop) -> str:
    n = np.cross(iop[:3], iop[3:])
    ax = int(np.argmax(np.abs(n)))
    return ["SAG", "COR", "AX"][ax]


def _bvalue(ds) -> float | None:
    b = _f(ds.DiffusionBValue)
    if b is None and ds.MRDiffusionSequence:
        b = _f(ds.MRDiffusionSequence[0].DiffusionBValue)
    if b is None:
        b = _f(ds.SiemensBValue)
    if b is None and ds.GEBValue is not None:
        v = _f(ds.GEBValue)
        if v is not None:
            b = v % 1e9 if v >= 1e9 else v  # GE stores b + 1e9 on some versions
    if b is None:
        b = _f(ds.PhilipsBValue)
    return b


PHI_KEYS = ["PatientName", "PatientID", "OtherPatientIDs", "OtherPatientNames", "AccessionNumber",
            "PatientBirthDate", "InstitutionName", "ReferringPhysicianName", "PerformingPhysicianName",
            "OperatorsName", "StationName", "DeviceSerialNumber", "StudyID"]


def safe_extract_zip(zpath: Path, dest: Path) -> list[Path]:
    """Extract a zip without path traversal and with size/count limits."""
    out = []
    dest = dest.resolve()
    with zipfile.ZipFile(zpath) as z:
        infos = z.infolist()
        if len(infos) > MAX_ZIP_FILES:
            raise IngestError("zip has too many files")
        if sum(i.file_size for i in infos) > MAX_ZIP_BYTES:
            raise IngestError("zip expands to more than 30 GB")
        for i, info in enumerate(infos):
            if info.is_dir():
                continue
            # flatten to a numbered name: original names may contain identifiers
            target = (dest / f"z{i:06d}_{Path(info.filename).suffix[:8]}").resolve()
            if not str(target).startswith(str(dest)):
                raise IngestError("unsafe path in zip")
            with z.open(info) as src, open(target, "wb") as dst:
                while True:
                    chunk = src.read(1 << 20)
                    if not chunk:
                        break
                    dst.write(chunk)
            out.append(target)
    return out


def iter_files(paths) -> list[Path]:
    files = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            for root, dirs, fs in os.walk(p):
                dirs[:] = [d for d in dirs if not d.startswith(".")]
                for f in fs:
                    if not f.startswith("."):
                        files.append(Path(root) / f)
        elif p.exists():
            files.append(p)
    return sorted(files)


# ---------------------------------------------------------------- scanning
def scan(paths, tmpdir: Path | None = None, log=None) -> list[SeriesInfo]:
    """Scan files/folders (zips are extracted into ``tmpdir``) and return series."""
    files = iter_files(paths)
    expanded = []
    for f in files:
        if f.suffix.lower() == ".zip":
            if tmpdir is None:
                raise IngestError("zip input needs a temporary directory")
            sub = Path(tmpdir) / f"zip{len(expanded):04d}"
            sub.mkdir(parents=True, exist_ok=True)
            expanded.extend(safe_extract_zip(f, sub))
        else:
            expanded.append(f)
    series: dict[tuple, SeriesInfo] = {}
    nifti = []
    skipped = 0
    for f in expanded:
        name = f.name.lower()
        if name.endswith(".nii") or name.endswith(".nii.gz"):
            nifti.append(f)
            continue
        if not dicom.is_dicom_file(f):
            skipped += 1
            continue
        try:
            ds = dicom.read(f, stop_before_pixels=False)
        except DicomError:
            skipped += 1
            continue
        _add_dicom(series, ds, str(f))
    out = list(series.values())
    for s in out:
        _finish(s)
    for f in nifti:
        out.append(_nifti_series(f))
    if log:
        log(f"scanned {len(expanded)} files: {len(out)} series, {skipped} non-image files skipped")
    return out


def _add_dicom(series: dict, ds: dicom.Dataset, path: str):
    sop_class = _s(ds.SOPClassUID or ds.MediaStorageSOPClassUID)
    modality = _s(ds.Modality)
    study = _s(ds.StudyInstanceUID) or "unknown-study"
    suid = _s(ds.SeriesInstanceUID) or "unknown-series"
    image_type = _lst(ds.ImageType)
    frames: list[Frame] = []
    reject = ""
    if modality and modality != "MR":
        reject = f"modality {modality}"
    if sop_class.startswith(SC_SOP_PREFIX):
        reject = "secondary capture (possible screenshot / burned-in text)"
    if _s(ds.BurnedInAnnotation).upper() == "YES":
        reject = "BurnedInAnnotation = YES"
    photo = _s(ds.PhotometricInterpretation).upper()
    if photo and not photo.startswith("MONOCHROME"):
        reject = f"colour pixels ({photo})"
    if dicom.PIXEL_TAG not in ds:
        reject = reject or "no pixel data"
    rows, cols = int(ds.Rows or 0), int(ds.Columns or 0)
    nfr = int(_f(ds.NumberOfFrames) or 1)
    b_series = _bvalue(ds)
    te_series = _f(ds.EchoTime)
    split_keys = []
    shared = (ds.SharedFunctionalGroupsSequence or [None])[0]
    perframe = ds.PerFrameFunctionalGroupsSequence or []
    if nfr > 1 or perframe:
        for k in range(nfr):
            pf = perframe[k] if k < len(perframe) else None

            def g(seq_kw, el_kw):
                for grp in (pf, shared):
                    if grp is None:
                        continue
                    sq = getattr(grp, seq_kw)
                    if sq:
                        v = getattr(sq[0], el_kw)
                        if v not in (None, ""):
                            return v
                return None
            ipp = g("PlanePositionSequence", "ImagePositionPatient")
            iop = g("PlaneOrientationSequence", "ImageOrientationPatient")
            sp = g("PixelMeasuresSequence", "PixelSpacing")
            th = g("PixelMeasuresSequence", "SliceThickness")
            if ipp is None or iop is None or sp is None:
                reject = reject or "multi-frame file without per-frame geometry"
                break
            b = _f(g("MRDiffusionSequence", "DiffusionBValue"))
            te = _f(g("MREchoSequence", "EffectiveEchoTime"))
            ftype = _lst(g("MRImageFrameTypeSequence", "FrameType"))
            frames.append(Frame(path, k, np.array(ipp, float), np.array(iop, float), tuple(map(float, sp)),
                                rows, cols, k, _f(th)))
            split_keys.append((b, te, tuple(ftype[2:3])))
    else:
        ipp, iop, sp = ds.ImagePositionPatient, ds.ImageOrientationPatient, ds.PixelSpacing
        if ipp is None or iop is None or sp is None or not isinstance(ipp, list) or len(ipp) != 3:
            reject = reject or "no image geometry"
        else:
            frames.append(Frame(path, 0, np.array(ipp, float), np.array(iop, float),
                                tuple(map(float, sp if isinstance(sp, list) else [sp, sp])),
                                rows, cols, int(_f(ds.InstanceNumber) or 0), _f(ds.SliceThickness),
                                _tm_seconds(ds.AcquisitionTime or ds.ContentTime)))
            split_keys.append((b_series, te_series if _f(ds.EchoNumbers) not in (None, 1.0) else None,
                               tuple(image_type[2:3])))
    if not frames:
        frames = [None]
        split_keys = [(None, None, ())]
    for fr, sk in zip(frames, split_keys):
        iop_key = tuple(np.round(fr.iop, 3)) if fr is not None else ()
        key = (study, suid, sk, iop_key, (fr.rows, fr.cols) if fr is not None else ())
        s = series.get(key)
        if s is None:
            s = SeriesInfo(uid=suid + (f"#{len(series)}" if any(k[:2] == (study, suid) for k in series) else ""),
                           study_uid=study, kind="dicom")
            s.description = _s(ds.SeriesDescription)
            s.protocol = _s(ds.ProtocolName)
            s.sequence = _s(ds.SequenceName or ds.PulseSequenceName)
            s.scanning = _s(ds.ScanningSequence)
            s.variant = _s(ds.SequenceVariant)
            s.acq_type = _s(ds.MRAcquisitionType)
            s.te = te_series if sk[1] is None else sk[1]
            if s.te is None and frames and frames[0] is not None:
                pass
            s.tr, s.ti, s.flip = _f(ds.RepetitionTime), _f(ds.InversionTime), _f(ds.FlipAngle)
            if s.ti is None and ds.InversionTimes is not None:
                s.ti = _f(ds.InversionTimes)
            s.field_t = _f(ds.MagneticFieldStrength)
            s.manufacturer = _s(ds.Manufacturer)
            s.model = _s(ds.ManufacturerModelName)
            agent = _s(ds.ContrastBolusAgent).strip()
            s.contrast_tag = bool(agent) or bool(ds.ContrastBolusAgentSequence) or bool(ds.ContrastBolusUsageSequence)
            s.image_type = image_type
            s.bvalue = sk[0]
            s.series_number = int(_f(ds.SeriesNumber) or 0)
            s.study_date = _s(ds.StudyDate or ds.SeriesDate or ds.AcquisitionDate)
            s.series_time = _s(ds.SeriesTime or ds.AcquisitionTime or ds.ContentTime)
            s.sop_class = sop_class
            s.phi = {k: _s(getattr(ds, k)) for k in PHI_KEYS if getattr(ds, k) not in (None, "")}
            if ds.has_overlays():
                s.warnings.append("overlay planes present (never imported)")
            s.reject = reject
            series[key] = s
        if fr is not None:
            s.frames.append(fr)
        if reject and not s.reject:
            s.reject = reject


def _finish(s: SeriesInfo):
    fr = [f for f in s.frames if f is not None]
    s.frames = fr
    if not fr:
        s.reject = s.reject or "no frames"
        return
    iop = fr[0].iop
    s.orientation = _orientation_label(iop)
    n = np.cross(iop[:3], iop[3:])
    pos = np.array([f.ipp @ n for f in fr])
    s.n_slices = len(np.unique(np.round(pos, 2)))
    s.in_plane = fr[0].spacing
    ths = [f.thickness for f in fr if f.thickness]
    s.thickness = float(np.median(ths)) if ths else None
    if s.n_slices < 10 and not s.reject:
        s.reject = f"only {s.n_slices} slices (localizer or single slice)"
    if "LOCALIZER" in s.image_type and not s.reject:
        s.reject = "localizer"


def _nifti_series(path: Path) -> SeriesInfo:
    s = SeriesInfo(uid="nifti:" + str(path), study_uid="nifti-study", kind="nifti", nifti_path=str(path))
    name = path.name
    for ext in (".nii.gz", ".nii"):
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
    s.description = name  # file name used only for classification; not stored after ingest
    try:
        v = read_nifti(path, dtype=None)
        s.n_slices = int(v.shape[2])
        sp = v.spacing
        s.in_plane = (float(sp[1]), float(sp[0]))
        s.thickness = float(sp[2])
        s.orientation = "AX"
        if v.meta.get("warning"):
            s.warnings.append(v.meta["warning"])
    except Exception as e:
        s.reject = f"unreadable NIfTI: {e}"
    return s


# ---------------------------------------------------------------- building
HASTE_HINTS = ("haste", "ssfse", "ss-fse", "ssh_tse", "sshtse", "ss_tse", "single shot", "singleshot", "h2d", "ssfp_t2",
               "t2_ss", "fiesta", "trufi", "bfffe", "ss-tse", "sstse")


def looks_like_haste(s: SeriesInfo) -> bool:
    t = (s.text() + " " + s.sequence).lower()
    return any(h in t for h in HASTE_HINTS)


def build_stack(s: SeriesInfo, name: str = "") -> Stack:
    """Return a float32 Stack for a series. DICOM slices keep their own positions; duplicate
    positions (repeated acquisitions within one series) keep the first instance."""
    if s.kind == "nifti":
        return nifti_stack(Path(s.nifti_path), name or s.description)
    fr = sorted(s.frames, key=lambda f: (f.instance, f.path, f.frame))
    iop = fr[0].iop
    row_dir, col_dir = iop[:3] / np.linalg.norm(iop[:3]), iop[3:] / np.linalg.norm(iop[3:])
    n = np.cross(row_dir, col_dir)
    pos = np.array([f.ipp @ n for f in fr])
    order = np.argsort(pos, kind="stable")
    fr = [fr[i] for i in order]
    pos = pos[order]
    keep = [0]
    for i in range(1, len(fr)):
        if abs(pos[i] - pos[keep[-1]]) > 1e-3:
            keep.append(i)
    if len(keep) < len(fr):
        s.warnings.append(f"{len(fr) - len(keep)} repeated slice position(s): first instance kept")
    fr = [fr[i] for i in keep]
    if len(fr) < 2:
        raise IngestError("fewer than two slices")
    if any(f.rows != fr[0].rows or f.cols != fr[0].cols for f in fr):
        raise IngestError("slices of different sizes in one series")
    cache: dict[str, tuple] = {}
    planes = []
    for f in fr:
        if f.path not in cache:
            ds = dicom.read(f.path)
            cache.clear()
            cache[f.path] = (ds, dicom.pixel_array(ds, f.path))
        arr = cache[f.path][1]
        planes.append(arr[f.frame])
    data = np.stack(planes, axis=-1).astype(np.float32)   # (rows, cols, k)
    data = np.transpose(data, (1, 0, 2))                   # (cols = i, rows = j, k)
    rs, cs = fr[0].spacing  # PixelSpacing = [row spacing (between rows), column spacing]
    flip = np.diag([-1.0, -1.0, 1.0, 1.0])
    sa = np.zeros((len(fr), 4, 4))
    for k, f in enumerate(fr):
        A = np.eye(4)
        A[:3, 0] = row_dir * cs
        A[:3, 1] = col_dir * rs
        A[:3, 2] = n
        A[:3, 3] = f.ipp
        sa[k] = flip @ A
    gaps = np.diff(pos[keep])
    thick = s.thickness or (float(np.median(gaps)) if len(gaps) else 1.0)
    if len(gaps) and np.max(np.abs(gaps - np.median(gaps))) > max(0.05, 0.1 * np.median(gaps)):
        s.warnings.append(f"uneven slice gaps ({gaps.min():.2f}-{gaps.max():.2f} mm): kept as acquired")
    meta = {"source": "dicom"}
    times = [f.acq_time for f in fr]
    if all(t is not None for t in times) and len(set(np.round(times, 3))) == len(times):
        meta["acq_order"] = [int(i) for i in np.argsort(np.argsort(times, kind="stable"))]  # rank in time per slice
    st = Stack(name or s.description or "series", np.ascontiguousarray(data), sa, float(thick), None,
               plane_label(flip[:3, :3] @ n), meta)
    return st


def nifti_stack(path: Path, name: str = "") -> Stack:
    """A NIfTI file as a stack. The slice thickness comes from a BIDS/dcm2niix JSON sidecar
    (SliceThickness) when one is next to the file; otherwise it is taken as the slice spacing."""
    v = read_nifti(path)
    thick = None
    base = path.name[:-7] if path.name.lower().endswith(".nii.gz") else path.stem
    side = path.with_name(base + ".json")
    meta = {"source": "nifti"}
    if side.exists():
        try:
            import json
            j = json.loads(side.read_text(encoding="utf-8", errors="replace"))
            if isinstance(j.get("SliceThickness"), (int, float)) and j["SliceThickness"] > 0:
                thick = float(j["SliceThickness"])
                meta["thickness_from"] = "JSON sidecar"
            # technical text used only to recognise single-shot series and fetal exams (not stored)
            meta["sidecar_text"] = " ".join(str(j.get(k, "")) for k in (
                "SeriesDescription", "ProtocolName", "SequenceName", "PulseSequenceDetails", "ProcedureStepDescription",
                "BodyPartExamined") if isinstance(j.get(k, ""), str))
        except (OSError, ValueError):
            pass
    if v.meta.get("warning"):
        meta["warning"] = v.meta["warning"]
    return Stack.from_volume(name or base, v.data, v.affine, thickness=thick, meta=meta)


def is_mask_name(name: str) -> bool:
    n = name.lower()
    return any(t in n for t in ("mask", "_seg", "label", "brainmask", "roi"))
