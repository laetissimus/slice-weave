"""Self-contained DICOM Part-10 reader (no third-party dependencies).

Supports implicit/explicit VR little endian, explicit VR big endian, deflated
explicit VR little endian and RLE Lossless. Other encapsulated (JPEG family)
transfer syntaxes are decoded through the optional pydicom hook
(:func:`sliceweave.hooks.pydicom_pixels`) when pydicom and a pixel-data plugin are
installed; otherwise a clear error is raised.

Only what the workbench needs is parsed into Python values; everything else is
kept as raw bytes and is never written anywhere.
"""
from __future__ import annotations

import struct
import zlib
from pathlib import Path

import numpy as np


class DicomError(Exception):
    pass


class UnsupportedTransferSyntax(DicomError):
    pass


TS_IMPLICIT_LE = "1.2.840.10008.1.2"
TS_EXPLICIT_LE = "1.2.840.10008.1.2.1"
TS_DEFLATED = "1.2.840.10008.1.2.1.99"
TS_EXPLICIT_BE = "1.2.840.10008.1.2.2"
TS_RLE = "1.2.840.10008.1.2.5"
NATIVE_TS = {TS_IMPLICIT_LE, TS_EXPLICIT_LE, TS_DEFLATED, TS_EXPLICIT_BE}

ITEM, ITEM_END, SEQ_END = 0xFFFEE000, 0xFFFEE00D, 0xFFFEE0DD
LONG_VRS = {b"OB", b"OD", b"OF", b"OL", b"OV", b"OW", b"SQ", b"SV", b"UC", b"UN", b"UR", b"UT", b"UV"}
TEXT_VRS = {"AE", "AS", "CS", "DA", "DS", "DT", "IS", "LO", "LT", "PN", "SH", "ST", "TM", "UC", "UI", "UR", "UT"}
NOSPLIT_VRS = {"LT", "ST", "UT", "UR"}
BIN_FMT = {"US": "H", "SS": "h", "UL": "I", "SL": "i", "FL": "f", "FD": "d", "SV": "q", "UV": "Q"}

# keyword -> (tag, VR). VR is used when the file is implicit VR.
TAGS = {
    "TransferSyntaxUID": (0x00020010, "UI"), "MediaStorageSOPClassUID": (0x00020002, "UI"),
    "SpecificCharacterSet": (0x00080005, "CS"), "ImageType": (0x00080008, "CS"),
    "SOPClassUID": (0x00080016, "UI"), "SOPInstanceUID": (0x00080018, "UI"),
    "StudyDate": (0x00080020, "DA"), "SeriesDate": (0x00080021, "DA"),
    "AcquisitionDate": (0x00080022, "DA"), "ContentDate": (0x00080023, "DA"),
    "AcquisitionDateTime": (0x0008002A, "DT"), "StudyTime": (0x00080030, "TM"),
    "SeriesTime": (0x00080031, "TM"), "AcquisitionTime": (0x00080032, "TM"),
    "ContentTime": (0x00080033, "TM"), "AccessionNumber": (0x00080050, "SH"),
    "Modality": (0x00080060, "CS"), "Manufacturer": (0x00080070, "LO"),
    "InstitutionName": (0x00080080, "LO"), "InstitutionAddress": (0x00080081, "ST"),
    "ReferringPhysicianName": (0x00080090, "PN"), "StationName": (0x00081010, "SH"),
    "StudyDescription": (0x00081030, "LO"), "SeriesDescription": (0x0008103E, "LO"),
    "PerformingPhysicianName": (0x00081050, "PN"), "OperatorsName": (0x00081070, "PN"),
    "ManufacturerModelName": (0x00081090, "LO"), "FrameType": (0x00089007, "CS"),
    "AcquisitionContrast": (0x00089209, "CS"),
    "PatientName": (0x00100010, "PN"), "PatientID": (0x00100020, "LO"),
    "PatientBirthDate": (0x00100030, "DA"), "PatientSex": (0x00100040, "CS"),
    "OtherPatientIDs": (0x00101000, "LO"), "OtherPatientNames": (0x00101001, "PN"),
    "PatientAge": (0x00101010, "AS"),
    "ContrastBolusAgent": (0x00180010, "LO"), "ContrastBolusAgentSequence": (0x00180012, "SQ"),
    "BodyPartExamined": (0x00180015, "CS"), "ScanningSequence": (0x00180020, "CS"),
    "SequenceVariant": (0x00180021, "CS"), "ScanOptions": (0x00180022, "CS"),
    "MRAcquisitionType": (0x00180023, "CS"), "SequenceName": (0x00180024, "SH"),
    "SliceThickness": (0x00180050, "DS"), "RepetitionTime": (0x00180080, "DS"),
    "EchoTime": (0x00180081, "DS"), "InversionTime": (0x00180082, "DS"),
    "EchoNumbers": (0x00180086, "IS"), "MagneticFieldStrength": (0x00180087, "DS"),
    "SpacingBetweenSlices": (0x00180088, "DS"), "DeviceSerialNumber": (0x00181000, "LO"),
    "ProtocolName": (0x00181030, "LO"), "ContrastBolusRoute": (0x00181040, "LO"),
    "TriggerTime": (0x00181060, "DS"), "FlipAngle": (0x00181314, "DS"),
    "PulseSequenceName": (0x00189005, "SH"), "InversionTimes": (0x00189079, "FD"),
    "EffectiveEchoTime": (0x00189082, "FD"), "DiffusionBValue": (0x00189087, "FD"),
    "MRTimingAndRelatedParametersSequence": (0x00189112, "SQ"),
    "MREchoSequence": (0x00189114, "SQ"), "MRModifierSequence": (0x00189115, "SQ"),
    "MRDiffusionSequence": (0x00189117, "SQ"), "MRImageFrameTypeSequence": (0x00189226, "SQ"),
    "ContrastBolusUsageSequence": (0x00189341, "SQ"),
    "StudyInstanceUID": (0x0020000D, "UI"), "SeriesInstanceUID": (0x0020000E, "UI"),
    "StudyID": (0x00200010, "SH"), "SeriesNumber": (0x00200011, "IS"),
    "AcquisitionNumber": (0x00200012, "IS"), "InstanceNumber": (0x00200013, "IS"),
    "ImagePositionPatient": (0x00200032, "DS"), "ImageOrientationPatient": (0x00200037, "DS"),
    "FrameOfReferenceUID": (0x00200052, "UI"), "TemporalPositionIdentifier": (0x00200100, "IS"),
    "SliceLocation": (0x00201041, "DS"), "ImageComments": (0x00204000, "LT"),
    "StackID": (0x00209056, "SH"), "InStackPositionNumber": (0x00209057, "UL"),
    "FrameContentSequence": (0x00209111, "SQ"), "PlanePositionSequence": (0x00209113, "SQ"),
    "PlaneOrientationSequence": (0x00209116, "SQ"), "TemporalPositionIndex": (0x00209128, "UL"),
    "DimensionIndexValues": (0x00209157, "UL"),
    "SamplesPerPixel": (0x00280002, "US"), "PhotometricInterpretation": (0x00280004, "CS"),
    "PlanarConfiguration": (0x00280006, "US"), "NumberOfFrames": (0x00280008, "IS"),
    "Rows": (0x00280010, "US"), "Columns": (0x00280011, "US"), "PixelSpacing": (0x00280030, "DS"),
    "BitsAllocated": (0x00280100, "US"), "BitsStored": (0x00280101, "US"),
    "HighBit": (0x00280102, "US"), "PixelRepresentation": (0x00280103, "US"),
    "BurnedInAnnotation": (0x00280301, "CS"), "WindowCenter": (0x00281050, "DS"),
    "WindowWidth": (0x00281051, "DS"), "RescaleIntercept": (0x00281052, "DS"),
    "RescaleSlope": (0x00281053, "DS"), "PixelMeasuresSequence": (0x00289110, "SQ"),
    "PixelValueTransformationSequence": (0x00289145, "SQ"),
    "SharedFunctionalGroupsSequence": (0x52009229, "SQ"),
    "PerFrameFunctionalGroupsSequence": (0x52009230, "SQ"),
    "PixelData": (0x7FE00010, "OW"),
    # private b-value tags read before de-identification (never stored)
    "SiemensBValue": (0x0019100C, "IS"), "GEBValue": (0x00431039, "IS"),
    "PhilipsBValue": (0x20011003, "FL"),
}
TAG2KW = {t: k for k, (t, _) in TAGS.items()}
TAG2VR = {t: vr for (t, vr) in TAGS.values()}
PIXEL_TAG = 0x7FE00010


class Dataset(dict):
    """dict tag -> (vr, value). Attribute access by keyword returns the parsed value."""

    def __getattr__(self, kw):
        if kw.startswith("_"):
            raise AttributeError(kw)
        tag = TAGS.get(kw, (None,))[0]
        if tag is None or tag not in self:
            return None
        return self[tag][1]

    def get_kw(self, kw, default=None):
        v = getattr(self, kw)
        return default if v is None else v

    def first(self, kw, default=None):
        v = getattr(self, kw)
        if v is None:
            return default
        if isinstance(v, (list, tuple)):
            return v[0] if len(v) else default
        return v

    def has_overlays(self) -> bool:
        return any(0x6000 <= (t >> 16) <= 0x601E and (t >> 16) % 2 == 0 for t in self)

    def has_private(self) -> bool:
        return any((t >> 16) % 2 == 1 for t in self)


# ------------------------------------------------------------------ parsing core
class _Reader:
    def __init__(self, buf: bytes, little: bool, explicit: bool, charset: str = "latin-1",
                 stop_before_pixels: bool = False):
        self.buf = buf
        self.little = little
        self.explicit = explicit
        self.e = "<" if little else ">"
        self.charset = charset
        self.stop = stop_before_pixels

    def tag_at(self, pos):
        g, el = struct.unpack_from(self.e + "HH", self.buf, pos)
        return (g << 16) | el

    def parse(self, pos: int, end: int, until_item_end: bool = False):
        ds = Dataset()
        buf = self.buf
        n = len(buf)
        while pos < end and pos + 8 <= n:
            tag = self.tag_at(pos)
            if tag == ITEM_END:
                return ds, pos + 8
            if tag in (ITEM, SEQ_END):
                # malformed; stop at this level
                return ds, pos
            if self.stop and tag >= PIXEL_TAG:
                return ds, end
            explicit = self.explicit and (tag >> 16) != 0xFFFE
            if explicit:
                vr = buf[pos + 4:pos + 6]
                if vr in LONG_VRS:
                    (length,) = struct.unpack_from(self.e + "I", buf, pos + 8)
                    hdr = 12
                elif vr.isalpha() and vr.isupper():
                    (length,) = struct.unpack_from(self.e + "H", buf, pos + 6)
                    hdr = 8
                else:  # VR missing even though explicit: fall back to implicit for this element
                    vr = TAG2VR.get(tag, "UN").encode()
                    (length,) = struct.unpack_from(self.e + "I", buf, pos + 4)
                    hdr = 8
                vr = vr.decode()
            else:
                vr = TAG2VR.get(tag, "UN")
                (length,) = struct.unpack_from(self.e + "I", buf, pos + 4)
                hdr = 8
            vstart = pos + hdr
            if length == 0xFFFFFFFF:
                if tag == PIXEL_TAG:
                    frags, pos = self._fragments(vstart)
                    ds[tag] = ("OB", frags)
                    continue
                # undefined length: SQ (or UN that must be read as implicit LE SQ)
                if vr == "UN" and self.explicit:
                    sub = _Reader(buf, True, False, self.charset)
                    items, pos = sub._sequence(vstart, None)
                else:
                    items, pos = self._sequence(vstart, None)
                ds[tag] = ("SQ", items)
                continue
            vend = vstart + length
            if vend > n:
                raise DicomError("truncated file")
            if vr == "SQ":
                items, _ = self._sequence(vstart, vend)
                ds[tag] = ("SQ", items)
            elif tag == PIXEL_TAG:
                ds[tag] = (vr, (vstart, length))
            elif tag in TAG2VR or (tag >> 16) == 0x0002:
                ds[tag] = (vr, self._value(vr if vr != "UN" else TAG2VR.get(tag, "UN"), buf[vstart:vend]))
            else:
                ds[tag] = (vr, None)  # not needed: kept only as a marker (private/overlay checks)
            pos = vend
        return ds, pos

    def _sequence(self, pos, end):
        items = []
        buf = self.buf
        while end is None or pos < end:
            if pos + 8 > len(buf):
                break
            tag = self.tag_at(pos)
            (length,) = struct.unpack_from(self.e + "I", buf, pos + 4)
            if tag == SEQ_END:
                return items, pos + 8
            if tag != ITEM:
                raise DicomError("bad sequence item")
            if length == 0xFFFFFFFF:
                item, pos = self.parse(pos + 8, len(buf), until_item_end=True)
            else:
                item, _ = self.parse(pos + 8, pos + 8 + length)
                pos = pos + 8 + length
            items.append(item)
        return items, pos

    def _fragments(self, pos):
        frags = []
        buf = self.buf
        while pos + 8 <= len(buf):
            tag = self.tag_at(pos)
            (length,) = struct.unpack_from(self.e + "I", buf, pos + 4)
            if tag == SEQ_END:
                return frags, pos + 8
            if tag != ITEM:
                raise DicomError("bad encapsulated pixel data")
            frags.append((pos + 8, length))
            pos += 8 + length
        return frags, pos

    def _value(self, vr, raw: bytes):
        if vr in TEXT_VRS:
            enc = "ascii" if vr in ("UI", "DA", "TM", "DT", "DS", "IS", "AS", "CS") else self.charset
            s = raw.decode(enc, "replace").rstrip("\x00 ").lstrip(" ") if vr != "UI" else raw.decode("ascii", "replace").rstrip("\x00 ")
            if vr in NOSPLIT_VRS:
                return s
            parts = [p.strip() for p in s.split("\\")] if s else []
            if vr in ("DS", "IS"):
                out = []
                for p in parts:
                    try:
                        out.append(float(p) if vr == "DS" else int(float(p)))
                    except ValueError:
                        out.append(None)
                parts = out
            if not parts:
                return ""
            if len(parts) == 1:
                return parts[0]
            return parts
        if vr in BIN_FMT:
            fmt = BIN_FMT[vr]
            k = struct.calcsize(fmt)
            cnt = len(raw) // k
            vals = list(struct.unpack(self.e + fmt * cnt, raw[:cnt * k]))
            return vals[0] if cnt == 1 else vals
        return raw


_CHARSETS = {"ISO_IR 192": "utf-8", "ISO_IR 100": "latin-1", "": "latin-1", "ISO_IR 6": "ascii",
             "GB18030": "gb18030", "ISO_IR 101": "iso8859_2", "ISO_IR 144": "iso8859_5"}


def read(path, stop_before_pixels: bool = False) -> Dataset:
    """Parse a DICOM file. Raises DicomError for non-DICOM files."""
    raw = Path(path).read_bytes()
    return read_bytes(raw, stop_before_pixels)


def read_bytes(raw: bytes, stop_before_pixels: bool = False) -> Dataset:
    if len(raw) >= 132 and raw[128:132] == b"DICM":
        pos = 132
    elif len(raw) > 8 and raw[0:2] in (b"\x02\x00", b"\x08\x00") :
        pos = 0  # no preamble
    else:
        raise DicomError("not a DICOM file")
    meta = Dataset()
    ts = TS_IMPLICIT_LE if pos == 0 else None
    if pos == 132 or raw[0:2] == b"\x02\x00":
        r = _Reader(raw, True, True)
        # file meta: group 2 explicit VR LE
        end = len(raw)
        while pos + 8 <= end:
            g = struct.unpack_from("<H", raw, pos)[0]
            if g != 0x0002:
                break
            m, pos = r.parse(pos, pos + 1)  # parse a single element
            meta.update(m)
            # parse() returns after reading one element because end=pos+1
        ts = meta.TransferSyntaxUID or TS_EXPLICIT_LE
    body = raw
    if ts == TS_DEFLATED:
        body = raw[:pos] + zlib.decompress(raw[pos:], -15)
    little = ts != TS_EXPLICIT_BE
    explicit = ts != TS_IMPLICIT_LE
    # heuristic: some files claim implicit but are explicit (or no meta at all)
    if pos + 6 <= len(body):
        vr = body[pos + 4:pos + 6]
        looks_explicit = vr.isalpha() and vr.isupper()
        if not explicit and looks_explicit and ts == TS_IMPLICIT_LE and "TransferSyntaxUID" not in meta:
            explicit = True
    cs_reader = _Reader(body, little, explicit)
    ds, _ = cs_reader.parse(pos, len(body))
    cs = ds.SpecificCharacterSet
    if isinstance(cs, list):
        cs = cs[-1] if cs else ""
    charset = _CHARSETS.get(cs or "", "latin-1")
    if charset != "latin-1":
        ds, _ = _Reader(body, little, explicit, charset, stop_before_pixels).parse(pos, len(body))
    elif stop_before_pixels:
        ds.pop(PIXEL_TAG, None)
    ds.update(meta)
    ds._buf = body  # type: ignore[attr-defined]
    ds._ts = ts  # type: ignore[attr-defined]
    ds._little = little  # type: ignore[attr-defined]
    return ds


# --------------------------------------------------------------------- pixels
def _rle_decode_frame(data: bytes, rows: int, cols: int, spp: int, bits: int) -> np.ndarray:
    nseg = struct.unpack_from("<I", data, 0)[0]
    offsets = list(struct.unpack_from("<15I", data, 4))[:nseg]
    offsets.append(len(data))
    bytes_per_sample = bits // 8
    if nseg != spp * bytes_per_sample:
        raise DicomError("RLE segment count mismatch")
    npix = rows * cols
    segs = []
    for s in range(nseg):
        src = data[offsets[s]:offsets[s + 1]]
        out = bytearray()
        i = 0
        L = len(src)
        while i < L and len(out) < npix:
            n = src[i]
            i += 1
            if n < 128:
                out += src[i:i + n + 1]
                i += n + 1
            elif n > 128:
                out += bytes([src[i]]) * (257 - n)
                i += 1
        segs.append(np.frombuffer(bytes(out[:npix]).ljust(npix, b"\x00"), dtype=np.uint8))
    planes = []
    for smp in range(spp):
        acc = np.zeros(npix, dtype=np.uint64)
        for b in range(bytes_per_sample):  # most significant byte first
            acc = (acc << np.uint64(8)) | segs[smp * bytes_per_sample + b].astype(np.uint64)
        planes.append(acc)
    arr = np.stack(planes, axis=-1).reshape(rows, cols, spp)
    return arr


def pixel_array(ds: Dataset, path=None, rescale: bool = True) -> np.ndarray:
    """Return pixel data as float32 array (frames, rows, cols) for grayscale data.

    Colour images raise DicomError: the workbench never imports colour pixels
    (screenshots, burned-in annotations)."""
    rows, cols = int(ds.Rows or 0), int(ds.Columns or 0)
    spp = int(ds.SamplesPerPixel or 1)
    bits = int(ds.BitsAllocated or 16)
    stored = int(ds.BitsStored or bits)
    signed = int(ds.PixelRepresentation or 0) == 1
    nframes = int(ds.NumberOfFrames or 1) if ds.NumberOfFrames not in (None, "") else 1
    photo = str(ds.PhotometricInterpretation or "MONOCHROME2").upper()
    if spp != 1 or not photo.startswith("MONOCHROME"):
        raise DicomError(f"colour image ({photo}) refused: possible screenshot or burned-in annotation")
    if PIXEL_TAG not in ds:
        raise DicomError("no pixel data")
    vr, val = ds[PIXEL_TAG]
    buf = ds._buf  # type: ignore[attr-defined]
    ts = ds._ts  # type: ignore[attr-defined]
    if isinstance(val, list):  # encapsulated
        if ts == TS_RLE:
            frags = val[1:] if len(val) > nframes else val  # first item is basic offset table
            if len(frags) < nframes:
                raise DicomError("RLE: fewer fragments than frames")
            frames = [_rle_decode_frame(buf[o:o + n], rows, cols, spp, bits)[..., 0] for (o, n) in frags[:nframes]]
            arr = np.stack(frames).astype(np.int64)
            if signed:
                arr = np.where(arr >= (1 << (bits - 1)), arr - (1 << bits), arr)
        else:
            from . import hooks
            arr = hooks.pydicom_pixels(path, ts)
            arr = np.asarray(arr)
            if arr.ndim == 2:
                arr = arr[None]
            rescale_needed = rescale
            return _apply_rescale(ds, arr.astype(np.float32)) if rescale_needed else arr.astype(np.float32)
    else:
        start, length = val
        dt = {8: np.uint8, 16: np.uint16, 32: np.uint32}.get(bits)
        if dt is None:
            raise DicomError(f"unsupported BitsAllocated {bits}")
        dt = np.dtype(dt).newbyteorder("<" if ds._little else ">")  # type: ignore[attr-defined]
        count = rows * cols * nframes
        if length < count * dt.itemsize:
            raise DicomError("pixel data shorter than expected")
        arr = np.frombuffer(buf, dtype=dt, count=count, offset=start).astype(np.int64)
        arr = arr.reshape(nframes, rows, cols)
        if stored < bits:
            arr = arr & ((1 << stored) - 1)
            if signed:
                arr = np.where(arr >= (1 << (stored - 1)), arr - (1 << stored), arr)
        elif signed:
            arr = np.where(arr >= (1 << (bits - 1)), arr - (1 << bits), arr)
    arr = arr.astype(np.float32)
    return _apply_rescale(ds, arr) if rescale else arr


def _apply_rescale(ds: Dataset, arr: np.ndarray) -> np.ndarray:
    # Enhanced multi-frame: per-frame or shared PixelValueTransformationSequence
    nfr = arr.shape[0]
    slopes = np.ones(nfr, np.float32)
    inters = np.zeros(nfr, np.float32)
    shared = (ds.SharedFunctionalGroupsSequence or [None])[0]
    per = ds.PerFrameFunctionalGroupsSequence or []
    found = False
    for f in range(nfr):
        for grp in ([per[f]] if f < len(per) else []) + ([shared] if shared else []):
            pv = grp.PixelValueTransformationSequence if grp else None
            if pv:
                s = pv[0].RescaleSlope
                i = pv[0].RescaleIntercept
                slopes[f] = float(s) if s not in (None, "") else 1.0
                inters[f] = float(i) if i not in (None, "") else 0.0
                found = True
                break
    if not found:
        s = ds.RescaleSlope
        i = ds.RescaleIntercept
        slopes[:] = float(s) if isinstance(s, (int, float)) else 1.0
        inters[:] = float(i) if isinstance(i, (int, float)) else 0.0
    if np.any(slopes != 1) or np.any(inters != 0):
        arr = arr * slopes[:, None, None] + inters[:, None, None]
    return arr


def is_dicom_file(path) -> bool:
    try:
        with open(path, "rb") as f:
            head = f.read(132)
    except OSError:
        return False
    if len(head) >= 132 and head[128:132] == b"DICM":
        return True
    # preamble-less files: group 0x0008 or 0x0002 first, plausible VR or length
    return len(head) >= 8 and head[0:2] in (b"\x08\x00", b"\x02\x00") and head[2:4] in (b"\x00\x00", b"\x05\x00", b"\x08\x00", b"\x16\x00", b"\x10\x00", b"\x01\x00")
