"""Minimal, dependency-free NIfTI-1/NIfTI-2 reader and NIfTI-1 writer.

Conventions used throughout SliceWeave
-----------------------------------------
* ``data`` is a numpy array indexed ``[i, j, k]`` (NIfTI order).
* ``affine`` maps voxel indices to world millimetres in RAS+ (NIfTI convention).
* Headers of input files are never copied to outputs: every file this module
  writes gets a fresh header whose ``descrip`` field holds only what the caller
  passes (never names, IDs or dates).
"""
from __future__ import annotations

import gzip
import io
import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_DTYPES = {
    2: np.uint8, 4: np.int16, 8: np.int32, 16: np.float32, 64: np.float64,
    256: np.int8, 512: np.uint16, 768: np.uint32, 1024: np.int64, 1280: np.uint64,
}
_CODES = {np.dtype(v).name: k for k, v in _DTYPES.items()}


@dataclass
class Volume:
    data: np.ndarray
    affine: np.ndarray
    meta: dict = field(default_factory=dict)

    @property
    def shape(self):
        return self.data.shape[:3]

    @property
    def spacing(self) -> np.ndarray:
        return np.sqrt((self.affine[:3, :3] ** 2).sum(axis=0))

    def like(self, data: np.ndarray, **meta) -> "Volume":
        return Volume(np.asarray(data), self.affine.copy(), dict(meta))


# --------------------------------------------------------------------------- read
def _open_bytes(path: Path) -> bytes:
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw


def _quat_to_affine(b, c, d, qfac, pix, offs) -> np.ndarray:
    a2 = 1.0 - (b * b + c * c + d * d)
    a = np.sqrt(a2) if a2 > 1e-7 else 0.0
    if a2 <= 1e-7:  # renormalise as in nifti1_io
        n = np.sqrt(b * b + c * c + d * d)
        b, c, d = b / n, c / n, d / n
    R = np.array([
        [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
        [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
        [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - c * c - b * b],
    ])
    z = np.array([pix[0], pix[1], pix[2] * qfac], dtype=float)
    A = np.eye(4)
    A[:3, :3] = R * z
    A[:3, 3] = offs
    return A


def read_nifti(path, dtype=np.float32) -> Volume:
    """Read a .nii or .nii.gz file. Only the first 3-D volume of 4-D data is kept."""
    raw = _open_bytes(Path(path))
    if len(raw) < 348:
        raise ValueError(f"{Path(path).name}: too short to be NIfTI")
    meta: dict = {}
    for endian in ("<", ">"):
        (sizeof_hdr,) = struct.unpack(endian + "i", raw[:4])
        if sizeof_hdr in (348, 540):
            break
    else:
        raise ValueError(f"{Path(path).name}: not a NIfTI file")
    e = endian
    if sizeof_hdr == 348:
        dim = struct.unpack(e + "8h", raw[40:56])
        datatype, = struct.unpack(e + "h", raw[70:72])
        pixdim = struct.unpack(e + "8f", raw[76:108])
        vox_offset, = struct.unpack(e + "f", raw[108:112])
        scl_slope, scl_inter = struct.unpack(e + "2f", raw[112:120])
        qform_code, sform_code = struct.unpack(e + "2h", raw[252:256])
        qb, qc, qd, qx, qy, qz = struct.unpack(e + "6f", raw[256:280])
        srow = np.array(struct.unpack(e + "12f", raw[280:328]), dtype=float).reshape(3, 4)
        magic = raw[344:348]
        if magic not in (b"n+1\x00", b"ni1\x00"):
            raise ValueError(f"{Path(path).name}: unsupported NIfTI magic (Analyze header?)")
        vox_offset = int(vox_offset)
    else:
        dim = struct.unpack(e + "8q", raw[16:80])
        datatype, = struct.unpack(e + "h", raw[12:14])
        pixdim = struct.unpack(e + "8d", raw[104:168])
        vox_offset, = struct.unpack(e + "q", raw[168:176])
        scl_slope, scl_inter = struct.unpack(e + "2d", raw[176:192])
        qform_code, sform_code = struct.unpack(e + "2i", raw[344:352])
        qb, qc, qd, qx, qy, qz = struct.unpack(e + "6d", raw[352:400])
        srow = np.array(struct.unpack(e + "12d", raw[400:496]), dtype=float).reshape(3, 4)
    if datatype not in _DTYPES:
        raise ValueError(f"{Path(path).name}: unsupported NIfTI datatype {datatype}")
    ndim = int(dim[0])
    shape = [int(x) for x in dim[1:1 + ndim]]
    while len(shape) < 3:
        shape.append(1)
    nvox3 = int(np.prod(shape[:3]))
    dt = np.dtype(_DTYPES[datatype]).newbyteorder(e)
    arr = np.frombuffer(raw, dtype=dt, count=nvox3, offset=vox_offset)
    arr = arr.reshape(shape[:3], order="F")
    if ndim > 3 and int(np.prod(shape[3:])) > 1:
        meta["warning"] = "4-D file: only the first volume was read"
        meta["n_volumes"] = int(np.prod(shape[3:]))
    # free-text description field; read only so known tools can be recognised, never stored
    d0, d1 = (148, 228) if sizeof_hdr == 348 else (240, 320)
    meta["descrip"] = raw[d0:d1].split(b"\x00", 1)[0].decode("ascii", "replace")
    if scl_slope not in (0.0, 1.0) or (scl_inter != 0.0 and np.isfinite(scl_inter)):
        slope = scl_slope if scl_slope != 0 else 1.0
        arr = arr.astype(np.float64) * slope + scl_inter
    if sform_code > 0 and np.all(np.isfinite(srow)) and abs(np.linalg.det(srow[:, :3])) > 1e-9:
        aff = np.eye(4)
        aff[:3, :] = srow
    elif qform_code > 0:
        qfac = -1.0 if pixdim[0] < 0 else 1.0
        aff = _quat_to_affine(qb, qc, qd, qfac, pixdim[1:4], (qx, qy, qz))
    else:
        aff = np.diag([pixdim[1] or 1.0, pixdim[2] or 1.0, pixdim[3] or 1.0, 1.0])
        meta["warning"] = "no qform/sform: orientation unknown, assumed RAS"
    out = np.ascontiguousarray(arr.astype(dtype)) if dtype is not None else np.array(arr)
    return Volume(out, aff, meta)


# -------------------------------------------------------------------------- write
def _affine_to_quat(aff: np.ndarray):
    M = aff[:3, :3].astype(float)
    zooms = np.sqrt((M ** 2).sum(axis=0))
    zooms[zooms == 0] = 1.0
    R = M / zooms
    qfac = 1.0
    if np.linalg.det(R) < 0:
        qfac = -1.0
        R[:, 2] = -R[:, 2]
    # nearest proper rotation (polar decomposition)
    u, _, vt = np.linalg.svd(R)
    R = u @ vt
    tr = np.trace(R)
    if tr > 0:
        s = np.sqrt(tr + 1.0) * 2
        a = 0.25 * s
        b = (R[2, 1] - R[1, 2]) / s
        c = (R[0, 2] - R[2, 0]) / s
        d = (R[1, 0] - R[0, 1]) / s
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2
        a = (R[2, 1] - R[1, 2]) / s
        b = 0.25 * s
        c = (R[0, 1] + R[1, 0]) / s
        d = (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2
        a = (R[0, 2] - R[2, 0]) / s
        b = (R[0, 1] + R[1, 0]) / s
        c = 0.25 * s
        d = (R[1, 2] + R[2, 1]) / s
    else:
        s = np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2
        a = (R[1, 0] - R[0, 1]) / s
        b = (R[0, 2] + R[2, 0]) / s
        c = (R[1, 2] + R[2, 1]) / s
        d = 0.25 * s
    if a < 0:
        a, b, c, d = -a, -b, -c, -d
    return (b, c, d), qfac, zooms


def nifti_bytes(data: np.ndarray, affine: np.ndarray, descrip: str = "",
                intent_name: str = "") -> bytes:
    data = np.asarray(data)
    if data.dtype == np.bool_:
        data = data.astype(np.uint8)
    if data.dtype.name not in _CODES:
        data = data.astype(np.float32)
    data = data.astype(data.dtype.newbyteorder("<"))
    shape = list(data.shape[:3]) + [1] * (3 - data.ndim)
    hdr = bytearray(352)
    struct.pack_into("<i", hdr, 0, 348)
    dim = [3] + shape + [1, 1, 1, 1]
    struct.pack_into("<8h", hdr, 40, *dim)
    code = _CODES[data.dtype.name]
    struct.pack_into("<hh", hdr, 70, code, data.dtype.itemsize * 8)
    (qb, qc, qd), qfac, zooms = _affine_to_quat(affine)
    struct.pack_into("<8f", hdr, 76, qfac, zooms[0], zooms[1], zooms[2], 0, 0, 0, 0)
    struct.pack_into("<f", hdr, 108, 352.0)
    struct.pack_into("<2f", hdr, 112, 1.0, 0.0)
    hdr[123] = 10  # xyzt_units: mm + sec
    d = descrip.encode("ascii", "replace")[:79]
    hdr[148:148 + len(d)] = d
    struct.pack_into("<2h", hdr, 252, 1, 1)  # qform_code, sform_code = scanner
    struct.pack_into("<6f", hdr, 256, qb, qc, qd, *affine[:3, 3])
    struct.pack_into("<12f", hdr, 280, *np.asarray(affine[:3, :], dtype=float).ravel())
    nm = intent_name.encode("ascii", "replace")[:15]
    hdr[328:328 + len(nm)] = nm
    hdr[344:348] = b"n+1\x00"
    return bytes(hdr) + data.tobytes(order="F")


def write_nifti(path, data: np.ndarray, affine: np.ndarray, descrip: str = "",
                compress: bool | None = None) -> Path:
    path = Path(path)
    b = nifti_bytes(data, affine, descrip)
    if compress is None:
        compress = path.name.endswith(".gz")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    if compress:
        with gzip.open(tmp, "wb", compresslevel=4) as f:
            f.write(b)
    else:
        tmp.write_bytes(b)
    tmp.replace(path)
    return path


def save(path, vol: Volume, descrip: str = "") -> Path:
    return write_nifti(path, vol.data, vol.affine, descrip)


def load(path, dtype=np.float32) -> Volume:
    return read_nifti(path, dtype=dtype)
