"""Display grids and viewer layers.

The viewer shows every layer on one isotropic grid (directions of the reference stack,
permuted to be closest to RAS), so the three planes are always perpendicular:

* the *overview* grid covers the field of view of all included stacks (to place the region);
* the *result* grid is the reconstruction grid (possibly reoriented).

Input stacks are shown "as acquired": in-plane linear interpolation and the nearest slice
through-plane, so the slice thickness and gaps remain visible.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from .geometry import Grid, apply, make_grid, nearest_ras_axes
from .stacks import Stack


def overview_grid(stacks: list[Stack], ref: Stack, max_voxels: int = 8_000_000) -> Grid:
    R = nearest_ras_axes(ref.axes)
    pts = np.concatenate([s.corners_world() for s in stacks])
    q = pts @ R
    lo, hi = q.min(axis=0), q.max(axis=0)
    h = max(min(min(s.pixel) for s in stacks), 0.8)
    half = (hi - lo) / 2
    while np.prod(2 * half / h + 1) > max_voxels:
        h *= 1.1
    return make_grid(R, R @ ((lo + hi) / 2), half, h, fast=False)


def stack_on_grid(st: Stack, grid: Grid) -> np.ndarray:
    """Resample a stack onto the grid: bilinear in-plane, nearest slice through-plane; zero
    outside the stack (or farther than half a slice spacing from any slice)."""
    n = st.normal
    A0 = st.slice_aff[0].copy()
    A0[:3, 2] = n
    inv = np.linalg.inv(A0)
    pos = np.array([(A[:3, 3] - A0[:3, 3]) @ n for A in st.slice_aff])
    order = np.argsort(pos)
    pos_s = pos[order]
    half = max(st.spacing, st.thickness) / 2 + 1e-3
    out = np.zeros(grid.shape, np.float32)
    nx, ny, nzz = grid.shape
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    base = np.stack([ii.ravel(), jj.ravel(), np.zeros(ii.size)], axis=1).astype(np.float64)
    M = inv @ grid.affine
    for k in range(nzz):
        base[:, 2] = k
        q = apply(M, base)
        d = q[:, 2]
        idx = np.clip(np.searchsorted(pos_s, d), 1, len(pos_s) - 1) if len(pos_s) > 1 else np.zeros(len(d), int)
        if len(pos_s) > 1:
            left = pos_s[idx - 1]
            right = pos_s[idx]
            idx = np.where(np.abs(d - left) <= np.abs(d - right), idx - 1, idx)
        near = np.abs(d - pos_s[idx]) <= half
        s_idx = order[idx]
        # in-plane coordinates relative to the chosen slice (slices may not be exactly parallel-offset)
        coords = np.stack([q[:, 0], q[:, 1], s_idx.astype(np.float64)])
        v = ndi.map_coordinates(st.data, coords, order=1, mode="constant", cval=0.0, prefilter=False)
        inside = (q[:, 0] >= -0.5) & (q[:, 0] <= st.shape[0] - 0.5) & (q[:, 1] >= -0.5) & (q[:, 1] <= st.shape[1] - 0.5)
        out[:, :, k] = np.where(near & inside, v, 0).reshape(nx, ny)
    return out


def resample(vol: np.ndarray, src: Grid, dst: Grid, order: int = 1) -> np.ndarray:
    M = src.inv @ dst.affine
    out = np.empty(dst.shape, np.float32)
    nx, ny, nz = dst.shape
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    base = np.stack([ii.ravel(), jj.ravel(), np.zeros(ii.size)], axis=1).astype(np.float64)
    data = ndi.spline_filter(vol, order=order) if order > 1 else vol
    for k in range(nz):
        base[:, 2] = k
        q = apply(M, base).T
        out[:, :, k] = ndi.map_coordinates(data, q, order=order, mode="constant", cval=0.0,
                                           prefilter=False).reshape(nx, ny)
    if order > 1:
        np.maximum(out, 0, out=out)
    return out


def brain_window(pos: np.ndarray, full_hi: float) -> list[float]:
    """Window on the brain tissue: centred on the most common foreground intensity (parenchyma), so
    that grey and white matter spread over the grey levels and CSF / fat are near white."""
    if pos.size < 1000 or full_hi <= 0:
        return [0.0, float(full_hi)]
    v = pos[pos > 0.12 * full_hi]
    if v.size < 1000:
        return [0.0, float(full_hi)]
    if v.size > 2_000_000:
        v = v[:: v.size // 2_000_000 + 1]
    h, e = np.histogram(v, bins=200, range=(0.0, float(full_hi)))
    h = np.convolve(h.astype(float), np.ones(7) / 7, mode="same")
    peak = 0.5 * (e[int(np.argmax(h))] + e[int(np.argmax(h)) + 1])
    lo, hi = 0.4 * peak, min(1.8 * peak, float(full_hi))
    return [float(lo), float(max(hi, lo + 1e-6))]


def quantize(vol: np.ndarray) -> tuple[bytes, dict]:
    """uint8 (x fastest, then y, then z) plus the offset/scale to recover values and a default window."""
    v = vol[np.isfinite(vol)]
    pos = v[v > 0]
    lo = 0.0
    hi = float(np.percentile(pos, 99.8)) if pos.size else 1.0
    scale = (hi - lo) / 255.0 or 1.0
    q = np.clip(np.round((vol - lo) / scale), 0, 255).astype(np.uint8)
    w_hi = float(np.percentile(pos, 99.0)) if pos.size else hi
    meta = {"offset": lo, "scale": scale, "window": [lo, w_hi], "window_brain": brain_window(pos, w_hi)}
    return gzip.compress(np.ascontiguousarray(q.transpose(2, 1, 0)).tobytes(), 4), meta


def write_layer(cache_dir: Path, key: str, vol: np.ndarray):
    cache_dir.mkdir(parents=True, exist_ok=True)
    data, meta = quantize(vol)
    (cache_dir / f"{key}.u8.gz").write_bytes(data)
    (cache_dir / f"{key}.json").write_text(json.dumps(meta), encoding="utf-8")
    return meta


def read_layer(cache_dir: Path, key: str) -> tuple[bytes, dict] | None:
    p = cache_dir / f"{key}.u8.gz"
    m = cache_dir / f"{key}.json"
    if not p.exists() or not m.exists():
        return None
    return gzip.decompress(p.read_bytes()), json.loads(m.read_text(encoding="utf-8"))
