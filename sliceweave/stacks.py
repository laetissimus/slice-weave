"""2-D multi-slice stacks (one per acquired HASTE series) with per-slice geometry.

A stack keeps every slice's own position (``slice_aff[s]`` maps pixel (i, j, 0, 1) of slice s
to world millimetres, RAS+), so uneven gaps and interleaved packets are never resampled
before reconstruction. ``thickness`` is the nominal slice thickness (it sets the
through-plane point-spread function), which can differ from the slice spacing (gap or overlap).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .geometry import apply


@dataclass
class Stack:
    name: str
    data: np.ndarray            # (nx, ny, nz) float32: i along rows of the image, j along columns
    slice_aff: np.ndarray       # (nz, 4, 4)
    thickness: float
    mask: np.ndarray | None = None  # optional (nx, ny, nz) bool, e.g. a brain mask from another tool
    plane: str = ""             # AX | COR | SAG (from the slice normal)
    meta: dict = field(default_factory=dict)

    # -------------------------------------------------------------- geometry
    @property
    def shape(self):
        return self.data.shape

    @property
    def nz(self) -> int:
        return self.data.shape[2]

    @property
    def pixel(self) -> tuple[float, float]:
        A = self.slice_aff[0]
        return float(np.linalg.norm(A[:3, 0])), float(np.linalg.norm(A[:3, 1]))

    @property
    def resolution(self) -> tuple[float, float]:
        """In-plane resolution along i and j (mm): the acquired resolution when the header
        gives it (images are often interpolated, e.g. a 256 matrix shown at 512), never finer
        than the pixel size."""
        dx, dy = self.pixel
        acq = self.meta.get("acq_res")
        if not acq:
            return dx, dy
        return max(dx, float(acq[0] or 0)), max(dy, float(acq[1] or 0))

    @property
    def normal(self) -> np.ndarray:
        A = self.slice_aff[0]
        n = np.cross(A[:3, 0], A[:3, 1])
        return n / np.linalg.norm(n)

    @property
    def spacing(self) -> float:
        if self.nz < 2:
            return self.thickness
        c = np.array([A[:3, 3] for A in self.slice_aff])
        d = np.abs(np.diff(c @ self.normal))
        return float(np.median(d)) if len(d) else self.thickness

    @property
    def axes(self) -> np.ndarray:
        """Unit row, column and normal directions (columns of a 3x3 matrix)."""
        A = self.slice_aff[0]
        u = A[:3, 0] / np.linalg.norm(A[:3, 0])
        v = A[:3, 1] / np.linalg.norm(A[:3, 1])
        return np.stack([u, v, self.normal], axis=1)

    def volume_affine(self) -> np.ndarray:
        """Affine of the stack as a regular volume (median slice spacing)."""
        A = self.slice_aff[0].copy()
        A[:3, 2] = self.normal * self.spacing
        # order slices along the normal
        return A

    def pixel_world(self, s: int, ii=None, jj=None) -> np.ndarray:
        nx, ny = self.data.shape[:2]
        if ii is None:
            jj, ii = np.meshgrid(np.arange(ny), np.arange(nx))
            ii, jj = ii.ravel(), jj.ravel()
        pts = np.stack([ii, jj, np.zeros_like(ii), ], axis=1).astype(np.float64)
        return apply(self.slice_aff[s], pts)

    def corners_world(self) -> np.ndarray:
        nx, ny, _ = self.data.shape
        out = []
        for s in (0, self.nz - 1):
            for i in (-0.5, nx - 0.5):
                for j in (-0.5, ny - 0.5):
                    out.append(apply(self.slice_aff[s], [[i, j, 0.0]])[0])
        h = self.thickness / 2.0
        n = self.normal
        pts = np.array(out)
        return np.concatenate([pts - n * h, pts + n * h])

    def summary(self) -> dict:
        dx, dy = self.pixel
        return {"name": self.name, "shape": [int(x) for x in self.shape], "pixel_mm": [round(dx, 3), round(dy, 3)],
                "thickness_mm": round(float(self.thickness), 3), "spacing_mm": round(self.spacing, 3),
                "plane": self.plane, "has_mask": self.mask is not None,
                "resolution_mm": [round(v, 3) for v in self.resolution]}

    # --------------------------------------------------------- constructors
    @classmethod
    def from_volume(cls, name: str, data: np.ndarray, affine: np.ndarray, thickness: float | None = None,
                    mask: np.ndarray | None = None, meta: dict | None = None) -> "Stack":
        """A regular 3-D volume (NIfTI). The slice axis is the axis with the largest spacing."""
        data = np.asarray(data, np.float32)
        affine = np.asarray(affine, float)
        meta = dict(meta or {})
        sp = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))
        ax = int(np.argmax(sp)) if np.max(sp) > 1.2 * np.min(sp) else 2
        order = [0, 1, 2] if ax == 2 else [a for a in range(3) if a != ax] + [ax]
        res_axes = meta.pop("acq_res_axes", None)   # acquired resolution per file axis (None = unknown)
        if res_axes is not None:
            ra = [res_axes[order[0]], res_axes[order[1]]]
            if any(r for r in ra):
                meta["acq_res"] = [float(r) if r else 0.0 for r in ra]
        if ax != 2:  # make the slice axis the last one
            data = np.transpose(data, order)
            if mask is not None:
                mask = np.transpose(mask, order)
            P = np.eye(4)
            P[:3, :3] = np.eye(3)[:, order]
            affine = affine @ P
            sp = sp[order]
        nz = data.shape[2]
        sa = np.repeat(affine[None], nz, axis=0)
        sa[:, :3, 3] = affine[:3, 3][None] + np.arange(nz)[:, None] * affine[:3, 2][None]
        sa[:, :3, 2] = 0.0
        n = np.cross(affine[:3, 0], affine[:3, 1])
        sa[:, :3, 2] = n / np.linalg.norm(n)
        th = float(thickness) if thickness else float(sp[2])
        st = cls(name, np.ascontiguousarray(data), sa, th, None if mask is None else np.asarray(mask) > 0,
                 plane_label(n), meta)
        return st

    # ------------------------------------------------------------ resizing
    def crop_to(self, inside_fn, margin_px: int = 2) -> "Stack":
        """Crop rows/columns/slices to the pixels for which ``inside_fn(world_pts)`` is true."""
        nx, ny, nz = self.shape
        keep = np.zeros((nx, ny), bool)
        ks = []
        jj, ii = np.meshgrid(np.arange(ny), np.arange(nx))
        ii, jj = ii.ravel(), jj.ravel()
        for s in range(nz):
            m = inside_fn(self.pixel_world(s, ii, jj)).reshape(nx, ny)
            if m.any():
                keep |= m
                ks.append(s)
        if not ks:
            raise ValueError(f"{self.name}: no pixel inside the reconstruction region")
        ri = np.where(keep.any(axis=1))[0]
        rj = np.where(keep.any(axis=0))[0]
        i0, i1 = max(ri[0] - margin_px, 0), min(ri[-1] + margin_px + 1, nx)
        j0, j1 = max(rj[0] - margin_px, 0), min(rj[-1] + margin_px + 1, ny)
        k0, k1 = ks[0], ks[-1] + 1
        T = np.eye(4)
        T[:3, 3] = [i0, j0, 0]
        sa = np.array([A @ T for A in self.slice_aff[k0:k1]])
        st = Stack(self.name, np.ascontiguousarray(self.data[i0:i1, j0:j1, k0:k1]), sa, self.thickness,
                   None if self.mask is None else self.mask[i0:i1, j0:j1, k0:k1], self.plane, dict(self.meta))
        st.meta["crop"] = [int(i0), int(j0), int(k0)]
        if self.meta.get("acq_order") is not None and len(self.meta["acq_order"]) == nz:
            st.meta["acq_order"] = list(self.meta["acq_order"][k0:k1])
        st.meta["slice_index"] = list(range(k0, k1)) if "slice_index" not in self.meta else \
            self.meta["slice_index"][k0:k1]
        return st

    def binned(self, b: int) -> "Stack":
        """Average b x b pixel blocks in-plane (used when pixels are much smaller than the output voxel)."""
        if b <= 1:
            return self
        nx, ny, nz = self.shape
        mx, my = nx // b, ny // b
        d = self.data[: mx * b, : my * b].reshape(mx, b, my, b, nz).mean(axis=(1, 3)).astype(np.float32)
        m = None
        if self.mask is not None:
            m = self.mask[: mx * b, : my * b].reshape(mx, b, my, b, nz).mean(axis=(1, 3)) >= 0.5
        S = np.array([[b, 0, 0, (b - 1) / 2.0], [0, b, 0, (b - 1) / 2.0], [0, 0, 1, 0], [0, 0, 0, 1]], float)
        sa = np.array([A @ S for A in self.slice_aff])
        st = Stack(self.name, d, sa, self.thickness, m, self.plane, dict(self.meta))
        st.meta["binned"] = int(b)
        return st


def plane_label(normal) -> str:
    ax = int(np.argmax(np.abs(np.asarray(normal, float))))
    return ["SAG", "COR", "AX"][ax]


def robust_max(x: np.ndarray, q: float = 99.5) -> float:
    v = x[x > 0]
    return float(np.percentile(v, q)) if v.size else 1.0


def otsu(values: np.ndarray, bins: int = 256) -> float:
    v = values[np.isfinite(values)]
    if v.size == 0:
        return 0.0
    hist, edges = np.histogram(v, bins=bins)
    p = hist.astype(float) / max(hist.sum(), 1)
    w = np.cumsum(p)
    mu = np.cumsum(p * (edges[:-1] + edges[1:]) / 2)
    mt = mu[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        sb = (mt * w - mu) ** 2 / (w * (1 - w))
    sb[~np.isfinite(sb)] = 0
    return float((edges[:-1] + edges[1:])[int(np.argmax(sb))] / 2)


def head_mask(stack: Stack) -> np.ndarray:
    """Rough foreground (head) mask of a stack: Otsu threshold, closing and hole filling per slice."""
    from scipy import ndimage as ndi
    d = stack.data
    t = otsu(np.clip(d, 0, robust_max(d)).ravel()) * 0.5
    m = d > t
    out = np.zeros_like(m)
    st = ndi.generate_binary_structure(2, 1)
    for s in range(d.shape[2]):
        sl = ndi.binary_closing(m[:, :, s], st, iterations=3)
        sl = ndi.binary_fill_holes(sl)
        sl = ndi.binary_opening(sl, st, iterations=2)
        lab, n = ndi.label(sl)
        if n > 1:
            sizes = ndi.sum(sl, lab, range(1, n + 1))
            sl = lab == (1 + int(np.argmax(sizes)))
        out[:, :, s] = sl
    return out


def packets(st: Stack, n_packets: int = 4) -> list[list[int]]:
    """Groups of slices acquired close together in time. Uses the acquisition order when it is
    known (DICOM acquisition times); otherwise assumes an interleaved acquisition (all even
    slices, then all odd slices), the usual order of single-shot sequences."""
    nz = st.nz
    order = st.meta.get("acq_order")
    if order is not None and len(order) == nz:
        by_time = [int(s) for s in np.argsort(order, kind="stable")]
    else:
        by_time = list(range(0, nz, 2)) + list(range(1, nz, 2))
    n = max(1, min(n_packets, nz // 3))
    return [list(map(int, c)) for c in np.array_split(np.array(by_time), n) if len(c)]
