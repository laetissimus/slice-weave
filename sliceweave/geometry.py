"""Rigid transforms and the reconstruction grid.

Conventions
-----------
* World coordinates are millimeters in RAS+ (NIfTI convention).
* A rigid transform is a 4x4 matrix acting on world points (column vectors).
* A :class:`Grid` maps voxel indices (i, j, k) to world points with an affine whose
  linear part is ``R * h`` (orthonormal directions ``R``, isotropic spacing ``h``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import fft as sfft


# ------------------------------------------------------------------ rotations
def skew(v) -> np.ndarray:
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def rotvec_to_matrix(w) -> np.ndarray:
    w = np.asarray(w, float)
    th = float(np.linalg.norm(w))
    if th < 1e-12:
        return np.eye(3) + skew(w)
    K = skew(w / th)
    return np.eye(3) + np.sin(th) * K + (1.0 - np.cos(th)) * (K @ K)


def matrix_to_rotvec(R) -> np.ndarray:
    R = np.asarray(R, float)
    c = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
    th = float(np.arccos(c))
    if th < 1e-9:
        return np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / 2.0
    if np.pi - th < 1e-6:
        # 180 degrees: axis from the largest diagonal element
        k = int(np.argmax(np.diag(R)))
        v = R[:, k] + np.eye(3)[k]
        return th * v / np.linalg.norm(v)
    v = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    return th * v / (2.0 * np.sin(th))


def euler_deg_to_matrix(rx: float, ry: float, rz: float) -> np.ndarray:
    """Rotation about x, then y, then z (degrees), in the grid's own axes."""
    a, b, c = np.deg2rad([rx, ry, rz])
    Rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    Rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def rigid(R=None, t=None) -> np.ndarray:
    T = np.eye(4)
    if R is not None:
        T[:3, :3] = R
    if t is not None:
        T[:3, 3] = t
    return T


def about_center(R, center, t=None) -> np.ndarray:
    """Rotation R about ``center`` followed by translation t: p -> R (p - c) + c + t."""
    c = np.asarray(center, float)
    tt = c - R @ c + (np.zeros(3) if t is None else np.asarray(t, float))
    return rigid(R, tt)


def apply(T, pts) -> np.ndarray:
    pts = np.asarray(pts, dtype=np.float64)
    return pts @ T[:3, :3].T + T[:3, 3]


def rigid_params(T, center) -> tuple[float, float]:
    """Rotation angle (degrees) and displacement of ``center`` (mm) of a rigid transform."""
    ang = float(np.rad2deg(np.linalg.norm(matrix_to_rotvec(T[:3, :3]))))
    c = np.asarray(center, float)
    disp = float(np.linalg.norm(apply(T, c[None])[0] - c))
    return ang, disp


def mean_rigid(Ts: list[np.ndarray], center) -> np.ndarray:
    """Average of rigid transforms (rotation vectors and center displacements averaged)."""
    if not Ts:
        return np.eye(4)
    c = np.asarray(center, float)
    w = np.mean([matrix_to_rotvec(T[:3, :3]) for T in Ts], axis=0)
    d = np.mean([apply(T, c[None])[0] - c for T in Ts], axis=0)
    return about_center(rotvec_to_matrix(w), c, d)


# ---------------------------------------------------------------------- grids
@dataclass
class Grid:
    shape: tuple
    affine: np.ndarray

    @property
    def spacing(self) -> float:
        return float(np.linalg.norm(self.affine[:3, 0]))

    @property
    def R(self) -> np.ndarray:
        return self.affine[:3, :3] / self.spacing

    @property
    def size(self) -> int:
        return int(np.prod(self.shape))

    @property
    def inv(self) -> np.ndarray:
        return np.linalg.inv(self.affine)

    def world_to_vox(self, pts) -> np.ndarray:
        return apply(self.inv, pts)

    def vox_to_world(self, ijk) -> np.ndarray:
        return apply(self.affine, ijk)

    @property
    def center(self) -> np.ndarray:
        return self.vox_to_world((np.array(self.shape, float) - 1) / 2.0)

    def to_dict(self) -> dict:
        return {"shape": [int(x) for x in self.shape], "affine": np.round(self.affine, 6).tolist(),
                "spacing": round(self.spacing, 4)}

    @classmethod
    def from_dict(cls, d) -> "Grid":
        return cls(tuple(int(x) for x in d["shape"]), np.array(d["affine"], float))

    def key(self) -> str:
        import hashlib
        return hashlib.sha1(repr((self.shape, np.round(self.affine, 4).tolist())).encode()).hexdigest()[:12]


def make_grid(R: np.ndarray, center, half_size_mm, spacing: float, margin_mm: float = 0.0,
              fast: bool = True, max_voxels: int | None = None) -> Grid:
    """Grid with directions R (columns), centered on ``center``, covering +/- half_size + margin."""
    R = np.asarray(R, float)
    half = np.asarray(half_size_mm, float) + margin_mm
    h = float(spacing)
    if max_voxels:
        while np.prod(np.ceil(2 * half / h) + 1) > max_voxels:
            h *= 1.05
    shape = []
    for n in (np.ceil(2 * half / h) + 1).astype(int):
        n = max(int(n), 8)
        shape.append(int(sfft.next_fast_len(n, real=True)) if fast else n)
    shape = tuple(shape)
    A = np.eye(4)
    A[:3, :3] = R * h
    A[:3, 3] = np.asarray(center, float) - (R * h) @ ((np.array(shape, float) - 1) / 2.0)
    return Grid(shape, A)


def nearest_ras_axes(dirs: np.ndarray) -> np.ndarray:
    """Given three orthonormal direction vectors (columns), permute and flip them so that
    column 0 is closest to +R(x), column 1 to +A(y) and column 2 to +S(z)."""
    dirs = np.asarray(dirs, float)
    best, best_score = None, -np.inf
    import itertools
    for perm in itertools.permutations(range(3)):
        M = dirs[:, perm]
        signs = np.sign(np.diag(M))
        signs[signs == 0] = 1
        M = M * signs
        score = float(np.trace(M))
        if score > best_score:
            best, best_score = M, score
    if np.linalg.det(best) < 0:  # keep a proper rotation (cannot happen for near-identity)
        best[:, 2] = -best[:, 2]
    return best
