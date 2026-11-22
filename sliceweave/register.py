"""Rigid stack-to-volume and slice-to-volume registration.

A set of pixels (one slice, or all slices of a stack) with world positions p and values y is
aligned to the current volume by finding a rigid transform D (rotation about the pixels' center
and translation) that minimizes

    sum_p w_p ( y_p - a * V(D p) - b )^2

with a, b re-fitted at every step (so the cost is equivalent to normalized cross-correlation)
and V the volume blurred with the slices' point-spread function (plus extra smoothing at coarse
levels). Gauss-Newton steps with Marquardt damping and step halving; Huber weights on the
residuals make the fit robust to partly corrupted slices.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi

from .geometry import Grid, about_center, apply, rotvec_to_matrix


class Sampler:
    """Trilinear sampling of a volume and of its spatial gradient (world units)."""

    def __init__(self, vol: np.ndarray, grid: Grid):
        self.vol = np.ascontiguousarray(vol, np.float32)
        self.grid = grid
        self.grads = [np.ascontiguousarray(g, np.float32) for g in np.gradient(self.vol)]
        self.M = grid.R / grid.spacing  # world gradient = M @ voxel gradient
        self.hi = np.array(grid.shape, float) - 1.0

    def sample(self, pts: np.ndarray, grad: bool = True):
        q = self.grid.world_to_vox(pts)
        ok = np.all((q >= 0) & (q <= self.hi), axis=1)
        qt = q.T
        v = ndi.map_coordinates(self.vol, qt, order=1, mode="nearest", prefilter=False)
        if not grad:
            return v, None, ok
        g = np.stack([ndi.map_coordinates(gg, qt, order=1, mode="nearest", prefilter=False) for gg in self.grads], axis=1)
        return v, g @ self.M.T, ok


@dataclass
class RegResult:
    D: np.ndarray
    ncc: float
    n: int
    iterations: int
    ok: bool
    ncc_start: float = 0.0


def _fit_ab(v, y, w):
    sw = w.sum()
    mv = (w * v).sum() / sw
    my = (w * y).sum() / sw
    var = (w * (v - mv) ** 2).sum()
    a = ((w * (v - mv) * (y - my)).sum() / var) if var > 0 else 0.0
    b = my - a * mv
    return a, b


def _ncc(v, y, w):
    sw = w.sum()
    mv = (w * v).sum() / sw
    my = (w * y).sum() / sw
    dv, dy = v - mv, y - my
    d = np.sqrt((w * dv * dv).sum() * (w * dy * dy).sum())
    return float((w * dv * dy).sum() / d) if d > 0 else 0.0


def register(sampler: Sampler, pts: np.ndarray, y: np.ndarray, center=None, iters: int = 25,
             huber_k: float | None = 1.5, min_valid: float = 0.4, max_rot_deg: float = 45.0,
             max_trans_mm: float = 40.0) -> RegResult:
    """Return the rigid correction D (world -> world) to apply on top of the pixels' positions."""
    pts = np.asarray(pts, np.float64)
    y = np.asarray(y, np.float64)
    n0 = len(y)
    center = pts.mean(axis=0) if center is None else np.asarray(center, float)
    D = np.eye(4)
    mu = 1e-3

    def evaluate(Dm, grad):
        p = apply(Dm, pts)
        v, g, ok = sampler.sample(p, grad)
        if ok.sum() < max(20, min_valid * n0):
            return None
        vv, yy = v[ok].astype(np.float64), y[ok]
        w = np.ones_like(yy)
        a, b = _fit_ab(vv, yy, w)
        r = yy - a * vv - b
        if huber_k:
            s = 1.4826 * np.median(np.abs(r)) + 1e-9
            w = 1.0 / np.maximum(1.0, np.abs(r) / (huber_k * s))
            a, b = _fit_ab(vv, yy, w)
            r = yy - a * vv - b
        nc = _ncc(vv, yy, w)
        # cost = 1 - (signed) weighted NCC: anti-correlated alignments are never preferred
        cost = 1.0 - nc
        return {"p": p[ok], "v": vv, "g": None if g is None else g[ok], "r": r, "w": w, "a": a, "b": b,
                "cost": cost, "ncc": nc, "n": int(ok.sum())}

    cur = evaluate(D, True)
    if cur is None:
        return RegResult(D, 0.0, 0, 0, False, 0.0)
    ncc0 = cur["ncc"]
    it = 0
    for it in range(1, iters + 1):
        c = apply(D, center[None])[0]
        rel = cur["p"] - c
        gw = cur["g"]
        a = cur["a"]
        J = np.empty((len(rel), 8))
        J[:, 0:3] = a * np.cross(rel, gw)
        J[:, 3:6] = a * gw
        J[:, 6] = cur["v"]
        J[:, 7] = 1.0
        w = cur["w"]
        JtW = J.T * w
        H = JtW @ J
        gvec = JtW @ cur["r"]
        Hd = H + mu * np.diag(np.diag(H) + 1e-12)
        try:
            delta = np.linalg.solve(Hd, gvec)
        except np.linalg.LinAlgError:
            break
        accepted = False
        step = 1.0
        for _ in range(5):
            om, tr = delta[:3] * step, delta[3:6] * step
            Dn = about_center(rotvec_to_matrix(om), c, tr) @ D
            new = evaluate(Dn, True)
            if new is not None and new["cost"] <= cur["cost"] + 1e-12:
                accepted = True
                break
            step *= 0.5
        if not accepted:
            mu *= 10.0
            if mu > 1e4:
                break
            continue
        D, cur = Dn, new
        mu = max(mu / 3.0, 1e-6)
        if np.linalg.norm(om) < 2e-5 and np.linalg.norm(tr) < 2e-3:
            break
    # sanity: refuse implausible results
    ang = np.rad2deg(np.linalg.norm(_rotvec(D)))
    disp = np.linalg.norm(apply(D, center[None])[0] - center)
    ok = ang <= max_rot_deg and disp <= max_trans_mm and cur["a"] > 0
    return RegResult(D if ok else np.eye(4), cur["ncc"] if ok else ncc0, cur["n"], it, bool(ok), ncc0)


def _rotvec(D):
    from .geometry import matrix_to_rotvec
    return matrix_to_rotvec(D[:3, :3])


def ncc_at(sampler: Sampler, pts: np.ndarray, y: np.ndarray, D=None) -> float:
    p = pts if D is None else apply(D, pts)
    v, _, ok = sampler.sample(p, grad=False)
    if ok.sum() < 20:
        return -1.0
    return _ncc(v[ok].astype(np.float64), np.asarray(y, np.float64)[ok], np.ones(int(ok.sum())))
