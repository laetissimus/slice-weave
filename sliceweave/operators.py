"""Slice acquisition model and super-resolution solver.

Forward model for one slice pixel p (Jiang et al. 2007; Kuklisova-Murgasova et al. 2012;
Ebner et al. 2020):

    y_p = s_slice * (A x)_p + noise,     (A x)_p = sum_v PSF_p(v) x_v

where x is the isotropic volume, the PSF is a 3-D Gaussian oriented with the slice
(FWHM = 1.2 x in-plane pixel size in-plane, FWHM = slice thickness through-plane) and the
pixel position includes the slice's rigid motion correction.

Implementation: slices are grouped by orientation (within ``group_deg`` degrees). For a group
the PSF is applied to the whole volume in the Fourier domain (exact for any orientation),
then sampled at the pixel centres with trilinear interpolation (a sparse matrix):

    A_g x = P_g F^-1[ G_g . F[x] ],      A_g^T r = F^-1[ G_g . F[P_g^T r] ]

The super-resolution volume minimises

    sum_p w_p (y_p / s_p - (A x)_p)^2 + lambda * ||grad x||^2      (first-order Tikhonov)

or, with ``huber`` on, an edge-preserving Huber penalty on the gradient (iteratively
reweighted), solved with preconditioned conjugate gradients.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from scipy import fft as sfft
from scipy import sparse

from .geometry import Grid, apply

FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))  # 2.3548
WORKERS = -1


@dataclass
class Group:
    stack: int
    slices: list                # slice indices in the stack
    n_pix: int
    P: sparse.csr_matrix        # (n_pix, n_vox)
    G: np.ndarray               # rfft-domain PSF transfer function (float32)
    kappa: float                # approx. sum of squared PSF weights per pixel (for the preconditioner)
    pix_slice: np.ndarray       # slice index of every pixel (int)
    pix_ij: np.ndarray          # (n_pix, 2) pixel indices (i, j) in the stack
    y: np.ndarray               # observed values (float32)
    cov_vox: np.ndarray         # 3x3 PSF covariance in voxel units


@dataclass
class Operator:
    grid: Grid
    groups: list = field(default_factory=list)

    # ------------------------------------------------------------------ build
    @classmethod
    def build(cls, grid: Grid, stacks, transforms, masks, include=None, group_deg: float = 12.0,
              psf_inplane: float = 1.2, psf_through: float = 1.0) -> "Operator":
        """``transforms[i]``: (nz, 4, 4) world rigid corrections for stack i;
        ``masks[i]``: (nx, ny, nz) bool pixels to use; ``include[i][s]``: slice used (bool)."""
        op = cls(grid)
        shape = grid.shape
        kfreq = _freq_grid(shape)
        for si, st in enumerate(stacks):
            dx, dy = st.pixel
            normals, rows = [], []
            used = []
            for s in range(st.nz):
                if include is not None and not include[si][s]:
                    continue
                if not masks[si][:, :, s].any():
                    continue
                M = transforms[si][s][:3, :3] @ st.slice_aff[s][:3, :3]
                u = M[:, 0] / np.linalg.norm(M[:, 0])
                n = np.cross(M[:, 0], M[:, 1])
                n /= np.linalg.norm(n)
                normals.append(n)
                rows.append(u)
                used.append(s)
            if not used:
                continue
            # greedy clustering of slice orientations
            clusters: list[list[int]] = []
            cdirs: list[np.ndarray] = []
            cos_t = np.cos(np.deg2rad(group_deg))
            for idx, n in enumerate(normals):
                for c, d in enumerate(cdirs):
                    if abs(float(n @ d)) >= cos_t:
                        clusters[c].append(idx)
                        m = np.sum([normals[q] * np.sign(normals[q] @ d) for q in clusters[c]], axis=0)
                        cdirs[c] = m / np.linalg.norm(m)
                        break
                else:
                    clusters.append([idx])
                    cdirs.append(n)
            for c, members in enumerate(clusters):
                n = cdirs[c]
                u = np.mean([rows[q] * np.sign(rows[q] @ rows[members[0]]) for q in members], axis=0)
                u = u - (u @ n) * n
                u /= np.linalg.norm(u)
                v = np.cross(n, u)
                sx = psf_inplane * dx / FWHM
                sy = psf_inplane * dy / FWHM
                sz = psf_through * st.thickness / FWHM
                Q = np.stack([u, v, n], axis=1)
                cov_w = Q @ np.diag([sx ** 2, sy ** 2, sz ** 2]) @ Q.T
                g = op._group(grid, st, si, [used[q] for q in members], transforms[si], masks[si], cov_w, kfreq)
                if g is not None:
                    op.groups.append(g)
        return op

    @staticmethod
    def _group(grid, st, si, slices, T, mask, cov_w, kfreq):
        shape = grid.shape
        h = grid.spacing
        Rg = grid.R
        cov_v = Rg.T @ cov_w @ Rg / h ** 2
        # trilinear interpolation already smooths by ~1/6 voxel^2 per axis
        cov_v = cov_v - np.eye(3) / 6.0
        w, V = np.linalg.eigh(cov_v)
        w = np.maximum(w, 0.02)
        cov_v = V @ np.diag(w) @ V.T
        G = _gauss_tf(kfreq, cov_v)
        kappa = float(min(1.0, 1.0 / ((4 * np.pi) ** 1.5 * np.sqrt(np.prod(w + 1 / 6.0)))))
        pts, ys, sl, ij = [], [], [], []
        for s in slices:
            ii, jj = np.nonzero(mask[:, :, s])
            P = apply(T[s] @ st.slice_aff[s], np.stack([ii, jj, np.zeros_like(ii)], axis=1))
            pts.append(P)
            ys.append(st.data[ii, jj, s])
            sl.append(np.full(ii.size, s, np.int32))
            ij.append(np.stack([ii, jj], axis=1).astype(np.int32))
        pts = np.concatenate(pts)
        q = grid.world_to_vox(pts)
        P, keep = interp_matrix(q, shape)
        if not keep.any():
            return None
        return Group(si, sorted(set(int(x) for x in np.concatenate(sl)[keep])), int(keep.sum()), P, G, kappa,
                     np.concatenate(sl)[keep], np.concatenate(ij)[keep], np.concatenate(ys)[keep].astype(np.float32),
                     cov_v)

    # ---------------------------------------------------------------- apply
    @property
    def n_pix(self) -> int:
        return sum(g.n_pix for g in self.groups)

    def blurred(self, x: np.ndarray, g: Group, extra_sigma_vox: float = 0.0, X=None) -> np.ndarray:
        X = sfft.rfftn(x, workers=WORKERS) if X is None else X
        G = g.G
        if extra_sigma_vox > 0:
            G = G * _gauss_tf(_freq_grid(self.grid.shape), np.eye(3) * extra_sigma_vox ** 2)
        return sfft.irfftn(X * G, s=self.grid.shape, workers=WORKERS).astype(np.float32, copy=False)

    def forward(self, x: np.ndarray) -> list[np.ndarray]:
        X = sfft.rfftn(x, workers=WORKERS)
        out = []
        for g in self.groups:
            xb = sfft.irfftn(X * g.G, s=self.grid.shape, workers=WORKERS)
            out.append((g.P @ xb.ravel()).astype(np.float32, copy=False))
        return out

    def adjoint(self, rs: list[np.ndarray]) -> np.ndarray:
        shape = self.grid.shape
        acc = None
        for g, r in zip(self.groups, rs):
            z = (g.P.T @ r.astype(np.float32, copy=False)).reshape(shape)
            Z = sfft.rfftn(z, workers=WORKERS) * g.G
            acc = Z if acc is None else acc + Z
        if acc is None:
            return np.zeros(shape, np.float32)
        return sfft.irfftn(acc, s=shape, workers=WORKERS).astype(np.float32, copy=False)

    def coverage(self, weights: list[np.ndarray]) -> np.ndarray:
        """A^T W 1: how much (weighted) data each voxel receives."""
        return self.adjoint(weights)

    def precond_diag(self, weights: list[np.ndarray]) -> np.ndarray:
        shape = self.grid.shape
        acc = None
        for g, w in zip(self.groups, weights):
            z = (g.P.T @ (w * g.kappa).astype(np.float32)).reshape(shape)
            Z = sfft.rfftn(z, workers=WORKERS) * g.G
            acc = Z if acc is None else acc + Z
        d = sfft.irfftn(acc, s=shape, workers=WORKERS).astype(np.float32)
        return np.maximum(d, 0)


# ---------------------------------------------------------------- helpers
def _freq_grid(shape):
    fx = sfft.fftfreq(shape[0]).astype(np.float32)
    fy = sfft.fftfreq(shape[1]).astype(np.float32)
    fz = sfft.rfftfreq(shape[2]).astype(np.float32)
    return fx, fy, fz


def _gauss_tf(kfreq, cov) -> np.ndarray:
    """Fourier transform of a normalised Gaussian with covariance ``cov`` (voxel units)."""
    fx, fy, fz = kfreq
    c = np.asarray(cov, np.float64)
    a = -2.0 * np.pi ** 2
    q = (c[0, 0] * fx[:, None, None] ** 2 + c[1, 1] * fy[None, :, None] ** 2 + c[2, 2] * fz[None, None, :] ** 2
         + 2 * c[0, 1] * fx[:, None, None] * fy[None, :, None]
         + 2 * c[0, 2] * fx[:, None, None] * fz[None, None, :]
         + 2 * c[1, 2] * fy[None, :, None] * fz[None, None, :])
    return np.exp(a * q).astype(np.float32)


def interp_matrix(q: np.ndarray, shape) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Trilinear interpolation matrix for points q (voxel coordinates). Points with less than
    half of their interpolation weight inside the grid are dropped (``keep`` mask)."""
    nx, ny, nz = shape
    q = np.asarray(q, np.float64)
    i0 = np.floor(q).astype(np.int64)
    f = (q - i0).astype(np.float32)
    rows, cols, vals = [], [], []
    wsum = np.zeros(len(q), np.float32)
    corner = []
    for dx in (0, 1):
        wx = f[:, 0] if dx else 1 - f[:, 0]
        ix = i0[:, 0] + dx
        for dy in (0, 1):
            wy = f[:, 1] if dy else 1 - f[:, 1]
            iy = i0[:, 1] + dy
            for dz in (0, 1):
                wz = f[:, 2] if dz else 1 - f[:, 2]
                iz = i0[:, 2] + dz
                w = wx * wy * wz
                ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & (iz >= 0) & (iz < nz)
                w = np.where(ok, w, 0).astype(np.float32)
                wsum += w
                corner.append((np.where(ok, (ix * ny + iy) * nz + iz, 0), w))
    keep = wsum >= 0.5
    idx = np.nonzero(keep)[0]
    remap = np.full(len(q), -1, np.int64)
    remap[idx] = np.arange(idx.size)
    for c, w in corner:
        sel = keep & (w > 0)
        rows.append(remap[sel])
        cols.append(c[sel])
        vals.append(w[sel] / wsum[sel])
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    vals = np.concatenate(vals).astype(np.float32)
    P = sparse.csr_matrix((vals, (rows, cols)), shape=(idx.size, nx * ny * nz), dtype=np.float32)
    return P, keep


# ------------------------------------------------------------------ regulariser
def grad_norm_op(x: np.ndarray, w=None) -> np.ndarray:
    """D^T W D x with forward differences and Neumann boundaries (W per axis or None)."""
    out = np.zeros_like(x)
    for ax in range(3):
        g = np.diff(x, axis=ax)
        if w is not None:
            g = g * w[ax]
        sl_hi = [slice(None)] * 3
        sl_lo = [slice(None)] * 3
        sl_hi[ax] = slice(1, None)
        sl_lo[ax] = slice(None, -1)
        out[tuple(sl_hi)] += g
        out[tuple(sl_lo)] -= g
    return out


def huber_weights(x: np.ndarray, delta: float) -> list[np.ndarray]:
    out = []
    for ax in range(3):
        g = np.abs(np.diff(x, axis=ax))
        out.append((1.0 / np.maximum(1.0, g / delta)).astype(np.float32))
    return out


# ------------------------------------------------------------------ solver
def solve(op: Operator, ys: list[np.ndarray], ws: list[np.ndarray], x0: np.ndarray | None = None,
          alpha: float = 0.03, iters: int = 20, huber: bool = False, huber_delta: float | None = None,
          irls: int = 3, tol: float = 1e-3, log=None) -> tuple[np.ndarray, dict]:
    """Preconditioned CG for the super-resolution problem. Returns the volume and statistics."""
    t0 = time.time()
    b = op.adjoint([w * y for w, y in zip(ws, ys)])
    cov = op.coverage(ws)
    diag = op.precond_diag(ws)
    has = cov > 0.05 * (float(cov.max()) if cov.size else 1.0)
    scale = float(np.mean(diag[has])) if has.any() else 1.0
    lam = alpha * scale
    if x0 is None:
        x0 = np.where(cov > 1e-6 * max(float(cov.max()), 1e-12), b / np.maximum(cov, 1e-12), 0).astype(np.float32)
    x = x0.astype(np.float32, copy=True)
    M = (diag + 6.0 * lam).astype(np.float32)
    Minv = 1.0 / np.maximum(M, 1e-12)
    rw = None
    stats = {"lambda": lam, "alpha": alpha, "cg_iterations": 0}
    outer = irls if huber else 1
    if huber and huber_delta is None:
        ref = x[has] if has.any() else x
        huber_delta = 0.05 * float(np.percentile(ref, 99)) if ref.size else 1.0
    for o in range(outer):
        if huber:
            rw = huber_weights(x, huber_delta)

        def H(v):
            return op.adjoint([w * a for w, a in zip(ws, op.forward(v))]) + lam * grad_norm_op(v, rw)

        r = b - H(x)
        z = Minv * r
        p = z.copy()
        rz = float(np.vdot(r, z))
        bn = float(np.sqrt(np.vdot(b, b))) or 1.0
        n_it = iters if o == outer - 1 or not huber else max(5, iters // 2)
        for it in range(n_it):
            Ap = H(p)
            pAp = float(np.vdot(p, Ap))
            if pAp <= 0:
                break
            a = rz / pAp
            x += a * p
            r -= a * Ap
            rel = float(np.sqrt(np.vdot(r, r))) / bn
            stats["cg_iterations"] += 1
            if rel < tol:
                break
            z = Minv * r
            rz_new = float(np.vdot(r, z))
            p = z + (rz_new / rz) * p
            rz = rz_new
        stats["residual_rel"] = rel if iters else None
    np.maximum(x, 0, out=x)
    stats["seconds"] = round(time.time() - t0, 1)
    if log:
        log(f"  SRR solve: {stats['cg_iterations']} CG iterations, relative residual {stats.get('residual_rel', 0):.3g}, "
            f"{stats['seconds']} s")
    return x, {"coverage": cov, **stats}
