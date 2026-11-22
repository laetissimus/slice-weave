"""Reconstruction pipeline: stacks -> motion-corrected, super-resolved isotropic volume.

Steps (each reported to the progress callback):

1. prepare     crop stacks to the region of interest, bin very small pixels, normalize
               intensities, estimate noise, build masks
2. stacks      register every stack to the reference stack (rigid, coarse to fine)
3. reconstruct first super-resolution volume from all stacks
4. motion      cycles of: slice-to-volume (or stack-to-volume) registration, slice
               intensity matching, bias-field estimation, outlier-slice rejection, SRR
5. final       final SRR with more iterations; coverage and region masking;
               "average" baseline volume (PSF-weighted mean, no deconvolution)

See docs/METHODS.md for the model and the references.
"""
from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field

import numpy as np
from scipy import ndimage as ndi

from . import denoise
from . import operators as ops
from .geometry import (Grid, about_center, apply, make_grid, matrix_to_rotvec, mean_rigid, nearest_ras_axes,
                       rigid, rigid_params, rotvec_to_matrix)
from .register import Sampler, ncc_at, register
from .stacks import Stack, head_mask, packets

STEPS = [("prepare", "Prepare stacks"), ("stacks", "Register stacks"), ("reconstruct", "First reconstruction"),
         ("motion", "Motion correction and outlier rejection"), ("final", "Final reconstruction")]


# ---------------------------------------------------------------------- settings
MODES = ("adult", "fetal", "flair")
MODE_LABELS = {"adult": "T2 single-shot", "fetal": "fetal T2 single-shot", "flair": "FLAIR"}
MODE_NAMES = {"adult": "adult/pediatric", "flair": "2-D FLAIR (adult/pediatric)", "fetal": "fetal"}
FLAIR_SLICE_MSG = ("2-D FLAIR is a multi-shot acquisition: each slice is built from echoes collected over the "
                   "whole series, so it has no single position to correct. Use motion correction 'between series' "
                   "(or none); exclude a series that is degraded by motion.")


@dataclass
class Params:
    mode: str = "adult"              # adult (single-shot T2) | fetal | flair (2-D multi-shot FLAIR, adult/pediatric)
    resolution: float = 1.0          # output voxel size (mm)
    reference: int = 0               # index into the stack list
    motion: str = "stack"            # none | stack | slice
    cycles: int = 2                  # motion-correction cycles
    alpha: float | None = None       # regularization weight; None = automatic from the noise level
    regularizer: str = "huber"       # huber (edge-preserving) | tikhonov
    outlier_rejection: bool = True
    bias_correction: bool = True
    ncc_threshold: float = 0.6       # slices whose NCC with the simulated slice is lower are excluded
    final_iterations: int = 30
    cycle_iterations: int = 12
    group_deg: float = 12.0
    psf_through: float = 1.0         # scales the through-plane PSF width (1 = FWHM equal to the slice thickness)
    max_voxels: int = 30_000_000
    mask_output: bool = True
    temporal_smoothing: bool = True  # median of each slice's motion with its neighbors in acquisition time
    denoise: bool = True  # non-local means on the final volume (with a weaker final weight)
    denoise_strength: float = 1.5

    @classmethod
    def preset(cls, mode: str, **kw) -> "Params":
        if mode == "fetal":
            # denoising off: on the fetal phantom it lowered SSIM (docs/VALIDATION.md, "Denoising")
            p = cls(mode="fetal", resolution=0.8, motion="slice", cycles=6, ncc_threshold=0.5,
                    cycle_iterations=10, final_iterations=30, denoise=False)
        elif mode == "flair":
            # multi-shot: a slice is not a snapshot, so only motion between series is corrected
            # denoising off: on the FLAIR phantom it lowered the contrast of 2-3 mm lesions
            p = cls(mode="flair", resolution=1.0, motion="stack", cycles=1, ncc_threshold=0.6,
                    temporal_smoothing=False, denoise=False)
        else:
            p = cls(mode="adult", resolution=1.0, motion="stack", cycles=1, ncc_threshold=0.6)
        for k, v in kw.items():
            if v is not None and hasattr(p, k):
                setattr(p, k, v)
        return p

    def to_dict(self):
        return asdict(self)


@dataclass
class ROI:
    """Reconstruction region: a box or ellipsoid centered on ``center`` with half-sizes ``half``
    (mm) along the grid axes ``R``. ``offsets[k]``: where the region is in stack k if the subject
    moved between stacks (world shift added to the center for that stack's initial alignment)."""
    center: np.ndarray
    half: np.ndarray
    R: np.ndarray
    shape: str = "box"
    offsets: dict = field(default_factory=dict)

    def inside(self, pts: np.ndarray, margin: float = 0.0) -> np.ndarray:
        q = (np.asarray(pts) - self.center) @ self.R
        h = self.half + margin
        if self.shape == "ellipsoid":
            return ((q / h) ** 2).sum(axis=1) <= 1.0
        return np.all(np.abs(q) <= h, axis=1)

    def to_dict(self):
        return {"center": np.round(self.center, 3).tolist(), "half": np.round(self.half, 3).tolist(),
                "R": np.round(self.R, 6).tolist(), "shape": self.shape,
                "offsets": {str(k): np.round(np.asarray(v), 3).tolist() for k, v in self.offsets.items()}}

    @classmethod
    def from_dict(cls, d):
        return cls(np.array(d["center"], float), np.array(d["half"], float), np.array(d["R"], float),
                   d.get("shape", "box"), {int(k): np.array(v, float) for k, v in (d.get("offsets") or {}).items()})


def default_roi(stacks: list[Stack], reference: int, mode: str) -> ROI:
    """Adult: bounding box of the head (union of the stacks' foreground). Fetal: an ellipsoid
    in the middle of the reference stack (the user is expected to place it on the brain)."""
    R = nearest_ras_axes(stacks[reference].axes)
    if mode == "fetal":
        ref = stacks[reference]
        if ref.mask is not None and ref.mask.any():
            pts = _mask_points(ref, ref.mask)
            q = pts @ R
            lo, hi = q.min(axis=0), q.max(axis=0)
            return ROI(R @ ((lo + hi) / 2), (hi - lo) / 2 + 4.0, R, "ellipsoid")
        c = ref.corners_world().mean(axis=0)
        return ROI(c, np.array([40.0, 46.0, 38.0]), R, "ellipsoid")
    pts = []
    for st in stacks:
        m = head_mask(st)
        if m.any():
            p = _mask_points(st, m)
            pts.append(p[:: max(1, len(p) // 20000)])
    pts = np.concatenate(pts) if pts else stacks[reference].corners_world()
    q = pts @ R
    lo, hi = np.percentile(q, 0.5, axis=0), np.percentile(q, 99.5, axis=0)
    return ROI(R @ ((lo + hi) / 2), (hi - lo) / 2 + 3.0, R, "box")


def _mask_points(st: Stack, m: np.ndarray) -> np.ndarray:
    ii, jj, ss = np.nonzero(m)
    A = st.slice_aff[ss]
    return np.einsum("nij,nj->ni", A[:, :3, :3], np.stack([ii, jj, np.zeros(len(ii))], 1)) + A[:, :3, 3]


def grid_for(roi: ROI, params: Params, stacks: list[Stack] | None = None) -> Grid:
    margin = 6.0
    if stacks:
        margin = max(margin, 1.2 * max(s.thickness for s in stacks))
    return make_grid(roi.R, roi.center, roi.half, params.resolution, margin_mm=margin, max_voxels=params.max_voxels)


# ---------------------------------------------------------------------- helpers
def noise_sigma(data: np.ndarray, mask: np.ndarray, step=(1.0, 1.0)) -> float:
    """Immerkaer-type estimate: median absolute 2-D Laplacian of the pixels in the mask.
    ``step``: acquired resolution / pixel size. Interpolated images (e.g. a 256 matrix shown at
    512) have correlated noise that the Laplacian under-reads, so they are first sampled at
    about the acquired spacing."""
    if max(step) > 1.05:
        ii = np.unique(np.round(np.arange(0, data.shape[0] - 0.5, max(step[0], 1.0))).astype(int))
        jj = np.unique(np.round(np.arange(0, data.shape[1] - 0.5, max(step[1], 1.0))).astype(int))
        data, mask = data[ii][:, jj], mask[ii][:, jj]
    k = np.array([[1, -2, 1], [-2, 4, -2], [1, -2, 1]], np.float32)
    vals = []
    for s in range(data.shape[2]):
        m = mask[:, :, s]
        if m.sum() < 50:
            continue
        lap = ndi.convolve(data[:, :, s], k, mode="nearest")
        vals.append(lap[ndi.binary_erosion(m, iterations=1)])
    if not vals:
        return 0.0
    v = np.concatenate(vals)
    return float(np.median(np.abs(v)) / (0.6745 * 6.0))


def initial_slice_scale(data: np.ndarray, fg: np.ndarray) -> np.ndarray:
    """Per-slice intensity scale relative to the neighboring slices (running median of the
    slice medians): catches slices with global signal loss before any volume exists."""
    nz = data.shape[2]
    med = np.ones(nz)
    for s in range(nz):
        v = data[:, :, s][fg[:, :, s]]
        v = v[v > 0.1 * (np.percentile(v, 99) if v.size else 1)] if v.size else v
        med[s] = float(np.median(v)) if v.size > 50 else np.nan
    if np.all(np.isnan(med)):
        return np.ones(nz)
    med = np.where(np.isnan(med), np.nanmedian(med), med)
    ref = ndi.median_filter(med, size=min(5, nz), mode="nearest")
    return np.clip(med / np.maximum(ref, 1e-6), 0.3, 3.0)


# Weight of the final volume relative to the noise-calibrated weight (which the registration volumes
# keep). The noise calibration was done on a formula phantom with smooth anatomy. FLAIR: on the
# FLAIR phantom about a third of it kept twice the contrast of 2-3 mm lesions. Adult/pediatric
# T2: on real anatomy (pseudo-T2 from a 3-D scan) a quarter of it kept the sulci and the cortical
# ribbon best at every noise level tested (docs/VALIDATION.md, "Gray-white matter and sulci").
FINAL_ALPHA_FACTOR = {"adult": 0.25, "flair": 0.35, "fetal": 1.0}
FLAIR_ALPHA_FACTOR = FINAL_ALPHA_FACTOR["flair"]
# With the denoising step (non-local means on the final volume) the final weight is lowered further:
# the penalty then smooths less across the cortex and the sulci, and non-local means, which averages
# only voxels with similar surroundings, removes the noise (docs/VALIDATION.md, "Denoising").
DENOISE_ALPHA_FACTOR = 0.4


def auto_alpha(rel_noise: float, regularizer: str, mode: str | None = None) -> float:
    """Regularization weight from the relative noise level (noise sigma / median tissue signal).
    Calibrated on simulated stacks with known ground truth (docs/VALIDATION.md, "Regularization"):
    the weight that maximized PSNR and SSIM grew about linearly with the noise level. With a mode,
    the weight of that mode's final volume (FINAL_ALPHA_FACTOR); without, the registration weight."""
    if regularizer == "huber":
        a = 0.1 + 65.0 * rel_noise
    else:
        a = 0.1 + 40.0 * rel_noise
    return float(a * FINAL_ALPHA_FACTOR.get(mode, 1.0))


def blur_volume(x: np.ndarray, grid: Grid, cov_world: np.ndarray, extra_sigma_mm: float = 0.0) -> np.ndarray:
    from scipy import fft as sfft
    h = grid.spacing
    cov_v = grid.R.T @ cov_world @ grid.R / h ** 2 + np.eye(3) * (extra_sigma_mm / h) ** 2 - np.eye(3) / 6.0
    w, V = np.linalg.eigh(cov_v)
    cov_v = V @ np.diag(np.maximum(w, 0.02)) @ V.T
    G = ops._gauss_tf(ops._freq_grid(grid.shape), cov_v)
    X = sfft.rfftn(x, workers=ops.WORKERS)
    return sfft.irfftn(X * G, s=grid.shape, workers=ops.WORKERS).astype(np.float32)


def psf_cov(st: Stack, Rmot: np.ndarray, psf_through: float = 1.0) -> np.ndarray:
    dx, dy = st.resolution
    Q = Rmot @ st.axes
    sx, sy = 1.2 * dx / ops.FWHM, 1.2 * dy / ops.FWHM
    sz = psf_through * st.thickness / ops.FWHM
    return Q @ np.diag([sx ** 2, sy ** 2, sz ** 2]) @ Q.T


# ---------------------------------------------------------------------- pipeline
@dataclass
class Result:
    volume: np.ndarray
    grid: Grid
    average: np.ndarray
    coverage: np.ndarray
    report: dict
    transforms: list
    include: list
    stacks: list


class Reconstruction:
    def __init__(self, stacks: list[Stack], roi: ROI, params: Params, include: list[bool] | None = None,
                 log=None, progress=None, cancel=None):
        self.all_stacks = stacks
        self.use = [i for i in range(len(stacks)) if include is None or include[i]]
        if params.reference not in self.use:
            raise ValueError("the reference stack must be included")
        if params.mode == "flair" and params.motion == "slice":
            raise ValueError(FLAIR_SLICE_MSG)
        self.roi = roi
        self.p = params
        self._log = log or (lambda m: None)
        self._progress = progress or (lambda step, frac, msg="": None)
        self.cancel = cancel or (lambda: False)
        self.timing: dict = {}
        self.events: list = []

    def log(self, msg):
        self._log(msg)

    def check_cancel(self):
        if self.cancel():
            raise RuntimeError("canceled")

    # -------------------------------------------------------------- prepare
    def prepare(self):
        p = self.p
        self._progress("prepare", 0.0, "cropping stacks to the region")
        self.grid = grid_for(self.roi, p, [self.all_stacks[i] for i in self.use])
        self.log(f"grid {self.grid.shape} at {self.grid.spacing:.3f} mm ({self.grid.size / 1e6:.1f} M voxels)")
        if abs(self.grid.spacing - p.resolution) > 1e-3:
            self.log(f"NOTE: voxel size raised to {self.grid.spacing:.3f} mm to stay under {p.max_voxels / 1e6:.0f} M voxels")
        h = self.grid.spacing
        self.stacks, self.norm, self.noise = [], [], []
        self.T, self.scale, self.bias, self.pixw, self.inc, self.fg = [], [], [], [], [], []
        ref_pos = self.use.index(p.reference)
        self.ref = ref_pos
        margin = max(1.5 * max(self.all_stacks[i].thickness for i in self.use), 8.0)
        for n, i in enumerate(self.use):
            st0 = self.all_stacks[i]
            off = np.asarray(self.roi.offsets.get(i, np.zeros(3)), float)
            roi_k = ROI(self.roi.center + off, self.roi.half, self.roi.R, self.roi.shape)
            st = st0.crop_to(lambda pts: roi_k.inside(pts, margin=margin), margin_px=2)
            b = max(1, int(h / min(st.pixel) + 0.1))
            st_full, raw_full = st, st.data   # unbinned pixels (raw values) for the noise estimate
            st = st.binned(b)
            if p.mode == "fetal":
                fg = st.mask.copy() if st.mask is not None else np.ones(st.shape, bool)
            else:
                fg = ndi.binary_dilation(head_mask(st), iterations=2)
                if st.mask is not None:
                    fg &= st.mask
            vals = st.data[fg & (st.data > 0)]
            f = float(np.median(vals[vals > np.percentile(vals, 20)])) if vals.size > 100 else float(np.median(st.data) or 1)
            st.data = (st.data / f).astype(np.float32)
            self.stacks.append(st)
            self.norm.append(f)
            self.fg.append(fg)
            # noise on the unbinned pixels, sampled at about the acquired resolution
            fg_full = np.zeros(st_full.shape, bool)
            up = np.repeat(np.repeat(fg, b, axis=0), b, axis=1)
            fg_full[:up.shape[0], :up.shape[1]] = up
            step = (st_full.resolution[0] / st_full.pixel[0], st_full.resolution[1] / st_full.pixel[1])
            self.noise.append(noise_sigma(raw_full / f, fg_full, step))
            T0 = rigid(t=-off) if np.any(off) else np.eye(4)
            self.T.append(np.repeat(T0[None], st.nz, axis=0))
            sc0 = initial_slice_scale(st.data, fg)
            self.scale.append(np.ones(st.nz))
            self.bias.append(np.ones(st.shape, np.float32))
            self.pixw.append(np.ones(st.shape, np.float32))
            # slices much darker or brighter than their neighbors start excluded (re-checked later)
            self.inc.append((sc0 > 0.6) & (sc0 < 1.6))
            self.log(f"stack {st.name}: {st.shape[0]}x{st.shape[1]}x{st.nz} after crop"
                     + (f", binned {b}x{b}" if b > 1 else "") + f", pixel {st.pixel[0]:.2f} mm"
                     + (f" (acquired {st.resolution[0]:.2f} x {st.resolution[1]:.2f} mm)"
                        if max(st.resolution) > 1.05 * max(st.pixel) else "") + ", thickness "
                     f"{st.thickness:.2f} mm, noise {self.noise[-1]:.3f} of median signal")
            self._progress("prepare", (n + 1) / len(self.use), st.name)
        rel = float(np.median(self.noise))
        self.rel_noise = rel
        self.alpha = p.alpha if p.alpha is not None else auto_alpha(rel, p.regularizer, p.mode)
        if p.alpha is None and p.denoise:
            self.alpha *= DENOISE_ALPHA_FACTOR
        # volumes used as registration targets keep the noise-calibrated weight (a noisier target
        # makes the motion estimates worse); only the final volume uses the mode's weight
        self.alpha_motion = p.alpha if p.alpha is not None else auto_alpha(rel, p.regularizer)
        self.log(f"relative noise {rel:.3f}; regularization {p.regularizer}, alpha {self.alpha:.3g}"
                 + (" (automatic" + (f"; {self.alpha_motion:.3g} during motion correction"
                                     if self.alpha_motion != self.alpha else "") + ")" if p.alpha is None else ""))

    # ---------------------------------------------------------------- masks
    def masks(self) -> list[np.ndarray]:
        out = []
        for k, st in enumerate(self.stacks):
            nx, ny, nz = st.shape
            m = np.zeros(st.shape, bool)
            jj, ii = np.meshgrid(np.arange(ny), np.arange(nx))
            base = np.stack([ii.ravel(), jj.ravel(), np.zeros(ii.size)], axis=1)
            for s in range(nz):
                P = apply(self.T[k][s] @ st.slice_aff[s], base)
                m[:, :, s] = self.roi.inside(P).reshape(nx, ny)
            out.append(m & self.fg[k])
        return out

    # ------------------------------------------------------------ operators
    def build(self, masks, include=None):
        inc = include if include is not None else self.inc
        return ops.Operator.build(self.grid, self.stacks, self.T, masks, include=inc, group_deg=self.p.group_deg,
                                  psf_through=self.p.psf_through)

    def data_and_weights(self, op) -> tuple[list, list]:
        ys, ws = [], []
        for g in op.groups:
            k = g.stack
            i, j = g.pix_ij[:, 0], g.pix_ij[:, 1]
            s = g.pix_slice
            f = self.scale[k][s] * self.bias[k][i, j, s]
            ys.append((g.y / f).astype(np.float32))
            ws.append((self.pixw[k][i, j, s] * f * f).astype(np.float32))
        return ys, ws

    def srr(self, op, x0=None, iters=20, final=False):
        ys, ws = self.data_and_weights(op)
        a = self.alpha if final else self.alpha_motion
        x, st = ops.solve(op, ys, ws, x0=x0, alpha=a, iters=iters, huber=self.p.regularizer == "huber",
                          log=self.log)
        return x, st

    # ---------------------------------------------------------- registration
    def smooth_data(self, k: int, sigma_mm: float) -> np.ndarray:
        """Stack data (bias- and scale-corrected) smoothed in-plane by ``sigma_mm``: at coarse
        registration levels the slices are smoothed like the volume."""
        st = self.stacks[k]
        d = st.data / (self.bias[k] * self.scale[k][None, None, :])
        if sigma_mm <= 0:
            return d
        dx, dy = st.pixel
        return ndi.gaussian_filter(d, (sigma_mm / dx, sigma_mm / dy, 0)).astype(np.float32)

    def register_stacks(self, x, levels, skip_reference=True, max_pts=80000):
        rng = np.random.default_rng(0)
        masks = self.masks()
        out = []
        for k, st in enumerate(self.stacks):
            if skip_reference and k == self.ref:
                continue
            self.check_cancel()
            ii, jj, ss = np.nonzero(masks[k])
            if ii.size < 200:
                self.log(f"  {st.name}: too few pixels in the region for registration")
                continue
            sel = rng.choice(ii.size, size=min(max_pts, ii.size), replace=False)
            ii, jj, ss = ii[sel], jj[sel], ss[sel]
            A = np.einsum("nab,nbc->nac", self.T[k][ss], st.slice_aff[ss])
            P = np.einsum("nij,nj->ni", A[:, :3, :3], np.stack([ii, jj, np.zeros(len(ii))], 1)) + A[:, :3, 3]
            Rm = mean_rigid(list(self.T[k]), self.roi.center)[:3, :3]
            cov = psf_cov(st, Rm, self.p.psf_through)
            Dacc = np.eye(4)
            for sig in levels:
                y = self.smooth_data(k, sig)[ii, jj, ss]
                smp = Sampler(blur_volume(x, self.grid, cov, sig), self.grid)
                r = register(smp, apply(Dacc, P), y, center=self.roi.center, iters=30, max_rot_deg=25,
                             max_trans_mm=25)
                if r.ok:
                    Dacc = r.D @ Dacc
            # accept only if the match at the finest level improved over no correction
            n0, n1 = ncc_at(smp, P, y), ncc_at(smp, P, y, Dacc)
            if n1 <= n0 + 1e-3:
                Dacc = np.eye(4)
                n1 = n0
            self.T[k] = np.einsum("ab,nbc->nac", Dacc, self.T[k])
            ang, disp = rigid_params(Dacc, self.roi.center)
            self.log(f"  {st.name}: rotated {ang:.2f} deg, moved {disp:.2f} mm (NCC {n0:.3f} -> {n1:.3f})")
            out.append({"stack": st.name, "rotation_deg": round(ang, 3), "translation_mm": round(disp, 3),
                        "ncc_before": round(n0, 4), "ncc_after": round(n1, 4)})
        return out

    def register_units(self, x, units, levels, max_rot=20.0, max_trans=15.0, min_gain=1e-3, label="slice",
                       max_pts=40000):
        """Rigid registration of groups of slices ("units": whole packets or single slices) to the
        volume, coarse to fine. An update is kept only if it improves the match at the finest
        level and stays within the limits."""
        rng = np.random.default_rng(1)
        masks = self.masks()
        before = [T.copy() for T in self.T]
        n_done, n_fail, n_kept = 0, 0, 0
        by_stack: dict[int, list] = {}
        for k, sl in units:
            by_stack.setdefault(k, []).append(sl)
        total = len(units)
        for k, unit_list in by_stack.items():
            self.check_cancel()
            st = self.stacks[k]
            Rm = mean_rigid(list(self.T[k]), self.roi.center)[:3, :3]
            cov = psf_cov(st, Rm, self.p.psf_through)
            smps = [(sig, Sampler(blur_volume(x, self.grid, cov, sig), self.grid), self.smooth_data(k, sig))
                    for sig in levels]
            npix = masks[k].sum(axis=(0, 1))
            min_pix = max(60, 0.35 * float(np.median(npix[npix > 0]))) if (npix > 0).any() else 60
            for sl in unit_list:
                # a single slice that only grazes the region does not constrain 6 parameters
                sel = [s for s in sl if npix[s] >= (min_pix if len(sl) == 1 else 60)]
                if not sel:
                    continue
                ii, jj, ss = [], [], []
                for s in sel:
                    a, b = np.nonzero(masks[k][:, :, s])
                    ii.append(a)
                    jj.append(b)
                    ss.append(np.full(a.size, s))
                ii, jj, ss = np.concatenate(ii), np.concatenate(jj), np.concatenate(ss)
                if ii.size > max_pts:
                    pick = rng.choice(ii.size, max_pts, replace=False)
                    ii, jj, ss = ii[pick], jj[pick], ss[pick]
                A = np.einsum("nab,nbc->nac", self.T[k][ss], st.slice_aff[ss])
                P = np.einsum("nij,nj->ni", A[:, :3, :3], np.stack([ii, jj, np.zeros(ii.size)], 1)) + A[:, :3, 3]
                ctr = P.mean(axis=0)
                D = np.eye(4)
                for sig, smp, d in smps:
                    r = register(smp, apply(D, P), d[ii, jj, ss], center=ctr, iters=20, max_rot_deg=max_rot,
                                 max_trans_mm=max_trans)
                    if r.ok:
                        D = r.D @ D
                sig, smp, d = smps[-1]
                y = d[ii, jj, ss]
                n0, n1 = ncc_at(smp, P, y), ncc_at(smp, P, y, D)
                ang, disp = rigid_params(D, ctr)
                if n1 > n0 + min_gain and ang <= max_rot and disp <= max_trans:
                    self.T[k][sel] = np.einsum("ab,nbc->nac", D, self.T[k][sel])
                    n_kept += 1
                else:
                    n_fail += 1
                n_done += 1
                if n_done % 25 == 0:
                    self._progress("motion", None, f"{label} {n_done}/{total}")
        if label == "slice":
            self.temporal_repair()
        # keep the reference stack (on average) where it was: the volume frame does not drift
        rel = [self.T[self.ref][s] @ np.linalg.inv(before[self.ref][s]) for s in range(self.stacks[self.ref].nz)]
        Mi = np.linalg.inv(mean_rigid(rel, self.roi.center))
        for k in range(len(self.T)):
            self.T[k] = np.einsum("ab,nbc->nac", Mi, self.T[k])
        moved = [rigid_params(self.T[k][s] @ np.linalg.inv(before[k][s]), self.roi.center)[1]
                 for k in range(len(self.T)) for s in range(len(self.T[k]))]
        self.log(f"  {label}-to-volume: {n_done} {label}s, {n_kept} updated, {n_fail} kept as they were; "
                 f"median slice displacement {np.median(moved) if moved else 0:.2f} mm")

    def temporal_repair(self, window: int = 2, max_dev_deg: float = 6.0, max_dev_mm: float = 4.0):
        """The fetus moves smoothly most of the time: a slice whose estimated pose differs a lot
        from the poses of the slices acquired just before and after it (median of +/- ``window``
        slices in acquisition order, robust to sudden jumps) is given that median pose."""
        c = self.roi.center
        n_fixed = 0
        for k, st in enumerate(self.stacks):
            order = st.meta.get("acq_order")
            by_time = (list(np.argsort(order, kind="stable")) if order is not None and len(order) == st.nz
                       else list(range(0, st.nz, 2)) + list(range(1, st.nz, 2)))
            par = np.array([np.concatenate([matrix_to_rotvec(self.T[k][s][:3, :3]),
                                            apply(self.T[k][s], c[None])[0] - c]) for s in by_time])
            new = par.copy()
            for t in range(len(by_time)):
                nb = [u for u in range(max(0, t - window), min(len(by_time), t + window + 1)) if u != t]
                if len(nb) < 2:
                    continue
                med = np.median(par[nb], axis=0)
                dev_r = np.rad2deg(np.linalg.norm(par[t, :3] - med[:3]))
                dev_t = np.linalg.norm(par[t, 3:] - med[3:])
                if dev_r > max_dev_deg or dev_t > max_dev_mm:
                    new[t] = med
                    n_fixed += 1
            for t, s in enumerate(by_time):
                if not np.allclose(new[t], par[t]):
                    self.T[k][s] = about_center(rotvec_to_matrix(new[t, :3]), c, new[t, 3:])
        if n_fixed:
            self.log(f"  {n_fixed} slice pose(s) inconsistent with the slices acquired just before and after: "
                     "replaced by their median")

    def smooth_in_time(self):
        """Fetal motion is mostly smooth in time: replace each slice's rigid motion by the median
        (per parameter) of itself and its two neighbors in acquisition order. Sudden jumps are
        kept (a median preserves steps); isolated registration failures are removed."""
        from .geometry import matrix_to_rotvec, rotvec_to_matrix
        for k, st in enumerate(self.stacks):
            order = st.meta.get("acq_order")
            by_time = [int(i) for i in np.argsort(order, kind="stable")] if order is not None and len(order) == st.nz \
                else list(range(0, st.nz, 2)) + list(range(1, st.nz, 2))
            if len(by_time) < 3:
                continue
            ctr = np.array([apply(st.slice_aff[s], [[st.shape[0] / 2, st.shape[1] / 2, 0]])[0] for s in range(st.nz)])
            prm = np.zeros((st.nz, 6))
            for s in range(st.nz):
                T = self.T[k][s]
                prm[s, :3] = matrix_to_rotvec(T[:3, :3])
                prm[s, 3:] = apply(T, ctr[s][None])[0] - ctr[s]
            seq = prm[by_time]
            sm = seq.copy()
            for i in range(1, len(seq) - 1):
                sm[i] = np.median(seq[i - 1:i + 2], axis=0)
            new = np.empty_like(prm)
            new[by_time] = sm
            for s in range(st.nz):
                self.T[k][s] = about_center(rotvec_to_matrix(new[s, :3]), ctr[s], new[s, 3:])

    # ---------------------------------------------------- intensity, outliers
    def update_weights(self, x, op, cycle: int, reject: bool, thr: float):
        sims = op.forward(x)
        ys, _ = self.data_and_weights(op)
        # robust pixel weights (Huber) from the residuals of all included pixels
        res = np.concatenate([y - s for y, s in zip(ys, sims)]) if sims else np.zeros(1)
        sig = 1.4826 * np.median(np.abs(res - np.median(res))) + 1e-6
        slice_stats = {}
        for g, y, sim in zip(op.groups, ys, sims):
            k = g.stack
            i, j, s = g.pix_ij[:, 0], g.pix_ij[:, 1], g.pix_slice
            r = y - sim
            self.pixw[k][i, j, s] = (1.0 / np.maximum(1.0, np.abs(r) / (2.0 * sig))).astype(np.float32)
            order = np.argsort(s, kind="stable")
            s_sorted = s[order]
            bounds = np.flatnonzero(np.diff(s_sorted)) + 1
            for idx in np.split(order, bounds):
                if idx.size == 0:
                    continue
                sl = int(s[idx[0]])
                yy = g.y[idx] / self.bias[k][i[idx], j[idx], sl]
                ss_ = sim[idx]
                den = float(ss_ @ ss_)
                sc = float(yy @ ss_) / den if den > 0 else 1.0
                slice_stats[(k, sl)] = (sc, _ncc(yy, ss_), idx.size, yy, ss_, i[idx], j[idx])
        # intensity matching: slice scale relative to the reference stack
        ref_sc = [v[0] for (k, s), v in slice_stats.items() if k == self.ref]
        norm = float(np.median(ref_sc)) if ref_sc else 1.0
        for (k, s), v in slice_stats.items():
            self.scale[k][s] = np.clip(v[0] / norm, 0.3, 3.0)
        # smooth bias field per stack (normalized convolution of the log intensity ratio)
        if self.p.bias_correction and cycle >= 0:
            for k, st in enumerate(self.stacks):
                num = np.zeros(st.shape, np.float32)
                den = np.zeros(st.shape, np.float32)
                for (kk, s), v in slice_stats.items():
                    if kk != k:
                        continue
                    sc, _, _, yy, ss_, ii, jj = v
                    ok = (yy > 0) & (ss_ > 1e-3)
                    lr = np.log(np.maximum(yy[ok], 1e-6) / np.maximum(self.scale[k][s] * norm * ss_[ok], 1e-6))
                    w = self.pixw[k][ii[ok], jj[ok], s]
                    num[ii[ok], jj[ok], s] = w * lr
                    den[ii[ok], jj[ok], s] = w
                if den.sum() < 100:
                    continue
                dx, dy = st.pixel
                sg = (20.0 / dx, 20.0 / dy, 20.0 / max(st.spacing, 1e-3))
                ns, ds = ndi.gaussian_filter(num, sg), ndi.gaussian_filter(den, sg)
                lb = np.where(ds > 1e-3, ns / np.maximum(ds, 1e-6), 0)
                wmean = float((lb * den).sum() / den.sum())
                self.bias[k] = np.clip(self.bias[k] * np.exp(lb - wmean), 0.5, 2.0).astype(np.float32)
        # outlier slices
        excluded = []
        if reject:
            nccs = np.array([v[1] for v in slice_stats.values()])
            med = float(np.median(nccs)) if nccs.size else 1.0
            mad = float(np.median(np.abs(nccs - med))) * 1.4826 if nccs.size else 0.0
            if self.p.mode != "fetal":  # only clear outliers: whole-brain slices of a still head correlate well
                thr = min(thr, med - 4.0 * mad)
            for (k, s), v in slice_stats.items():
                bad = v[1] < thr or not (0.4 <= self.scale[k][s] <= 2.5)
                self.inc[k][s] = not bad
                if bad:
                    excluded.append((k, s))
        self.slice_stats = {key: (v[0], v[1], v[2]) for key, v in slice_stats.items()}
        if reject:
            self.log(f"  outlier slices: {len(excluded)} excluded")
        return excluded

    # ------------------------------------------------------------------ run
    def run(self) -> Result:
        t_all = time.time()
        p = self.p
        t = time.time()
        self.prepare()
        self.timing["prepare"] = round(time.time() - t, 1)
        stack_reg = []
        # -- initial volume from the reference stack, then stack-to-volume registration
        t = time.time()
        if len(self.stacks) > 1 and p.motion != "none":
            masks = self.masks()
            if p.mode == "fetal":
                # every stack moved during its own acquisition: start from all stacks together
                self._progress("stacks", 0.0, "initial volume from all stacks")
                op = self.build(masks)
            else:
                self._progress("stacks", 0.0, "initial volume from the reference stack")
                only_ref = [np.zeros(st.nz, bool) if k != self.ref else self.inc[k] for k, st in enumerate(self.stacks)]
                op = self.build(masks, include=only_ref)
            ys, ws = self.data_and_weights(op)
            x_ref, _ = ops.solve(op, ys, ws, alpha=max(self.alpha_motion, 0.3), iters=8, huber=False)
            levels = [4.0, 2.0, 1.0] if p.mode == "fetal" else [3.0, 1.5]
            self._progress("stacks", 0.3, "registering stacks")
            stack_reg = self.register_stacks(x_ref, levels)
        self.timing["stacks"] = round(time.time() - t, 1)
        # -- first SRR
        t = time.time()
        self._progress("reconstruct", 0.0, "super-resolution")
        masks = self.masks()
        op = self.build(masks)
        x, _ = self.srr(op, iters=p.cycle_iterations)
        self.timing["reconstruct"] = round(time.time() - t, 1)
        # -- motion-correction cycles
        t = time.time()
        excluded = []
        # each cycle: motion correction against the current volume, then intensity correction and the
        # outlier test, then SRR again. Without motion correction (or with 0 cycles) intensity
        # correction and the outlier test still run once when they are switched on.
        n_motion = max(int(p.cycles), 0) if p.motion in ("stack", "slice") else 0
        n_cycles = max(n_motion, 1 if (p.outlier_rejection or p.bias_correction) else 0)
        slice_levels = {1: [2.0, 1.0], 2: [1.0, 0.5]}
        for c in range(n_cycles):
            self.check_cancel()
            self._progress("motion", c / max(n_cycles, 1), f"cycle {c + 1} of {n_cycles}")
            self.log(f"cycle {c + 1}/{n_cycles}" + ("" if c < n_motion else " (intensity and outlier check only)"))
            if c >= n_motion:
                pass
            elif p.motion == "slice":
                if c == 0:
                    # coarse to fine in the slices too: whole stacks, then packets of slices acquired together
                    if len(self.stacks) > 1:
                        stack_reg = self.register_stacks(x, [2.0, 1.0])
                    for n_pk in (2, 4):
                        units = [(k, pk) for k, st in enumerate(self.stacks) for pk in packets(st, n_pk)]
                        self.register_units(x, units, [2.0, 1.0], max_rot=20, max_trans=15, label="packet")
                elif c == 1:
                    units = [(k, pk) for k, st in enumerate(self.stacks) for pk in packets(st, 8)]
                    self.register_units(x, units, [1.5, 1.0], max_rot=15, max_trans=10, label="packet")
                    units = [(k, [s]) for k, st in enumerate(self.stacks) for s in range(st.nz)]
                    self.register_units(x, units, slice_levels[1], max_rot=10, max_trans=8, min_gain=2e-3,
                                        label="slice")
                else:
                    units = [(k, [s]) for k, st in enumerate(self.stacks) for s in range(st.nz)]
                    self.register_units(x, units, slice_levels.get(c, [0.0]), max_rot=10, max_trans=8,
                                        min_gain=2e-3, label="slice")
                if c >= 1 and p.temporal_smoothing:
                    self.smooth_in_time()
            elif p.motion == "stack" and len(self.stacks) > 1:
                stack_reg = self.register_stacks(x, [1.0, 0.0])
            masks = self.masks()
            op = self.build(masks, include=[np.ones(st.nz, bool) for st in self.stacks])
            thr = p.ncc_threshold
            if p.mode == "fetal":  # stricter as the motion estimates improve (as in NiftyMIC)
                thr = p.ncc_threshold + 0.1 * min(c, 2) * (c > 0)
            excluded = self.update_weights(x, op, c, p.outlier_rejection, thr)
            op = self.build(masks)
            x, _ = self.srr(op, x0=x, iters=p.cycle_iterations)
        self.timing["motion"] = round(time.time() - t, 1)
        # -- final
        t = time.time()
        self._progress("final", 0.0, "final super-resolution")
        op = self.build(self.masks())
        x, st = self.srr(op, x0=x, iters=p.final_iterations, final=True)
        self.denoise_info = None
        if p.denoise:
            t1 = time.time()
            self._progress("final", 0.5, "denoising (non-local means)")
            fg = _support(st["coverage"]) > 0.5
            if fg.sum() > 1000:
                fg &= x > 0.2 * float(np.percentile(x[fg], 95))
            sig = noise_sigma(x, fg)
            x = denoise.nlm3d(x, sig, p.denoise_strength, cancel=self.check_cancel)
            self.denoise_info = {"method": "non-local means, 3x3x3 patches, 5x5x5 search", "noise": round(sig, 5),
                                 "strength": p.denoise_strength, "seconds": round(time.time() - t1, 1)}
            self.log(f"denoising: non-local means, noise {sig:.4f}, strength {p.denoise_strength:g}, "
                     f"{time.time() - t1:.1f} s")
        ys, ws = self.data_and_weights(op)
        cov = st["coverage"]
        b = op.adjoint([w * y for w, y in zip(ws, ys)])
        avg = np.where(cov > 1e-6, b / np.maximum(cov, 1e-6), 0).astype(np.float32)
        support = _support(cov)
        mask = support
        if p.mask_output and p.mode == "fetal":
            ijk = np.indices(self.grid.shape).reshape(3, -1).T
            inside = self.roi.inside(self.grid.vox_to_world(ijk), margin=2.0).reshape(self.grid.shape)
            soft = ndi.gaussian_filter(inside.astype(np.float32), 1.0)
            mask = support * soft
        f = self.norm[self.ref]
        vol = (x * mask * f).astype(np.float32)
        avg = (avg * mask * f).astype(np.float32)
        self.timing["final"] = round(time.time() - t, 1)
        self.timing["total"] = round(time.time() - t_all, 1)
        self._progress("final", 1.0, "done")
        report = self.make_report(stack_reg, excluded, st)
        return Result(vol, self.grid, avg, cov, report, [T.copy() for T in self.T], [i.copy() for i in self.inc],
                      self.stacks)

    def make_report(self, stack_reg, excluded, st) -> dict:
        stacks = []
        for k, s in enumerate(self.stacks):
            Ts = self.T[k]
            Tm = mean_rigid(list(Ts), self.roi.center)
            ang, disp = rigid_params(Tm, self.roi.center)
            per = []
            for sl in range(s.nz):
                a, d = rigid_params(Ts[sl] @ np.linalg.inv(Tm), self.roi.center)
                stt = getattr(self, "slice_stats", {}).get((k, sl))
                idx = s.meta.get("slice_index", list(range(s.nz)))
                per.append({"slice": int(idx[sl]) + 1, "included": bool(self.inc[k][sl]),
                            "ncc": None if stt is None else round(float(stt[1]), 4),
                            "scale": round(float(self.scale[k][sl]), 4),
                            "pixels": None if stt is None else int(stt[2]),
                            "rotation_deg": round(a, 3), "displacement_mm": round(d, 3)})
            n_inc = int(sum(x["included"] for x in per))
            stacks.append({"name": s.name, "reference": k == self.ref, "slices": s.nz, "slices_used": n_inc,
                           "thickness_mm": round(s.thickness, 3), "pixel_mm": [round(v, 3) for v in s.pixel],
                           "inplane_resolution_mm": [round(v, 3) for v in s.resolution],
                           "binned": int(s.meta.get("binned", 1)), "noise_rel": round(self.noise[k], 4),
                           "stack_rotation_deg": round(ang, 3), "stack_translation_mm": round(disp, 3),
                           "slices_detail": per})
        return {"params": self.p.to_dict(), "alpha_used": round(float(self.alpha), 5),
                "denoise": getattr(self, "denoise_info", None),
                "alpha_motion": round(float(self.alpha_motion), 5), "grid": self.grid.to_dict(),
                "roi": self.roi.to_dict(), "stacks": stacks, "stack_registration": stack_reg,
                "excluded_slices": len(excluded), "timing_s": self.timing,
                "solver": {k: v for k, v in st.items() if k != "coverage"}}


def _ncc(a, b):
    a = a - a.mean()
    b = b - b.mean()
    d = np.sqrt((a @ a) * (b @ b))
    return float((a @ b) / d) if d > 0 else 0.0


def _support(cov: np.ndarray) -> np.ndarray:
    """Soft mask of voxels with enough data (coverage relative to the typical in-region coverage)."""
    pos = cov[cov > 0]
    if pos.size == 0:
        return np.zeros_like(cov)
    ref = float(np.percentile(pos, 75))
    return np.clip(cov / (0.15 * ref), 0, 1).astype(np.float32)
