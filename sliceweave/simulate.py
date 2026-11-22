"""Synthetic test data: procedural head phantoms and simulated 2-D HASTE-like stacks.

Nothing here comes from a patient. The phantoms are generated from formulas, so the self-test
and the demo run without any downloaded data. Simulated stacks are labeled SIMULATED.

Stacks are simulated from a high-resolution "truth" volume by sampling the volume at each pixel
with a slice profile that is deliberately *not* the reconstruction model (a box-shaped slice
profile of width = slice thickness through-plane, Gauss-Hermite quadrature in-plane), so the
reconstruction is not tested with its own forward model.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from .geometry import about_center, apply, rotvec_to_matrix
from .stacks import Stack, plane_label

# T2-weighted (HASTE-like) relative intensities
T2 = {"bg": 0.0, "scalp": 0.55, "skull": 0.05, "csf": 1.0, "gm": 0.52, "wm": 0.38, "deep": 0.42,
      "fetal_wm": 0.62, "fetal_gm": 0.36, "fluid": 0.95, "uterus": 0.45}
# FLAIR: CSF suppressed, gray matter brighter than white matter, lesions bright
FLAIR = {"bg": 0.0, "scalp": 0.55, "skull": 0.05, "csf": 0.06, "gm": 0.62, "wm": 0.45, "deep": 0.55, "lesion": 0.95}
LESION_RADII_MM = (1.0, 1.0, 1.25, 1.5, 1.5, 2.0, 2.0, 2.5, 2.5, 3.0, 3.0, 4.0)


# ----------------------------------------------------------------- phantoms
def _grid(fov_mm, spacing):
    n = int(round(fov_mm / spacing))
    c = (np.arange(n) - (n - 1) / 2.0) * spacing
    aff = np.diag([spacing, spacing, spacing, 1.0])
    aff[:3, 3] = c[0]
    return n, c, aff


def _field(rng, shape, sigma_vox):
    f = ndi.gaussian_filter(rng.normal(0, 1, shape).astype(np.float32), sigma_vox)
    return f / (f.std() + 1e-9)


def phantom_head(kind: str = "adult", spacing: float = 0.6, seed: int = 7) -> tuple[np.ndarray, np.ndarray, dict]:
    """Return (data, affine, info): a T2-weighted-like head phantom generated from formulas.

    Adult: 150 mm head (scalp, skull, CSF, folded cortex with winding sulci, white matter,
    deep gray nuclei, ventricles, cerebellum, brainstem, interhemispheric fissure).
    Fetal (about 30 weeks): 80 mm brain with shallow folding, a thin dark cortical ribbon,
    bright unmyelinated white matter, asymmetric ventricles, cerebellum and brainstem, inside
    amniotic fluid and a textured uterine background.
    FLAIR (kind "flair"): the adult anatomy with FLAIR contrast and 12 bright white-matter
    lesions of 2-8 mm diameter (ellipsoids, listed in info["lesions"]).
    All include resolution targets (thin plates in three orientations, small spheres)."""
    rng = np.random.default_rng(seed)
    flair = kind == "flair"
    adult = kind in ("adult", "flair")
    tab = FLAIR if flair else T2
    fov = 176.0 if adult else 136.0
    n, c, aff = _grid(fov, spacing)
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    vol = np.zeros((n, n, n), np.float32)
    sc = 1.0 if adult else 0.52                     # anatomical scale
    a, b, cc = (66.0, 80.0, 64.0) if adult else (34.0, 41.0, 31.0)
    R = np.sqrt((X / a) ** 2 + (Y / b) ** 2 + (Z / cc) ** 2)   # 1 at the cortical surface
    f1 = _field(rng, (n, n, n), (4.0 * sc) / spacing)          # folding depth
    f2 = _field(rng, (n, n, n), (2.6 * sc) / spacing)          # sulcal pattern
    gm, wm, csf = (tab["gm"], tab["wm"], tab["csf"]) if adult else (T2["fetal_gm"], T2["fetal_wm"], T2["csf"])
    if adult:
        vol[R < 1.18] = tab["scalp"]
        vol[R < 1.10] = tab["skull"]
        vol[R < 1.04] = csf
    else:
        tex = _field(rng, (n, n, n), 3.0 / spacing)
        vol[:] = np.clip(T2["uterus"] + 0.12 * tex, 0.05, 1.0)
        Rf = np.sqrt(((X - 4) / 62) ** 2 + (Y / 60) ** 2 + ((Z + 3) / 58) ** 2)
        vol[Rf < 1.0] = T2["fluid"]
        vol[R < 1.22] = T2["scalp"] * 0.8
        vol[R < 1.15] = T2["skull"] * 2
        vol[R < 1.09] = csf
    surf = 1.0 + (0.012 if adult else 0.008) * f1
    brain = R < surf
    vol[brain] = gm
    depth = 0.14 if adult else 0.05                         # cortical thickness relative to the radius
    vol[brain & (R < surf - depth + (0.05 if adult else 0.02) * f1)] = wm
    sulci = brain & (np.abs(f2) < (0.16 if adult else 0.10)) & (R > (0.72 if adult else 0.86))
    vol[sulci] = csf
    vol[brain & (np.abs(X) < 0.9 * sc + 0.4) & (Z > -18 * sc)] = csf     # interhemispheric fissure
    # deep gray nuclei, ventricles (asymmetric), cerebellum, brainstem
    for sx in (-1, 1):
        dg = (((X - sx * 16 * sc) / (11 * sc)) ** 2 + ((Y + 2 * sc) / (17 * sc)) ** 2 + ((Z + 2 * sc) / (11 * sc)) ** 2) < 1
        vol[dg] = tab["deep"] if adult else gm * 1.15
        big = 1.25 if sx < 0 else 1.0
        vent = (((np.abs(X) - 8 * sc) / (5 * sc * big)) ** 2 + ((Y - 6 * sc) / (24 * sc)) ** 2
                + ((Z - 10 * sc) / (8 * sc * big)) ** 2 < 1) & (np.sign(X) == sx)
        horn = (((X - sx * 14 * sc) / (4 * sc)) ** 2 + ((Y + 30 * sc) / (12 * sc)) ** 2
                + ((Z - 4 * sc) / (6 * sc * big)) ** 2) < 1
        vol[vent | horn] = csf
    cereb = (X / (26 * sc)) ** 2 + ((Y + 44 * sc) / (15 * sc)) ** 2 + ((Z + 30 * sc) / (13 * sc)) ** 2 < 1
    vol[cereb] = gm if adult else T2["fetal_gm"] * 1.3
    folia = cereb & (np.sin((Z + 30 * sc) * 2 * np.pi / (2.4 * sc + 0.6) + 2 * f1) > 0.55)
    vol[folia] = wm if adult else T2["fetal_wm"]
    stem = ((X / (8 * sc)) ** 2 + ((Y + 22 * sc) / (9 * sc)) ** 2 < 1) & (Z < -6 * sc) & (Z > -52 * sc)
    vol[stem] = wm if adult else T2["fetal_gm"] * 1.4
    # resolution targets inside the brain: thin plates in three orientations and small spheres
    plate = 0.9 if adult else 0.7
    targets = [
        (np.abs(Z - 32 * sc) < plate / 2) & (X ** 2 + (Y - 10 * sc) ** 2 < (14 * sc) ** 2),
        (np.abs(Y + 50 * sc) < plate / 2) & ((X - 22 * sc) ** 2 + (Z - 6 * sc) ** 2 < (11 * sc) ** 2),
        (np.abs(X + 36 * sc) < plate / 2) & ((Y - 22 * sc) ** 2 + (Z - 12 * sc) ** 2 < (11 * sc) ** 2),
    ]
    for t in targets:
        vol[t] = csf
    for k, r in enumerate([1.0, 1.5, 2.0, 3.0]):
        ctr = np.array([-24 + 14 * k, 34, -14]) * sc
        sph = (X - ctr[0]) ** 2 + (Y - ctr[1]) ** 2 + (Z - ctr[2]) ** 2 < (r * (1 if adult else 0.8)) ** 2
        vol[sph] = csf
    lesions = _add_lesions(vol, c, spacing, wm, rng, R) if flair else []
    vol = ndi.gaussian_filter(vol, 0.35).astype(np.float32)
    info = {"kind": kind, "spacing": spacing, "brain_half_axes_mm": [a, b, cc], "simulated": True}
    if flair:
        info["lesions"] = lesions
    return vol, aff, info


def _add_lesions(vol, c, spacing, wm, rng, R) -> list[dict]:
    """Bright ellipsoidal white-matter lesions (axis ratio up to 1.6, random orientation), at
    least 12 mm apart, centered in white matter."""
    from .geometry import rotvec_to_matrix
    n = len(c)
    out = []
    for r in LESION_RADII_MM:
        for _ in range(2000):
            idx = rng.integers(0, n, 3)
            if not (0.3 < R[tuple(idx)] < 0.78) or abs(vol[tuple(idx)] - wm) > 1e-6:
                continue
            ctr = c[idx]
            if any(np.linalg.norm(ctr - np.array(L["center"])) < 12.0 for L in out):
                continue
            radii = np.array([r * rng.uniform(1.0, 1.6), r, r])
            Rm = rotvec_to_matrix(rng.normal(0, 1, 3))
            h = int(np.ceil((radii.max() + 1.5) / spacing))
            sl = tuple(slice(max(i - h, 0), min(i + h + 1, n)) for i in idx)
            gx, gy, gz = np.meshgrid(c[sl[0]], c[sl[1]], c[sl[2]], indexing="ij")
            q = (np.stack([gx, gy, gz], -1) - ctr) @ Rm
            inside = ((q / radii) ** 2).sum(-1) <= 1.0
            vol[sl][inside] = FLAIR["lesion"]
            out.append({"center": [round(float(v), 3) for v in ctr], "radii": [round(float(v), 3) for v in radii],
                        "R": np.round(Rm, 6).tolist(), "diameter_mm": round(2 * r, 2)})
            break
    return out


# -------------------------------------------------------------- stack simulation
_GH3 = (np.array([-1.224744871, 0.0, 1.224744871]), np.array([0.2954089752, 1.181635901, 0.2954089752]) / np.sqrt(np.pi))


def stack_geometry(axes: np.ndarray, center, n_slices: int, pixel, fov_mm: tuple, spacing: float,
                   n_rows: int | None = None) -> np.ndarray:
    """slice_aff (nz, 4, 4) for a stack with in-plane directions axes[:,0], axes[:,1].
    ``pixel``: one size or (along i, along j)."""
    u, v, n = axes[:, 0], axes[:, 1], axes[:, 2]
    px, py = (pixel, pixel) if np.isscalar(pixel) else pixel
    nx = int(round(fov_mm[0] / px))
    ny = int(round(fov_mm[1] / py))
    out = np.zeros((n_slices, 4, 4))
    for s in range(n_slices):
        off = (s - (n_slices - 1) / 2.0) * spacing
        o = np.asarray(center, float) - u * px * (nx - 1) / 2 - v * py * (ny - 1) / 2 + n * off
        A = np.eye(4)
        A[:3, 0] = u * px
        A[:3, 1] = v * py
        A[:3, 2] = n
        A[:3, 3] = o
        out[s] = A
    return out, (nx, ny)


def sample_slice(truth, truth_aff, A, T, shape2d, pixel, thickness, profile="box", static=None, fetal_mask_fn=None):
    """Simulate one slice: pixel grid A (4x4), object motion T (pixel world -> object world)."""
    nx, ny = shape2d
    jj, ii = np.meshgrid(np.arange(ny), np.arange(nx))
    base = np.stack([ii.ravel(), jj.ravel(), np.zeros(ii.size)], axis=1).astype(np.float64)
    u = A[:3, 0] / np.linalg.norm(A[:3, 0])
    v = A[:3, 1] / np.linalg.norm(A[:3, 1])
    n = np.cross(u, v)
    px, py = (pixel, pixel) if np.isscalar(pixel) else pixel
    sx, sy = 1.2 * px / 2.3548, 1.2 * py / 2.3548
    if profile == "box":
        zs = (np.arange(7) - 3) / 7.0 * thickness  # 7 equally spaced points over the slice width
        zw = np.ones(7) / 7.0
    else:
        g = np.polynomial.hermite.hermgauss(5)
        zs, zw = g[0] * np.sqrt(2) * thickness / 2.3548, g[1] / np.sqrt(np.pi)
    inv = np.linalg.inv(truth_aff)
    P0 = apply(A, base)
    out = np.zeros(len(base), np.float64)
    for a, wa in zip(*_GH3):
        for b, wb in zip(*_GH3):
            for z, wz in zip(zs, zw):
                P = P0 + u * (a * np.sqrt(2) * sx) + v * (b * np.sqrt(2) * sy) + n * z
                Q = apply(T, P)
                val = ndi.map_coordinates(truth, apply(inv, Q).T, order=1, mode="constant", cval=0.0)
                if static is not None:
                    sv = ndi.map_coordinates(static, apply(inv, P).T, order=1, mode="constant", cval=0.0)
                    inside = fetal_mask_fn(Q)
                    val = np.where(inside, val, sv)
                out += wa * wb * wz * val
    return out.reshape(nx, ny).astype(np.float32)


def random_motion(n_slices: int, rng, max_rot_deg: float, max_trans_mm: float, center, order: str = "interleaved",
                  jumps: int = 1, jump_scale: float = 0.5):
    """Smooth random rigid motion over acquisition time (slice order) plus a few sudden jumps."""
    if order == "interleaved":
        acq = list(range(0, n_slices, 2)) + list(range(1, n_slices, 2))
    else:
        acq = list(range(n_slices))
    t = np.linspace(0, 1, n_slices)
    def smooth(scale):
        k = rng.normal(0, 1, (4, 3))
        traj = sum(np.outer(np.sin((i + 1) * np.pi * t + rng.uniform(0, 6.3)), k[i]) / (i + 1) for i in range(4))
        traj = traj / (np.abs(traj).max() + 1e-9) * scale
        return traj
    rot = smooth(np.deg2rad(max_rot_deg))
    tra = smooth(max_trans_mm)
    for _ in range(jumps):
        j = rng.integers(n_slices // 3, n_slices)
        rot[j:] += rng.normal(0, np.deg2rad(max_rot_deg) * jump_scale, 3)
        tra[j:] += rng.normal(0, max_trans_mm * jump_scale, 3)
    Ts = np.zeros((n_slices, 4, 4))
    for time_idx, s in enumerate(acq):
        Ts[s] = about_center(rotvec_to_matrix(rot[time_idx]), center, tra[time_idx])
    return Ts, acq


def simulate_stacks(truth, truth_aff, kind: str = "adult", seed: int = 3, motion: bool | None = None,
                    noise: float = 0.02, bias: float = 0.15, corrupt: int | None = None, log=None,
                    static=None, fetal_radius=None, motion_deg: float = 6.0, motion_mm: float = 3.0,
                    tilt_deg: float | None = None, center=None, cover_mm: float | None = None,
                    fov_mm: float | None = None, pixel_mm: float | None = None, thick_mm: float | None = None,
                    gap_mm: float | None = None) -> tuple[list[Stack], list[np.ndarray]]:
    """Simulate 3-plane stacks (adult: AX, COR, SAG once each; fetal: each orientation twice,
    oblique, with slice motion). Returns the stacks and the true per-slice transforms
    (pixel world -> truth world)."""
    rng = np.random.default_rng(seed)
    flair = kind == "flair"
    adult = kind in ("adult", "flair")
    motion = (not adult) if motion is None else motion
    recon_matrix = None
    if flair:
        # 2-D FLAIR: 256 (readout) x 192 (phase) acquired over 230 mm, zero-filled to 512 x 512,
        # 5 mm slices with a 1.5 mm gap; phase encoding along the first image axis
        pixel, thick, gap, fov = (230.0 / 192, 230.0 / 256), 5.0, 1.5, (230.0, 230.0)
        recon_matrix = (512, 512)
        planes = [("AX", 0), ("COR", 0), ("SAG", 0)]
        cover = 150.0
    elif adult:
        pixel, thick, gap, fov = 0.9, 4.0, 0.4, (176.0, 176.0)
        planes = [("AX", 0), ("COR", 0), ("SAG", 0)]
        cover = 150.0
    else:
        pixel, thick, gap, fov = 1.1, 3.0, 0.3, (130.0, 130.0)
        planes = [("AX", 0), ("COR", 0), ("SAG", 0), ("AX", 1), ("COR", 1), ("SAG", 1)]
        cover = 96.0
    corrupt = (0 if adult else 2) if corrupt is None else corrupt
    center = np.zeros(3) if center is None else np.asarray(center, float)
    if cover_mm:
        cover = float(cover_mm)
    if fov_mm:
        fov = (float(fov_mm), float(fov_mm))
    if pixel_mm and not flair:
        pixel = float(pixel_mm)
    if thick_mm:
        thick = float(thick_mm)
    if gap_mm is not None:
        gap = float(gap_mm)
    base_axes = {"AX": np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]], float).T,   # rows x, columns -y
                 "COR": np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]], float).T,
                 "SAG": np.array([[0, 1, 0], [0, 0, -1], [-1, 0, 0]], float).T}
    stacks, truths = [], []
    for name, rep in planes:
        ax = base_axes[name].copy()
        # small obliquity (adult: planned on the head; fetal: planned on the moving fetus)
        tilt = rng.normal(0, (2.0 if adult else 6.0) if tilt_deg is None else tilt_deg, 3)
        Rt = rotvec_to_matrix(np.deg2rad(tilt))
        ax = Rt @ ax
        spacing = thick + gap
        nz = int(np.ceil(cover / spacing))
        off = rng.uniform(-spacing / 2, spacing / 2) if rep else 0.0
        sa, (nx, ny) = stack_geometry(ax, center + ax[:, 2] * off, nz, pixel, fov, spacing)
        if motion:
            Ts, acq = random_motion(nz, rng, motion_deg, motion_mm, center)
        else:
            Ts, acq = np.repeat(np.eye(4)[None], nz, axis=0), list(range(nz))
            if adult:  # patient moved a little between series (stack-level motion)
                Tstack = about_center(rotvec_to_matrix(np.deg2rad(rng.normal(0, 2.0, 3))), center,
                                      rng.normal(0, 1.5, 3))
                if stacks:
                    Ts = np.repeat(Tstack[None], nz, axis=0)
        data = np.zeros((nx, ny, nz), np.float32)
        fm = None
        if static is not None:
            def fm(Q, r=fetal_radius):
                q = np.asarray(Q)
                return (q[:, 0] / r[0]) ** 2 + (q[:, 1] / r[1]) ** 2 + (q[:, 2] / r[2]) ** 2 < 1
        for s in range(nz):
            data[:, :, s] = sample_slice(truth, truth_aff, sa[s], Ts[s], (nx, ny), pixel, thick,
                                         static=static, fetal_mask_fn=fm)
        # smooth multiplicative bias, Rician-like noise
        if bias:
            g = [np.linspace(-1, 1, k) for k in (nx, ny, nz)]
            Gx, Gy, Gz = np.meshgrid(*g, indexing="ij")
            cvec = rng.normal(0, 1, 3)
            field = 1.0 + bias * (cvec[0] * Gx + cvec[1] * Gy + 0.5 * cvec[2] * Gz) / np.sqrt((cvec ** 2).sum())
            data *= field.astype(np.float32)
        meta_extra = {}
        if recon_matrix is not None:
            # complex noise at the acquired resolution, then zero-filled interpolation (as the scanner does)
            cplx = data + (rng.normal(0, noise, data.shape) + 1j * rng.normal(0, noise, data.shape) if noise else 0)
            data = np.abs(zero_fill(cplx, recon_matrix)).astype(np.float32)
            fx, fy = nx / recon_matrix[0], ny / recon_matrix[1]
            sa = np.array([A @ np.diag([fx, fy, 1.0, 1.0]) for A in sa])
            meta_extra = {"acq_res": [round(float(pixel[0]), 4), round(float(pixel[1]), 4)],
                          "acquired_matrix": [nx, ny], "recon_matrix": list(recon_matrix)}
            nx, ny = recon_matrix
        elif noise:
            data = np.sqrt((data + rng.normal(0, noise, data.shape)) ** 2 + rng.normal(0, noise, data.shape) ** 2)
        data = data.astype(np.float32)
        bad = []
        for _ in range(corrupt):
            # partial signal loss over a band of the slice (motion during the echo train)
            s = int(rng.integers(2, nz - 2))
            r0 = int(rng.integers(0, ny // 2))
            data[:, r0:r0 + ny // 2, s] *= np.float32(rng.uniform(0.2, 0.5))
            bad.append(s)
        st = Stack(f"SIM {name}{'' if not rep else ' 2'}", data * 1000.0, sa, thick, None, plane_label(ax[:, 2]),
                   {"simulated": True, "acq_order": [int(v) for v in np.argsort(acq)], "corrupted_slices": bad,
                    **meta_extra})
        stacks.append(st)
        truths.append(Ts)
        if log:
            log(f"simulated {st.name}: {nx}x{ny}x{nz}, pixel {st.pixel[0]:.2f} mm (acquired {st.resolution[0]:.2f} x "
                f"{st.resolution[1]:.2f}), thickness {thick} mm, gap {gap} mm")
    return stacks, truths


def zero_fill(img: np.ndarray, shape2d: tuple) -> np.ndarray:
    """Zero-filled (sinc) interpolation of each slice to ``shape2d``, keeping the first sample's
    position: k-space is padded symmetrically, as scanners do for an interpolated matrix."""
    nx, ny = img.shape[:2]
    mx, my = shape2d
    K = np.fft.fftshift(np.fft.fft2(img, axes=(0, 1)), axes=(0, 1))
    P = np.zeros((mx, my) + img.shape[2:], complex)
    ox, oy = (mx - nx) // 2, (my - ny) // 2
    P[ox:ox + nx, oy:oy + ny] = K
    return np.fft.ifft2(np.fft.ifftshift(P, axes=(0, 1)), axes=(0, 1)) * (mx * my) / (nx * ny)
