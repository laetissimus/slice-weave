"""Reproduce the numbers and figures in docs/VALIDATION.md.

    python tools/validate.py --out docs/validation
        [--adult-t1 PATH/chris_t1.nii.gz]          # optional: a real 3-D head scan as ground truth
        [--fetal-dir PATH/mialsuperresolutiontoolkit/data]   # optional: public fetal HASTE stacks

Cases
  phantom-adult  formula phantom, 3 orthogonal stacks (0.9 x 0.9 x 4 mm, 0.4 mm gap), noise,
                 bias field, patient motion between series. Truth known.
  t1-adult       the same protocol simulated from a real 0.88 mm 3-D head scan (structural truth).
  phantom-flair  adult formula phantom with FLAIR contrast and 12 white-matter lesions (2-8 mm),
                 3 orthogonal 2-D FLAIR stacks: 256 x 192 acquired over 230 mm, zero-filled to 512,
                 5 mm slices with a 1.5 mm gap. Also lesion contrast recovery and a pixel-size-PSF ablation.
  phantom-fetal  fetal formula phantom, 6 stacks (1.1 x 1.1 x 3 mm) with slice-by-slice motion,
                 2 partly corrupted slices per stack. Truth and true motion known.
  real-fetal     6 real fetal HASTE stacks (1.5 T, 1.125 x 1.125 x 3 mm) with manual brain masks.
                 No truth: slice-to-volume consistency before/after motion correction.
  cortex         (only when named; needs --adult-t1) gray-white matter and sulci: pseudo-T2 contrast
                 (adult and infant) made from the real 3-D head scan, child protocol 0.6 x 0.6 x 3 mm
                 without gap. Final regularization weight at three noise levels, output voxel size,
                 and slice motion within the series (between-series vs slice-to-volume correction).
  cycles         (only when named) number of motion-correction cycles, fetal and adult phantoms.

Stacks are simulated with a box-shaped slice profile and Gauss-Hermite in-plane sampling, not
with the reconstruction's Gaussian model. Everything runs offline on the CPU.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sliceweave import metrics, simulate  # noqa: E402
from sliceweave.geometry import Grid, apply  # noqa: E402
from sliceweave.nifti import read_nifti  # noqa: E402
from sliceweave.recon import Params, Reconstruction, default_roi  # noqa: E402
from sliceweave.register import Sampler, register  # noqa: E402
from sliceweave.selftest import _pose_error  # noqa: E402
from sliceweave.stacks import Stack  # noqa: E402


def on_grid(vol, aff, grid, D=None, order=1):
    W = grid.vox_to_world(np.indices(grid.shape).reshape(3, -1).T)
    M = np.linalg.inv(aff) @ (np.eye(4) if D is None else D)
    return ndi.map_coordinates(vol, apply(M, W).T, order=order, mode="constant").reshape(grid.shape).astype(np.float32)


def stack_on(st, grid):
    """A single stack interpolated (trilinear) onto the grid: what a reformat of one series shows."""
    return on_grid(st.data, st.volume_affine(), grid)


def figure(path, rows, titles, center=None, window=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(rows)
    fig, axs = plt.subplots(n, 3, figsize=(10.5, 3.5 * n), squeeze=False, facecolor="black")
    ref = rows[-1]
    c = [s // 2 for s in ref.shape] if center is None else center
    for r, (v, t) in enumerate(zip(rows, titles)):
        hi = window or float(np.percentile(v[v > 0], 99.5))  # each row its own window
        sl = [v[:, :, c[2]].T[::-1, ::-1], v[:, c[1], :].T[::-1, ::-1], v[c[0], :, :].T[::-1, ::-1]]
        for k, im in enumerate(sl):
            ax = axs[r, k]
            ax.imshow(im, cmap="gray", vmin=0, vmax=hi, interpolation="nearest")
            ax.set_axis_off()
        axs[r, 0].set_title(t, color="white", fontsize=10, loc="left")
    plt.tight_layout()
    plt.savefig(path, dpi=90, facecolor="black")
    plt.close(fig)


def adult_case(name, truth, taff, mask_fn, out, log, center=None, cover=None, fov=None):
    t0 = time.time()
    stacks, Ts = simulate.simulate_stacks(truth, taff, "adult", motion=False, noise=0.02, bias=0.15,
                                          center=center, cover_mm=cover, fov_mm=fov)
    p = Params.preset("adult")
    roi = default_roi(stacks, 0, "adult")
    rec = Reconstruction(stacks, roi, p, log=lambda m: None)
    res = rec.run()
    g = res.grid
    tg = on_grid(truth, taff, g)
    m = mask_fn(g, tg)
    single = stack_on(stacks[0], g)
    out_m = {"single_axial_stack": metrics.compare(single, tg, m), "average": metrics.compare(res.average, tg, m),
             "srr": metrics.compare(res.volume, tg, m)}
    errs = []
    from sliceweave.geometry import rigid_params
    for k in (1, 2):
        errs.append(rigid_params(np.linalg.inv(Ts[k][0]) @ res.transforms[k][0], roi.center))
        errs[-1] = (errs[-1], rigid_params(Ts[k][0], roi.center))
    out_m["motion_between_series"] = [{"true_deg": round(t[0], 2), "true_mm": round(t[1], 2),
                                       "residual_deg": round(e[0], 3), "residual_mm": round(e[1], 3)} for e, t in errs]
    out_m["seconds"] = round(time.time() - t0, 1)
    out_m["recon_seconds"] = res.report["timing_s"]["total"]
    out_m["grid"] = res.grid.to_dict()["shape"]
    figure(out / f"{name}.png", [tg, single, res.average, res.volume],
           ["truth", "single axial stack (reformatted)", "average of stacks (no SRR)", "SliceWeave SRR"])
    log(f"{name}: {json.dumps(out_m)}")
    return out_m


def lesion_contrast(vol, grid, lesions, truth_grid, wm=0.45, lesion=0.95):
    """Contrast recovery per lesion: (value at the lesion center - median of the surrounding
    white matter) / (true lesion - white-matter difference)."""
    W = grid.vox_to_world(np.indices(grid.shape).reshape(3, -1).T)
    out = []
    for L in lesions:
        c, radii, Rm = np.array(L["center"]), np.array(L["radii"]), np.array(L["R"])
        near = np.all(np.abs(W - c) < 3.2 * radii.max() + 2, axis=1)
        q = (W[near] - c) @ Rm
        d = np.sqrt(((q / radii) ** 2).sum(axis=1))
        shell = (d >= 2.0) & (d <= 3.5) & (np.abs(truth_grid.reshape(-1)[near] - wm) < 0.05)
        if shell.sum() < 10:
            continue
        v_c = ndi.map_coordinates(vol, grid.world_to_vox(c[None]).T, order=1)[0]
        bg = float(np.median(vol.reshape(-1)[near][shell]))
        out.append((L["diameter_mm"], (v_c - bg) / (lesion - wm)))
    bins = {"2-3 mm": (0, 3.01), "4-5 mm": (3.5, 5.01), "6-8 mm": (5.5, 9)}
    return {k: round(float(np.median([cr for dmm, cr in out if lo <= dmm < hi])), 3)
            for k, (lo, hi) in bins.items() if any(lo <= dmm < hi for dmm, _ in out)}


def flair_case(out, log):
    t0 = time.time()
    truth, taff, info = simulate.phantom_head("flair", spacing=0.6)
    ax = info["brain_half_axes_mm"]
    stacks, Ts = simulate.simulate_stacks(truth, taff, "flair", motion=False, noise=0.02, bias=0.15)
    roi = default_roi(stacks, 0, "flair")
    rec = Reconstruction(stacks, roi, Params.preset("flair"), log=lambda m: None)
    res = rec.run()
    g = res.grid
    tg = on_grid(truth, taff, g)
    W = g.vox_to_world(np.indices(g.shape).reshape(3, -1).T)
    m = (((W[:, 0] / ax[0]) ** 2 + (W[:, 1] / ax[1]) ** 2 + (W[:, 2] / ax[2]) ** 2) < 1.0).reshape(g.shape)
    # ablation: the same reconstruction assuming the pixel size is the resolution (headers ignored)
    plain = []
    for st in stacks:
        st2 = Stack(st.name, st.data, st.slice_aff, st.thickness, None, st.plane, {k: v for k, v in st.meta.items()
                                                                                   if k != "acq_res"})
        plain.append(st2)
    # same regularization weights as the main run, so only the point-spread function differs
    a_fin, a_mot = res.report["alpha_used"], res.report["alpha_motion"]

    class SameAlpha(Reconstruction):
        def prepare(self):
            super().prepare()
            self.alpha, self.alpha_motion = a_fin, a_mot

    res_px = SameAlpha(plain, roi, Params.preset("flair"), log=lambda m: None).run()
    vol_px = on_grid(res_px.volume, res_px.grid.affine, g)
    # ablation: the noise-calibrated weight of the other modes also for the final volume
    res_cal = Reconstruction(stacks, roi, Params.preset("flair", alpha=a_mot, denoise=False), log=lambda m: None).run()
    # with the denoising step (off by default for FLAIR)
    res_nod = Reconstruction(stacks, roi, Params.preset("flair", denoise=True), log=lambda m: None).run()
    vol_nod = on_grid(res_nod.volume, res_nod.grid.affine, g)
    vol_cal = on_grid(res_cal.volume, res_cal.grid.affine, g)
    single = stack_on(stacks[0], g)
    vols = {"single_axial_stack": single, "average": res.average, "srr": res.volume, "srr_pixel_psf": vol_px,
            "srr_noise_calibrated_alpha": vol_cal, "srr_with_denoising": vol_nod}
    r = {}
    for k, v in vols.items():
        s_ = metrics.fit_scale(v, tg, m)
        r[k] = metrics.compare(v, tg, m)
        r[k]["lesion_contrast_recovery"] = lesion_contrast(v * s_, g, info["lesions"], tg)
    r["truth_on_output_grid"] = {"lesion_contrast_recovery": lesion_contrast(tg, g, info["lesions"], tg)}
    r["noise_rel_estimated"] = [s["noise_rel"] for s in res.report["stacks"]]
    r["alpha"] = res.report["alpha_used"]
    r["alpha_motion_correction"] = res.report["alpha_motion"]
    r["motion_between_series"] = []
    from sliceweave.geometry import rigid_params
    for k in (1, 2):
        tr = rigid_params(Ts[k][0], roi.center)
        er = rigid_params(np.linalg.inv(Ts[k][0]) @ res.transforms[k][0], roi.center)
        r["motion_between_series"].append({"true_deg": round(tr[0], 2), "true_mm": round(tr[1], 2),
                                           "residual_deg": round(er[0], 3), "residual_mm": round(er[1], 3)})
    r["recon_seconds"] = res.report["timing_s"]["total"]
    r["seconds"] = round(time.time() - t0, 1)
    L = sorted(info["lesions"], key=lambda x: abs(x["diameter_mm"] - 5))[0]
    ctr = [int(round(v)) for v in g.world_to_vox(np.array(L["center"])[None])[0]]
    figure(out / "phantom-flair.png", [tg, single, res.average, res.volume],
           ["truth (FLAIR phantom, lesions 2-8 mm)", "single axial stack (reformatted)", "average of stacks (no SRR)",
            "SliceWeave SRR (2-D FLAIR preset)"], center=ctr)
    log(f"phantom-flair: {json.dumps(r)}")
    return r


def fetal_phantom(out, log):
    t0 = time.time()
    truth, taff, info = simulate.phantom_head("fetal", spacing=0.6)
    a, b, c = info["brain_half_axes_mm"]
    stacks, Ts = simulate.simulate_stacks(truth, taff, "fetal", static=truth, fetal_radius=(a * 1.3, b * 1.3, c * 1.3),
                                          noise=0.02, motion_deg=6.0, motion_mm=3.0)
    roi = default_roi(stacks, 0, "fetal")
    roi.center, roi.half = np.zeros(3), np.array([a + 4, b + 4, c + 4])
    results = {}

    def brain(g):
        W = g.vox_to_world(np.indices(g.shape).reshape(3, -1).T)
        return W, (((W[:, 0] / a) ** 2 + (W[:, 1] / b) ** 2 + (W[:, 2] / c) ** 2) < 1.0)

    def aligned_truth(res):
        """Recon frame -> truth frame (rigid), so that the volumes can be compared."""
        W, inb = brain(res.grid)
        r = register(Sampler(truth, Grid(truth.shape, taff)), W[inb][::4], res.volume.reshape(-1)[inb][::4],
                     center=[0, 0, 0], iters=40, huber_k=None)
        return on_grid(truth, taff, res.grid, r.D), inb.reshape(res.grid.shape)

    vols = {}
    for label, kw in [("no_motion_correction", {"motion": "none"}), ("sliceweave", {}),
                      ("sliceweave_with_denoising", {"denoise": True})]:
        rec = Reconstruction(stacks, roi, Params.preset("fetal", **kw), log=lambda m: None)
        rec.prepare()
        e0 = _pose_error(rec, Ts)
        res = rec.run()
        e1 = _pose_error(rec, Ts)
        tg, m = aligned_truth(res)
        results[label] = {"srr": metrics.compare(res.volume, tg, m), "average": metrics.compare(res.average, tg, m),
                          "slice_pose_error_mm_before": round(e0, 2), "slice_pose_error_mm_after": round(e1, 2),
                          "slices_used": int(sum(s["slices_used"] for s in res.report["stacks"])),
                          "slices": int(sum(s["slices"] for s in res.report["stacks"])),
                          "seconds": res.report["timing_s"]["total"]}
        vols[label] = (res, tg)
        log(f"phantom-fetal {label}: {json.dumps(results[label])}")
    # oracle: the true motion given to the same super-resolution step (upper bound)
    rec = Reconstruction(stacks, roi, Params.preset("fetal", motion="none", outlier_rejection=False), log=lambda m: None)
    rec.prepare()
    for k, st in enumerate(rec.stacks):
        rec.T[k] = Ts[k][np.array(st.meta["slice_index"])].copy()
    op = rec.build(rec.masks())
    x, _ = rec.srr(op, iters=30)
    g = rec.grid
    W, inb = brain(g)
    tg0 = on_grid(truth, taff, g)
    results["true_motion_upper_bound"] = {"srr": metrics.compare(x, tg0, inb.reshape(g.shape))}
    log(f"phantom-fetal oracle: {json.dumps(results['true_motion_upper_bound'])}")
    res, tg = vols["sliceweave"]
    res0, _ = vols["no_motion_correction"]
    figure(out / "phantom-fetal.png", [tg, res0.volume, res.volume],
           ["truth (aligned)", "SRR without motion correction", "SliceWeave SRR with slice-to-volume correction"])
    results["seconds"] = round(time.time() - t0, 1)
    return results


# ---------------------------------------------------------------- cortex (gray-white matter, sulci)
def _largest(m):
    lab, n = ndi.label(m)
    if n == 0:
        return m
    return lab == (1 + int(np.argmax(ndi.sum(m, lab, range(1, n + 1)))))


def pseudo_t2(t1_path, contrast="adult"):
    """T2-like truth from a 3-D T1 head scan: CSF bright; adult: gray matter brighter than white
    matter; infant (unmyelinated): white matter brighter than gray matter, lower contrast. Returns
    the truth, its affine, the normalized T1 (for tissue labels), the brain region and its center."""
    v = read_nifti(t1_path)
    t = (v.data / np.percentile(v.data, 99.5)).astype(np.float32)
    head = _largest(ndi.binary_fill_holes(ndi.binary_closing(t > 0.12, iterations=3)))
    brain = _largest(ndi.binary_opening((t > 0.45) & (t < 0.95), iterations=2))
    brain_ext = ndi.binary_fill_holes(ndi.binary_closing(brain, iterations=4)) & head
    knots = [0.0, 0.2, 0.45, 0.62, 0.80, 0.95, 1.3]
    vals = [2.4, 2.4, 1.55, 1.15, 0.80, 0.80, 0.9] if contrast == "adult" else [2.4, 2.4, 1.6, 1.0, 1.2, 1.2, 1.2]
    outside = np.interp(t, [0, 0.15, 0.5, 0.8, 1.0, 1.3], [0, 0.15, 0.6, 1.2, 1.6, 1.6])
    T2 = np.where(brain_ext, np.interp(t, knots, vals), np.where(head, outside, 0)).astype(np.float32)
    ctr = apply(v.affine, np.argwhere(brain)[::50]).mean(axis=0)
    return T2, v.affine, t, brain_ext, ctr


def cortex_metrics(vol, g, cov, truth, taff, t, brain_ext, D=None):
    """Brain PSNR / SSIM, and on tissue labels from the T1: gray-white matter contrast recovery,
    sulcal CSF - gray matter contrast recovery, and the correlation of band-pass detail (0.5-2 mm)
    in the cortical ribbon. Only voxels with enough data are evaluated. D: rigid map from the
    reconstruction frame to the truth frame (when the reference stack itself moved)."""
    tr = on_grid(truth, taff, g, D)
    tg = on_grid(t, taff, g, D)
    pos = cov[cov > 0]
    be = (on_grid(brain_ext.astype(np.float32), taff, g, D) > 0.5) & (cov >= 0.15 * float(np.percentile(pos, 75)))
    gm = be & (tg > 0.55) & (tg < 0.68)
    wm = ndi.binary_erosion(be & (tg > 0.76) & (tg < 0.9), iterations=1)
    csf = be & (tg < 0.35)
    vent = ndi.binary_erosion(csf, iterations=2)
    sulcal = csf & ndi.binary_dilation(gm, iterations=2) & ~ndi.binary_dilation(vent, iterations=2)
    cortex = ndi.binary_dilation(gm, iterations=2) & be
    v = vol * metrics.fit_scale(vol, tr, be)
    out = metrics.compare(vol, tr, be)
    out.pop("ncc", None)

    def recovery(a, b):
        return round(float((np.median(v[a]) - np.median(v[b])) / (np.median(tr[a]) - np.median(tr[b]))), 3)
    out["gm_wm_contrast"] = recovery(gm, wm)
    out["sulcal_csf_contrast"] = recovery(sulcal, gm)
    h = g.spacing
    dv = ndi.gaussian_filter(v, 0.5 / h) - ndi.gaussian_filter(v, 2.0 / h)
    dt = ndi.gaussian_filter(tr, 0.5 / h) - ndi.gaussian_filter(tr, 2.0 / h)
    out["cortex_detail_corr"] = round(float(np.corrcoef(dv[cortex], dt[cortex])[0, 1]), 3)
    return out


def cortex_case(t1_path, out, log):
    from sliceweave import operators as ops
    from sliceweave.denoise import nlm3d
    from sliceweave.recon import FINAL_ALPHA_FACTOR, _support, noise_sigma
    res_all = {"protocol": "AX, COR, SAG; 0.6 x 0.6 mm in-plane, 3 mm slices, no gap; bias 15 %; motion between "
                           "series; output 1.0 mm unless stated", "final_alpha_factor": FINAL_ALPHA_FACTOR}
    sim_kw = dict(bias=0.15, cover_mm=126, fov_mm=200, pixel_mm=0.6, thick_mm=3.0, gap_mm=0.0)

    class Keep(Reconstruction):
        def srr(self, op, x0=None, iters=20, final=False):
            if final:
                self.x_before_final = None if x0 is None else x0.copy()
            return super().srr(op, x0=x0, iters=iters, final=final)

    # 1. weight of the final volume (the registration volumes keep the noise-calibrated weight)
    calib = []
    for contrast, noise in (("adult", 0.015), ("adult", 0.03), ("adult", 0.06), ("infant", 0.03)):
        truth, taff, t, be, ctr = pseudo_t2(t1_path, contrast)
        stacks, _ = simulate.simulate_stacks(truth, taff, "adult", noise=noise, center=ctr, **sim_kw)
        rec = Keep(stacks, default_roi(stacks, 0, "adult"), Params.preset("adult"), log=lambda m: None)
        res = rec.run()
        op = rec.build(rec.masks())
        ys, ws = rec.data_and_weights(op)
        f = rec.norm[rec.ref]
        row = {"contrast": contrast, "noise": noise, "noise_estimated": round(rec.rel_noise, 4),
               "alpha_noise_calibrated": round(rec.alpha_motion, 3), "by_factor": {}}
        for fa in (1.0, 0.5, 0.35, 0.25, 0.15, 0.1):
            x, st = ops.solve(op, ys, ws, x0=rec.x_before_final, alpha=rec.alpha_motion * fa,
                              iters=rec.p.final_iterations, huber=True)
            row["by_factor"][str(fa)] = cortex_metrics(x * f * _support(st["coverage"]), rec.grid, st["coverage"],
                                                       truth, taff, t, be)
        # the default: 0.25 x 0.4 = 0.1 x the weight, then non-local means (as in Reconstruction.run)
        x, st = ops.solve(op, ys, ws, x0=rec.x_before_final, alpha=rec.alpha, iters=rec.p.final_iterations, huber=True)
        sup = _support(st["coverage"])
        fg = sup > 0.5
        fg &= x > 0.2 * float(np.percentile(x[fg], 95))
        x = nlm3d(x, noise_sigma(x, fg), rec.p.denoise_strength)
        row["denoised_default"] = {"factor": round(rec.alpha / rec.alpha_motion, 3),
                                   **cortex_metrics(x * f * sup, rec.grid, st["coverage"], truth, taff, t, be)}
        calib.append(row)
        log(f"cortex weight: {json.dumps(row)}")
    res_all["final_weight"] = calib
    # 2. output voxel size (defaults), adult contrast, noise 0.03
    truth, taff, t, be, ctr = pseudo_t2(t1_path, "adult")
    stacks, _ = simulate.simulate_stacks(truth, taff, "adult", noise=0.03, center=ctr, **sim_kw)
    sizes, vols = [], {}
    for r_mm in (1.5, 1.0, 0.8):
        t0 = time.time()
        res = Reconstruction(stacks, default_roi(stacks, 0, "adult"), Params.preset("adult", resolution=r_mm),
                             log=lambda m: None).run()
        row = {"voxel_mm": r_mm, **cortex_metrics(res.volume, res.grid, res.coverage, truth, taff, t, be),
               "seconds": round(time.time() - t0)}
        sizes.append(row)
        vols[r_mm] = res
        log(f"cortex voxel: {json.dumps(row)}")
    res_all["voxel_size"] = sizes
    # 3. slice motion within the series (a child moving during each series)
    stacks_m, Ts_m = simulate.simulate_stacks(truth, taff, "adult", noise=0.03, center=ctr, motion=True,
                                              motion_deg=3.0, motion_mm=2.0, corrupt=0, **sim_kw)
    mot, Ds = [], {}
    for label, kw in (("between series, 1 cycle", {}), ("slice-to-volume, 6 cycles", {"motion": "slice", "cycles": 6})):
        t0 = time.time()
        rec = Reconstruction(stacks_m, default_roi(stacks_m, 0, "adult"), Params.preset("adult", **kw),
                             log=lambda m: None)
        rec.prepare()
        e0 = _pose_error(rec, Ts_m)
        res = rec.run()
        # the reference stack moved too: compare after the best rigid alignment to the truth
        W = res.grid.vox_to_world(np.indices(res.grid.shape).reshape(3, -1).T)
        inb = on_grid(be.astype(np.float32), taff, res.grid).reshape(-1) > 0.5
        D = register(Sampler(truth, Grid(truth.shape, taff)), W[inb][::4], res.volume.reshape(-1)[inb][::4],
                     center=ctr, iters=40, huber_k=None).D
        row = {"correction": label, **cortex_metrics(res.volume, res.grid, res.coverage, truth, taff, t, be, D),
               "slice_pose_error_mm_before": round(e0, 2), "slice_pose_error_mm_after": round(_pose_error(rec, Ts_m), 2),
               "slices_used": int(sum(s["slices_used"] for s in res.report["stacks"])),
               "slices": int(sum(s["slices"] for s in res.report["stacks"])), "seconds": round(time.time() - t0)}
        mot.append(row)
        vols[label], Ds[label] = res, D
        log(f"cortex motion: {json.dumps(row)}")
    res_all["slice_motion"] = mot
    # figure: truth, 1.5 mm, 1.0 mm (default) on the 1.0 mm grid; moving child without / with slice correction
    g = vols[1.0].grid
    def moved(label):  # into the truth frame: x_truth = D x_rec, so sample the result at D^-1 x
        r = vols[label]
        return on_grid(r.volume, r.grid.affine, g, np.linalg.inv(Ds[label]))
    rows = [on_grid(truth, taff, g), on_grid(vols[1.5].volume, vols[1.5].grid.affine, g), vols[1.0].volume,
            moved("between series, 1 cycle"), moved("slice-to-volume, 6 cycles")]
    rows = [r * (metrics.fit_scale(r, rows[0], rows[0] > 0.05) if i else 1.0) for i, r in enumerate(rows)]
    nx, ny, nz = g.shape
    figure(out / "cortex-t2.png", rows, ["truth (pseudo-T2 from a 3-D scan)", "SRR 1.5 mm voxels",
                                          "SRR 1.0 mm (default)", "child moving: between-series correction",
                                          "child moving: slice-to-volume, 6 cycles"],
           center=[int(nx * 0.38), ny // 2, int(nz * 0.62)], window=1.8)
    return res_all


def cycles_sweep(out, log):
    """Effect of the number of motion-correction cycles (fetal phantom, slice-to-volume; adult phantom,
    between series). Same stacks for every setting."""
    res_all = {"fetal_slice_to_volume": [], "adult_between_series": []}
    truth, taff, info = simulate.phantom_head("fetal", spacing=0.6)
    a, b, c = info["brain_half_axes_mm"]
    stacks, Ts = simulate.simulate_stacks(truth, taff, "fetal", static=truth, fetal_radius=(a * 1.3, b * 1.3, c * 1.3),
                                          noise=0.02, motion_deg=6.0, motion_mm=3.0)
    roi = default_roi(stacks, 0, "fetal")
    roi.center, roi.half = np.zeros(3), np.array([a + 4, b + 4, c + 4])
    for n in (1, 2, 3, 4, 6, 8):
        rec = Reconstruction(stacks, roi, Params.preset("fetal", cycles=n), log=lambda m: None)
        res = rec.run()
        W = res.grid.vox_to_world(np.indices(res.grid.shape).reshape(3, -1).T)
        inb = ((W[:, 0] / a) ** 2 + (W[:, 1] / b) ** 2 + (W[:, 2] / c) ** 2) < 1.0
        r = register(Sampler(truth, Grid(truth.shape, taff)), W[inb][::4], res.volume.reshape(-1)[inb][::4],
                     center=[0, 0, 0], iters=40, huber_k=None)
        tg = on_grid(truth, taff, res.grid, r.D)
        row = {"cycles": n, "srr": metrics.compare(res.volume, tg, inb.reshape(res.grid.shape)),
               "slice_pose_error_mm": round(_pose_error(rec, Ts), 2),
               "slices_used": int(sum(s["slices_used"] for s in res.report["stacks"])),
               "slices": int(sum(s["slices"] for s in res.report["stacks"])),
               "seconds": res.report["timing_s"]["total"]}
        res_all["fetal_slice_to_volume"].append(row)
        log(f"cycles fetal: {json.dumps(row)}")
    truth, taff, info = simulate.phantom_head("adult", spacing=0.6)
    ax = info["brain_half_axes_mm"]
    stacks, Ts = simulate.simulate_stacks(truth, taff, "adult", motion=False, noise=0.02, bias=0.15)
    roi = default_roi(stacks, 0, "adult")
    from sliceweave.geometry import rigid_params
    for n in (0, 1, 2, 3):
        res = Reconstruction(stacks, roi, Params.preset("adult", cycles=n), log=lambda m: None).run()
        g = res.grid
        W = g.vox_to_world(np.indices(g.shape).reshape(3, -1).T)
        m = (((W[:, 0] / ax[0]) ** 2 + (W[:, 1] / ax[1]) ** 2 + (W[:, 2] / ax[2]) ** 2) < 1.0).reshape(g.shape)
        errs = [rigid_params(np.linalg.inv(Ts[k][0]) @ res.transforms[k][0], roi.center) for k in (1, 2)]
        row = {"cycles": n, "srr": metrics.compare(res.volume, on_grid(truth, taff, g), m),
               "residual_deg_max": round(max(e[0] for e in errs), 3),
               "residual_mm_max": round(max(e[1] for e in errs), 3),
               "seconds": res.report["timing_s"]["total"]}
        res_all["adult_between_series"].append(row)
        log(f"cycles adult: {json.dumps(row)}")
    return res_all


def real_fetal(data_dir: Path, out, log):
    from sliceweave.series import nifti_stack
    from sliceweave.stacks import Stack
    anat = data_dir / "sub-01" / "anat"
    mdir = data_dir / "derivatives" / "manual_masks" / "sub-01" / "anat"
    stacks = []
    for f in sorted(anat.glob("*_T2w.nii.gz")):
        st = nifti_stack(f)
        mf = mdir / f.name.replace("_T2w.nii.gz", "_T2w_desc-brain_mask.nii.gz")
        if mf.exists():
            mv = read_nifti(mf)
            st.mask = Stack.from_volume("m", mv.data, mv.affine).data > 0.5
        st.name = f"{st.plane} run-{len(stacks) + 1}"
        stacks.append(st)
    from sliceweave.ingest import suggest_reference
    entries = [{"id": str(i), "include": True, "plane": s.plane} for i, s in enumerate(stacks)]
    ref = int(suggest_reference(stacks, entries, "fetal"))
    results = {"stacks": [s.summary() for s in stacks], "reference": stacks[ref].name}
    vols = {}
    for label, kw in [("no_motion_correction", {"motion": "none"}), ("sliceweave", {}),
                      ("sliceweave_with_denoising", {"denoise": True})]:
        roi = default_roi(stacks, ref, "fetal")
        rec = Reconstruction(stacks, roi, Params.preset("fetal", reference=ref, **kw), log=lambda m: None)
        res = rec.run()
        nccs = [v[1] for v in rec.slice_stats.values()]
        results[label] = {"median_slice_ncc": round(float(np.median(nccs)), 4),
                          "p10_slice_ncc": round(float(np.percentile(nccs, 10)), 4),
                          "slices_used": int(sum(s["slices_used"] for s in res.report["stacks"])),
                          "slices": int(sum(s["slices"] for s in res.report["stacks"])),
                          "stack_motion_deg_mm": [(s["name"], s["stack_rotation_deg"], s["stack_translation_mm"])
                                                  for s in res.report["stacks"]],
                          "seconds": res.report["timing_s"]["total"]}
        vols[label] = res
        log(f"real-fetal {label}: {json.dumps(results[label])}")
    res = vols["sliceweave"]
    single = on_grid(stacks[ref].data, stacks[ref].volume_affine(), res.grid) * \
        (res.volume.max() / max(np.percentile(stacks[ref].data, 99.9), 1e-6))
    single *= (res.volume > 0)
    figure(out / "real-fetal.png", [single, vols["no_motion_correction"].volume, res.volume],
           [f"reference stack {stacks[ref].name} (reformatted)", "SRR without motion correction",
            "SliceWeave SRR with slice-to-volume correction"])
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="docs/validation")
    ap.add_argument("--adult-t1")
    ap.add_argument("--fetal-dir")
    ap.add_argument("--only", help="comma-separated subset: phantom-adult,t1-adult,phantom-flair,phantom-fetal,real-fetal; "
                                     "cycles (number-of-cycles sweep, only when named)")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    only = set(a.only.split(",")) if a.only else None
    resf = out / "results.json"
    allres = json.loads(resf.read_text()) if resf.exists() else {}
    log = lambda m: print(m, flush=True)  # noqa: E731

    def want(n):
        return only is None or n in only

    if want("phantom-adult"):
        truth, taff, info = simulate.phantom_head("adult", spacing=0.6)
        ax = info["brain_half_axes_mm"]

        def brain(g, tg):
            W = g.vox_to_world(np.indices(g.shape).reshape(3, -1).T)
            return (((W[:, 0] / ax[0]) ** 2 + (W[:, 1] / ax[1]) ** 2 + (W[:, 2] / ax[2]) ** 2) < 1.0).reshape(g.shape)
        allres["phantom-adult"] = adult_case("phantom-adult", truth, taff, brain, out, log)
    if a.adult_t1 and want("t1-adult"):
        v = read_nifti(a.adult_t1)
        truth = (v.data / np.percentile(v.data, 99.5)).astype(np.float32)
        fg = ndi.binary_fill_holes(ndi.binary_closing(truth > 0.15, iterations=3))
        W = apply(v.affine, np.argwhere(fg)[::50])
        ctr = W.mean(axis=0)

        def head(g, tg):
            m = ndi.binary_erosion(ndi.binary_fill_holes(tg > 0.15), iterations=3)
            return m
        allres["t1-adult"] = adult_case("t1-adult", truth, v.affine, head, out, log, center=ctr, cover=176.0,
                                        fov=240.0)
        allres["t1-adult"]["source"] = Path(a.adult_t1).name
    if want("phantom-flair"):
        allres["phantom-flair"] = flair_case(out, log)
    if want("phantom-fetal"):
        allres["phantom-fetal"] = fetal_phantom(out, log)
    if a.fetal_dir and want("real-fetal"):
        allres["real-fetal"] = real_fetal(Path(a.fetal_dir), out, log)
    if a.adult_t1 and only is not None and "cortex" in only:  # not part of the default run (about 45 min)
        allres["cortex"] = cortex_case(a.adult_t1, out, log)
    if only is not None and "cycles" in only:  # not part of the default run (about 25 min)
        allres["cycles"] = cycles_sweep(out, log)
    resf.write_text(json.dumps(allres, indent=1), encoding="utf-8")
    print(f"wrote {resf}")


if __name__ == "__main__":
    main()
