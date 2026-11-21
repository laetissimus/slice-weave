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
    """Contrast recovery per lesion: (value at the lesion centre - median of the surrounding
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
    # same regularisation weights as the main run, so only the point-spread function differs
    a_fin, a_mot = res.report["alpha_used"], res.report["alpha_motion"]

    class SameAlpha(Reconstruction):
        def prepare(self):
            super().prepare()
            self.alpha, self.alpha_motion = a_fin, a_mot

    res_px = SameAlpha(plain, roi, Params.preset("flair"), log=lambda m: None).run()
    vol_px = on_grid(res_px.volume, res_px.grid.affine, g)
    # ablation: the noise-calibrated weight of the other modes also for the final volume
    res_cal = Reconstruction(stacks, roi, Params.preset("flair", alpha=a_mot), log=lambda m: None).run()
    vol_cal = on_grid(res_cal.volume, res_cal.grid.affine, g)
    single = stack_on(stacks[0], g)
    vols = {"single_axial_stack": single, "average": res.average, "srr": res.volume, "srr_pixel_psf": vol_px,
            "srr_noise_calibrated_alpha": vol_cal}
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
    for label, kw in [("no_motion_correction", {"motion": "none"}), ("sliceweave", {})]:
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
    for label, kw in [("no_motion_correction", {"motion": "none"}), ("sliceweave", {})]:
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
    if only is not None and "cycles" in only:  # not part of the default run (about 25 min)
        allres["cycles"] = cycles_sweep(out, log)
    resf.write_text(json.dumps(allres, indent=1), encoding="utf-8")
    print(f"wrote {resf}")


if __name__ == "__main__":
    main()
