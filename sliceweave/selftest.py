"""Self-test and demo cases on synthetic phantoms (no patient data, nothing downloaded).

The self-test checks, end to end on simulated stacks with known ground truth:
* DICOM write/read round trip of a simulated stack (geometry preserved);
* adult mode: stack-to-stack motion is recovered and the super-resolution volume is closer to
  the truth than the plain average of the stacks (PSNR, SSIM);
* fetal mode (``--full`` only): slice-to-volume motion correction reduces the slice pose error.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from . import metrics, simulate
from .geometry import apply, rigid_params
from .recon import Params, Reconstruction, default_roi


def make_demo(store, cid: str, kind: str = "adult", log=None):
    """Create a simulated case (stacks saved like imported ones)."""
    log = log or (lambda m: None)
    truth, aff, info = simulate.phantom_head(kind, spacing=0.7 if kind == "fetal" else 1.0)
    a, b, c = info["brain_half_axes_mm"]
    kw = {"static": truth, "fetal_radius": (a * 1.3, b * 1.3, c * 1.3)} if kind == "fetal" else {}
    stacks, _ = simulate.simulate_stacks(truth, aff, kind, log=log, **kw)
    case = store.load(cid)
    entries = []
    for n, st in enumerate(stacks):
        sid = f"s{n + 1:02d}"
        store.save_stack(cid, sid, st)
        dx, dy = st.pixel
        rx, ry = st.resolution
        fl = kind == "flair"
        entries.append({"id": sid, "label": f"{st.plane} {n + 1} (SIM)",
                        "description": st.name + (" 2D FLAIR" if fl else " HASTE") + " - SIMULATED",
                        "plane": st.plane, "shape": [int(v) for v in st.shape], "pixel_mm": [round(dx, 3), round(dy, 3)],
                        "resolution_mm": [round(rx, 3), round(ry, 3)],
                        "resolution_source": "acquisition matrix" if st.meta.get("acq_res") else "pixel size",
                        "thickness_mm": st.thickness, "spacing_mm": round(st.spacing, 3),
                        "thickness_source": "simulation", "source": "simulated", "haste_like": not fl,
                        "flair_like": fl, "sequence": "flair" if fl else "single-shot", "include": True,
                        "has_mask": False, "acq_order_known": True, "series_number": n + 1, "warnings": []})
    case["stacks"] = entries
    case["simulated"] = True
    case["status"] = "imported"
    case["mode_hint"] = kind
    store.save(cid, case)
    store.log(cid, f"demo: {len(stacks)} simulated {kind} stacks (phantom, no patient data)")


def _truth_on(grid, truth, aff, D=None):
    W = grid.vox_to_world(np.indices(grid.shape).reshape(3, -1).T)
    M = np.linalg.inv(aff) if D is None else np.linalg.inv(aff) @ D
    return ndi.map_coordinates(truth, apply(M, W).T, order=1).reshape(grid.shape), W


def run(full: bool = False, out: Path | None = None, log=print) -> dict:
    t0 = time.time()
    results = {"checks": []}

    def check(name, passed, **detail):
        results["checks"].append({"name": name, "passed": bool(passed), **detail})
        log(f"[{'PASS' if passed else 'FAIL'}] {name} {json.dumps(detail)}")

    # 1. DICOM round trip
    import tempfile
    from . import series
    from .dicom_write import write_series
    with tempfile.TemporaryDirectory() as td:
        rng = np.random.default_rng(0)
        vol = (rng.random((40, 36, 12)) * 1000).astype(np.float32)
        aff = np.array([[0.9, 0, 0, -20], [0, 0, -4.4, 30], [0, 0.9, 0, -15], [0, 0, 0, 1]], float)
        write_series(vol, aff, Path(td) / "s", {"SeriesDescription": "T2 HASTE COR TEST", "SliceThickness": 4.0,
                                                 "SeriesInstanceUID": "2.25.1", "StudyInstanceUID": "2.25.2"},
                     slope=1.0)
        found = series.scan([Path(td) / "s"])
        st = series.build_stack(found[0])
        same = np.allclose(np.sort(st.data.ravel()), np.sort(np.round(vol).ravel()), atol=0.6)
        A = np.diag([-1.0, -1.0, 1.0, 1.0]) @ aff
        geo = np.allclose(np.sort([s[:3, 3] @ st.normal for s in st.slice_aff]),
                          np.sort([(aff @ [0, 0, k, 1])[:3] @ st.normal for k in range(12)]), atol=1e-3)
        check("DICOM write and read back", same and geo and st.plane == "COR", plane=st.plane, slices=st.nz,
              haste_like=series.looks_like_haste(found[0]))
        del A
        # 2-D FLAIR: recognised, acquired (not displayed) resolution read from the header
        v2 = np.zeros((128, 128, 6), np.float32)
        v2[40:90, 40:90] = 500
        a2 = np.array([[1.8, 0, 0, -115], [0, 1.8, 0, -115], [0, 0, 6.5, 0], [0, 0, 0, 1]], float)
        write_series(v2, a2, Path(td) / "f", {"SeriesDescription": "AX T2 FLAIR TEST", "SliceThickness": 5.0,
                                               "AcquisitionMatrix": [0, 64, 48, 0], "InPlanePhaseEncodingDirection": "ROW",
                                               "InversionTime": 2500.0, "EchoTime": 86.0,
                                               "SeriesInstanceUID": "2.25.3", "StudyInstanceUID": "2.25.2"}, slope=1.0)
        f = series.scan([Path(td) / "f"])[0]
        stf = series.build_stack(f)
        ok_res = np.allclose(stf.resolution, (128 * 1.8 / 48, 128 * 1.8 / 64), atol=1e-3)
        try:
            Reconstruction([stf], default_roi([stf], 0, "flair"), Params.preset("flair", motion="slice"))
            guard = False
        except ValueError:
            guard = True
        check("2-D FLAIR: recognised, acquired resolution read, slice-by-slice motion refused",
              series.looks_like_flair(f) and not series.looks_like_haste(f) and ok_res and guard,
              resolution_mm=[round(v, 3) for v in stf.resolution])
    # 2. adult mode on the phantom
    truth, taff, info = simulate.phantom_head("adult", spacing=1.0 if not full else 0.7)
    stacks, Ts = simulate.simulate_stacks(truth, taff, "adult", motion=False, noise=0.02)
    p = Params.preset("adult", resolution=1.0 if not full else 0.9)
    if not full:
        p.final_iterations, p.cycle_iterations = 20, 10
    roi = default_roi(stacks, 0, "adult")
    res = Reconstruction(stacks, roi, p, log=lambda m: None).run()
    tg, W = _truth_on(res.grid, truth, taff)
    a, b, c = info["brain_half_axes_mm"]
    brain = (((W[:, 0] / a) ** 2 + (W[:, 1] / b) ** 2 + (W[:, 2] / c) ** 2) < 1.0).reshape(res.grid.shape)
    m_srr = metrics.compare(res.volume, tg, brain)
    m_avg = metrics.compare(res.average, tg, brain)
    errs = []
    for k in (1, 2):
        est = res.transforms[k][0]
        errs.append(rigid_params(np.linalg.inv(Ts[k][0]) @ est, roi.center))
    max_deg = max(e[0] for e in errs)
    max_mm = max(e[1] for e in errs)
    check("adult: motion between series recovered", max_deg < 0.5 and max_mm < 0.5,
          residual_deg=round(max_deg, 3), residual_mm=round(max_mm, 3))
    check("adult: super-resolution closer to the truth than the average of the stacks",
          m_srr["psnr_db"] > m_avg["psnr_db"] and m_srr["ssim"] > m_avg["ssim"], srr=m_srr, average=m_avg)
    results["adult"] = {"srr": m_srr, "average": m_avg, "seconds": res.report["timing_s"]}
    # 3. fetal mode (slow)
    if full:
        truth, taff, info = simulate.phantom_head("fetal", spacing=0.7)
        a, b, c = info["brain_half_axes_mm"]
        stacks, Ts = simulate.simulate_stacks(truth, taff, "fetal", static=truth,
                                              fetal_radius=(a * 1.3, b * 1.3, c * 1.3))
        p = Params.preset("fetal")
        roi = default_roi(stacks, 0, "fetal")
        roi.center = np.zeros(3)
        roi.half = np.array([a + 4, b + 4, c + 4])
        rec = Reconstruction(stacks, roi, p, log=lambda m: None)
        rec.prepare()
        e0 = _pose_error(rec, Ts)
        res = rec.run()
        e1 = _pose_error(rec, Ts)
        check("fetal: slice motion correction reduces the slice pose error", e1 < 0.75 * e0,
              before_mm=round(e0, 2), after_mm=round(e1, 2))
        results["fetal"] = {"pose_error_before_mm": e0, "pose_error_after_mm": e1}
    results["seconds"] = round(time.time() - t0, 1)
    results["passed"] = all(c["passed"] for c in results["checks"])
    if out:
        Path(out).write_text(json.dumps(results, indent=1), encoding="utf-8")
    return results


def _pose_error(rec, Ts) -> float:
    """Median slice pose error (mm, mean over points +/-30 mm around the slice centre) after the
    best global rigid alignment of the estimated to the true poses."""
    off = np.array([[0, 0, 0], [30, 0, 0], [-30, 0, 0], [0, 30, 0], [0, -30, 0], [0, 0, 30], [0, 0, -30]], float)
    E, T = [], []
    for k, st in enumerate(rec.stacks):
        idx = np.array(st.meta["slice_index"])
        for s in range(st.nz):
            pc = apply(st.slice_aff[s], [[st.shape[0] / 2, st.shape[1] / 2, 0]])[0]
            E.append(apply(rec.T[k][s], pc + off))
            T.append(apply(Ts[k][idx[s]], pc + off))
    E, T = np.concatenate(E), np.concatenate(T)
    me, mt = E.mean(0), T.mean(0)
    U, _, Vt = np.linalg.svd((E - me).T @ (T - mt))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    res = np.linalg.norm((E - me) @ R.T + mt - T, axis=1).reshape(-1, len(off)).mean(1)
    return float(np.median(res))
