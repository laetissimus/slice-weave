"""Settings and display that affect cortical detail: resolution warning, cycles following the motion
correction, the brain window of the viewer, the weaker final regularisation of adult / paediatric T2."""
import numpy as np

from sliceweave.layers import brain_window
from sliceweave.recon import FINAL_ALPHA_FACTOR, auto_alpha
from sliceweave.workflow import _cycles_follow_motion, plan_warnings


def _entries(res=(0.6, 0.6)):
    return [{"id": f"s{i}", "include": True, "plane": p, "sequence": "single-shot", "resolution_mm": list(res)}
            for i, p in enumerate(("AX", "COR", "SAG"))]


def test_coarse_output_voxel_is_flagged():
    assert any("coarser than the acquired" in w for w in plan_warnings(_entries(), "adult", 1.5))
    assert not any("coarser" in w for w in plan_warnings(_entries(), "adult", 1.0))
    # FLAIR 1.2 x 0.9 mm acquired: 1.2 mm output is fine
    assert not any("coarser" in w for w in plan_warnings(_entries((1.2, 0.9)), "flair", 1.2))


def test_cycles_follow_motion():
    plan = {"mode": "adult", "params": {"motion": "slice", "cycles": 1}}
    _cycles_follow_motion(plan, {"motion": "stack", "cycles": 1})
    assert plan["params"]["cycles"] == 6
    plan = {"mode": "adult", "params": {"motion": "stack", "cycles": 6}}
    _cycles_follow_motion(plan, {"motion": "slice", "cycles": 6})
    assert plan["params"]["cycles"] == 1
    # a cycle count the user set in the same change is kept
    plan = {"mode": "adult", "params": {"motion": "slice", "cycles": 3}}
    _cycles_follow_motion(plan, {"motion": "stack", "cycles": 1})
    assert plan["params"]["cycles"] == 3


def test_brain_window_spreads_parenchyma():
    rng = np.random.default_rng(0)
    wm, gm, csf = rng.normal(0.8, 0.03, 40000), rng.normal(1.15, 0.03, 30000), rng.normal(2.4, 0.05, 15000)
    v = np.concatenate([wm, gm, csf]).astype(np.float32)
    full_hi = float(np.percentile(v, 99.0))
    lo, hi = brain_window(v, full_hi)
    assert lo < 0.8 and 1.15 < hi < 2.4          # grey and white matter inside, CSF saturated
    assert (1.15 - 0.8) / (hi - lo) > 1.5 * (1.15 - 0.8) / full_hi


def test_adult_final_weight_is_weaker():
    assert FINAL_ALPHA_FACTOR["adult"] < 0.5
    assert auto_alpha(0.03, "huber", "adult") < auto_alpha(0.03, "huber")


def test_dicom_export_carries_brain_and_full_windows(tmp_path):
    from sliceweave import dicom
    from sliceweave.dicom_write import write_series
    from sliceweave.export import dicom_windows
    rng = np.random.default_rng(1)
    vol = np.zeros((40, 40, 10), np.float32)
    vol[5:35, 5:35] = rng.normal(800, 30, (30, 30, 10))      # white matter
    vol[10:30, 10:20] = rng.normal(1150, 30, (20, 10, 10))   # grey matter
    vol[18:22, 22:30] = rng.normal(2400, 50, (4, 8, 10))     # CSF
    win = dicom_windows(vol)
    (bc, fc), (bw, fw) = win["WindowCenter"], win["WindowWidth"]
    assert bw < 0.8 * fw and bc - bw / 2 < 800 < 1150 < bc + bw / 2 < 2400
    write_series(vol, np.eye(4), tmp_path / "s", {"SeriesDescription": "SRR TEST", **win}, slope=1.0)
    ds = dicom.read(sorted((tmp_path / "s").iterdir())[0])
    assert len(ds.WindowCenter) == 2 and np.isclose(float(ds.WindowWidth[0]), bw, rtol=1e-3)
