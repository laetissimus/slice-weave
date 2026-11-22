"""Settings and display that affect cortical detail: resolution warning, cycles following the motion
correction, the brain window of the viewer, the weaker final regularization of adult/pediatric T2."""
import numpy as np
import pytest

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
    assert lo < 0.8 and 1.15 < hi < 2.4          # gray and white matter inside, CSF saturated
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
    vol[10:30, 10:20] = rng.normal(1150, 30, (20, 10, 10))   # gray matter
    vol[18:22, 22:30] = rng.normal(2400, 50, (4, 8, 10))     # CSF
    win = dicom_windows(vol)
    (bc, fc), (bw, fw) = win["WindowCenter"], win["WindowWidth"]
    assert bw < 0.8 * fw and bc - bw / 2 < 800 < 1150 < bc + bw / 2 < 2400
    write_series(vol, np.eye(4), tmp_path / "s", {"SeriesDescription": "SRR TEST", **win}, slope=1.0)
    ds = dicom.read(sorted((tmp_path / "s").iterdir())[0])
    assert len(ds.WindowCenter) == 2 and np.isclose(float(ds.WindowWidth[0]), bw, rtol=1e-3)


def test_nlm_removes_noise_and_keeps_edges():
    from sliceweave.denoise import nlm2d, nlm3d
    rng = np.random.default_rng(3)
    clean = np.zeros((40, 40, 24), np.float32)
    clean[:, 20:] = 1.0                      # a sharp edge
    clean[10:14, 5:9, 8:12] = 2.0            # a small bright structure
    noisy = clean + rng.normal(0, 0.1, clean.shape).astype(np.float32)
    for f in (lambda x: nlm3d(x, 0.1, 1.5), lambda x: nlm3d(x, 0.1, 1.5, chunk=8), lambda x: nlm2d(x, 0.1, 1.5)):
        d = f(noisy)
        assert np.sqrt(np.mean((d - clean) ** 2)) < 0.6 * np.sqrt(np.mean((noisy - clean) ** 2))
        assert abs(float(d[:, 21:].mean() - d[:, :19].mean()) - 1.0) < 0.05      # edge height kept
        assert float(d[11:13, 6:8, 9:11].mean()) > 1.8                            # small structure kept
    assert np.array_equal(nlm3d(noisy, 0.0), noisy)


def test_chunked_nlm_equals_whole_volume():
    from sliceweave.denoise import nlm3d
    x = np.random.default_rng(4).random((20, 18, 30)).astype(np.float32)
    assert np.allclose(nlm3d(x, 0.2, chunk=200), nlm3d(x, 0.2, chunk=7), atol=1e-6)


def test_denoise_setting_in_plan_and_params():
    from sliceweave.workflow import _coerce, EDITABLE
    from sliceweave.recon import Params, DENOISE_ALPHA_FACTOR
    assert "denoise" in EDITABLE and _coerce("denoise", 0) is False and _coerce("denoise_strength", 99) == 4.0
    assert Params.preset("adult").denoise and not Params.preset("adult", denoise=False).denoise
    assert not Params.preset("flair").denoise and not Params.preset("fetal").denoise
    assert 0 < DENOISE_ALPHA_FACTOR < 1


def test_version_follows_the_release_tag(monkeypatch):
    import importlib, sys as _sys
    from pathlib import Path
    _sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "packaging"))
    bw = importlib.import_module("build_windows")
    monkeypatch.setenv("SLICEWEAVE_VERSION", "v0.1.2")
    assert bw.resolve_version() == "0.1.2"
    monkeypatch.setenv("SLICEWEAVE_VERSION", "nightly")
    with pytest.raises(SystemExit):
        bw.resolve_version()
    from sliceweave.config import BASE_VERSION, VERSION
    assert VERSION.startswith(BASE_VERSION)
