"""2-D FLAIR preset: recognition, acquired resolution, planning rules and the motion guard."""
import json
from types import SimpleNamespace

import numpy as np
import pytest

from sliceweave import series
from sliceweave.dicom_write import write_series
from sliceweave.nifti import write_nifti
from sliceweave.recon import FLAIR_SLICE_MSG, Params, Reconstruction, noise_sigma, psf_cov
from sliceweave.simulate import zero_fill
from sliceweave.stacks import Stack
from sliceweave.workflow import _select_for_mode, plan_warnings


def _s(desc="", seq="", ti=None, te=None):
    return SimpleNamespace(text=lambda: desc, sequence=seq, ti=ti, te=te)


def test_flair_recognised_by_name_sequence_or_timing():
    assert series.looks_like_flair(_s("AX T2 FLAIR"))
    assert series.looks_like_flair(_s("t2_tirm_tra_dark-fluid"))
    assert series.looks_like_flair(_s("ax brain", "*tir2d1_21", ti=2500, te=86))
    assert series.looks_like_flair(_s("Ax T2 IR", ti=2400, te=120))        # long TI and long TE
    assert not series.looks_like_flair(_s("t2_haste_tra", "*h2d1_256", te=90))
    assert not series.looks_like_flair(_s("AX STIR", ti=160, te=40))        # short-TI inversion recovery
    # T1-FLAIR is also multi-shot inversion recovery: recognised by name and handled by the same preset
    assert series.looks_like_flair(_s("T1 FLAIR", ti=800, te=10))


def test_haste_is_not_flair():
    s = series.SeriesInfo(uid="1", study_uid="2", kind="dicom", description="t2_haste_cor", sequence="*h2d1_320", te=90)
    assert series.looks_like_haste(s) and not series.looks_like_flair(s)
    f = series.SeriesInfo(uid="1", study_uid="2", kind="dicom", description="AX FLAIR single shot", te=90, ti=2500)
    assert series.looks_like_flair(f) and not series.looks_like_haste(f)


def test_acquired_resolution_from_dicom(tmp_path):
    """256 (readout) x 192 (phase, along the row) acquired over 230 mm, shown at 512 x 512."""
    vol = np.zeros((512, 512, 6), np.float32)
    vol[200:300, 200:300] = 100
    sp = 230.0 / 512
    aff = np.array([[sp, 0, 0, -115], [0, sp, 0, -115], [0, 0, 6.5, 0], [0, 0, 0, 1]], float)
    write_series(vol, aff, tmp_path / "s", {"SeriesDescription": "AX T2 FLAIR", "SliceThickness": 5.0,
                                             "AcquisitionMatrix": [0, 256, 192, 0],
                                             "InPlanePhaseEncodingDirection": "ROW", "InversionTime": 2500.0,
                                             "EchoTime": 86.0, "RepetitionTime": 9000.0,
                                             "SeriesInstanceUID": "2.25.7", "StudyInstanceUID": "2.25.8"}, slope=1.0)
    s = series.scan([tmp_path / "s"])[0]
    assert series.looks_like_flair(s)
    st = series.build_stack(s)
    assert np.allclose(st.pixel, (sp, sp), atol=1e-4)
    assert np.allclose(st.resolution, (230.0 / 192, 230.0 / 256), atol=1e-3)


def test_acquired_resolution_from_sidecar(tmp_path):
    d = np.zeros((320, 320, 8), np.float32)
    write_nifti(tmp_path / "ax.nii.gz", d, np.diag([1.125, 1.125, 3.3, 1.0]), "t")
    (tmp_path / "ax.json").write_text(json.dumps({"SliceThickness": 3, "AcquisitionMatrixPE": 224, "ReconMatrixPE": 320,
                                                  "PhaseEncodingDirection": "i", "BaseResolution": 320}))
    st = series.nifti_stack(tmp_path / "ax.nii.gz")
    assert np.allclose(st.resolution, (1.125 * 320 / 224, 1.125), atol=1e-3)


def test_psf_follows_the_acquired_resolution():
    sa = np.repeat(np.eye(4)[None], 4, 0)
    sa[:, :3, :3] = np.diag([0.45, 0.45, 1.0])
    st = Stack("s", np.zeros((8, 8, 4), np.float32), sa, 5.0)
    c0 = psf_cov(st, np.eye(3))
    st.meta["acq_res"] = [1.2, 0.9]
    c1 = psf_cov(st, np.eye(3))
    assert c1[0, 0] > 6 * c0[0, 0] and c1[1, 1] > 3.9 * c0[1, 1] and np.isclose(c1[2, 2], c0[2, 2])


def test_noise_estimate_on_interpolated_images():
    rng = np.random.default_rng(0)
    img = rng.normal(0, 1, (96, 96, 4)) + 1j * rng.normal(0, 1, (96, 96, 4))
    up = zero_fill(img, (192, 192)).real.astype(np.float32)
    m = np.ones(up.shape, bool)
    white = noise_sigma(up, m)
    acq = noise_sigma(up, m, (2.0, 2.0))
    assert white < 0.3 and 0.85 < acq < 1.15


def test_flair_preset_and_slice_motion_guard():
    p = Params.preset("flair")
    assert p.mode == "flair" and p.motion == "stack" and p.resolution == 1.0
    sa = np.repeat(np.eye(4)[None], 6, 0)
    st = Stack("s", np.ones((20, 20, 6), np.float32), sa, 4.0)
    from sliceweave.recon import ROI
    roi = ROI(np.zeros(3), np.array([10.0, 10, 10]), np.eye(3))
    with pytest.raises(ValueError, match="multi-shot"):
        Reconstruction([st], roi, Params.preset("flair", motion="slice"))
    assert "multi-shot" in FLAIR_SLICE_MSG


def test_mode_switch_selects_one_contrast_and_warns():
    def e(i, plane, fl, te, ti=None):
        return {"id": f"s{i}", "plane": plane, "include": True, "flair_like": fl, "haste_like": not fl,
                "sequence": "flair" if fl else "single-shot", "te_ms": te, "ti_ms": ti}
    entries = [e(1, "AX", False, 90), e(2, "COR", False, 90), e(3, "SAG", False, 90),
               e(4, "AX", True, 86, 2500), e(5, "COR", True, 86, 2500), e(6, "SAG", True, 86, 2500)]
    w = plan_warnings(entries, "adult")
    assert any("FLAIR and non-FLAIR" in x for x in w)
    _select_for_mode(entries, "flair")
    assert [x["include"] for x in entries] == [False] * 3 + [True] * 3
    assert plan_warnings(entries, "flair") == []
    _select_for_mode(entries, "adult")
    assert [x["include"] for x in entries] == [True] * 3 + [False] * 3
    entries[4]["te_ms"] = 140
    _select_for_mode(entries, "flair")
    assert any("echo time" in x for x in plan_warnings(entries, "flair"))
    entries[5]["include"] = False
    assert any("planes" in x for x in plan_warnings(entries, "flair"))
