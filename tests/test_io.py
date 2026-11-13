"""DICOM and NIfTI input/output."""
import numpy as np

from sliceweave import series
from sliceweave.dicom_write import write_series
from sliceweave.nifti import read_nifti, write_nifti

from conftest import DATA


def test_dicom_series_round_trip(tmp_path):
    rng = np.random.default_rng(0)
    vol = (rng.random((40, 36, 12)) * 1000).astype(np.float32)
    aff = np.array([[0.9, 0, 0, -20], [0, 0, -4.4, 30], [0, 0.9, 0, -15], [0, 0, 0, 1]], float)
    write_series(vol, aff, tmp_path / "s", {"SeriesDescription": "T2 HASTE COR TEST", "SliceThickness": 4.0,
                                             "SeriesInstanceUID": "2.25.1", "StudyInstanceUID": "2.25.2"}, slope=1.0)
    found = series.scan([tmp_path / "s"])
    assert len(found) == 1 and series.looks_like_haste(found[0])
    st = series.build_stack(found[0])
    assert st.nz == 12 and st.plane == "COR" and np.isclose(st.thickness, 4.0) and np.isclose(st.spacing, 4.4)
    assert np.allclose(np.sort(st.data.ravel()), np.sort(np.round(vol).ravel()), atol=0.6)


def test_nifti_stack_reads_sidecar_thickness(tmp_path):
    d = np.zeros((20, 20, 6), np.float32)
    write_nifti(tmp_path / "ax.nii.gz", d, np.diag([1.0, 1.0, 3.3, 1.0]), "t")
    (tmp_path / "ax.json").write_text('{"SliceThickness": 3.0}')
    st = series.nifti_stack(tmp_path / "ax.nii.gz")
    assert np.isclose(st.thickness, 3.0) and np.isclose(st.spacing, 3.3)


def test_nifti_round_trip(tmp_path):
    d = np.random.default_rng(1).random((5, 6, 7)).astype(np.float32)
    A = np.array([[-0.8, 0, 0, 10], [0, 0.8, 0, -3], [0, 0, 0.8, 2], [0, 0, 0, 1]], float)
    write_nifti(tmp_path / "v.nii.gz", d, A, "t")
    v = read_nifti(tmp_path / "v.nii.gz")
    assert np.allclose(v.data, d) and np.allclose(v.affine, A)


def test_uncompressed_dicom_samples_decode():
    from sliceweave import dicom
    for name in ("MR_small.dcm", "MR_small_implicit.dcm", "MR_small_bigendian.dcm", "MR_small_RLE.dcm"):
        ds = dicom.read(DATA / "dicom" / name)
        px = dicom.pixel_array(ds, DATA / "dicom" / name)
        assert px.shape[-2:] == (64, 64)
