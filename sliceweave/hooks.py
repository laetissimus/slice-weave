"""Optional decoder for JPEG-family compressed DICOM (pydicom + a pixel-data plugin such as GDCM).

Switched on automatically when installed; the Windows executable bundles both.
"""
from __future__ import annotations

from functools import lru_cache

from .dicom import DicomError


@lru_cache(maxsize=1)
def pydicom_version() -> str | None:
    try:
        import pydicom  # noqa: F401
        return pydicom.__version__
    except Exception:  # noqa: BLE001
        return None


@lru_cache(maxsize=1)
def gdcm_version() -> str | None:
    try:
        import gdcm  # noqa: F401
        return str(gdcm.Version.GetVersion())
    except Exception:  # noqa: BLE001
        return None


def pydicom_pixels(path, transfer_syntax: str):
    if path is None:
        raise DicomError(f"compressed transfer syntax {transfer_syntax}: file path needed for decoding")
    if pydicom_version() is None:
        raise DicomError(
            f"compressed DICOM (transfer syntax {transfer_syntax}). Install the optional decoder: "
            "pip install pydicom python-gdcm  (the Windows executable includes them)")
    import pydicom
    try:
        ds = pydicom.dcmread(str(path))
        return ds.pixel_array
    except Exception as e:  # missing plugin or corrupt stream
        raise DicomError(f"pydicom could not decode {transfer_syntax}: {e}") from e


def status() -> dict:
    return {"pydicom": pydicom_version(), "gdcm": gdcm_version()}
