# Third-party components

## Libraries in the executable

| Component | Licence | Use |
|---|---|---|
| Python 3 | PSF | runtime |
| NumPy | BSD-3-Clause | arrays |
| SciPy | BSD-3-Clause | FFT, sparse matrices, interpolation, filters |
| Starlette | BSD-3-Clause | local web server |
| Uvicorn | BSD-3-Clause | local web server |
| python-multipart | Apache-2.0 | file upload parsing |
| pydicom (optional) | MIT | JPEG-family compressed DICOM |
| GDCM / python-gdcm (optional) | BSD-3-Clause | JPEG-family compressed DICOM |
| PyInstaller bootloader | GPL-2.0 with bootloader exception | packaging the executable |

## Bundled files

* `sliceweave/web/fonts/`: Barlow, Barlow Semi Condensed and IBM Plex Mono (SIL Open Font
  License 1.1, `OFL.txt`).
* `tests/data/dicom/`: small DICOM test files from pydicom's test data (licence in
  `LICENSE-pydicom.txt`). Used by the tests only.

No patient data, atlases or trained models are bundled. The demo and the self-test use
phantoms generated from formulas (`sliceweave/simulate.py`).

## Data used for validation (not bundled)

* Fetal HASTE sample (sub-01, 6 stacks with manual brain masks) from the MIAL Super-Resolution
  Toolkit repository (Medical-Image-Analysis-Laboratory/mialsuperresolutiontoolkit, `data/`),
  BSD-3-Clause.
* `chris_t1.nii.gz` (0.88 mm 3-D T1-weighted head scan) from the NiiVue repository
  (niivue/niivue), BSD-2-Clause.

## Methods

The algorithms follow the published SVR/SRR literature (see `docs/METHODS.md`); no code was
copied from NiftyMIC, SVRTK, MIALSRTK, BTK or NeSVoR.
