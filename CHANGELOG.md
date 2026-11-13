# Changelog

## 1.0.0 (2026-11)

First version.

* Reconstruction engine: Gaussian slice PSF applied in the Fourier domain, edge-preserving
  (Huber) or Tikhonov super-resolution with automatic regularisation from the noise level,
  stack-to-volume and hierarchical slice-to-volume registration (stacks, packets, slices), temporal
  repair and smoothing of slice motion, per-slice intensity matching, bias-field estimation and
  outlier-slice rejection.
* Adult / paediatric and fetal presets; region drawn on the viewer or brain masks from another tool.
* Local web app: case management, DICOM/NIfTI import with de-identification, stack table
  (include, reference, thickness), 3-plane viewer of inputs and result with compare mode,
  reorientation, report and slice-quality plots.
* Exports: NIfTI (SRR and average), DICOM (de-identified, or attached to the original study),
  research bundle, slice table.
* Command line (`recon`, `simulate`, `selftest`, `doctor`, `dicom-info`, `wipe`), Windows
  executable built by GitHub Actions, offline network guard.
