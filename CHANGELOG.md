# Changelog

## 1.1.0 (2026-11)

* **2-D FLAIR preset** (adult / paediatric): multi-shot FLAIR stacks in three planes. Motion between
  series is corrected; slice-to-volume correction is refused (a multi-shot slice has no single
  position). FLAIR series are recognised by name (FLAIR, TIRM, dark fluid, ...), Siemens
  `*tir2d`, or TI > 1500 ms with TE > 60 ms. Import, the viewer and the command line keep one
  contrast per reconstruction (switching mode selects the FLAIR or the single-shot stacks).
  The final FLAIR volume is regularised less (0.35 x the noise-calibrated weight; registration
  targets keep the full weight) to keep the contrast of small lesions.
* **Acquired in-plane resolution** (all modes): read from DICOM AcquisitionMatrix +
  InPlanePhaseEncodingDirection or from a dcm2niix sidecar (AcquisitionMatrixPE / ReconMatrixPE /
  BaseResolution). The slice PSF uses it instead of the (often interpolated) pixel size; the noise
  is estimated at the acquired spacing (zero-filled images have correlated noise that the
  pixel-level estimate under-reads by up to 10x).
* Plan warnings: mixed contrasts, different TE/TR/TI between the selected series, fewer than three
  planes.
* FLAIR phantom (12 white-matter lesions, 2-8 mm) and simulation of a zero-filled 256 x 192
  acquisition; validation case `phantom-flair` with lesion contrast recovery.
* CLI: `--mode flair`; `recon` now keeps only the series of the chosen kind when a whole study is
  given; `simulate --kind flair`. Viewer: FLAIR demo.
* Cycles: with 0 cycles or without motion correction, intensity correction and the outlier check
  now run once when switched on (before, they were skipped unless outlier rejection was on). The
  viewer explains the setting for the chosen motion correction and disables it without motion
  correction. README section "Cycles"; `tools/validate.py --only cycles` measures its effect.
* **Grey-white matter and sulci (adult / paediatric T2):** the final volume is regularised at
  0.25 x the noise-calibrated weight (the motion-correction volumes keep the full weight). On
  pseudo-T2 from a real head scan (0.6 x 0.6 x 3 mm child protocol) sulcal CSF contrast
  0.90 -> 1.01, cortex detail 0.955 -> 0.973, PSNR +1.0 dB, SSIM 0.93 -> 0.96; infant grey-white
  contrast 0.50 -> 0.64. Adult phantom SRR 22.5 -> 23.2 dB.
* Viewer **Brain window** (default; centred on the parenchyma, CSF near white) next to Full range;
  DICOM exports carry the brain window first and the full range second.
* Warning when the output voxel is much coarser than the acquired in-plane resolution; switching
  to slice-to-volume raises an unchanged low cycle count to 6 (and back); the slice check warns
  when many slices were excluded (suggests slice-to-volume, or another reference stack when the
  reference lost many slices).
* Validation case `cortex` (grey-white contrast, sulcal CSF contrast, cortical detail; final
  weight, voxel size, a child moving during each series). README section "Grey-white matter and
  sulci".
* Windows build also runs on pull requests.
* Published releases (or a manual run with `release_tag`) get the tested executables attached;
  the build adds THIRD-PARTY-NOTICES.txt with the licences of the bundled Python and libraries.

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
