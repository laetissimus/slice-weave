# SliceWeave methods

SliceWeave estimates one isotropic volume from several 2-D multi-slice stacks (single-shot T2
HASTE / SSFSE / SSh-TSE, or 2-D FLAIR, in axial, coronal and sagittal planes). It follows the
slice-to-volume reconstruction (SVR) and super-resolution reconstruction (SRR) framework of
Rousseau et al. (2006), Jiang et al. (2007), Gholipour et al. (2010), Kuklisova-Murgasova et al.
(2012) and Ebner et al. (2020, NiftyMIC). It is an independent implementation in NumPy/SciPy
(no code from those packages), CPU only, designed to run offline in one executable.

## 1. Data model

* Each stack keeps every slice's own position (`slice_aff[s]`: pixel (i, j) -> world mm, RAS+).
  Uneven gaps and interleaved packets are never resampled before reconstruction.
* Slice thickness comes from DICOM `SliceThickness`, a dcm2niix/BIDS JSON sidecar, or the
  user; without one it is taken as the slice spacing. It sets the through-plane point-spread
  function (PSF), so a wrong value blurs or sharpens the result. Check it in the stack table.
* In-plane resolution: images are often interpolated (a 256 x 192 acquisition shown at 512 x 512),
  so the pixel size overstates the resolution. The acquired resolution along each image axis is
  read from DICOM `AcquisitionMatrix` and `InPlanePhaseEncodingDirection` (FOV / acquired
  samples; `PercentPhaseFieldOfView` for square images), or from a dcm2niix sidecar
  (`AcquisitionMatrixPE`, `ReconMatrixPE`, `PhaseEncodingDirection`, `BaseResolution`). It is
  never taken finer than the pixel, and an implausible value (more than 4 pixels) is ignored.
  Without these fields the pixel size is used.
* Acquisition order (for packets and temporal regularisation) comes from DICOM
  `AcquisitionTime`; without it an interleaved order (even slices, then odd slices) is assumed.
* Each slice s of stack k has a rigid correction `T[k][s]` (world -> world), a scalar intensity
  scale, and each stack a smooth multiplicative bias field.

## 2. Acquisition model

    y_p = s_slice * b(p) * (A x)_p + noise,     (A x)_p = sum_v PSF_p(v) x_v

The PSF is a 3-D Gaussian aligned with the (motion-corrected) slice: FWHM = 1.2 x acquired
in-plane resolution along each in-plane axis (the pixel size when the resolution is unknown) and
FWHM = slice thickness through-plane (Kuklisova-Murgasova 2012; Ebner 2020).

Implementation: slices are grouped by orientation (within 12 degrees). For each group the PSF
is applied to the whole volume in the Fourier domain (exact for any orientation; the trilinear
interpolation's own smoothing is subtracted), then the volume is sampled at the pixel centres with
trilinear interpolation (one sparse matrix per group). The adjoint applies the same steps in
reverse. A dot-product test of the operator and its adjoint is part of the test suite.

Pixels much smaller than the output voxel are averaged in b x b blocks first (for example 0.45 mm
pixels for a 1 mm output), with the PSF widened accordingly.

## 3. Super-resolution

    x = argmin  sum_p w_p (y_p / (s_p b_p) - (A x)_p)^2 + lambda R(x)

* `R` is an edge-preserving Huber penalty on the forward differences (iteratively reweighted,
  3 reweightings) or first-order Tikhonov. Default: Huber.
* `lambda = alpha * mean(diag(A^T W A))`, so alpha does not depend on the image scale or the
  number of stacks. With "Automatic" regularisation alpha = 0.1 + 65 x (relative noise), where the
  relative noise is estimated from each stack (median absolute 2-D Laplacian, Immerkaer 1996)
  divided by the median tissue signal. For interpolated images the Laplacian is computed on
  pixels sampled at about the acquired spacing: zero-filled interpolation correlates neighbouring
  pixels and the pixel-level estimate would read the noise up to 10 times too low. The coefficient
  was chosen on simulated stacks with known truth (docs/VALIDATION.md).
* The volumes that serve as registration targets use this weight; the final volume uses a
  fraction of it: adult / paediatric T2 0.25 x (keeps the sulci and the cortical ribbon; docs/
  VALIDATION.md, "Grey-white matter and sulci"), FLAIR 0.35 x (keeps the contrast of small
  lesions; "2-D FLAIR"), fetal 1.0 x. A fixed alpha applies to all volumes.
* Solver: preconditioned conjugate gradients (Jacobi preconditioner), non-negativity by
  clipping at the end.
* Output masking: voxels with too little data (coverage `A^T W 1` below 15 % of the typical
  in-region value) fade to zero; in fetal mode the result is also limited to the region.

The "average" volume exported alongside is `A^T W y / A^T W 1`: the PSF-weighted mean of the
same corrected slices, without deconvolution. It is the natural baseline for comparison.

## 4. Registration

A set of pixels (a stack, a packet or one slice) with positions p and values y is aligned to
the current volume by a rigid D minimising `sum w (y - a V(D p) - b)^2`, with a, b re-fitted at
every step (equivalent to normalised cross-correlation; anti-correlated solutions are refused),
V the volume blurred with the slices' PSF, and Huber weights on the residuals. Gauss-Newton steps
with Marquardt damping and step halving; coarse-to-fine (extra Gaussian smoothing of both the
volume and the slices, e.g. 4 -> 2 -> 1 mm). An update is kept only if it improves the NCC at
the finest level and stays within limits (stacks 25 deg / 25 mm, slices 10 deg / 8 mm per
cycle).

## 5. Pipeline

1. **Prepare.** Crop each stack to the region (+ margin), bin small pixels, normalise each
   stack by its median tissue signal, estimate noise, flag slices much darker or brighter than
   their neighbours (global signal loss) as initially excluded.
2. **Stacks.** Initial volume (adult: reference stack only; fetal: all stacks), then rigid
   stack-to-volume registration of every other stack. The reference stack defines the frame.
3. **First SRR** from all stacks.
4. **Motion cycles** (adult default 1, fetal default 6):
   * "Between series" (adult): stack-to-volume registration again, against the SRR.
   * "Slice-to-volume" (fetal; optional for a moving child): cycle 1 registers whole stacks,
     then packets (2, then 4 per stack) of slices acquired together; cycle 2 smaller packets
     (8 per stack) then single slices; later cycles single slices, coarse to fine.
     After each slice pass: (a) temporal repair - a slice whose pose differs from the median of its
     +/- 2 neighbours in acquisition time by more than 6 deg or 4 mm gets that median;
     (b) temporal smoothing - each slice's motion parameters are replaced by the median of itself
     and its two time neighbours (a median keeps sudden jumps). The reference stack is then
     re-anchored so the volume frame does not drift.
   * Per-slice intensity scale (least squares against the simulated slice) and a smooth bias
     field per stack (normalised Gaussian convolution of the log ratio, 20 mm).
   * Outlier slices: excluded if the NCC between the slice and the simulated slice is below a
     threshold (fetal 0.5, rising by 0.1 per cycle to 0.7; adult: only robust outliers, below
     min(0.6, median - 4 MAD)), or if the slice scale is outside 0.4-2.5. Excluded slices are
     re-tested every cycle. Pixel weights: Huber on the residuals.
   * SRR (warm start).
   With 0 cycles, or without motion correction, one pass of intensity correction and the outlier
   test still runs (when switched on). Effect of the number of cycles: docs/VALIDATION.md.
5. **Final SRR** with more iterations, masking, and the average baseline.

## 6. Modes

| | Adult / paediatric T2 | 2-D FLAIR | Fetal |
|---|---|---|---|
| Series | single-shot T2 | 2-D FLAIR (multi-shot inversion recovery) | single-shot T2 |
| Region | head, automatic | head, automatic | user-placed ellipsoid or brain masks |
| Motion | between series (slice-to-volume optional) | between series only | slice-to-volume |
| Output voxel | 1.0 mm | 1.0 mm | 0.8 mm |
| Initial volume | reference stack | reference stack | all stacks |
| Outlier slices | only robust outliers | only robust outliers | NCC 0.5 -> 0.7 |
| Final regularisation | 0.25 x noise-calibrated | 0.35 x noise-calibrated | noise-calibrated |

**Why FLAIR is handled differently.** A single-shot slice is acquired in under a second, so it is
a snapshot with one position, and slices can be repositioned one by one. A multi-shot 2-D FLAIR
slice is filled over the whole series (several echo trains per slice, interleaved with the other
slices, minutes in total); motion during the series produces ghosting and blurring within the
slice, which no rigid repositioning can undo. SliceWeave therefore corrects only the motion
between the three FLAIR series and refuses slice-to-volume correction in this mode (viewer,
command line and engine). A FLAIR series degraded by motion should be excluded.

**Recognition.** A series is FLAIR-like when its description, protocol or sequence name contains
FLAIR, TIRM, dark fluid, T2_FL, IR-FSE (and similar), when the Siemens sequence name is `*tir2d`,
or when TI > 1500 ms with TE > 60 ms. T1-FLAIR series named as such are also FLAIR-like (same
multi-shot inversion-recovery handling). One reconstruction combines one contrast: import and the
command line select the FLAIR series in FLAIR mode and the single-shot series otherwise; switching
mode in the viewer switches the selection.

**Plan warnings** (shown before running): FLAIR and non-FLAIR series selected together; selected
series differing in TE (> 15 %), TR (> 25 %) or TI (> 10 %); fewer than three stacks or fewer than
three planes.

## 7. Orientation of the output

The grid axes are the reference stack's axes, permuted and flipped to be closest to R, A and S
(so an AC-PC-angled axial reference gives an AC-PC-angled volume). The optional reorientation
step rotates the result about its centre in its own axes and resamples it with cubic B-splines;
the affine (NIfTI) and ImageOrientationPatient (DICOM) always describe the true position in the
scanner frame. For fetal studies this lets the volume be shown in fetal standard planes.

## 8. Display window

The viewer opens every image in a *Brain* window: the most common intensity inside the head
above 12 % of the bright end (the parenchyma peak of the histogram) is found, and the window runs
from 0.4 x to 1.8 x that value. Grey and white matter then spread over the grey levels and CSF
is near white, instead of the parenchyma being squeezed into the lower half by bright CSF.
*Full range* (0 to the 99th percentile) and right-drag (width, level) remain. DICOM exports carry
both windows (BRAIN first, then FULL; WindowCenterWidthExplanation), so PACS opens the brain
window. The window changes the display only, not the stored values.

## 9. What it does not do (yet)

* No automatic fetal brain localisation/segmentation: the region is placed by the user, or
  brain masks from another tool (e.g. a NIfTI per stack named `*mask*`) are used.
* No deformable (non-rigid) motion, no slice-profile estimation, no learned priors.
* No correction of motion within a multi-shot (FLAIR) series, and no modelling of inflow or
  CSF-pulsation artefacts, which differ between planes (they are down-weighted as outliers at best).
* No atlas-based standard-plane reorientation (manual rotation only).
* CPU only; fetal cases take minutes.

## References

* Rousseau F, et al. Registration-based approach for reconstruction of high-resolution in utero
  fetal MR brain images. Acad Radiol 2006;13:1072-1081.
* Jiang S, et al. MRI of moving subjects using multislice snapshot images with volume
  reconstruction (SVR). IEEE TMI 2007;26:967-980.
* Gholipour A, Estroff JA, Warfield SK. Robust super-resolution volume reconstruction from slice
  acquisitions. IEEE TMI 2010;29:1739-1758.
* Kuklisova-Murgasova M, et al. Reconstruction of fetal brain MRI with intensity matching and
  complete outlier removal. Med Image Anal 2012;16:1550-1564. (SVRTK lineage)
* Tourbier S, et al. An efficient total variation algorithm for super-resolution in fetal brain
  MRI with adaptive regularization. NeuroImage 2015;118:584-597. (MIALSRTK)
* Ebner M, et al. An automated framework for localization, segmentation and super-resolution
  reconstruction of fetal brain MRI. NeuroImage 2020;206:116324. (NiftyMIC)
* Xu J, et al. NeSVoR: implicit neural representation for slice-to-volume reconstruction in MRI.
  IEEE TMI 2023;42:1707-1719.
* Immerkaer J. Fast noise variance estimation. CVIU 1996;64:300-302.
