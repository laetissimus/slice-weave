# SliceWeave methods

SliceWeave estimates one isotropic volume from several 2-D multi-slice stacks (typically
single-shot T2 HASTE / SSFSE / SSh-TSE in axial, coronal and sagittal planes). It follows the
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
* Acquisition order (for packets and temporal regularisation) comes from DICOM
  `AcquisitionTime`; without it an interleaved order (even slices, then odd slices) is assumed.
* Each slice s of stack k has a rigid correction `T[k][s]` (world -> world), a scalar intensity
  scale, and each stack a smooth multiplicative bias field.

## 2. Acquisition model

    y_p = s_slice * b(p) * (A x)_p + noise,     (A x)_p = sum_v PSF_p(v) x_v

The PSF is a 3-D Gaussian aligned with the (motion-corrected) slice: FWHM = 1.2 x pixel size
in-plane and FWHM = slice thickness through-plane (Kuklisova-Murgasova 2012; Ebner 2020).

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
  divided by the median tissue signal. The coefficient was chosen on simulated stacks with
  known truth (docs/VALIDATION.md).
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
5. **Final SRR** with more iterations, masking, and the average baseline.

## 6. Orientation of the output

The grid axes are the reference stack's axes, permuted and flipped to be closest to R, A and S
(so an AC-PC-angled axial reference gives an AC-PC-angled volume). The optional reorientation
step rotates the result about its centre in its own axes and resamples it with cubic B-splines;
the affine (NIfTI) and ImageOrientationPatient (DICOM) always describe the true position in the
scanner frame. For fetal studies this lets the volume be shown in fetal standard planes.

## 7. What it does not do (yet)

* No automatic fetal brain localisation/segmentation: the region is placed by the user, or
  brain masks from another tool (e.g. a NIfTI per stack named `*mask*`) are used.
* No deformable (non-rigid) motion, no slice-profile estimation, no learned priors.
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
