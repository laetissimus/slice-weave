"""Image-similarity metrics used by the self-test and the validation scripts."""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def fit_scale(x: np.ndarray, ref: np.ndarray, mask: np.ndarray) -> float:
    a, b = x[mask].astype(np.float64), ref[mask].astype(np.float64)
    return float((a @ b) / max(a @ a, 1e-12))


def psnr(x: np.ndarray, ref: np.ndarray, mask: np.ndarray, data_range: float | None = None) -> float:
    dr = data_range or float(np.percentile(ref[mask], 99.9))
    mse = float(np.mean((x[mask].astype(np.float64) - ref[mask]) ** 2))
    return float(10 * np.log10(dr ** 2 / max(mse, 1e-20)))


def ssim(x: np.ndarray, ref: np.ndarray, mask: np.ndarray, sigma: float = 1.5, data_range: float | None = None) -> float:
    """Mean SSIM (Gaussian window) over the mask (Wang et al. 2004)."""
    dr = data_range or float(np.percentile(ref[mask], 99.9))
    c1, c2 = (0.01 * dr) ** 2, (0.03 * dr) ** 2
    x = x.astype(np.float64)
    y = ref.astype(np.float64)
    f = lambda a: ndi.gaussian_filter(a, sigma)  # noqa: E731
    mx, my = f(x), f(y)
    sxx = f(x * x) - mx * mx
    syy = f(y * y) - my * my
    sxy = f(x * y) - mx * my
    s = ((2 * mx * my + c1) * (2 * sxy + c2)) / ((mx ** 2 + my ** 2 + c1) * (sxx + syy + c2))
    return float(np.mean(s[mask]))


def ncc(a: np.ndarray, b: np.ndarray, mask: np.ndarray | None = None) -> float:
    if mask is not None:
        a, b = a[mask], b[mask]
    a = a.astype(np.float64) - a.mean()
    b = b.astype(np.float64) - b.mean()
    d = np.sqrt((a @ a) * (b @ b))
    return float((a @ b) / d) if d > 0 else 0.0


def compare(x: np.ndarray, truth: np.ndarray, mask: np.ndarray) -> dict:
    s = fit_scale(x, truth, mask)
    xs = x * s
    return {"psnr_db": round(psnr(xs, truth, mask), 2), "ssim": round(ssim(xs, truth, mask), 4),
            "ncc": round(ncc(xs, truth, mask), 4)}
