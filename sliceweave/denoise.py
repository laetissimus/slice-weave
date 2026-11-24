"""Non-local means denoising (Buades, Coll and Morel 2005), vectorized over the image.

Each pixel becomes a weighted mean of the pixels in a search window around it; the weight of a
neighbor falls with the mean squared difference between the patches around the two pixels, after
subtracting what noise alone gives (2 sigma^2). Edges and fine structures are kept because only
pixels with similar surroundings are averaged.

``nlm2d`` works on each slice of a stack separately (axes 0 and 1); ``nlm3d`` on a volume.
"""
from __future__ import annotations

import itertools

import numpy as np
from scipy import ndimage as ndi


def _nlm(x: np.ndarray, sigma: float, h: float, patch: int, search: int, axes: tuple, mask=None,
         cancel=None) -> np.ndarray:
    x = np.asarray(x, np.float32)
    if sigma <= 0 or h <= 0:
        return x.copy()
    pad = [(search, search) if a in axes else (0, 0) for a in range(x.ndim)]
    xp = np.pad(x, pad, mode="reflect")
    size = [2 * patch + 1 if a in axes else 1 for a in range(x.ndim)]
    acc = np.zeros_like(x)
    wsum = np.zeros_like(x)
    wmax = np.zeros_like(x)
    var2 = np.float32(2.0 * sigma * sigma)
    inv_h2 = np.float32(1.0 / (h * h))
    rng = range(-search, search + 1)
    for off in itertools.product(*[rng if a in axes else (0,) for a in range(x.ndim)]):
        if not any(off):
            continue
        if cancel is not None:
            cancel()
        sl = tuple(slice(pad[a][0] + o, pad[a][0] + o + x.shape[a]) for a, o in enumerate(off))
        nb = xp[sl]
        d2 = ndi.uniform_filter((x - nb) ** 2, size=size, mode="reflect")
        w = np.exp(-np.maximum(d2 - var2, 0) * inv_h2)
        acc += w * nb
        wsum += w
        np.maximum(wmax, w, out=wmax)
    wmax = np.where(wsum > 0, wmax, 1.0).astype(np.float32)  # the pixel itself: the largest neighbor weight
    out = (acc + wmax * x) / (wsum + wmax)
    if mask is not None:
        out = np.where(mask, out, x)
    return out.astype(np.float32)


def nlm2d(stack: np.ndarray, sigma: float, strength: float = 1.0, patch: int = 1, search: int = 3,
          cancel=None) -> np.ndarray:
    """Slice by slice (in-plane axes 0, 1): 3 x 3 patches, 7 x 7 search window. h = 0.4 sigma x
    strength (Buades et al., IPOL 2011, for 3 x 3 patches)."""
    return _nlm(stack, sigma, 0.4 * sigma * strength, patch, search, (0, 1), cancel=cancel)


def nlm3d(vol: np.ndarray, sigma: float, strength: float = 1.0, patch: int = 1, search: int = 2, mask=None,
          cancel=None, chunk: int = 48) -> np.ndarray:
    """Volume: 3 x 3 x 3 patches, 5 x 5 x 5 search window. Works on slabs of ``chunk`` planes (with
    overlap) so memory stays at a few copies of one slab."""
    vol = np.asarray(vol, np.float32)
    nz = vol.shape[2]
    halo = search + patch
    if nz <= chunk + 2 * halo:
        return _nlm(vol, sigma, 0.4 * sigma * strength, patch, search, (0, 1, 2), mask=mask, cancel=cancel)
    out = np.empty_like(vol)
    for z0 in range(0, nz, chunk):
        z1 = min(z0 + chunk, nz)
        a, b = max(z0 - halo, 0), min(z1 + halo, nz)
        part = _nlm(vol[:, :, a:b], sigma, 0.4 * sigma * strength, patch, search, (0, 1, 2),
                    mask=None if mask is None else mask[:, :, a:b], cancel=cancel)
        out[:, :, z0:z1] = part[:, :, z0 - a:z0 - a + (z1 - z0)]
    return out
