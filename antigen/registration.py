"""
Reusable image registration and polynomial surface helpers.

These utilities provide small, well-tested building blocks for tasks such as:
- Building a 2D registration image from a 1D profile
- Estimating subpixel shifts between image tiles (robust to small gaps/NaNs)
- Fitting low-order 2D polynomials for smooth shift surfaces
- Matching width of 1D profiles to a target size

They are intentionally lightweight and have no dependencies on Antigen internals,
so they can be imported from other modules without circular references.
"""
from __future__ import annotations

from typing import Tuple, Callable
import numpy as np

Array = np.ndarray


def match_width(profile_1d: Array, width: int) -> Array:
    """Pad or truncate a 1D array to the requested width.

    - If profile is longer than width, return the leading slice.
    - If profile is shorter, zero-pad the tail.
    - If equal, return the input view.
    """
    profile_1d = np.asarray(profile_1d, dtype=float)
    if profile_1d.ndim != 1:
        raise ValueError("profile must be 1D")
    n = int(profile_1d.shape[0])
    if width < 0:
        raise ValueError("width must be non-negative")
    if n == width:
        return profile_1d
    if n > width:
        return profile_1d[:width]
    out = np.zeros(width, dtype=float)
    out[:n] = profile_1d
    return out


def build_registration_image(profile_1d: Array, height: int) -> Array:
    """Tile a 1D profile into a 2D registration image of shape (height, width)."""
    profile_1d = np.asarray(profile_1d, dtype=float)
    if profile_1d.ndim != 1:
        raise ValueError("profile must be 1D")
    if height <= 0:
        raise ValueError("height must be positive")
    return np.tile(profile_1d[np.newaxis, :], (int(height), 1))


def estimate_subpixel_shift_1d_xcorr(ref: Array, src: Array, *, max_lag: int | None = None) -> Tuple[float, float]:
    """Estimate subpixel (dy, dx) shift mapping src -> ref using 1D cross-correlation.

    Approach (robust for small shifts, insensitive to NaNs):
    - Collapse each tile to a 1D profile along x using the median across rows.
    - Mean-center and standardize each profile; replace NaNs with 0 after centering.
    - Compute discrete cross-correlation within a limited lag window.
    - Fit a quadratic to the best peak and its neighbors to get subpixel dx.
    - dy is 0 in this 1D scheme by construction.
    """
    ref = np.asarray(ref, dtype=float)
    src = np.asarray(src, dtype=float)
    if ref.shape != src.shape:
        raise ValueError("ref and src tiles must have same shape")
    H, W = ref.shape
    if W < 3:
        raise ValueError("tile width too small for correlation (need >=3)")

    def _standardize(a: Array) -> Array:
        m = np.isfinite(a)
        if not np.any(m):
            raise ValueError("profile is all-NaN; cannot estimate shift")
        mu = np.nanmean(a[m])
        b = a - mu
        b = np.where(np.isfinite(b), b, 0.0)  # fill remaining NaNs post-centering
        sd = float(np.sqrt(np.mean(b ** 2)))  # RMS scale
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError("profile has zero or invalid variance; cannot estimate shift")
        return b / sd

    ref_1d = _standardize(np.nanmedian(ref, axis=0))
    src_1d = _standardize(np.nanmedian(src, axis=0))

    if max_lag is None:
        max_lag = int(min(12, W // 4))
    lags = np.arange(-max_lag, max_lag + 1)
    c = np.empty(lags.size, dtype=float)
    for i, lag in enumerate(lags):
        if lag >= 0:
            a = ref_1d[lag:]
            b = src_1d[: W - lag]
        else:
            a = ref_1d[: W + lag]
            b = src_1d[-lag:]
        c[i] = -np.inf if a.size < 3 else float(np.dot(a, b) / max(a.size, 1))

    k = int(np.argmax(c))
    lag0 = int(lags[k])

    def _quad_subpixel(y_minus: float, y0: float, y_plus: float) -> float:
        denom = (y_minus - 2.0 * y0 + y_plus)
        if not np.isfinite(denom) or abs(denom) < 1e-12:
            return 0.0
        return 0.5 * (y_minus - y_plus) / denom

    if 0 < k < len(c) - 1 and np.isfinite(c[k - 1]) and np.isfinite(c[k]) and np.isfinite(c[k + 1]):
        delta = np.clip(_quad_subpixel(c[k - 1], c[k], c[k + 1]), -0.5, 0.5)
    else:
        delta = 0.0
    dy = 0.0
    dx = float(lag0) + float(delta)
    return dy, dx


def fit_poly2d_least_squares(x: Array, y: Array, z: Array, order: int = 3) -> Callable[[Array, Array], Array]:
    """Fit a 2D polynomial surface z = f(x, y) up to total degree 'order'.

    Returns an evaluator function eval_fn(X, Y) broadcasting over inputs.
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    if any(arr.ndim != 1 for arr in (x, y, z)):
        raise ValueError("x, y, z must be 1D")
    if not (len(x) == len(y) == len(z)):
        raise ValueError("x, y, z lengths must match")
    terms: list[tuple[int, int]] = []
    for i in range(order + 1):
        for j in range(order + 1 - i):
            terms.append((i, j))
    A = np.column_stack([(x ** i) * (y ** j) for (i, j) in terms])
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)

    def eval_fn(X, Y):
        X = np.asarray(X, dtype=float)
        Y = np.asarray(Y, dtype=float)
        out = np.zeros_like(X, dtype=float)
        for c, (i, j) in zip(coef, terms):
            out += c * (X ** i) * (Y ** j)
        return out

    return eval_fn
