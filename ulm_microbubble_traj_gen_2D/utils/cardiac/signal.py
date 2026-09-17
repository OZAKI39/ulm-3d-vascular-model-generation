"""Numerical helpers for sampling, smoothing, and synthetic ECG pulsatility."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import maximum_filter1d


def percentile_linear(values: np.ndarray, percentile: float) -> float:
    """Return a linearly interpolated percentile, matching MATLAB-style behavior."""

    clean = np.asarray(values, dtype=float).ravel()
    clean = clean[np.isfinite(clean)]
    if clean.size == 0:
        return float("nan")

    clean.sort()
    pct = min(max(float(percentile), 0.0), 100.0)
    if clean.size == 1:
        return float(clean[0])

    position = (clean.size - 1) * pct / 100.0
    lower = int(np.floor(position))
    upper = int(np.ceil(position))
    if lower == upper:
        return float(clean[lower])

    alpha = position - lower
    return float((1.0 - alpha) * clean[lower] + alpha * clean[upper])


def sample_from_pdf(
    domain: np.ndarray,
    pdf: np.ndarray,
    sample_size: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw samples from a discrete probability density by inverse-CDF sampling."""

    x = np.asarray(domain, dtype=float).ravel()
    p = np.asarray(pdf, dtype=float).ravel()
    if x.size != p.size:
        raise ValueError("domain and pdf must have the same length.")

    # Negative probability mass is not meaningful, so clamp it to zero.  If the
    # whole PDF is invalid or empty after clamping, use a uniform distribution.
    p = np.maximum(p, 0.0)
    cdf = np.cumsum(p)
    if cdf[-1] <= 0 or not np.isfinite(cdf[-1]):
        p = np.ones_like(p)
        cdf = np.cumsum(p)

    cdf = cdf / cdf[-1]
    rnd = rng.random(int(sample_size))
    return np.interp(rnd, cdf, x, left=0.0, right=x[-1])


def moving_average_rows(values: np.ndarray, fraction: float = 0.02) -> np.ndarray:
    """Smooth each column along the row axis with a short moving average window."""

    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2 or arr.shape[0] < 3:
        return arr

    window = max(3, int(round(arr.shape[0] * fraction)))
    if window % 2 == 0:
        window += 1

    pad = window // 2
    padded = np.pad(arr, ((pad, pad), (0, 0)), mode="edge")
    kernel = np.ones(window, dtype=float) / window
    out = np.empty_like(arr)
    for col in range(arr.shape[1]):
        out[:, col] = np.convolve(padded[:, col], kernel, mode="valid")
    return out


def generate_ecg_normalized(bpm: float, dt: float, duration_s: float) -> np.ndarray:
    """Create the ECG-like multiplicative velocity waveform used for pulsatility.

    The waveform follows the same P/QRS/T/U component idea as the original
    MATLAB simulator.  The raw synthetic ECG is shifted to be positive, scaled
    to roughly [0.5, 1.5], and passed through a short maximum filter to mimic the
    peak-envelope operation used by the MATLAB implementation.
    """

    x = np.arange(0.0, float(duration_s), float(dt), dtype=float)
    li = 30.0 / float(bpm)

    raw = (
        _p_wav(x, 0.25, 0.009, 0.016, li)
        + _q_wav(x, 0.025, 0.066, 0.166, li)
        + _qrs_wav(x, 1.0, 0.08, li)
        + _s_wav(x, 0.6, 0.066, 0.2, li)
        + _t_wav(x, 0.2, 0.142, 0.2, li)
        + _u_wav(x, 0.035, 0.0476, 0.433, li)
    )

    shifted = raw - np.min(raw)
    if np.max(shifted) > 0:
        shifted = shifted / np.max(shifted)
    shifted = shifted + 0.5
    envelope = maximum_filter1d(shifted, size=3, mode="nearest")
    return np.minimum(envelope, 1.5)


def _cos_series(x: np.ndarray, n: int, li: float, coefficient_fn) -> np.ndarray:
    """Evaluate a finite cosine series used by the synthetic ECG components."""

    out = np.zeros_like(x, dtype=float)
    for i in range(1, n + 1):
        out += coefficient_fn(i) * np.cos((i * np.pi * x) / li)
    return out


def _p_wav(x: np.ndarray, a: float, d: float, t: float, li: float) -> np.ndarray:
    """P-wave component of the synthetic ECG."""

    x = x + t
    b = (2.0 * li) / d

    def coeff(i: int) -> float:
        return (
            (
                np.sin((np.pi / (2.0 * b)) * (b - 2.0 * i)) / (b - 2.0 * i)
                + np.sin((np.pi / (2.0 * b)) * (b + 2.0 * i)) / (b + 2.0 * i)
            )
            * (2.0 / np.pi)
        )

    return a * ((1.0 / li) + _cos_series(x, 100, li, coeff))


def _q_wav(x: np.ndarray, a: float, d: float, t: float, li: float) -> np.ndarray:
    """Q-wave component of the synthetic ECG."""

    x = x + t
    b = (2.0 * li) / d
    q1 = (a / (2.0 * b)) * (2.0 - b)

    def coeff(i: int) -> float:
        return ((2.0 * b * a) / (i * i * np.pi * np.pi)) * (1.0 - np.cos((i * np.pi) / b))

    return -(q1 + _cos_series(x, 100, li, coeff))


def _qrs_wav(x: np.ndarray, a: float, d: float, li: float) -> np.ndarray:
    """QRS complex component of the synthetic ECG."""

    b = (2.0 * li) / d
    qrs1 = (a / (2.0 * b)) * (2.0 - b)

    def coeff(i: int) -> float:
        return ((2.0 * b * a) / (i * i * np.pi * np.pi)) * (1.0 - np.cos((i * np.pi) / b))

    return qrs1 + _cos_series(x, 100, li, coeff)


def _s_wav(x: np.ndarray, a: float, d: float, t: float, li: float) -> np.ndarray:
    """S-wave component of the synthetic ECG."""

    x = x - t
    b = (2.0 * li) / d
    s1 = (a / (2.0 * b)) * (2.0 - b)

    def coeff(i: int) -> float:
        return ((2.0 * b * a) / (i * i * np.pi * np.pi)) * (1.0 - np.cos((i * np.pi) / b))

    return -(s1 + _cos_series(x, 100, li, coeff))


def _t_wav(x: np.ndarray, a: float, d: float, t: float, li: float) -> np.ndarray:
    """T-wave component of the synthetic ECG."""

    x = x - t - 0.045
    b = (2.0 * li) / d

    def coeff(i: int) -> float:
        return (
            (
                np.sin((np.pi / (2.0 * b)) * (b - 2.0 * i)) / (b - 2.0 * i)
                + np.sin((np.pi / (2.0 * b)) * (b + 2.0 * i)) / (b + 2.0 * i)
            )
            * (2.0 / np.pi)
        )

    return a * ((1.0 / li) + _cos_series(x, 100, li, coeff))


def _u_wav(x: np.ndarray, a: float, d: float, t: float, li: float) -> np.ndarray:
    """U-wave component of the synthetic ECG."""

    x = x - t
    b = (2.0 * li) / d

    def coeff(i: int) -> float:
        return (
            (
                np.sin((np.pi / (2.0 * b)) * (b - 2.0 * i)) / (b - 2.0 * i)
                + np.sin((np.pi / (2.0 * b)) * (b + 2.0 * i)) / (b + 2.0 * i)
            )
            * (2.0 / np.pi)
        )

    return a * ((1.0 / li) + _cos_series(x, 100, li, coeff))
