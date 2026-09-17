"""Dimension-independent vector primitives used by the XYZ generator."""

from __future__ import annotations

import math

import numpy as np


def _norm(vector: np.ndarray) -> float:
    return math.sqrt(float(np.dot(vector, vector)))


def _distance(first: np.ndarray, second: np.ndarray) -> float:
    return _norm(np.asarray(first, dtype=float) - np.asarray(second, dtype=float))


def _angle_deg(first: np.ndarray, second: np.ndarray) -> float:
    first_norm = _norm(first)
    second_norm = _norm(second)
    if first_norm < 1.0e-12 or second_norm < 1.0e-12:
        return 0.0
    cosine = float(np.dot(first, second)) / (first_norm * second_norm)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _angles_deg(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    if first.size == 0:
        return np.empty(0, dtype=float)
    first_norms = np.linalg.norm(first, axis=1)
    if second.ndim == 1:
        second_norms = float(_norm(second))
        dots = first @ second
        denominator = first_norms * second_norms
    else:
        second_norms = np.linalg.norm(second, axis=1)
        dots = np.einsum("ij,ij->i", first, second)
        denominator = first_norms * second_norms
    cosines = np.ones_like(first_norms)
    valid = denominator >= 1.0e-12
    cosines[valid] = dots[valid] / denominator[valid]
    return np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
