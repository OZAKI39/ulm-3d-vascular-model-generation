"""
Initial 2D Poiseuille-like velocity field for planar lumen flow.
"""

from __future__ import annotations
import numpy as np
from scipy import ndimage
from ..core.types import RasterizedVessels


def build_initial_velocity_field(raster: RasterizedVessels, *, direction_smoothing_cells: float = 1.0) -> np.ndarray:
    """
    Build a 2D Poiseuille profile from local segment radius, flux, and centerline distance.
    """
    lumen           = np.asarray(raster.lumen_mask, dtype=bool)
    radius          = np.asarray(raster.radius_um, dtype=float)
    q2d             = np.asarray(raster.q2d_flow_um2_s, dtype=float)
    center_distance = np.asarray(raster.distance_to_centerline_um, dtype=float)
    direction       = _smoothed_direction(raster.direction_xz, lumen, float(direction_smoothing_cells))

    profile         = np.zeros(radius.shape, dtype=float)
    valid           = lumen & (radius > 0.0) & np.isfinite(center_distance)

    # Compute the Poiseuille profile: v = vmax * (1 - (r/R)^2)
    relative        = np.zeros(radius.shape, dtype=float)
    relative[valid] = np.clip(center_distance[valid] / np.maximum(radius[valid], np.finfo(float).eps), 0.0, 1.0)
    profile[valid]  = np.maximum(0.0, 1.0 - relative[valid] ** 2)

    # Calculate the maximum velocity from the 2D flow rate and the local radius
    speed           = np.zeros(radius.shape, dtype=float)
    speed[valid]    = 3.0 * q2d[valid] / (4.0 * np.maximum(radius[valid], np.finfo(float).eps)) * profile[valid]

    velocity            = direction * speed[..., None]
    velocity[~lumen]    = 0.0

    return velocity.astype(np.float32)


def _smoothed_direction(direction_xz: np.ndarray, lumen_mask: np.ndarray, sigma_cells: float) -> np.ndarray:
    """
    Our pipes are pixelated (rasterized), so the pipe routing may appear jagged at bends.
    Uses a Gaussian filter to process the direction vector.
    """
    direction = np.asarray(direction_xz, dtype=float)
    if float(sigma_cells) <= 0.0:
        return _normalize(direction)

    weight          = lumen_mask.astype(float)
    smoothed        = np.zeros_like(direction, dtype=float)
    weight_smooth   = ndimage.gaussian_filter(weight, sigma=float(sigma_cells), mode="nearest")
    for channel in range(2):
        smoothed[..., channel] = ndimage.gaussian_filter(direction[..., channel] * weight, sigma=float(sigma_cells), mode="nearest")
        smoothed[..., channel] = np.divide(
            smoothed[..., channel],
            np.maximum(weight_smooth, np.finfo(float).eps),
            out=np.zeros_like(smoothed[..., channel]),
            where=weight_smooth > 0.0,
        )
    smoothed[~lumen_mask] = 0.0
    return _normalize(smoothed)


def _normalize(vectors: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vectors, axis=-1, keepdims=True)
    return np.divide(vectors, np.maximum(norm, np.finfo(float).eps), out=np.zeros_like(vectors, dtype=float), where=norm > 0.0)
