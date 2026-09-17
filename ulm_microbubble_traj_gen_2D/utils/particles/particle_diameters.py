"""Per-bubble diameter configuration and reproducible sampling."""

from __future__ import annotations

from ..core.config import ParticleConfig


def particle_diameter_bounds_um(cfg: ParticleConfig) -> tuple[float, float]:
    """Return the validated inclusive diameter range in micrometers."""

    return float(cfg.bubble_diameter_min_um), float(cfg.bubble_diameter_max_um)
