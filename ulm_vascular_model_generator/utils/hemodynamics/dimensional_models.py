"""Dimension-specific haemodynamic models for vascular generation."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class Planar2DHaemodynamics:
    """True two-dimensional strip-flow model, per unit out-of-plane depth."""

    inlet_flux_um2_s: float
    root_half_width_um: float
    minimum_half_width_um: float
    murray_gamma: float
    velocity_min_mm_s: float
    velocity_max_mm_s: float

    geometry_mode: str = "planar_2d"
    flow_quantity: str = "planar_flux_per_unit_depth"
    flow_unit: str = "um2/s"
    size_name: str = "half-width"

    @classmethod
    def from_config(cls, cfg):
        return cls(
            inlet_flux_um2_s=float(cfg.planar_inlet_total_flux_um2_s),
            root_half_width_um=float(cfg.planar_root_half_width_um),
            minimum_half_width_um=float(cfg.planar_min_half_width_um),
            murray_gamma=float(cfg.planar_murray_gamma),
            velocity_min_mm_s=float(cfg.planar_velocity_min_mm_s),
            velocity_max_mm_s=float(cfg.planar_velocity_max_mm_s),
        )

    @property
    def inlet_value(self) -> float:
        return self.inlet_flux_um2_s

    @property
    def root_size_um(self) -> float:
        return self.root_half_width_um

    @property
    def minimum_size_um(self) -> float:
        return self.minimum_half_width_um

    def mean_velocity_um_s(self, flux_um2_s: float, half_width_um: float) -> float:
        return float(flux_um2_s) / (2.0 * float(half_width_um))

    def lumen_measure_cost(self, vessels) -> float:
        """Return total lumen area per unit out-of-plane depth [um^2]."""
        return 2.0 * math.fsum(
            vessel.length() * vessel.radius for vessel in vessels
        )

    def reference_lumen_measure(self, length_um: float) -> float:
        return 2.0 * self.root_half_width_um * float(length_um)

    def transport_metadata(self) -> dict[str, object]:
        return {
            "geometry_semantics": "planar_xz_lumen_strips_per_unit_depth",
            "flow_model": "planar_2d_fixed_total_flux_equal_terminal_shares",
            "size_model": "planar_murray_fixed_root_half_width",
            "flow_quantity": self.flow_quantity,
            "flow_unit": self.flow_unit,
            "planar_inlet_total_flux_um2_s": self.inlet_flux_um2_s,
            "planar_root_half_width_um": self.root_half_width_um,
            "planar_murray_gamma": self.murray_gamma,
        }


@dataclass(frozen=True)
class Volumetric3DHaemodynamics:
    """Three-dimensional circular-tube volume-flow model."""

    inlet_flow_um3_s: float
    root_radius_um: float
    minimum_radius_um: float
    murray_gamma: float
    velocity_min_mm_s: float
    velocity_max_mm_s: float

    geometry_mode: str = "volumetric_3d"
    flow_quantity: str = "volume_flow"
    flow_unit: str = "um3/s"
    size_name: str = "radius"

    @classmethod
    def from_config(cls, cfg):
        return cls(
            inlet_flow_um3_s=float(cfg.volumetric_inlet_total_flow_um3_s),
            root_radius_um=float(cfg.volumetric_root_radius_um),
            minimum_radius_um=float(cfg.volumetric_min_radius_um),
            murray_gamma=float(cfg.volumetric_murray_gamma),
            velocity_min_mm_s=float(cfg.volumetric_velocity_min_mm_s),
            velocity_max_mm_s=float(cfg.volumetric_velocity_max_mm_s),
        )

    @property
    def inlet_value(self) -> float:
        return self.inlet_flow_um3_s

    @property
    def root_size_um(self) -> float:
        return self.root_radius_um

    @property
    def minimum_size_um(self) -> float:
        return self.minimum_radius_um

    def mean_velocity_um_s(self, flow_um3_s: float, radius_um: float) -> float:
        radius = float(radius_um)
        return float(flow_um3_s) / (math.pi * radius * radius)

    def lumen_measure_cost(self, vessels) -> float:
        """Return total circular-cylinder lumen volume [um^3]."""
        return math.pi * math.fsum(
            vessel.length() * vessel.radius * vessel.radius for vessel in vessels
        )

    def reference_lumen_measure(self, length_um: float) -> float:
        return math.pi * self.root_radius_um**2 * float(length_um)

    def transport_metadata(self) -> dict[str, object]:
        return {
            "geometry_semantics": "three_dimensional_circular_cylinders",
            "flow_model": "volumetric_3d_fixed_total_flow_equal_terminal_shares",
            "size_model": "volumetric_murray_fixed_root_radius",
            "flow_quantity": self.flow_quantity,
            "flow_unit": self.flow_unit,
            "volumetric_inlet_total_flow_um3_s": self.inlet_flow_um3_s,
            "volumetric_root_radius_um": self.root_radius_um,
            "volumetric_murray_gamma": self.murray_gamma,
        }


def haemodynamics_for_config(cfg):
    """Construct exactly one dimensional model from an already validated config."""
    if cfg.geometry_mode == "planar_2d":
        return Planar2DHaemodynamics.from_config(cfg)
    if cfg.geometry_mode == "volumetric_3d":
        return Volumetric3DHaemodynamics.from_config(cfg)
    raise ValueError(f"Unsupported geometry_mode: {cfg.geometry_mode!r}")
