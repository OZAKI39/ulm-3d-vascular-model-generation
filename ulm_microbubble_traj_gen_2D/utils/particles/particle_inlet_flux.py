"""Physical root-inlet cross sections and finite-size microbubble number flux."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel

from ..core.config import ParticleConfig
from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from .particle_diameters import particle_diameter_bounds_um
from .particle_field_sampling import sample_bilinear
from .particle_geometry_tolerance import (
    project_roundoff_negative_wall_gap_xz_um,
)
from ..core.types import FlowField, GridDomain


UM3_TO_M3 = 1.0e-18
MB_PER_ML_TO_MB_PER_M3 = 1.0e6


@dataclass(frozen=True)
class InjectionSection:
    """One fixed line section normal to a root-vessel flow direction."""

    vessel_id: int
    center_xz_um: np.ndarray
    inward_normal_xz: np.ndarray
    transverse_coordinate_um: np.ndarray
    positions_xz_um: np.ndarray
    positions_grid: np.ndarray
    inward_speed_um_s: np.ndarray
    wall_distance_um: np.ndarray
    inlet_connected_clearance_um: np.ndarray
    inside_lumen: np.ndarray


@dataclass(frozen=True)
class InletFluxModel:
    """Normalized finite-size joint flux and its physical injection rate."""

    boundary_geometry: ContinuousVesselGeometry
    sections: tuple[InjectionSection, ...]
    radius_grid_um: np.ndarray
    radius_cdf: np.ndarray
    radius_min_um: float
    radius_max_um: float
    number_concentration_mb_per_ml: float
    number_concentration_mb_per_m3: float
    effective_thickness_um: float
    injection_rate_per_s: float
    mean_injection_interval_s: float
    raw_section_flow_um2_s: float
    size_accessible_section_flow_um2_s: float
    reference_inlet_flow_um2_s: float
    section_flux_relative_error: float

    def sample_radius_um(self, quantile: float) -> float:
        """Invert the flux-weighted marginal radius CDF deterministically."""

        if self.radius_min_um == self.radius_max_um:
            return float(self.radius_min_um)
        value = min(max(float(quantile), np.finfo(float).eps), 1.0 - np.finfo(float).eps)
        cdf, indices = np.unique(self.radius_cdf, return_index=True)
        radii = self.radius_grid_um[indices]
        return float(np.interp(value, cdf, radii))

    def sample_position_grid(self, radius_um: float, quantile: float) -> np.ndarray:
        """Invert the conditional flux CDF along all admissible inlet sections."""

        return self.sample_position_and_vessel_id(radius_um, quantile)[0]

    def sample_position_and_vessel_id(
        self, radius_um: float, quantile: float
    ) -> tuple[np.ndarray, int]:
        """Return one inlet position and its one-based persistent root owner."""

        section_intervals: list[tuple[InjectionSection, np.ndarray]] = []
        totals: list[float] = []
        for section in self.sections:
            interval_flux = _valid_interval_flux(section, float(radius_um))
            section_intervals.append((section, interval_flux))
            totals.append(float(np.sum(interval_flux)))
        total = float(np.sum(totals))
        if not math.isfinite(total) or total <= 0.0:
            raise ValueError(f"No inlet location can admit a bubble with radius {radius_um:.6g} um.")

        target = min(max(float(quantile), 0.0), 1.0 - np.finfo(float).eps) * total
        section_offset = 0.0
        selected_section = self.sections[0]
        selected_flux = section_intervals[0][1]
        for (section, interval_flux), section_total in zip(section_intervals, totals):
            if target < section_offset + section_total or section is self.sections[-1]:
                selected_section = section
                selected_flux = interval_flux
                target -= section_offset
                break
            section_offset += section_total

        cumulative = np.cumsum(selected_flux)
        interval_index = int(np.searchsorted(cumulative, target, side="right"))
        interval_index = min(interval_index, selected_flux.size - 1)
        previous = float(cumulative[interval_index - 1]) if interval_index > 0 else 0.0
        local_target = max(0.0, target - previous)
        s = selected_section.transverse_coordinate_um
        ds = float(s[interval_index + 1] - s[interval_index])
        u0 = float(selected_section.inward_speed_um_s[interval_index])
        u1 = float(selected_section.inward_speed_um_s[interval_index + 1])
        fraction = _linear_density_fraction(u0, u1, ds, local_target)
        transverse_position = float(s[interval_index]) + fraction * ds
        position_grid = np.asarray(
            [
                np.interp(
                    transverse_position,
                    selected_section.transverse_coordinate_um,
                    selected_section.positions_grid[:, 0],
                ),
                np.interp(
                    transverse_position,
                    selected_section.transverse_coordinate_um,
                    selected_section.positions_grid[:, 1],
                ),
            ],
            dtype=np.float64,
        )
        position_xz_um = self.boundary_geometry.grid_to_world_xz(position_grid)
        position_xz_um, true_gap_um, _ = project_roundoff_negative_wall_gap_xz_um(
            position_xz_um,
            float(radius_um),
            self.boundary_geometry,
        )
        position_grid = np.asarray(
            self.boundary_geometry.world_xz_to_grid(position_xz_um),
            dtype=np.float64,
        )
        admitted = (
            bool(self.boundary_geometry.contains_xz_um(position_xz_um))
            and bool(self.boundary_geometry.is_accessible_grid(position_grid, radius_um))
            and true_gap_um >= 0.0
        )
        if not admitted:
            raise RuntimeError(
                "The inverse inlet-flux sample left the canonical finite-radius inlet domain. "
                "Refine the inlet-section quadrature before running particle transport."
            )
        return position_grid, int(selected_section.vessel_id) + 1


def build_inlet_flux_model(
    domain: GridDomain,
    flow_field: FlowField,
    vessels: list[Vessel],
    boundary_geometry: ContinuousVesselGeometry,
    cfg: ParticleConfig,
    *,
    effective_thickness_um: float,
    boundary_depth_cells: float,
) -> InletFluxModel:
    """Build ``j(s,R)`` and its deterministic inverse distributions."""

    if boundary_geometry.shape != tuple(domain.shape):
        raise ValueError("boundary_geometry must describe the same grid as domain.")
    concentration_ml = float(cfg.inlet_number_concentration_mb_per_ml)
    if not math.isfinite(concentration_ml) or concentration_ml <= 0.0:
        raise ValueError("Continuous perfusion requires a positive inlet MB/mL concentration.")
    thickness = float(effective_thickness_um)
    if not math.isfinite(thickness) or thickness <= 0.0:
        raise ValueError("The effective flow thickness must be finite and positive.")

    roots = [vessel for vessel in vessels if int(vessel.parent_id) < 0]
    if not roots:
        raise ValueError("Continuous perfusion requires at least one root vessel.")
    sections = tuple(
        _build_root_section(
            domain,
            flow_field,
            boundary_geometry,
            vessel,
            boundary_depth_cells,
        )
        for vessel in roots
    )

    diameter_min_um, diameter_max_um = particle_diameter_bounds_um(cfg)
    radius_min = 0.5 * diameter_min_um
    radius_max = 0.5 * diameter_max_um
    if radius_min == radius_max:
        radius_grid = np.asarray([radius_min], dtype=np.float64)
        radius_cdf = np.asarray([1.0], dtype=np.float64)
        accessible_q2d = _section_flux_for_radius(sections, radius_min)
    else:
        radius_grid = np.linspace(radius_min, radius_max, 2049, dtype=np.float64)
        flux_by_radius = np.asarray(
            [_section_flux_for_radius(sections, radius) for radius in radius_grid],
            dtype=np.float64,
        )
        dr = np.diff(radius_grid)
        cumulative = np.concatenate(
            ([0.0], np.cumsum(0.5 * (flux_by_radius[:-1] + flux_by_radius[1:]) * dr))
        )
        if cumulative[-1] <= 0.0:
            raise ValueError("The configured finite-size bubble distribution cannot pass the inlet section.")
        radius_cdf = cumulative / float(cumulative[-1])
        accessible_q2d = float(cumulative[-1] / (radius_max - radius_min))

    raw_q2d = float(sum(_section_flux_for_radius((section,), 0.0) for section in sections))
    concentration_m3 = concentration_ml * MB_PER_ML_TO_MB_PER_M3
    injection_rate = concentration_m3 * thickness * accessible_q2d * UM3_TO_M3
    if not math.isfinite(injection_rate) or injection_rate <= 0.0:
        raise ValueError("The accepted flow and finite-size inlet geometry produce zero MB/s flux.")

    reference_q2d = _reference_inlet_flux_um2_s(flow_field, roots, thickness)
    relative_error = (
        abs(raw_q2d - reference_q2d) / max(abs(reference_q2d), np.finfo(float).eps)
        if reference_q2d > 0.0
        else float("nan")
    )
    return InletFluxModel(
        boundary_geometry=boundary_geometry,
        sections=sections,
        radius_grid_um=radius_grid,
        radius_cdf=radius_cdf,
        radius_min_um=radius_min,
        radius_max_um=radius_max,
        number_concentration_mb_per_ml=concentration_ml,
        number_concentration_mb_per_m3=concentration_m3,
        effective_thickness_um=thickness,
        injection_rate_per_s=injection_rate,
        mean_injection_interval_s=1.0 / injection_rate,
        raw_section_flow_um2_s=raw_q2d,
        size_accessible_section_flow_um2_s=accessible_q2d,
        reference_inlet_flow_um2_s=reference_q2d,
        section_flux_relative_error=relative_error,
    )


def _build_root_section(
    domain: GridDomain,
    flow_field: FlowField,
    boundary_geometry: ContinuousVesselGeometry,
    vessel: Vessel,
    boundary_depth_cells: float,
) -> InjectionSection:
    p0 = np.asarray([vessel.x_p[0], vessel.x_p[2]], dtype=np.float64)
    p1 = np.asarray([vessel.x_d[0], vessel.x_d[2]], dtype=np.float64)
    direction = p1 - p0
    length = float(np.linalg.norm(direction))
    spacing = float(domain.spacing_um)
    if length <= 2.0 * spacing:
        raise ValueError(f"Root vessel {int(vessel.vid)} is too short for an interior injection section.")
    inward = direction / length
    transverse = np.asarray([-inward[1], inward[0]], dtype=np.float64)
    requested_offset = max((float(boundary_depth_cells) + 1.0) * spacing, 2.0 * spacing)
    offset = min(requested_offset, 0.5 * length)
    if offset >= length - spacing:
        raise ValueError(f"Root vessel {int(vessel.vid)} has no section before its first distal junction.")
    center = p0 + inward * offset
    radius = float(vessel.radius)
    sample_count = max(9, int(math.ceil(4.0 * radius / spacing)) + 1)
    s = np.linspace(-radius, radius, sample_count, dtype=np.float64)
    positions_um = center[None, :] + s[:, None] * transverse[None, :]
    positions_grid = np.column_stack(
        (
            (positions_um[:, 0] - float(domain.origin_um[0])) / spacing,
            (positions_um[:, 1] - float(domain.origin_um[2])) / spacing,
        )
    )
    velocity = sample_bilinear(flow_field.velocity_xz_um_s, positions_grid)
    inside_lumen = np.asarray(
        boundary_geometry.contains_xz_um(positions_um),
        dtype=bool,
    )
    clearance = np.asarray(
        boundary_geometry.solid_wall_distance_at_xz_um(positions_um),
        dtype=np.float64,
    )
    connected_clearance = _nearest_inlet_connected_clearance_um(
        boundary_geometry,
        positions_grid,
    )
    speed = np.maximum(velocity @ inward, 0.0)
    admissible_point = (
        inside_lumen
        & np.isfinite(clearance)
        & (clearance > 0.0)
        & np.isfinite(connected_clearance)
        & (connected_clearance > 0.0)
    )
    speed = np.where(admissible_point, speed, 0.0)
    clearance = np.where(admissible_point, clearance, 0.0)
    connected_clearance = np.where(admissible_point, connected_clearance, 0.0)
    if not np.any(speed > 0.0):
        raise ValueError(f"Root injection section {int(vessel.vid)} has no inward flow.")
    return InjectionSection(
        vessel_id=int(vessel.vid),
        center_xz_um=center,
        inward_normal_xz=inward,
        transverse_coordinate_um=s,
        positions_xz_um=positions_um,
        positions_grid=positions_grid,
        inward_speed_um_s=np.asarray(speed, dtype=np.float64),
        wall_distance_um=np.asarray(clearance, dtype=np.float64),
        inlet_connected_clearance_um=np.asarray(connected_clearance, dtype=np.float64),
        inside_lumen=np.asarray(inside_lumen, dtype=bool),
    )


def _valid_interval_flux(section: InjectionSection, radius_um: float) -> np.ndarray:
    radius = float(radius_um)
    valid = (
        section.inside_lumen
        & (section.wall_distance_um >= radius)
        & (section.inlet_connected_clearance_um >= radius)
        & (section.inlet_connected_clearance_um > 0.0)
    )
    ds = np.diff(section.transverse_coordinate_um)
    usable = valid[:-1] & valid[1:]
    return np.where(
        usable,
        0.5 * (section.inward_speed_um_s[:-1] + section.inward_speed_um_s[1:]) * ds,
        0.0,
    )


def _section_flux_for_radius(sections: tuple[InjectionSection, ...], radius_um: float) -> float:
    return float(sum(np.sum(_valid_interval_flux(section, radius_um)) for section in sections))


def _nearest_inlet_connected_clearance_um(
    boundary_geometry: ContinuousVesselGeometry,
    positions_grid: np.ndarray,
) -> np.ndarray:
    """Return continuous-wall clearance at root-section samples."""

    points = np.asarray(positions_grid, dtype=np.float64)
    world = boundary_geometry.grid_to_world_xz(points)
    inside = np.asarray(
        boundary_geometry.contains_xz_um(world), dtype=bool
    )
    distance = np.asarray(
        boundary_geometry.solid_wall_distance_at_xz_um(world),
        dtype=np.float64,
    )
    return np.where(inside, distance, 0.0)


def _linear_density_fraction(u0: float, u1: float, ds: float, target: float) -> float:
    total = 0.5 * (u0 + u1) * ds
    if total <= 0.0:
        return 0.0
    target = min(max(float(target), 0.0), total)
    delta = u1 - u0
    if abs(delta) <= 1.0e-14 * max(abs(u0), abs(u1), 1.0):
        return target / total
    normalized = target / ds
    discriminant = max(u0 * u0 + 2.0 * delta * normalized, 0.0)
    fraction = (-u0 + math.sqrt(discriminant)) / delta
    return min(max(float(fraction), 0.0), 1.0)


def _reference_inlet_flux_um2_s(
    flow_field: FlowField,
    roots: list[Vessel],
    effective_thickness_um: float,
) -> float:
    if flow_field.inlet_actual_by_label_um2_s is not None:
        return float(np.sum(np.asarray(flow_field.inlet_actual_by_label_um2_s, dtype=float)))
    metadata_value = flow_field.solver_metadata.get("actual_inlet_flux_um2_s")
    if metadata_value is not None:
        return float(metadata_value)
    return float(sum(max(float(vessel.flow_rate), 0.0) for vessel in roots) / effective_thickness_um)
