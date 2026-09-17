"""Revised-v16 continuous vessel geometry for grid-carried flow fields.

The vascular lumen is constructed directly from exported vessel centre-line
segments and physical radii before any CFD rasterization.  This object is the
only geometry queried by finite-size particles.  The regular grid remains a
sampled carrier for velocity, pressure, viscosity, and visualization fields.

Shapely is deliberately confined to geometry construction and coarse batch
membership checks.  Runtime wall queries operate on an adaptively subdivided,
grid-independent boundary representation with a cKDTree broad phase.  Adjacent
elements of one continuous boundary chain are treated as one local wall so a
curve tessellation vertex cannot recreate the old raster-corner two-wall
failure.
"""

from dataclasses import dataclass, field
from hashlib import sha256
import math

import numpy as np
from scipy.spatial import cKDTree
from shapely import intersects_xy, prepare
from shapely.geometry import GeometryCollection, LineString, MultiPolygon, Point, Polygon
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

try:
    from shapely import make_valid
except ImportError:  # pragma: no cover - Shapely < 2 fallback
    make_valid = None

from ..core.types import GridDomain
from .continuous_vessel_numba import (
    exact_continuous_wall_states,
    inspect_swept_continuous_wall_paths,
    inspect_swept_continuous_wall_paths_with_end_distance,
    inspect_swept_continuous_wall_paths_with_end_state,
)
from ..particles.particle_boundary_numba import (
    directed_open_section_crossing_mask,
    first_directed_open_section_crossing_arrays,
    first_outlet_crossing_arrays,
)


_GEOMETRY_SCHEMA = "v16_continuous_swept_vessel_boundary"


@dataclass(frozen=True, slots=True)
class ExactContinuousWallState:
    """Nearest point, distance, normal, and stable continuous-wall owner."""

    distance_um: np.ndarray
    inward_normal_xz: np.ndarray
    primary_face_index: np.ndarray
    nearest_point_xz_um: np.ndarray
    unique_nearest_wall: np.ndarray


@dataclass(frozen=True, slots=True)
class SweptContinuousWallContact:
    """Complete straight-centre-chord audit against the continuous wall."""

    minimum_gap_um: float
    first_contact_fraction: float | None
    contact_position_xz_um: np.ndarray | None
    primary_face_index: int
    multiple_wall_contact: bool


@dataclass(frozen=True, slots=True)
class ContinuousVesselGeometry:
    """Immutable X-Z lumen boundary built before CFD rasterization."""

    domain: GridDomain
    lumen_polygon: Polygon = field(repr=False, compare=False)
    full_boundary_start_xz_um: np.ndarray
    full_boundary_end_xz_um: np.ndarray
    full_boundary_inward_normal_xz: np.ndarray
    solid_face_start_xz_um: np.ndarray
    solid_face_end_xz_um: np.ndarray
    solid_face_center_xz_um: np.ndarray
    solid_face_length_um: np.ndarray
    solid_face_outward_normal_xz: np.ndarray
    solid_face_inward_normal_xz: np.ndarray
    solid_face_ring_index: np.ndarray
    solid_face_arclength_start_um: np.ndarray
    solid_face_arclength_end_um: np.ndarray
    boundary_ring_length_um: float
    full_boundary_segment_count: int
    open_section_point_xz_um: np.ndarray
    open_section_outward_normal_xz: np.ndarray
    open_section_tangent_xz: np.ndarray
    open_section_half_width_um: np.ndarray
    open_section_label: np.ndarray
    open_section_kind: np.ndarray
    open_section_vessel_id: np.ndarray
    outlet_section_point_xz_um: np.ndarray = field(repr=False, compare=False)
    outlet_section_outward_normal_xz: np.ndarray = field(repr=False, compare=False)
    outlet_section_tangent_xz: np.ndarray = field(repr=False, compare=False)
    outlet_section_half_width_um: np.ndarray = field(repr=False, compare=False)
    outlet_section_label: np.ndarray = field(repr=False, compare=False)
    _inlet_section_point_xz_um: np.ndarray = field(repr=False, compare=False)
    _inlet_section_outward_normal_xz: np.ndarray = field(repr=False, compare=False)
    _inlet_section_tangent_xz: np.ndarray = field(repr=False, compare=False)
    _inlet_section_half_width_um: np.ndarray = field(repr=False, compare=False)
    _outlet_query_bin_edge_origin_xz_um: np.ndarray = field(repr=False, compare=False)
    _outlet_query_bin_size_um: float = field(repr=False, compare=False)
    _outlet_query_bin_shape: np.ndarray = field(repr=False, compare=False)
    _outlet_query_bin_offsets: np.ndarray = field(repr=False, compare=False)
    _outlet_query_section_indices: np.ndarray = field(repr=False, compare=False)
    curve_quad_segs: int
    maximum_boundary_element_length_um: float
    geometry_hash_sha256: str
    _solid_tree: cKDTree = field(repr=False, compare=False)
    _full_tree: cKDTree = field(repr=False, compare=False)
    _exact_bin_edge_origin_xz_um: np.ndarray = field(repr=False, compare=False)
    _exact_bin_size_um: float = field(repr=False, compare=False)
    _exact_bin_shape: np.ndarray = field(repr=False, compare=False)
    _exact_bin_offsets: np.ndarray = field(repr=False, compare=False)
    _exact_bin_segment_indices: np.ndarray = field(repr=False, compare=False)
    _sweep_bin_edge_origin_xz_um: np.ndarray = field(repr=False, compare=False)
    _sweep_bin_size_um: float = field(repr=False, compare=False)
    _sweep_bin_shape: np.ndarray = field(repr=False, compare=False)
    _sweep_bin_offsets: np.ndarray = field(repr=False, compare=False)
    _sweep_bin_segment_indices: np.ndarray = field(repr=False, compare=False)

    @property
    def shape(self):
        return tuple(int(value) for value in self.domain.shape)

    @property
    def geometry_schema(self):
        return _GEOMETRY_SCHEMA

    def world_xz_to_grid(self, positions_xz_um):
        points = np.asarray(positions_xz_um, dtype=np.float64)
        origin = np.asarray([self.domain.origin_um[0], self.domain.origin_um[2]])
        return (points - origin) / float(self.domain.spacing_um)

    def grid_to_world_xz(self, positions_grid):
        positions = np.asarray(positions_grid, dtype=np.float64)
        origin = np.asarray([self.domain.origin_um[0], self.domain.origin_um[2]])
        return origin + float(self.domain.spacing_um) * positions

    def contains_xz_um(self, positions_xz_um, *, tolerance_um=0.0):
        points, output_shape = _flatten_points(positions_xz_um)
        tolerance = float(tolerance_um)
        if not math.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("tolerance_um must be finite and non-negative.")
        polygon = self.lumen_polygon if tolerance == 0.0 else self.lumen_polygon.buffer(tolerance)
        prepare(polygon)
        result = np.asarray(
            intersects_xy(polygon, points[:, 0], points[:, 1]), dtype=bool
        )
        return _restore(result, output_shape)

    def contains_grid(self, positions_grid, *, tolerance_um=0.0):
        return self.contains_xz_um(
            self.grid_to_world_xz(positions_grid), tolerance_um=tolerance_um
        )

    def signed_lumen_distance_at_xz_um(self, positions_xz_um):
        points, output_shape = _flatten_points(positions_xz_um)
        distance, _, _, _ = _nearest_segment_rows(
            points,
            self.full_boundary_start_xz_um,
            self.full_boundary_end_xz_um,
            self.full_boundary_inward_normal_xz,
            self._full_tree,
        )
        prepare(self.lumen_polygon)
        inside = np.asarray(
            intersects_xy(self.lumen_polygon, points[:, 0], points[:, 1]),
            dtype=bool,
        )
        signed = np.where(inside, distance, -distance)
        return _restore(signed, output_shape)

    def exact_solid_wall_state_xz_um(self, positions_xz_um):
        points, output_shape = _flatten_points(positions_xz_um)
        distance, normal, primary, projected = _nearest_segment_rows(
            points,
            self.solid_face_start_xz_um,
            self.solid_face_end_xz_um,
            self.solid_face_inward_normal_xz,
            self._solid_tree,
            ring_indices=self.solid_face_ring_index,
            ring_size=int(self.full_boundary_segment_count),
        )
        unique = _unique_nearest_wall_rows(
            points,
            distance,
            primary,
            self.solid_face_start_xz_um,
            self.solid_face_end_xz_um,
            self.solid_face_ring_index,
            int(self.full_boundary_segment_count),
            self._solid_tree,
        )
        return ExactContinuousWallState(
            distance_um=_restore(distance, output_shape),
            inward_normal_xz=_restore_vectors(normal, output_shape),
            primary_face_index=_restore(primary, output_shape),
            nearest_point_xz_um=_restore_vectors(projected, output_shape),
            unique_nearest_wall=_restore(unique, output_shape),
        )

    def exact_solid_wall_state_xz_um_accelerated(self, positions_xz_um):
        points, output_shape = _flatten_points(positions_xz_um)
        distance, normal, primary, projected, unique = exact_continuous_wall_states(
            np.ascontiguousarray(points, dtype=np.float64),
            self.solid_face_start_xz_um,
            self.solid_face_end_xz_um,
            self.solid_face_inward_normal_xz,
            self.solid_face_ring_index,
            int(self.full_boundary_segment_count),
            self._exact_bin_edge_origin_xz_um,
            float(self._exact_bin_size_um),
            self._exact_bin_shape,
            self._exact_bin_offsets,
            self._exact_bin_segment_indices,
            use_numba=True,
        )
        return ExactContinuousWallState(
            distance_um=_restore(distance, output_shape),
            inward_normal_xz=_restore_vectors(normal, output_shape),
            primary_face_index=_restore(primary, output_shape),
            nearest_point_xz_um=_restore_vectors(projected, output_shape),
            unique_nearest_wall=_restore(unique, output_shape),
        )

    def solid_wall_distance_at_xz_um(self, positions_xz_um):
        return self.exact_solid_wall_state_xz_um(positions_xz_um).distance_um

    def exact_true_gap_at_xz_um(self, positions_xz_um, radius_um):
        return np.asarray(self.solid_wall_distance_at_xz_um(positions_xz_um)) - np.asarray(
            radius_um, dtype=np.float64
        )

    def true_gap_at_xz_um(self, positions_xz_um, radius_um):
        return self.exact_true_gap_at_xz_um(positions_xz_um, radius_um)

    def nearest_solid_face_projection_xz_um(self, positions_xz_um):
        points, _ = _flatten_points(positions_xz_um)
        distance, _, primary, projected = _nearest_segment_rows(
            points,
            self.solid_face_start_xz_um,
            self.solid_face_end_xz_um,
            self.solid_face_inward_normal_xz,
            self._solid_tree,
            ring_indices=self.solid_face_ring_index,
            ring_size=int(self.full_boundary_segment_count),
        )
        tied = []
        maximum_half = 0.5 * float(np.max(self.solid_face_length_um, initial=0.0))
        for lane, point in enumerate(points):
            radius = float(distance[lane]) + maximum_half + _scale_epsilon(point)
            candidates = np.asarray(self._solid_tree.query_ball_point(point, radius), dtype=np.int64)
            if candidates.size == 0:
                candidates = np.asarray([primary[lane]], dtype=np.int64)
            candidate_distance, candidate_projection = _point_segment_rows(
                point,
                self.solid_face_start_xz_um[candidates],
                self.solid_face_end_xz_um[candidates],
            )
            tolerance = 256.0 * np.finfo(np.float64).eps * max(
                float(np.max(np.abs(point))), float(distance[lane]), 1.0
            )
            local = np.flatnonzero(np.abs(candidate_distance - distance[lane]) <= tolerance)
            indices = np.unique(candidates[local])
            tied.append(indices if indices.size else np.asarray([primary[lane]], dtype=np.int64))
            selected = int(np.flatnonzero(candidates == primary[lane])[0])
            projected[lane] = candidate_projection[selected]
        return projected, primary, tuple(tied)

    def inspect_swept_solid_wall_path_xz_um(
        self,
        start_xz_um,
        end_xz_um,
        radius_um,
        *,
        tolerance_um=0.0,
    ):
        start = np.asarray(start_xz_um, dtype=np.float64).reshape(2)
        end = np.asarray(end_xz_um, dtype=np.float64).reshape(2)
        radius = float(radius_um)
        tolerance = float(tolerance_um)
        if radius < 0.0 or not math.isfinite(radius):
            raise ValueError("radius_um must be finite and non-negative.")
        midpoint = 0.5 * (start + end)
        half_chord = 0.5 * float(np.linalg.norm(end - start))
        maximum_half = 0.5 * float(np.max(self.solid_face_length_um, initial=0.0))
        midpoint_distance = float(
            self.exact_solid_wall_state_xz_um(midpoint).distance_um
        )
        candidates = np.asarray(
            self._solid_tree.query_ball_point(
                midpoint,
                midpoint_distance
                + 2.0 * half_chord
                + maximum_half
                + tolerance,
            ),
            dtype=np.int64,
        )
        if candidates.size == 0:
            nearest = self.exact_solid_wall_state_xz_um(midpoint)
            candidates = np.asarray([int(np.asarray(nearest.primary_face_index))], dtype=np.int64)

        minimum_distance = math.inf
        best_fraction = math.inf
        selected = -1
        simultaneous = []
        fraction_tolerance = 256.0 * np.finfo(np.float64).eps
        for raw_index in np.unique(candidates):
            index = int(raw_index)
            wall_start = self.solid_face_start_xz_um[index]
            wall_end = self.solid_face_end_xz_um[index]
            minimum_distance = min(
                minimum_distance,
                _segment_segment_distance(start, end, wall_start, wall_end),
            )
            fraction = _first_moving_point_segment_capsule_contact(
                start, end, wall_start, wall_end, radius, tolerance
            )
            if fraction is None:
                continue
            if fraction < best_fraction - fraction_tolerance:
                best_fraction = fraction
                selected = index
                simultaneous = [index]
            elif abs(fraction - best_fraction) <= fraction_tolerance:
                simultaneous.append(index)
                if index < selected:
                    selected = index

        multiple = _contains_distinct_wall_indices(
            np.asarray(simultaneous, dtype=np.int64),
            self.solid_face_ring_index,
            int(self.full_boundary_segment_count),
            self.solid_face_inward_normal_xz,
        )
        fraction_value = None if not math.isfinite(best_fraction) else float(
            np.clip(best_fraction, 0.0, 1.0)
        )
        contact = None if fraction_value is None else start + fraction_value * (end - start)
        return SweptContinuousWallContact(
            minimum_gap_um=float(minimum_distance - radius),
            first_contact_fraction=fraction_value,
            contact_position_xz_um=contact,
            primary_face_index=int(selected),
            multiple_wall_contact=bool(multiple),
        )

    def inspect_swept_solid_wall_paths_xz_um_accelerated(
        self,
        starts_xz_um,
        ends_xz_um,
        radii_um,
        *,
        tolerance_um=0.0,
    ):
        starts = np.asarray(starts_xz_um, dtype=np.float64).reshape(-1, 2)
        ends = np.asarray(ends_xz_um, dtype=np.float64).reshape(-1, 2)
        radii = np.asarray(radii_um, dtype=np.float64).reshape(-1)
        if starts.shape != ends.shape or radii.shape != (starts.shape[0],):
            raise ValueError("Swept wall-query batches must have matching lengths.")
        start_distance = np.asarray(
            self.exact_solid_wall_state_xz_um_accelerated(starts).distance_um,
            dtype=np.float64,
        ).reshape(-1)
        end_distance = np.asarray(
            self.exact_solid_wall_state_xz_um_accelerated(ends).distance_um,
            dtype=np.float64,
        ).reshape(-1)
        broad_phase_radius = np.maximum(
            radii, np.minimum(start_distance, end_distance)
        )
        return inspect_swept_continuous_wall_paths(
            np.ascontiguousarray(starts),
            np.ascontiguousarray(ends),
            np.ascontiguousarray(radii),
            np.ascontiguousarray(broad_phase_radius),
            self.solid_face_start_xz_um,
            self.solid_face_end_xz_um,
            self.solid_face_inward_normal_xz,
            self.solid_face_ring_index,
            int(self.full_boundary_segment_count),
            self.solid_face_arclength_start_um,
            self.solid_face_arclength_end_um,
            float(self.boundary_ring_length_um),
            self._sweep_bin_edge_origin_xz_um,
            float(self._sweep_bin_size_um),
            self._sweep_bin_shape,
            self._sweep_bin_offsets,
            self._sweep_bin_segment_indices,
            float(tolerance_um),
            use_numba=True,
        )

    def inspect_swept_solid_wall_paths_xz_um_precomputed(
        self,
        starts_xz_um,
        ends_xz_um,
        radii_um,
        start_distance_um,
        end_distance_um,
        *,
        tolerance_um=0.0,
    ):
        """Run the compiled sweep with exact endpoint distances already available."""

        starts = np.asarray(starts_xz_um, dtype=np.float64).reshape(-1, 2)
        ends = np.asarray(ends_xz_um, dtype=np.float64).reshape(-1, 2)
        radii = np.asarray(radii_um, dtype=np.float64).reshape(-1)
        start_distance = np.asarray(start_distance_um, dtype=np.float64).reshape(-1)
        end_distance = np.asarray(end_distance_um, dtype=np.float64).reshape(-1)
        count = int(starts.shape[0])
        if (
            starts.shape != ends.shape
            or radii.shape != (count,)
            or start_distance.shape != (count,)
            or end_distance.shape != (count,)
        ):
            raise ValueError(
                "Swept wall-query batches and endpoint distances must have matching lengths."
            )
        broad_phase_radius = np.maximum(
            radii, np.minimum(start_distance, end_distance)
        )
        return inspect_swept_continuous_wall_paths(
            np.ascontiguousarray(starts),
            np.ascontiguousarray(ends),
            np.ascontiguousarray(radii),
            np.ascontiguousarray(broad_phase_radius),
            self.solid_face_start_xz_um,
            self.solid_face_end_xz_um,
            self.solid_face_inward_normal_xz,
            self.solid_face_ring_index,
            int(self.full_boundary_segment_count),
            self.solid_face_arclength_start_um,
            self.solid_face_arclength_end_um,
            float(self.boundary_ring_length_um),
            self._sweep_bin_edge_origin_xz_um,
            float(self._sweep_bin_size_um),
            self._sweep_bin_shape,
            self._sweep_bin_offsets,
            self._sweep_bin_segment_indices,
            float(tolerance_um),
            use_numba=True,
        )

    def inspect_swept_solid_wall_paths_with_end_state_xz_um_precomputed(
        self,
        starts_xz_um,
        ends_xz_um,
        radii_um,
        start_distance_um,
        *,
        tolerance_um=0.0,
    ):
        """Evaluate exact endpoints and their shared swept paths in one dispatch."""

        starts = np.asarray(starts_xz_um, dtype=np.float64).reshape(-1, 2)
        ends = np.asarray(ends_xz_um, dtype=np.float64).reshape(-1, 2)
        radii = np.asarray(radii_um, dtype=np.float64).reshape(-1)
        start_distance = np.asarray(start_distance_um, dtype=np.float64).reshape(-1)
        count = int(starts.shape[0])
        if (
            starts.shape != ends.shape
            or radii.shape != (count,)
            or start_distance.shape != (count,)
        ):
            raise ValueError(
                "Swept wall-query batches and start distances must have matching lengths."
            )
        return inspect_swept_continuous_wall_paths_with_end_state(
            np.ascontiguousarray(starts),
            np.ascontiguousarray(ends),
            np.ascontiguousarray(radii),
            np.ascontiguousarray(start_distance),
            self.solid_face_start_xz_um,
            self.solid_face_end_xz_um,
            self.solid_face_inward_normal_xz,
            self.solid_face_ring_index,
            int(self.full_boundary_segment_count),
            self.solid_face_arclength_start_um,
            self.solid_face_arclength_end_um,
            float(self.boundary_ring_length_um),
            self._exact_bin_edge_origin_xz_um,
            float(self._exact_bin_size_um),
            self._exact_bin_shape,
            self._exact_bin_offsets,
            self._exact_bin_segment_indices,
            self._sweep_bin_edge_origin_xz_um,
            float(self._sweep_bin_size_um),
            self._sweep_bin_shape,
            self._sweep_bin_offsets,
            self._sweep_bin_segment_indices,
            float(tolerance_um),
            use_numba=True,
        )

    def inspect_swept_solid_wall_paths_with_end_distance_xz_um_precomputed(
        self,
        starts_xz_um,
        ends_xz_um,
        radii_um,
        start_distance_um,
        *,
        tolerance_um=0.0,
    ):
        """Evaluate only the needed endpoint distance and shared swept paths."""

        starts = np.asarray(starts_xz_um, dtype=np.float64).reshape(-1, 2)
        ends = np.asarray(ends_xz_um, dtype=np.float64).reshape(-1, 2)
        radii = np.asarray(radii_um, dtype=np.float64).reshape(-1)
        start_distance = np.asarray(start_distance_um, dtype=np.float64).reshape(-1)
        count = int(starts.shape[0])
        if (
            starts.shape != ends.shape
            or radii.shape != (count,)
            or start_distance.shape != (count,)
        ):
            raise ValueError(
                "Swept wall-query batches and start distances must have matching lengths."
            )
        return inspect_swept_continuous_wall_paths_with_end_distance(
            np.ascontiguousarray(starts),
            np.ascontiguousarray(ends),
            np.ascontiguousarray(radii),
            np.ascontiguousarray(start_distance),
            self.solid_face_start_xz_um,
            self.solid_face_end_xz_um,
            self.solid_face_inward_normal_xz,
            self.solid_face_ring_index,
            int(self.full_boundary_segment_count),
            self.solid_face_arclength_start_um,
            self.solid_face_arclength_end_um,
            float(self.boundary_ring_length_um),
            self._exact_bin_edge_origin_xz_um,
            float(self._exact_bin_size_um),
            self._exact_bin_shape,
            self._exact_bin_offsets,
            self._exact_bin_segment_indices,
            self._sweep_bin_edge_origin_xz_um,
            float(self._sweep_bin_size_um),
            self._sweep_bin_shape,
            self._sweep_bin_offsets,
            self._sweep_bin_segment_indices,
            float(tolerance_um),
            use_numba=True,
        )

    def multiple_wall_contact_mask_xz_um(
        self,
        positions_xz_um,
        radius_um,
        *,
        tolerance_um=1.0e-9,
    ):
        points, output_shape = _flatten_points(positions_xz_um)
        radii = np.broadcast_to(np.asarray(radius_um, dtype=np.float64), output_shape or (1,)).reshape(-1)
        result = np.zeros(points.shape[0], dtype=bool)
        maximum_half = 0.5 * float(np.max(self.solid_face_length_um, initial=0.0))
        for lane, point in enumerate(points):
            candidates = np.asarray(
                self._solid_tree.query_ball_point(
                    point, float(radii[lane]) + maximum_half + float(tolerance_um)
                ),
                dtype=np.int64,
            )
            if candidates.size < 2:
                continue
            distances, _ = _point_segment_rows(
                point,
                self.solid_face_start_xz_um[candidates],
                self.solid_face_end_xz_um[candidates],
            )
            contacting = candidates[distances <= radii[lane] + float(tolerance_um)]
            result[lane] = _contains_distinct_wall_indices(
                contacting,
                self.solid_face_ring_index,
                int(self.full_boundary_segment_count),
                self.solid_face_inward_normal_xz,
            )
        return _restore(result, output_shape)

    def outlet_signed_distances_xz_um(self, positions_xz_um):
        points, output_shape = _flatten_points(positions_xz_um)
        section_points = self.outlet_section_point_xz_um
        normals = self.outlet_section_outward_normal_xz
        values = np.einsum(
            "pki,ki->pk", section_points[None, :, :] - points[:, None, :], normals
        )
        if output_shape:
            return values.reshape(*output_shape, section_points.shape[0])
        return values[0]

    def first_outlet_crossing_xz_um(self, start_xz_um, end_xz_um):
        starts = np.asarray(start_xz_um, dtype=np.float64).reshape(1, 2)
        ends = np.asarray(end_xz_um, dtype=np.float64).reshape(1, 2)
        fraction, position, section_index, signed_start, signed_end = (
            self._first_outlet_crossing_rows(starts, ends)
        )
        if not np.isfinite(fraction[0]):
            return None
        index = int(section_index[0])
        return _Crossing(
            fraction=float(fraction[0]),
            outlet_index=index,
            label=int(self.outlet_section_label[index]),
            position_xz_um=position[0],
            signed_distance_start_um=float(signed_start[0]),
            signed_distance_end_um=float(signed_end[0]),
        )

    def first_outlet_crossing_grid(self, start_grid, end_grid):
        return self.first_outlet_crossing_xz_um(
            self.grid_to_world_xz(start_grid), self.grid_to_world_xz(end_grid)
        )

    def first_outlet_crossings_grid_accelerated(self, starts_grid, ends_grid):
        return self.first_outlet_crossings_grid(starts_grid, ends_grid)

    def first_outlet_crossings_grid(self, starts_grid, ends_grid):
        starts = self.grid_to_world_xz(starts_grid)
        ends = self.grid_to_world_xz(ends_grid)
        fractions, positions, indices, signed_start, signed_end = (
            self._first_outlet_crossing_rows(starts, ends)
        )
        result = []
        labels = self.outlet_section_label
        for lane, raw_index in enumerate(indices):
            index = int(raw_index)
            if index < 0:
                result.append(None)
                continue
            result.append(
                _Crossing(
                    fraction=float(fractions[lane]),
                    outlet_index=index,
                    label=int(labels[index]),
                    position_xz_um=np.asarray(positions[lane], dtype=np.float64),
                    signed_distance_start_um=float(signed_start[lane]),
                    signed_distance_end_um=float(signed_end[lane]),
                )
            )
        return tuple(result)

    def first_outlet_crossing_arrays_grid_accelerated(self, starts_grid, ends_grid):
        starts = self.grid_to_world_xz(starts_grid)
        ends = self.grid_to_world_xz(ends_grid)
        coordinate_scale = max(
            float(np.max(np.abs(starts), initial=0.0)),
            float(np.max(np.abs(ends), initial=0.0)),
            float(
                np.max(
                    np.abs(self.outlet_section_point_xz_um), initial=0.0
                )
            ),
            1.0,
        )
        plane_tolerance = (
            4096.0 * np.finfo(np.float64).eps * coordinate_scale
        )
        fractions, section_indices, positions, signed_start, signed_end = (
            first_outlet_crossing_arrays(
                np.ascontiguousarray(starts),
                np.ascontiguousarray(ends),
                self.outlet_section_point_xz_um,
                self.outlet_section_outward_normal_xz,
                self.outlet_section_tangent_xz,
                self.outlet_section_half_width_um,
                self._outlet_query_bin_edge_origin_xz_um,
                float(self._outlet_query_bin_size_um),
                self._outlet_query_bin_shape,
                self._outlet_query_bin_offsets,
                self._outlet_query_section_indices,
                plane_tolerance,
                use_numba=True,
            )
        )
        spacing = float(self.domain.spacing_um)
        return (
            fractions,
            section_indices,
            self.world_xz_to_grid(positions),
            signed_start / spacing,
            signed_end / spacing,
        )

    def directed_inlet_crossing_mask_grid_accelerated(self, starts_grid, ends_grid):
        """Identify reverse exits through an anatomical inlet cap.

        A path that starts in the accepted lumen, crosses neither a solid wall
        nor any open cap, must end inside the same continuous polygon.  The
        normal production path already queries outlets; this compact inlet
        guard completes that proof without an expensive GEOS membership call.
        """

        starts = np.ascontiguousarray(
            self.grid_to_world_xz(starts_grid), dtype=np.float64
        )
        ends = np.ascontiguousarray(
            self.grid_to_world_xz(ends_grid), dtype=np.float64
        )
        coordinate_scale = max(
            float(np.max(np.abs(starts), initial=0.0)),
            float(np.max(np.abs(ends), initial=0.0)),
            float(
                np.max(
                    np.abs(self._inlet_section_point_xz_um), initial=0.0
                )
            ),
            1.0,
        )
        plane_tolerance = (
            4096.0 * np.finfo(np.float64).eps * coordinate_scale
        )
        return directed_open_section_crossing_mask(
            starts,
            ends,
            self._inlet_section_point_xz_um,
            self._inlet_section_outward_normal_xz,
            self._inlet_section_tangent_xz,
            self._inlet_section_half_width_um,
            plane_tolerance,
            use_numba=True,
        )

    def first_inlet_crossing_arrays_grid_accelerated(self, starts_grid, ends_grid):
        """Return chronological directed reverse-inlet crossing data."""

        starts = np.ascontiguousarray(
            self.grid_to_world_xz(starts_grid), dtype=np.float64
        )
        ends = np.ascontiguousarray(
            self.grid_to_world_xz(ends_grid), dtype=np.float64
        )
        coordinate_scale = max(
            float(np.max(np.abs(starts), initial=0.0)),
            float(np.max(np.abs(ends), initial=0.0)),
            float(
                np.max(
                    np.abs(self._inlet_section_point_xz_um), initial=0.0
                )
            ),
            1.0,
        )
        plane_tolerance = (
            4096.0 * np.finfo(np.float64).eps * coordinate_scale
        )
        fractions, indices, positions = (
            first_directed_open_section_crossing_arrays(
                starts,
                ends,
                self._inlet_section_point_xz_um,
                self._inlet_section_outward_normal_xz,
                self._inlet_section_tangent_xz,
                self._inlet_section_half_width_um,
                plane_tolerance,
                use_numba=True,
            )
        )
        return fractions, indices, self.world_xz_to_grid(positions)

    def _first_outlet_crossing_rows(self, starts_xz_um, ends_xz_um):
        starts = np.asarray(starts_xz_um, dtype=np.float64).reshape(-1, 2)
        ends = np.asarray(ends_xz_um, dtype=np.float64).reshape(-1, 2)
        points = self.outlet_section_point_xz_um
        normals = self.outlet_section_outward_normal_xz
        tangents = self.outlet_section_tangent_xz
        widths = self.outlet_section_half_width_um
        count = starts.shape[0]
        fractions = np.full(count, np.nan, dtype=np.float64)
        positions = np.full((count, 2), np.nan, dtype=np.float64)
        result_indices = np.full(count, -1, dtype=np.int32)
        result_signed_start = np.full(count, np.nan, dtype=np.float64)
        result_signed_end = np.full(count, np.nan, dtype=np.float64)
        for lane in range(count):
            delta = ends[lane] - starts[lane]
            coordinate_scale = max(
                float(np.max(np.abs(starts[lane]))),
                float(np.max(np.abs(ends[lane]))),
                1.0,
            )
            plane_tolerance = (
                4096.0 * np.finfo(np.float64).eps * coordinate_scale
            )
            best = math.inf
            best_index = -1
            best_signed_start = math.nan
            best_signed_end = math.nan
            for section in range(points.shape[0]):
                start_signed = float((points[section] - starts[lane]) @ normals[section])
                end_signed = float((points[section] - ends[lane]) @ normals[section])
                effective_start = (
                    0.0 if abs(start_signed) <= plane_tolerance else start_signed
                )
                effective_end = (
                    0.0 if abs(end_signed) <= plane_tolerance else end_signed
                )
                denominator = effective_start - effective_end
                if (
                    start_signed < -plane_tolerance
                    or end_signed > plane_tolerance
                    or end_signed >= start_signed - plane_tolerance
                    or denominator <= 0.0
                ):
                    continue
                fraction = effective_start / denominator
                if fraction < -plane_tolerance or fraction > 1.0 + plane_tolerance:
                    continue
                fraction = float(np.clip(fraction, 0.0, 1.0))
                crossing = starts[lane] + fraction * delta
                lateral = abs(float((crossing - points[section]) @ tangents[section]))
                if lateral > widths[section] + max(
                    plane_tolerance, _scale_epsilon(crossing)
                ):
                    continue
                if fraction < best or (fraction == best and section < best_index):
                    best = fraction
                    best_index = section
                    best_signed_start = start_signed
                    best_signed_end = end_signed
            if math.isfinite(best):
                fractions[lane] = best
                positions[lane] = starts[lane] + best * delta
                result_indices[lane] = best_index
                result_signed_start[lane] = best_signed_start
                result_signed_end[lane] = best_signed_end
        return (
            fractions,
            positions,
            result_indices,
            result_signed_start,
            result_signed_end,
        )

    def accessible_mask(self, radius_um):
        grid_x, grid_z = np.meshgrid(
            self.domain.x_coordinates_um,
            self.domain.z_coordinates_um,
            indexing="ij",
        )
        points = np.column_stack((grid_x.ravel(), grid_z.ravel()))
        result = self.is_accessible_xz_um(points, radius_um)
        return np.asarray(result, dtype=bool).reshape(self.shape)

    def is_accessible_xz_um(self, positions_xz_um, radius_um):
        inside = np.asarray(self.contains_xz_um(positions_xz_um), dtype=bool)
        gap = np.asarray(self.exact_true_gap_at_xz_um(positions_xz_um, radius_um))
        return inside & (gap >= 0.0)

    def is_accessible_grid(self, positions_grid, radius_um):
        return self.is_accessible_xz_um(self.grid_to_world_xz(positions_grid), radius_um)


@dataclass(frozen=True, slots=True)
class _Crossing:
    fraction: float
    outlet_index: int
    label: int
    position_xz_um: np.ndarray
    signed_distance_start_um: float
    signed_distance_end_um: float


def build_continuous_vessel_geometry(vessels, domain, *, curve_quad_segs=64, maximum_boundary_element_length_um=1.0):
    """
    Build the grid-independent v16 lumen and its runtime boundary arrays.
    """

    vessel_list     = list(vessels)
    raw_quad_segs   = float(curve_quad_segs)
    if (
        isinstance(curve_quad_segs, bool)
        or not math.isfinite(raw_quad_segs)
        or not raw_quad_segs.is_integer()
        or raw_quad_segs < 8
    ):
        raise ValueError("curve_quad_segs must be an integer of at least 8.")
    quad_segs = int(raw_quad_segs)
    maximum_length = float(maximum_boundary_element_length_um)
    if (
        isinstance(maximum_boundary_element_length_um, bool)
        or not math.isfinite(maximum_length)
        or maximum_length <= 0.0
    ):
        raise ValueError("maximum_boundary_element_length_um must be positive.")
    if not vessel_list:
        raise ValueError("At least one vessel is required.")
    for vessel in vessel_list:
        endpoints = np.asarray(
            [
                vessel.x_p[0],
                vessel.x_p[2],
                vessel.x_d[0],
                vessel.x_d[2],
            ],
            dtype=np.float64,
        )
        radius = float(vessel.radius)
        if not np.all(np.isfinite(endpoints)) or not math.isfinite(radius):
            raise ValueError(
                f"Vessel {vessel.vid} coordinates and radius must be finite."
            )
        if radius <= 0.0:
            raise ValueError(f"Vessel {vessel.vid} radius must be positive.")
        if np.linalg.norm(endpoints[2:] - endpoints[:2]) <= np.finfo(float).eps:
            raise ValueError(
                f"Vessel {vessel.vid} has zero projected X-Z length."
            )
    polygon = _continuous_lumen_polygon(vessel_list, quad_segs=quad_segs)
    polygon = orient(polygon, sign=1.0)
    prepare(polygon)

    sections = _continuous_open_sections(vessel_list)
    outlet_mask = np.asarray(sections["kind"] > 0, dtype=bool)
    outlet_points = np.ascontiguousarray(sections["point_xz_um"][outlet_mask])
    outlet_normals = np.ascontiguousarray(sections["normal_xz"][outlet_mask])
    outlet_tangents = np.ascontiguousarray(sections["tangent_xz"][outlet_mask])
    outlet_widths = np.ascontiguousarray(sections["half_width_um"][outlet_mask])
    outlet_labels = np.ascontiguousarray(sections["label"][outlet_mask])
    inlet_mask = np.asarray(sections["kind"] < 0, dtype=bool)
    inlet_points = np.ascontiguousarray(sections["point_xz_um"][inlet_mask])
    inlet_normals = np.ascontiguousarray(sections["normal_xz"][inlet_mask])
    inlet_tangents = np.ascontiguousarray(sections["tangent_xz"][inlet_mask])
    inlet_widths = np.ascontiguousarray(sections["half_width_um"][inlet_mask])
    outlet_count = int(outlet_points.shape[0])
    outlet_query_origin = np.asarray(
        [domain.origin_um[0], domain.origin_um[2]], dtype=np.float64
    ) - float(domain.spacing_um)
    outlet_query_extent = np.asarray(
        [
            domain.shape[0] * domain.spacing_um,
            domain.shape[1] * domain.spacing_um,
        ],
        dtype=np.float64,
    )
    outlet_query_bin_size = float(
        max(outlet_query_extent) + 2.0 * domain.spacing_um
    )
    outlet_query_shape = np.asarray([1, 1], dtype=np.int32)
    outlet_query_offsets = np.asarray([0, outlet_count], dtype=np.int64)
    outlet_query_indices = np.arange(outlet_count, dtype=np.int32)
    coordinates = np.asarray(polygon.exterior.coords, dtype=np.float64)
    full_start, full_end = _subdivide_ring(coordinates, maximum_length)
    full_tangent = full_end - full_start
    full_length = np.linalg.norm(full_tangent, axis=1)
    if np.any(full_length <= 0.0):
        raise ValueError("Continuous boundary contains a zero-length element.")
    full_tangent /= full_length[:, None]
    full_inward = np.column_stack((-full_tangent[:, 1], full_tangent[:, 0]))

    open_mask = _open_cap_segment_mask(full_start, full_end, sections)
    solid = ~open_mask
    if not np.any(solid):
        raise ValueError("Continuous vessel geometry contains no closed solid wall.")
    solid_start = full_start[solid]
    solid_end = full_end[solid]
    solid_length = full_length[solid]
    solid_center = 0.5 * (solid_start + solid_end)
    solid_inward = full_inward[solid]
    solid_ring_index = np.flatnonzero(solid).astype(np.int64)
    full_arclength_start = np.concatenate(
        (np.asarray([0.0]), np.cumsum(full_length[:-1], dtype=np.float64))
    )
    solid_arclength_start = full_arclength_start[solid]
    solid_arclength_end = solid_arclength_start + solid_length

    solid_tree = cKDTree(solid_center)
    full_tree = cKDTree(0.5 * (full_start + full_end))
    runtime_bins = _build_continuous_runtime_bins(
        domain,
        solid_start,
        solid_end,
        solid_center,
        solid_tree,
    )
    digest = _geometry_hash(
        full_start,
        full_end,
        sections,
        quad_segs,
        maximum_length,
    )

    return ContinuousVesselGeometry(
        domain=domain,
        lumen_polygon=polygon,
        full_boundary_start_xz_um=np.ascontiguousarray(full_start),
        full_boundary_end_xz_um=np.ascontiguousarray(full_end),
        full_boundary_inward_normal_xz=np.ascontiguousarray(full_inward),
        solid_face_start_xz_um=np.ascontiguousarray(solid_start),
        solid_face_end_xz_um=np.ascontiguousarray(solid_end),
        solid_face_center_xz_um=np.ascontiguousarray(solid_center),
        solid_face_length_um=np.ascontiguousarray(solid_length),
        solid_face_outward_normal_xz=np.ascontiguousarray(-solid_inward),
        solid_face_inward_normal_xz=np.ascontiguousarray(solid_inward),
        solid_face_ring_index=np.ascontiguousarray(solid_ring_index),
        solid_face_arclength_start_um=np.ascontiguousarray(solid_arclength_start),
        solid_face_arclength_end_um=np.ascontiguousarray(solid_arclength_end),
        boundary_ring_length_um=float(np.sum(full_length)),
        full_boundary_segment_count=int(full_start.shape[0]),
        open_section_point_xz_um=sections["point_xz_um"],
        open_section_outward_normal_xz=sections["normal_xz"],
        open_section_tangent_xz=sections["tangent_xz"],
        open_section_half_width_um=sections["half_width_um"],
        open_section_label=sections["label"],
        open_section_kind=sections["kind"],
        open_section_vessel_id=sections["vessel_id"],
        outlet_section_point_xz_um=outlet_points,
        outlet_section_outward_normal_xz=outlet_normals,
        outlet_section_tangent_xz=outlet_tangents,
        outlet_section_half_width_um=outlet_widths,
        outlet_section_label=outlet_labels,
        _inlet_section_point_xz_um=inlet_points,
        _inlet_section_outward_normal_xz=inlet_normals,
        _inlet_section_tangent_xz=inlet_tangents,
        _inlet_section_half_width_um=inlet_widths,
        _outlet_query_bin_edge_origin_xz_um=np.ascontiguousarray(
            outlet_query_origin
        ),
        _outlet_query_bin_size_um=outlet_query_bin_size,
        _outlet_query_bin_shape=np.ascontiguousarray(outlet_query_shape),
        _outlet_query_bin_offsets=np.ascontiguousarray(outlet_query_offsets),
        _outlet_query_section_indices=np.ascontiguousarray(outlet_query_indices),
        curve_quad_segs=quad_segs,
        maximum_boundary_element_length_um=maximum_length,
        geometry_hash_sha256=digest,
        _solid_tree=solid_tree,
        _full_tree=full_tree,
        _exact_bin_edge_origin_xz_um=runtime_bins["edge_origin_xz_um"],
        _exact_bin_size_um=float(runtime_bins["bin_size_um"]),
        _exact_bin_shape=runtime_bins["shape"],
        _exact_bin_offsets=runtime_bins["exact_offsets"],
        _exact_bin_segment_indices=runtime_bins["exact_segment_indices"],
        _sweep_bin_edge_origin_xz_um=runtime_bins["edge_origin_xz_um"],
        _sweep_bin_size_um=float(runtime_bins["bin_size_um"]),
        _sweep_bin_shape=runtime_bins["shape"],
        _sweep_bin_offsets=runtime_bins["sweep_offsets"],
        _sweep_bin_segment_indices=runtime_bins["sweep_segment_indices"],
    )


def _build_continuous_runtime_bins(
    domain,
    starts,
    ends,
    centers,
    tree,
):
    """Construct conservative exact-nearest and swept-path CSR broad phases."""

    bin_width_cells = 4
    bin_size = float(bin_width_cells * domain.spacing_um)
    shape = np.asarray(
        [
            (int(domain.shape[0]) + bin_width_cells - 1) // bin_width_cells,
            (int(domain.shape[1]) + bin_width_cells - 1) // bin_width_cells,
        ],
        dtype=np.int32,
    )
    edge_origin = np.asarray(
        [domain.origin_um[0], domain.origin_um[2]], dtype=np.float64
    ) - 0.5 * float(domain.spacing_um)
    bin_x = np.repeat(
        np.arange(int(shape[0]), dtype=np.float64), int(shape[1])
    )
    bin_z = np.tile(
        np.arange(int(shape[1]), dtype=np.float64), int(shape[0])
    )
    bin_centers = np.empty((bin_x.size, 2), dtype=np.float64)
    bin_centers[:, 0] = edge_origin[0] + bin_size * (bin_x + 0.5)
    bin_centers[:, 1] = edge_origin[1] + bin_size * (bin_z + 0.5)

    nearest_distance = _nearest_segment_rows(
        bin_centers,
        starts,
        ends,
        np.zeros_like(starts),
        tree,
    )[0]
    half_diagonal = math.sqrt(0.5) * bin_size
    maximum_half_length = 0.5 * float(
        np.max(np.linalg.norm(ends - starts, axis=1), initial=0.0)
    )
    guard = 1024.0 * np.finfo(np.float64).eps * np.maximum(
        np.max(np.abs(bin_centers), axis=1), 1.0
    )
    radii = nearest_distance + 2.0 * half_diagonal + maximum_half_length + guard
    exact_lists = tree.query_ball_point(
        bin_centers, radii, return_sorted=True
    )
    exact_lists = [
        values
        if len(values)
        else [int(tree.query(bin_centers[index], k=1)[1])]
        for index, values in enumerate(exact_lists)
    ]

    sweep_lists = [
        [] for _ in range(int(shape[0]) * int(shape[1]))
    ]
    for segment in range(starts.shape[0]):
        scale = max(
            float(np.max(np.abs(starts[segment]))),
            float(np.max(np.abs(ends[segment]))),
            1.0,
        )
        padding = 512.0 * np.finfo(np.float64).eps * scale
        minimum = np.minimum(starts[segment], ends[segment]) - padding
        maximum = np.maximum(starts[segment], ends[segment]) + padding
        first = np.floor((minimum - edge_origin) / bin_size).astype(np.int64)
        last = np.floor((maximum - edge_origin) / bin_size).astype(np.int64)
        first = np.maximum(first, 0)
        last = np.minimum(last, shape.astype(np.int64) - 1)
        if np.any(first > last):
            continue
        for bin_x in range(int(first[0]), int(last[0]) + 1):
            for bin_z in range(int(first[1]), int(last[1]) + 1):
                sweep_lists[bin_x * int(shape[1]) + bin_z].append(segment)

    exact_offsets, exact_indices = _lists_to_csr(exact_lists)
    sweep_offsets, sweep_indices = _lists_to_csr(sweep_lists)
    return {
        "edge_origin_xz_um": np.ascontiguousarray(edge_origin),
        "bin_size_um": bin_size,
        "shape": np.ascontiguousarray(shape),
        "exact_offsets": exact_offsets,
        "exact_segment_indices": exact_indices,
        "sweep_offsets": sweep_offsets,
        "sweep_segment_indices": sweep_indices,
    }


def _lists_to_csr(values):
    offsets = np.zeros(len(values) + 1, dtype=np.int64)
    for index, row in enumerate(values):
        offsets[index + 1] = offsets[index] + len(row)
    indices = np.empty(int(offsets[-1]), dtype=np.int32)
    for index, row in enumerate(values):
        indices[offsets[index] : offsets[index + 1]] = row
    return np.ascontiguousarray(offsets), np.ascontiguousarray(indices)


def _continuous_lumen_polygon(vessels, *, quad_segs):
    polygons = []
    junctions = {}
    open_keys = _open_endpoint_keys(vessels)
    for vessel in vessels:
        p0 = np.asarray([vessel.x_p[0], vessel.x_p[2]], dtype=np.float64)
        p1 = np.asarray([vessel.x_d[0], vessel.x_d[2]], dtype=np.float64)
        delta = p1 - p0
        length = float(np.linalg.norm(delta))
        radius = float(vessel.radius)
        polygons.append(
            LineString([tuple(p0), tuple(p1)]).buffer(
                radius,
                cap_style=2,
                join_style=1,
                quad_segs=quad_segs,
            )
        )
        unit = delta / length
        junctions.setdefault(_endpoint_key(p0), []).append((unit, radius, length))
        junctions.setdefault(_endpoint_key(p1), []).append((-unit, radius, length))

    for key, connections in junctions.items():
        if key in open_keys or len(connections) < 2:
            continue
        # Join unequal-radius branches with external common tangents between a
        # central physical disk and one disk inside each incident corridor.
        # Unlike a bare disk/rectangle union, this produces a smooth taper and
        # does not replace a raster corner with a polygon corner.  Every length
        # below depends only on vessel radii/lengths, never CFD grid spacing.
        junction_point = np.asarray(key, dtype=np.float64)
        junction_radius = max(value for _, value, _ in connections)
        central_disk = Point(key).buffer(junction_radius, quad_segs=quad_segs)
        polygons.append(central_disk)
        for direction, radius, vessel_length in connections:
            transition_length = min(
                0.45 * float(vessel_length),
                max(
                    2.0 * junction_radius,
                    junction_radius + float(radius),
                    4.0 * abs(junction_radius - float(radius)),
                ),
            )
            if transition_length <= np.finfo(np.float64).eps:
                continue
            transition_center = junction_point + transition_length * direction
            branch_disk = Point(tuple(transition_center)).buffer(
                float(radius), quad_segs=quad_segs
            )
            polygons.append(unary_union((central_disk, branch_disk)).convex_hull)

    merged = _repair_polygon(unary_union(polygons))
    merged = _remove_polygon_holes(merged)
    merged = _repair_polygon(merged)
    merged = _smooth_continuous_junctions(
        merged, vessels, quad_segs=quad_segs
    )
    if isinstance(merged, MultiPolygon):
        raise ValueError("Continuous vessel geometry is disconnected.")
    if not isinstance(merged, Polygon):
        raise ValueError(f"Unsupported continuous lumen type: {type(merged).__name__}.")
    return merged


def _smooth_continuous_junctions(
    geometry,
    vessels,
    *,
    quad_segs,
):
    """Round residual union intersections on a physical, grid-free scale.

    Common-tangent transition hulls remove the large radius-mismatch corner.
    Their union can still contain tiny convex/concave intersections where
    several branches overlap.  A symmetric close/open fillet removes both
    kinds without choosing a CFD-dependent length.  Original flat collars are
    then restored at anatomical openings so outlet planes stay exact.
    """

    positive_radii = [float(vessel.radius) for vessel in vessels if vessel.radius > 0.0]
    if not positive_radii or geometry.is_empty:
        return geometry
    endpoint_counts = {}
    for vessel in vessels:
        for point in (vessel.x_p, vessel.x_d):
            key = _endpoint_key((float(point[0]), float(point[2])))
            endpoint_counts[key] = endpoint_counts.get(key, 0) + 1
    open_keys = _open_endpoint_keys(vessels)
    if not any(count >= 2 and key not in open_keys for key, count in endpoint_counts.items()):
        return geometry
    fillet_radius = 0.05 * min(positive_radii)
    if fillet_radius <= np.finfo(np.float64).eps:
        return geometry
    smoothed = geometry.buffer(
        fillet_radius, join_style=1, quad_segs=quad_segs
    ).buffer(
        -2.0 * fillet_radius, join_style=1, quad_segs=quad_segs
    ).buffer(
        fillet_radius, join_style=1, quad_segs=quad_segs
    )

    collars = []
    for vessel in vessels:
        p0, p1, unit = _vessel_xz(vessel)
        length = float(np.linalg.norm(p1 - p0))
        collar_length = min(0.45 * length, 2.0 * float(vessel.radius))
        if collar_length <= 0.0:
            continue
        if int(vessel.parent_id) < 0:
            collars.append(
                LineString((p0, p0 + collar_length * unit)).buffer(
                    float(vessel.radius),
                    cap_style=2,
                    join_style=1,
                    quad_segs=quad_segs,
                )
            )
        if len(vessel.children) == 0:
            collars.append(
                LineString((p1, p1 - collar_length * unit)).buffer(
                    float(vessel.radius),
                    cap_style=2,
                    join_style=1,
                    quad_segs=quad_segs,
                )
            )
    if collars:
        smoothed = unary_union((smoothed, *collars))
    smoothed = _repair_polygon(_remove_polygon_holes(smoothed))
    return smoothed


def _continuous_open_sections(vessels):
    rows = []
    roots = [vessel for vessel in vessels if int(vessel.parent_id) < 0]
    terminals = [vessel for vessel in vessels if len(vessel.children) == 0]
    for label, vessel in enumerate(roots, start=1):
        p0, p1, unit = _vessel_xz(vessel)
        normal = -unit
        tangent = np.asarray([-normal[1], normal[0]], dtype=np.float64)
        rows.append((p0, normal, tangent, float(vessel.radius), label, -1, int(vessel.vid)))
    for label, vessel in enumerate(terminals, start=1):
        p0, p1, unit = _vessel_xz(vessel)
        normal = unit
        tangent = np.asarray([-normal[1], normal[0]], dtype=np.float64)
        rows.append((p1, normal, tangent, float(vessel.radius), label, 1, int(vessel.vid)))
    if not rows:
        raise ValueError("Continuous geometry requires at least one open section.")
    return {
        "point_xz_um": np.ascontiguousarray(np.vstack([row[0] for row in rows])),
        "normal_xz": np.ascontiguousarray(np.vstack([row[1] for row in rows])),
        "tangent_xz": np.ascontiguousarray(np.vstack([row[2] for row in rows])),
        "half_width_um": np.asarray([row[3] for row in rows], dtype=np.float64),
        "label": np.asarray([row[4] for row in rows], dtype=np.int32),
        "kind": np.asarray([row[5] for row in rows], dtype=np.int8),
        "vessel_id": np.asarray([row[6] for row in rows], dtype=np.int32),
    }


def _nearest_segment_rows(
    points,
    starts,
    ends,
    inward,
    tree,
    *,
    ring_indices=None,
    ring_size=0,
):
    points = np.asarray(points, dtype=np.float64).reshape(-1, 2)
    segment_count = starts.shape[0]
    k = min(16, segment_count)
    centre_distance, candidates = tree.query(points, k=k)
    if k == 1:
        centre_distance = np.asarray(centre_distance)[:, None]
        candidates = np.asarray(candidates)[:, None]
    candidate_start = starts[candidates]
    candidate_end = ends[candidates]
    delta = candidate_end - candidate_start
    length_squared = np.einsum("nki,nki->nk", delta, delta)
    relative = points[:, None, :] - candidate_start
    numerator = np.einsum("nki,nki->nk", relative, delta)
    # Treat a degenerate segment as its start point.  The continuous builder
    # normally removes such elements, but this guard keeps reference queries
    # warning-free for imported/legacy geometry and matches the compiled
    # point-segment kernel exactly.
    fraction = np.zeros_like(length_squared)
    np.divide(
        numerator,
        length_squared,
        out=fraction,
        where=length_squared > 0.0,
    )
    np.clip(fraction, 0.0, 1.0, out=fraction)
    projected_all = candidate_start + fraction[..., None] * delta
    distance_all = np.linalg.norm(points[:, None, :] - projected_all, axis=2)
    local = np.argmin(distance_all, axis=1)
    rows = np.arange(points.shape[0])
    distance = distance_all[rows, local]
    primary = candidates[rows, local].astype(np.int64)
    projected = projected_all[rows, local]
    normal = inward[primary].copy()

    # All runtime elements are at most one micrometre long.  If a candidate
    # centre just outside the initial k-neighbour set could still beat the best
    # exact segment distance, expand only those uncertain rows.
    maximum_half = 0.5 * float(np.max(np.linalg.norm(ends - starts, axis=1), initial=0.0))
    uncertain = distance > np.asarray(centre_distance)[:, -1] - maximum_half
    if np.any(uncertain) and k < segment_count:
        expanded_k = min(max(64, 4 * k), segment_count)
        for lane in np.flatnonzero(uncertain):
            _, expanded = tree.query(points[lane], k=expanded_k)
            expanded = np.atleast_1d(expanded).astype(np.int64)
            expanded_distance, expanded_projection = _point_segment_rows(
                points[lane], starts[expanded], ends[expanded]
            )
            best_local = int(np.argmin(expanded_distance))
            distance[lane] = expanded_distance[best_local]
            primary[lane] = expanded[best_local]
            projected[lane] = expanded_projection[best_local]
            normal[lane] = inward[primary[lane]]

    if ring_indices is not None:
        # Average only normals of adjacent elements meeting at the same smooth
        # boundary vertex.  Non-adjacent equal-distance walls remain distinct.
        for lane, point in enumerate(points):
            radius = distance[lane] + maximum_half + _scale_epsilon(point)
            nearby = np.asarray(tree.query_ball_point(point, radius), dtype=np.int64)
            if nearby.size < 2:
                continue
            nearby_distance, _ = _point_segment_rows(point, starts[nearby], ends[nearby])
            tolerance = 256.0 * np.finfo(np.float64).eps * max(
                float(np.max(np.abs(point))), float(distance[lane]), 1.0
            )
            tied = nearby[np.abs(nearby_distance - distance[lane]) <= tolerance]
            adjacent = [int(primary[lane])]
            for candidate in tied:
                if _ring_adjacent(
                    int(ring_indices[primary[lane]]),
                    int(ring_indices[candidate]),
                    ring_size,
                ):
                    adjacent.append(int(candidate))
            if len(adjacent) > 1:
                averaged = np.sum(inward[np.unique(adjacent)], axis=0)
                norm = float(np.linalg.norm(averaged))
                if norm > 0.0:
                    normal[lane] = averaged / norm
    return distance, normal, primary, projected


def _unique_nearest_wall_rows(
    points,
    distance,
    primary,
    starts,
    ends,
    ring_indices,
    ring_size,
    tree,
):
    """Distinguish a smooth tessellation vertex from two opposing nearest walls."""

    unique = np.ones(points.shape[0], dtype=bool)
    maximum_half = 0.5 * float(
        np.max(np.linalg.norm(ends - starts, axis=1), initial=0.0)
    )
    for lane, point in enumerate(points):
        radius = float(distance[lane]) + maximum_half + _scale_epsilon(point)
        nearby = np.asarray(tree.query_ball_point(point, radius), dtype=np.int64)
        if nearby.size < 2:
            continue
        nearby_distance, _ = _point_segment_rows(point, starts[nearby], ends[nearby])
        tolerance = 256.0 * np.finfo(np.float64).eps * max(
            float(np.max(np.abs(point))), float(distance[lane]), 1.0
        )
        tied = nearby[np.abs(nearby_distance - distance[lane]) <= tolerance]
        selected_ring_index = int(ring_indices[int(primary[lane])])
        for candidate in tied:
            if not _ring_adjacent(
                selected_ring_index,
                int(ring_indices[int(candidate)]),
                ring_size,
            ):
                unique[lane] = False
                break
    return unique


def _point_segment_rows(point, starts, ends):
    delta = ends - starts
    length_squared = np.einsum("ni,ni->n", delta, delta)
    numerator = np.einsum("ni,ni->n", point - starts, delta)
    # A zero-length segment is a point.  Projecting onto its start is both the
    # exact point-segment result and the scalar equivalent of the compiled
    # kernel's guarded projection.  Do not perturb valid short segments by
    # adding an arbitrary epsilon to their denominator.
    fraction = np.zeros_like(length_squared)
    np.divide(
        numerator,
        length_squared,
        out=fraction,
        where=length_squared > 0.0,
    )
    np.clip(fraction, 0.0, 1.0, out=fraction)
    projected = starts + fraction[:, None] * delta
    return np.linalg.norm(point - projected, axis=1), projected


def _segment_segment_distance(
    first_start,
    first_end,
    second_start,
    second_end,
):
    first_delta = first_end - first_start
    second_delta = second_end - second_start
    first_length_squared = float(first_delta @ first_delta)
    second_length_squared = float(second_delta @ second_delta)
    if first_length_squared <= 0.0:
        return float(
            _point_segment_rows(
                first_start,
                second_start[None, :],
                second_end[None, :],
            )[0][0]
        )
    if second_length_squared <= 0.0:
        return float(
            _point_segment_rows(
                second_start,
                first_start[None, :],
                first_end[None, :],
            )[0][0]
        )
    if _segments_intersect(first_start, first_end, second_start, second_end):
        return 0.0
    return min(
        float(
            _point_segment_rows(
                first_start, second_start[None, :], second_end[None, :]
            )[0][0]
        ),
        float(
            _point_segment_rows(first_end, second_start[None, :], second_end[None, :])[
                0
            ][0]
        ),
        float(
            _point_segment_rows(second_start, first_start[None, :], first_end[None, :])[
                0
            ][0]
        ),
        float(
            _point_segment_rows(second_end, first_start[None, :], first_end[None, :])[
                0
            ][0]
        ),
    )


def _segments_intersect(a, b, c, d):
    r = b - a
    s = d - c
    r_length_squared = float(r @ r)
    s_length_squared = float(s @ s)
    coordinate_scale = max(float(np.max(np.abs([a, b, c, d]))), 1.0)
    coordinate_tolerance = 128.0 * np.finfo(np.float64).eps * coordinate_scale
    if r_length_squared <= 0.0:
        distance = float(_point_segment_rows(a, c[None, :], d[None, :])[0][0])
        return distance <= coordinate_tolerance
    if s_length_squared <= 0.0:
        distance = float(_point_segment_rows(c, a[None, :], b[None, :])[0][0])
        return distance <= coordinate_tolerance

    denominator = _cross(r, s)
    denominator_scale = abs(float(r[0] * s[1])) + abs(float(r[1] * s[0]))
    denominator_tolerance = (
        128.0 * np.finfo(np.float64).eps * denominator_scale
    )
    offset = c - a
    if abs(denominator) > denominator_tolerance:
        t = _cross(offset, s) / denominator
        u = _cross(offset, r) / denominator
        parameter_tolerance = 128.0 * np.finfo(np.float64).eps
        return (
            -parameter_tolerance <= t <= 1.0 + parameter_tolerance
            and -parameter_tolerance <= u <= 1.0 + parameter_tolerance
        )

    collinearity = _cross(offset, r)
    collinearity_scale = abs(float(offset[0] * r[1])) + abs(
        float(offset[1] * r[0])
    )
    collinearity_tolerance = (
        128.0 * np.finfo(np.float64).eps * collinearity_scale
    )
    if abs(collinearity) > collinearity_tolerance:
        return False
    axis = int(np.argmax(np.abs(r)))
    return max(min(a[axis], b[axis]), min(c[axis], d[axis])) <= min(
        max(a[axis], b[axis]), max(c[axis], d[axis])
    ) + coordinate_tolerance


def _first_moving_point_segment_capsule_contact(
    start,
    end,
    wall_start,
    wall_end,
    radius,
    tolerance,
):
    start_distance = float(
        _point_segment_rows(start, wall_start[None, :], wall_end[None, :])[0][0]
    )
    if start_distance <= radius + tolerance:
        return 0.0
    direction = end - start
    if float(direction @ direction) == 0.0:
        return None
    tangent = wall_end - wall_start
    length = float(np.linalg.norm(tangent))
    tangent /= length
    normal = np.asarray([-tangent[1], tangent[0]])
    centre = 0.5 * (wall_start + wall_end)
    relative = start - centre
    u0 = float(relative @ tangent)
    v0 = float(relative @ normal)
    du = float(direction @ tangent)
    dv = float(direction @ normal)
    half = 0.5 * length
    best = math.inf

    if dv != 0.0:
        for target_v in (-radius, radius):
            fraction = (target_v - v0) / dv
            if 0.0 <= fraction <= 1.0:
                u = u0 + fraction * du
                if -half - tolerance <= u <= half + tolerance:
                    best = min(best, fraction)

    quadratic_a = float(direction @ direction)
    for endpoint in (wall_start, wall_end):
        offset = start - endpoint
        quadratic_b = 2.0 * float(offset @ direction)
        quadratic_c = float(offset @ offset) - radius * radius
        discriminant = quadratic_b * quadratic_b - 4.0 * quadratic_a * quadratic_c
        if discriminant < 0.0:
            continue
        root = math.sqrt(max(discriminant, 0.0))
        for fraction in (
            (-quadratic_b - root) / (2.0 * quadratic_a),
            (-quadratic_b + root) / (2.0 * quadratic_a),
        ):
            if 0.0 <= fraction <= 1.0:
                best = min(best, fraction)
    return None if not math.isfinite(best) else float(best)


def _contains_distinct_wall_indices(
    indices,
    ring_indices,
    ring_size,
    inward,
):
    unique = np.unique(np.asarray(indices, dtype=np.int64))
    if unique.size < 2:
        return False
    ordered = unique[np.argsort(ring_indices[unique], kind="stable")]
    clusters = [[int(ordered[0])]]
    for raw_index in ordered[1:]:
        index = int(raw_index)
        previous = clusters[-1][-1]
        if int(ring_indices[index]) - int(ring_indices[previous]) <= 1:
            clusters[-1].append(index)
        else:
            clusters.append([index])
    if (
        len(clusters) > 1
        and int(ring_indices[clusters[0][0]]) == 0
        and int(ring_indices[clusters[-1][-1]]) == ring_size - 1
    ):
        clusters[0] = clusters[-1] + clusters[0]
        clusters.pop()
    if len(clusters) < 2:
        return False
    representatives = [
        np.mean(inward[np.asarray(cluster, dtype=np.int64)], axis=0)
        for cluster in clusters
    ]
    for first_offset, first_normal in enumerate(representatives):
        first_norm = float(np.linalg.norm(first_normal))
        if first_norm <= 0.0:
            continue
        first_normal = first_normal / first_norm
        for second_normal in representatives[first_offset + 1 :]:
            second_norm = float(np.linalg.norm(second_normal))
            if second_norm <= 0.0:
                continue
            if float(first_normal @ (second_normal / second_norm)) < 1.0 - 1.0e-12:
                return True
    return False


def _ring_adjacent(first, second, ring_size):
    difference = abs(int(first) - int(second))
    return difference <= 1 or (ring_size > 1 and difference == ring_size - 1)


def _subdivide_ring(coordinates, maximum_length):
    starts = []
    ends = []
    for raw_start, raw_end in zip(coordinates[:-1], coordinates[1:]):
        delta = raw_end - raw_start
        length = float(np.linalg.norm(delta))
        pieces = max(1, int(math.ceil(length / maximum_length)))
        values = np.linspace(0.0, 1.0, pieces + 1)
        samples = raw_start[None, :] + values[:, None] * delta[None, :]
        starts.extend(samples[:-1])
        ends.extend(samples[1:])
    return np.ascontiguousarray(np.vstack(starts)), np.ascontiguousarray(np.vstack(ends))


def _open_cap_segment_mask(starts, ends, sections):
    selected = np.zeros(starts.shape[0], dtype=bool)
    for point, normal, tangent, half_width in zip(
        sections["point_xz_um"],
        sections["normal_xz"],
        sections["tangent_xz"],
        sections["half_width_um"],
    ):
        scale = max(float(np.max(np.abs(point))), float(half_width), 1.0)
        tolerance = 4096.0 * np.finfo(np.float64).eps * scale
        start_relative = starts - point
        end_relative = ends - point
        on_plane = (
            np.abs(start_relative @ normal) <= tolerance
        ) & (np.abs(end_relative @ normal) <= tolerance)
        within = (
            np.abs(start_relative @ tangent) <= float(half_width) + tolerance
        ) & (np.abs(end_relative @ tangent) <= float(half_width) + tolerance)
        selected |= on_plane & within
    return selected


def _geometry_hash(
    starts,
    ends,
    sections,
    quad_segs,
    maximum_length,
):
    digest = sha256()
    digest.update(_GEOMETRY_SCHEMA.encode("utf-8"))
    digest.update(np.asarray([quad_segs, maximum_length], dtype=np.float64).tobytes())
    for array in (starts, ends):
        digest.update(np.ascontiguousarray(array, dtype=np.float64).tobytes())
    for key in (
        "point_xz_um",
        "normal_xz",
        "tangent_xz",
        "half_width_um",
        "label",
        "kind",
        "vessel_id",
    ):
        digest.update(key.encode("utf-8"))
        digest.update(np.ascontiguousarray(sections[key]).tobytes())
    return digest.hexdigest().upper()


def _vessel_xz(vessel):
    p0 = np.asarray([vessel.x_p[0], vessel.x_p[2]], dtype=np.float64)
    p1 = np.asarray([vessel.x_d[0], vessel.x_d[2]], dtype=np.float64)
    delta = p1 - p0
    length = float(np.linalg.norm(delta))
    return p0, p1, delta / length


def _open_endpoint_keys(vessels):
    result = set()
    for vessel in vessels:
        if int(vessel.parent_id) < 0:
            result.add(_endpoint_key(np.asarray([vessel.x_p[0], vessel.x_p[2]])))
        if len(vessel.children) == 0:
            result.add(_endpoint_key(np.asarray([vessel.x_d[0], vessel.x_d[2]])))
    return result


def _endpoint_key(point):
    return round(float(point[0]), 6), round(float(point[1]), 6)


def _repair_polygon(geometry):
    repaired = make_valid(geometry) if make_valid is not None and not geometry.is_valid else geometry
    if not repaired.is_valid:
        repaired = repaired.buffer(0)
    return repaired


def _remove_polygon_holes(geometry):
    if geometry.is_empty:
        return geometry
    if isinstance(geometry, Polygon):
        return Polygon(geometry.exterior)
    if isinstance(geometry, MultiPolygon):
        return unary_union([Polygon(part.exterior) for part in geometry.geoms if not part.is_empty])
    if isinstance(geometry, GeometryCollection):
        parts = [
            _remove_polygon_holes(part)
            for part in geometry.geoms
            if isinstance(part, (Polygon, MultiPolygon))
        ]
        return unary_union(parts) if parts else geometry
    return geometry


def _flatten_points(values):
    array = np.asarray(values, dtype=np.float64)
    if array.shape == (2,):
        return array.reshape(1, 2), ()
    if array.ndim < 2 or array.shape[-1] != 2:
        raise ValueError("X-Z coordinates must have final dimension 2.")
    return array.reshape(-1, 2), array.shape[:-1]


def _restore(values, shape):
    reshaped = np.asarray(values).reshape(shape or (1,))
    return reshaped if shape else reshaped[0]


def _restore_vectors(values, shape):
    reshaped = np.asarray(values).reshape(*(shape or (1,)), 2)
    return reshaped if shape else reshaped[0]


def _scale_epsilon(point):
    return 512.0 * np.finfo(np.float64).eps * max(
        float(np.max(np.abs(point))), 1.0
    )


def _cross(first, second):
    return float(first[0] * second[1] - first[1] * second[0])
