"""Static, open-boundary-aware molecular target fields on solid vessel walls.

The molecular binding model treats endothelial targets as a continuous surface
density rather than as randomly placed point receptors.  This module performs
the geometry-only part of that model:

* identify solid wall sites already accepted by the particle hydrodynamic
  fields while excluding inlet and outlet openings;
* intersect those wall sites with a physical-space Boolean mask selected from
  topology-based candidate vessel beds;
* evaluate the target indicator along a locally planar wall tangent; and
* integrate the target-positive part of a projected reaction disk.

All public coordinates are physical X-Z coordinates in micrometres.  Target
density is stored in molecules per square metre, while reaction areas are
returned in square micrometres so that the caller must perform the unit
conversion explicitly when forming molecule counts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy.spatial import cKDTree

from .molecular_binding import target_positive_reaction_area_um2
from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from ..particles.particle_hydrodynamic_fields import ParticleHydrodynamicFields
from ..core.types import GridDomain


_MISSING = object()
_SUPPORTED_REGION_MODES = frozenset(
    {"disabled", "mask_npz", "continuous_wall_npz"}
)


@dataclass(frozen=True)
class MolecularTargetField:
    """Immutable spatial inputs used by deterministic molecular binding.

    ``solid_wall_mask`` and ``target_wall_mask`` are cell-centred overlays whose
    shape matches the CFD domain.  Molecular reaction geometry itself uses the
    exact authoritative solid-face centres, axes, lengths, and normals shared
    with particle boundary construction; an opening is never turned into a
    target-bearing end cap.
    """

    enabled: bool
    region_mode: str
    target_density_molecules_per_m2: float
    spacing_um: float
    x_coordinates_um: np.ndarray
    z_coordinates_um: np.ndarray
    solid_wall_mask: np.ndarray
    target_wall_mask: np.ndarray
    target_density_field_molecules_per_m2: np.ndarray
    wall_coordinates_xz_um: np.ndarray
    wall_normal_xz: np.ndarray
    wall_axis: np.ndarray
    wall_length_um: np.ndarray
    wall_start_xz_um: np.ndarray
    wall_end_xz_um: np.ndarray
    wall_tangent_xz: np.ndarray
    wall_ring_index: np.ndarray
    wall_arclength_start_um: np.ndarray
    wall_arclength_end_um: np.ndarray
    boundary_ring_length_um: float
    wall_target_positive: np.ndarray
    open_boundary_mask: np.ndarray
    _maximum_wall_half_length_um: float = field(repr=False, compare=False)
    _boundary_geometry: ContinuousVesselGeometry = field(repr=False, compare=False)
    region_x_coordinates_um: np.ndarray | None = None
    region_z_coordinates_um: np.ndarray | None = None
    region_mask: np.ndarray | None = None
    source_mask_npz_path: Path | None = None
    _wall_tree: cKDTree | None = field(default=None, repr=False, compare=False)
    _target_positive_wall_tree: cKDTree | None = field(
        default=None, repr=False, compare=False
    )
    _maximum_target_face_half_length_um: float = field(
        default=float("nan"), repr=False, compare=False
    )

    @property
    def shape(self) -> tuple[int, int]:
        """Return the cell-centred X-Z field shape."""

        return (int(self.x_coordinates_um.size), int(self.z_coordinates_um.size))

    def polyline_target_candidate_mask(
        self,
        path_points_xz_um: np.ndarray,
        path_point_offsets: np.ndarray,
        bubble_radius_um: np.ndarray | float,
        capture_distance_um: float,
    ) -> np.ndarray:
        """Conservatively flag CSR polylines that may reach target-positive wall.

        Lane ``i`` owns ``points[offsets[i]:offsets[i + 1]]``.  Every chord is
        enclosed by a ball about its midpoint whose radius is its half-length.
        Expanding that ball by ``R + c + sqrt(2 R c)`` and the maximum
        target-positive face half-length makes the face-centre query
        conservative even when target support continues around a boundary-ring
        corner.  A rejected lane therefore cannot reach any target-positive
        wall element.  A lane containing one point is treated as a stationary
        zero-length chord.

        Historical or test-constructed fields may not carry the optional
        target-positive spatial index.  In that case every lane is retained so
        this broad phase can never alter molecular exposure results.
        """

        points = _as_point_rows(path_points_xz_um, "path_points_xz_um")
        offsets = np.asarray(path_point_offsets, dtype=np.int64)
        if offsets.ndim != 1 or offsets.size < 1:
            raise ValueError("path_point_offsets must be a non-empty one-dimensional array.")
        lane_count = int(offsets.size - 1)
        if (
            int(offsets[0]) != 0
            or int(offsets[-1]) != int(points.shape[0])
            or np.any(np.diff(offsets) < 1)
        ):
            raise ValueError(
                "path_point_offsets must start at zero, end at point_count, and "
                "assign at least one point to every lane."
            )

        radii = _broadcast_scalars(
            bubble_radius_um, lane_count, "bubble_radius_um"
        )
        if np.any(~np.isfinite(radii)) or np.any(radii < 0.0):
            raise ValueError(
                "bubble_radius_um must contain finite, non-negative values."
            )
        capture = float(capture_distance_um)
        if not np.isfinite(capture) or capture < 0.0:
            raise ValueError(
                "capture_distance_um must be finite and non-negative."
            )
        if lane_count == 0:
            return np.empty(0, dtype=bool)

        target_tree = self._target_positive_wall_tree
        maximum_face_half_length = float(
            self._maximum_target_face_half_length_um
        )
        if (
            target_tree is None
            or int(getattr(target_tree, "n", 0)) <= 0
            or not np.isfinite(maximum_face_half_length)
            or maximum_face_half_length < 0.0
        ):
            return np.ones(lane_count, dtype=bool)

        point_counts = np.diff(offsets)
        segment_counts = point_counts - 1
        segment_lane = np.repeat(
            np.arange(lane_count, dtype=np.int64), segment_counts
        )
        possible_segment_start = np.arange(
            max(0, points.shape[0] - 1), dtype=np.int64
        )
        if possible_segment_start.size:
            lane_end = np.zeros(possible_segment_start.size, dtype=bool)
            boundary_start = offsets[1:-1] - 1
            lane_end[boundary_start] = True
            segment_start_index = possible_segment_start[~lane_end]
        else:
            segment_start_index = np.empty(0, dtype=np.int64)

        if segment_start_index.size:
            segment_start = points[segment_start_index]
            segment_end = points[segment_start_index + 1]
            segment_delta = segment_end - segment_start
            midpoint = segment_start + 0.5 * segment_delta
            half_chord = 0.5 * np.linalg.norm(segment_delta, axis=1)
            coordinate_scale = np.maximum(
                np.max(np.abs(segment_start), axis=1),
                np.max(np.abs(segment_end), axis=1),
            )
        else:
            midpoint = np.empty((0, 2), dtype=np.float64)
            half_chord = np.empty(0, dtype=np.float64)
            coordinate_scale = np.empty(0, dtype=np.float64)

        stationary_lane = np.flatnonzero(point_counts == 1)
        if stationary_lane.size:
            stationary_point = points[offsets[stationary_lane]]
            midpoint = np.vstack((midpoint, stationary_point))
            half_chord = np.concatenate(
                (half_chord, np.zeros(stationary_lane.size, dtype=np.float64))
            )
            coordinate_scale = np.concatenate(
                (
                    coordinate_scale,
                    np.max(np.abs(stationary_point), axis=1),
                )
            )
            query_lane = np.concatenate((segment_lane, stationary_lane))
        else:
            query_lane = segment_lane

        if (
            np.any(~np.isfinite(midpoint))
            or np.any(~np.isfinite(half_chord))
        ):
            raise ValueError(
                "path_points_xz_um produced a non-finite chord midpoint or length."
            )
        lane_radius = radii[query_lane]
        maximum_reaction_disk_radius = np.sqrt(
            2.0 * lane_radius * capture
        )
        roundoff_scale = np.maximum.reduce(
            (
                np.ones(query_lane.size, dtype=np.float64),
                coordinate_scale,
                lane_radius,
                np.full(query_lane.size, capture, dtype=np.float64),
                maximum_reaction_disk_radius,
                half_chord,
                np.full(
                    query_lane.size,
                    maximum_face_half_length,
                    dtype=np.float64,
                ),
                np.full(query_lane.size, abs(float(self.spacing_um)), dtype=np.float64),
            )
        )
        roundoff_um = 256.0 * np.finfo(np.float64).eps * roundoff_scale
        search_radius = (
            lane_radius
            + capture
            + maximum_reaction_disk_radius
            + half_chord
            + maximum_face_half_length
            + roundoff_um
        )
        if np.any(~np.isfinite(search_radius)):
            raise ValueError("The target broad-phase search radius is not finite.")

        candidate_count = np.asarray(
            target_tree.query_ball_point(
                midpoint,
                search_radius,
                return_length=True,
            ),
            dtype=np.int64,
        ).reshape(-1)
        candidate = np.zeros(lane_count, dtype=bool)
        np.logical_or.at(candidate, query_lane, candidate_count > 0)
        return candidate

    def reaction_area_um2(
        self,
        points_xz_um: np.ndarray,
        tangents_xz: np.ndarray,
        reaction_radius_um: np.ndarray | float,
    ) -> np.ndarray:
        r"""Integrate the target-positive projected reaction area.

        The Revised-v7 local-plane expression is

        .. math::

           A_T = 2\int_{-a}^{a}\chi_T(s+\eta)
                 \sqrt{a^2-\eta^2}\,\mathrm d\eta.

        The target ROI and eligible closed-wall support are converted into
        exact signed tangent intervals.  The circular-segment antiderivative is
        then evaluated analytically, so target-edge crossings vary continuously
        and a fully target-positive disk gives exactly ``pi * a**2``.
        """

        points = _as_point_rows(points_xz_um, "points_xz_um")
        tangents = _broadcast_point_rows(tangents_xz, points.shape[0], "tangents_xz")
        radii = _broadcast_scalars(reaction_radius_um, points.shape[0], "reaction_radius_um")
        if np.any(~np.isfinite(radii)) or np.any(radii < 0.0):
            raise ValueError("reaction_radius_um must contain finite, non-negative values.")

        areas = np.zeros(points.shape[0], dtype=np.float64)
        active = radii > 0.0
        if not self.enabled or not np.any(active):
            return areas

        tangent_norm = np.linalg.norm(tangents, axis=1)
        if np.any(~np.isfinite(tangent_norm)) or np.any(tangent_norm <= 0.0):
            raise ValueError("tangents_xz must contain finite, non-zero vectors.")
        unit_tangents = tangents / tangent_norm[:, None]

        # Project the particle centre itself onto the authoritative solid-face
        # union.  Reconstructing a wall point from an interpolated distance and
        # normal can choose the wrong wall at a raster corner.  The returned tie
        # list keeps equally close, distinct faces separate; their normals are
        # never averaged into a direction that does not belong to either wall.
        _projected, _face_indices, tied_faces = self._nearest_solid_face_projections(
            points
        )
        for index in np.flatnonzero(active):
            radius = float(radii[index])
            maximum_area = 0.0
            for raw_face_index in tied_faces[index]:
                face_index = int(raw_face_index)
                point = _project_point_to_wall_segment(
                    points[index],
                    self.wall_start_xz_um[face_index],
                    self.wall_end_xz_um[face_index],
                )
                normal = np.asarray(self.wall_normal_xz[face_index], dtype=np.float64)
                tangent = np.asarray([-normal[1], normal[0]], dtype=np.float64)
                if float(np.dot(tangent, unit_tangents[index])) < 0.0:
                    tangent *= -1.0
                target_starts, target_ends = self._target_wall_intervals_along_tangent(
                    point,
                    tangent,
                    radius,
                    face_index,
                )
                candidate_area = target_positive_reaction_area_um2(
                    radius,
                    target_starts,
                    target_ends,
                )
                # Distinct tied planes describe alternative local wall frames,
                # not additive copies of the same reaction disk.  Taking the
                # largest target-positive area avoids a false zero at a mixed
                # target corner while keeping area <= pi * radius**2.
                maximum_area = max(maximum_area, float(candidate_area))
            areas[index] = maximum_area
        return areas

    def to_npz_payload(self, prefix: str = "") -> dict[str, np.ndarray]:
        """Return plain NumPy values ready to merge into an ``np.savez`` call."""

        if not isinstance(prefix, str):
            raise TypeError("prefix must be a string.")
        payload: dict[str, np.ndarray] = {
            f"{prefix}enabled": np.asarray(self.enabled, dtype=bool),
            f"{prefix}region_mode": np.asarray(self.region_mode),
            f"{prefix}target_density_molecules_per_m2": np.asarray(
                self.target_density_molecules_per_m2, dtype=np.float64
            ),
            f"{prefix}spacing_um": np.asarray(self.spacing_um, dtype=np.float64),
            f"{prefix}x_coordinates_um": np.asarray(self.x_coordinates_um, dtype=np.float64),
            f"{prefix}z_coordinates_um": np.asarray(self.z_coordinates_um, dtype=np.float64),
            f"{prefix}solid_wall_mask": np.asarray(self.solid_wall_mask, dtype=bool),
            f"{prefix}target_wall_mask": np.asarray(self.target_wall_mask, dtype=bool),
            f"{prefix}target_density_field_molecules_per_m2": np.asarray(
                self.target_density_field_molecules_per_m2, dtype=np.float64
            ),
            f"{prefix}open_boundary_mask": np.asarray(
                self.open_boundary_mask, dtype=bool
            ),
            f"{prefix}wall_coordinates_xz_um": np.asarray(self.wall_coordinates_xz_um, dtype=np.float64),
            f"{prefix}wall_normal_xz": np.asarray(self.wall_normal_xz, dtype=np.float64),
            f"{prefix}wall_axis": np.asarray(self.wall_axis, dtype=np.int8),
            f"{prefix}wall_length_um": np.asarray(self.wall_length_um, dtype=np.float64),
            f"{prefix}wall_start_xz_um": np.asarray(self.wall_start_xz_um, dtype=np.float64),
            f"{prefix}wall_end_xz_um": np.asarray(self.wall_end_xz_um, dtype=np.float64),
            f"{prefix}wall_tangent_xz": np.asarray(self.wall_tangent_xz, dtype=np.float64),
            f"{prefix}wall_ring_index": np.asarray(self.wall_ring_index, dtype=np.int64),
            f"{prefix}wall_arclength_start_um": np.asarray(
                self.wall_arclength_start_um, dtype=np.float64
            ),
            f"{prefix}wall_arclength_end_um": np.asarray(
                self.wall_arclength_end_um, dtype=np.float64
            ),
            f"{prefix}boundary_ring_length_um": np.asarray(
                self.boundary_ring_length_um, dtype=np.float64
            ),
            f"{prefix}wall_target_positive": np.asarray(
                self.wall_target_positive, dtype=bool
            ),
        }
        if self.region_mask is not None:
            payload[f"{prefix}region_x_coordinates_um"] = np.asarray(
                self.region_x_coordinates_um, dtype=np.float64
            )
            payload[f"{prefix}region_z_coordinates_um"] = np.asarray(
                self.region_z_coordinates_um, dtype=np.float64
            )
            payload[f"{prefix}region_mask"] = np.asarray(self.region_mask, dtype=bool)
        if self.source_mask_npz_path is not None:
            payload[f"{prefix}source_mask_npz_path"] = np.asarray(
                str(self.source_mask_npz_path.resolve())
            )
        return payload

    def nearest_solid_face_frame_xz_um(
        self,
        points_xz_um: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Project points onto ``Gamma_w`` and return its exact local wall frame.

        The returned tuple contains projected X-Z points, solid-face indices,
        inward normals, and tangents.  Equal-distance corner ties are resolved by
        the stable face-array order; normals from distinct walls are never
        averaged.
        """

        points = _as_point_rows(points_xz_um, "points_xz_um")
        projected, indices, _tied_faces = self._nearest_solid_face_projections(points)
        normals = np.asarray(self.wall_normal_xz[indices], dtype=np.float64).copy()
        tangents = np.column_stack((-normals[:, 1], normals[:, 0]))
        return projected, indices, normals, tangents

    def _nearest_solid_face_projections(
        self,
        points_xz_um: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, tuple[np.ndarray, ...]]:
        """Delegate exact nearest-face ownership to the canonical geometry."""

        projected, indices, tied_faces = (
            self._boundary_geometry.nearest_solid_face_projection_xz_um(
                points_xz_um
            )
        )
        return (
            np.asarray(projected, dtype=np.float64).reshape(-1, 2),
            np.asarray(indices, dtype=np.int64).reshape(-1),
            tuple(np.asarray(faces, dtype=np.int64) for faces in tied_faces),
        )

    def _target_wall_intervals_along_tangent(
        self,
        point_xz_um: np.ndarray,
        tangent_xz: np.ndarray,
        radius_um: float,
        face_index: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        if (
            self.boundary_ring_length_um > 0.0
            and self.wall_arclength_start_um.shape == self.wall_length_um.shape
            and np.all(np.isfinite(self.wall_arclength_start_um))
        ):
            return self._target_wall_intervals_along_arclength(
                point_xz_um, radius_um, face_index
            )
        if self._wall_tree is None or self.wall_coordinates_xz_um.shape[0] == 0:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        tolerance_um = 256.0 * np.finfo(np.float64).eps * max(
            1.0,
            float(self.spacing_um),
            float(radius_um),
        )
        candidate_indices = self._wall_tree.query_ball_point(
            point_xz_um,
            radius_um + float(self._maximum_wall_half_length_um) + tolerance_um,
        )
        starts: list[float] = []
        ends: list[float] = []
        for wall_index in candidate_indices:
            index = int(wall_index)
            if not bool(self.wall_target_positive[index]):
                continue
            face_tangent = np.asarray(self.wall_tangent_xz[index], dtype=np.float64)
            alignment = float(np.dot(face_tangent, tangent_xz))
            if abs(alignment) < 1.0 - 1.0e-12:
                continue
            relative = self.wall_coordinates_xz_um[index] - point_xz_um
            perpendicular = abs(
                float(relative[0] * tangent_xz[1] - relative[1] * tangent_xz[0])
            )
            if perpendicular > tolerance_um:
                continue
            centre = float(np.dot(relative, tangent_xz))
            half_length = 0.5 * float(self.wall_length_um[index])
            projected_half_length = abs(alignment) * half_length
            starts.append(max(-radius_um, centre - projected_half_length))
            ends.append(min(radius_um, centre + projected_half_length))
        return _merge_intervals(starts, ends, -radius_um, radius_um)

    def _target_wall_intervals_along_arclength(
        self,
        point_xz_um: np.ndarray,
        radius_um: float,
        face_index: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Map target support by continuous boundary arclength near one wall."""

        index = int(face_index)
        segment = self.wall_end_xz_um[index] - self.wall_start_xz_um[index]
        length_squared = float(segment @ segment)
        fraction = 0.0 if length_squared <= 0.0 else float(
            np.clip(
                ((point_xz_um - self.wall_start_xz_um[index]) @ segment)
                / length_squared,
                0.0,
                1.0,
            )
        )
        centre_arclength = float(self.wall_arclength_start_um[index]) + (
            fraction * float(self.wall_length_um[index])
        )
        total = float(self.boundary_ring_length_um)
        starts: list[float] = []
        ends: list[float] = []
        arc_start = self.wall_arclength_start_um
        arc_end = self.wall_arclength_end_um
        for shift in (-total, 0.0, total):
            lower = centre_arclength - radius_um - shift
            upper = centre_arclength + radius_um - shift
            first = int(np.searchsorted(arc_end, lower, side="left"))
            last = int(np.searchsorted(arc_start, upper, side="right"))
            for candidate in range(first, min(last, arc_start.size)):
                if not bool(self.wall_target_positive[candidate]):
                    continue
                relative_start = float(arc_start[candidate] + shift - centre_arclength)
                relative_end = float(arc_end[candidate] + shift - centre_arclength)
                starts.append(max(-radius_um, relative_start))
                ends.append(min(radius_um, relative_end))
        return _merge_intervals(starts, ends, -radius_um, radius_um)


def build_molecular_target_field(
    domain: GridDomain,
    hydrodynamic_fields: ParticleHydrodynamicFields,
    config: object | Mapping[str, Any] | None,
) -> MolecularTargetField:
    """
    Build a fixed target-density field from a configuration object or mapping.
    """
    # =============================================================================
    # ==== Read and check the basic grid and wall information.  
    # =============================================================================
    # Read the grid size and the width of one grid cell.
    shape       = tuple(int(value) for value in domain.shape)
    spacing_um  = float(domain.spacing_um)

    # Copy the physical X and Z positions of the grid cells into NumPy arrays.
    # These positions let us change a grid address, such as row 10, into a real position in micrometres.
    x_um = np.asarray(domain.x_coordinates_um, dtype=np.float64)
    z_um = np.asarray(domain.z_coordinates_um, dtype=np.float64)

    # Check the basic grid information before using it.
    # A bad grid would place targets at the wrong positions, so it is safer to stop with a clear message.
    if len(shape) != 2 or min(shape) < 2:
        raise ValueError("Molecular targets require a two-dimensional grid of at least 2 x 2 cells.")
    if not np.isfinite(spacing_um) or spacing_um <= 0.0:
        raise ValueError("domain.spacing_um must be finite and positive.")
    if x_um.shape != (shape[0],) or z_um.shape != (shape[1],):
        raise ValueError("Domain coordinate arrays must match domain.shape.")
    if np.any(~np.isfinite(x_um)) or np.any(~np.isfinite(z_um)):
        raise ValueError("Domain coordinate arrays must contain only finite values.")
    if np.any(np.diff(x_um) <= 0.0) or np.any(np.diff(z_um) <= 0.0):
        raise ValueError("Domain coordinate arrays must be strictly increasing.")

    # Read the on/off switch before requiring target geometry.  A disabled
    # target model is a strict no-op, so it must not fail only because an older
    # caller did not prepare molecular wall support that will never be used.
    enabled = bool(_config_value(config, ("enabled",), default=False))

    # =============================================================================
    # Read the user's three wall-related maps and check that they match the grid.
    # =============================================================================
    # Read the one boundary object shared by CFD openings, particle contact, and molecular targets.
    # Its solid faces are the real target-support surface; no shifted wall-cell approximation is built here.
    boundary_geometry = hydrodynamic_fields.boundary_geometry
    solid_sites   = np.asarray(hydrodynamic_fields.solid_site_mask, dtype=bool)
    open_boundary = np.asarray(hydrodynamic_fields.open_boundary_mask, dtype=bool)

    # Make sure all three maps fit the same grid.
    if solid_sites.shape != shape:
        raise ValueError("hydrodynamic_fields.solid_site_mask must match domain.shape.")
    if open_boundary.shape != shape:
        raise ValueError("hydrodynamic_fields.open_boundary_mask must match domain.shape.")
    # Open faces were already removed when Gamma_w was built.  Do not dilate the
    # opening mask: doing so would erase legitimate side-wall targets near an outlet.
    solid_wall_mask = np.ascontiguousarray(solid_sites, dtype=bool)

    # =============================================================================
    # Extract wall cell coordinates and normals.
    # =============================================================================
    # Keep wall-cell addresses for the saved Boolean overlay used by the selector and visualizer.
    wall_indices    = np.argwhere(solid_wall_mask)
    wall_coordinates = np.asarray(
        boundary_geometry.solid_face_center_xz_um,
        dtype=np.float64,
    ).copy()
    wall_normals = np.asarray(
        boundary_geometry.solid_face_inward_normal_xz,
        dtype=np.float64,
    )
    wall_length_um = np.asarray(
        boundary_geometry.solid_face_length_um,
        dtype=np.float64,
    )
    wall_axis = np.full(wall_coordinates.shape[0], -1, dtype=np.int8)
    wall_start = np.asarray(
        boundary_geometry.solid_face_start_xz_um, dtype=np.float64
    )
    wall_end = np.asarray(
        boundary_geometry.solid_face_end_xz_um, dtype=np.float64
    )
    wall_tangent = wall_end - wall_start
    wall_tangent /= np.linalg.norm(wall_tangent, axis=1)[:, None]
    wall_selection_coordinates = (
        wall_coordinates + 0.5 * spacing_um * wall_normals
    )
    wall_ring_index = np.asarray(
        boundary_geometry.solid_face_ring_index, dtype=np.int64
    )
    wall_arclength_start = np.asarray(
        boundary_geometry.solid_face_arclength_start_um,
        dtype=np.float64,
    )
    wall_arclength_end = np.asarray(
        boundary_geometry.solid_face_arclength_end_um,
        dtype=np.float64,
    )
    boundary_ring_length = float(boundary_geometry.boundary_ring_length_um)

    # =============================================================================
    # Use a KD-tree to find the nearest wall point for any given X-Z position.
    # =============================================================================
    # Build a fast face-centre lookup table; exact segment distances are checked after the lookup.
    wall_tree = cKDTree(wall_coordinates) if wall_coordinates.shape[0] else None
    
    # =============================================================================
    # Use the configuration to determine whether molecular targets are enabled and how they are defined.
    # =============================================================================
    # Read the user's main on/off switch for molecular targets.
    # Keeping this choice here gives the rest of the program one clear result to use.
    # When targets are disabled, prepare empty values with the correct shapes.
    # Later code can then use the same result object without needing many special cases.
    if not enabled:
        region_mode             = "disabled"
        density                 = 0.0
        region_x                = None
        region_z                = None
        region_mask             = None
        source_mask_npz_path    = None
        target_wall_mask        = np.zeros(shape, dtype=bool)
        loaded_continuous_wall_flags = None
    else:
        # Only the supported mask choice is accepted so the target cannot be placed by a typo.
        region_mode = str(_config_value(config, ("region_mode", "mode"))).strip().lower()
        if region_mode not in _SUPPORTED_REGION_MODES or region_mode == "disabled":
            raise ValueError(
                "Enabled molecular targets require region_mode 'mask_npz' "
                "or 'continuous_wall_npz'."
            )
        
        # Read how many target molecules cover one square metre of selected wall.
        # A negative or invalid amount has no physical meaning, so it must be rejected.
        density = float(_config_value(config, ("target_density_molecules_per_m2",)))
        if not np.isfinite(density) or density < 0.0:
            raise ValueError("target_density_molecules_per_m2 must be finite and non-negative when targets are enabled.")

        # Start all optional region data as empty.
        region_x                = None
        region_z                = None
        region_mask             = None
        source_mask_npz_path    = None
        loaded_continuous_wall_flags = None

        # The manual candidate union or automatic positive wall patches are stored as an NPZ True/False target map.
        # We also read the map's own X and Z positions so it can be matched to this simulation grid.
        mask_path               = Path(_config_value(config, ("mask_npz_path", "mask_path", "region_mask_path", "target_region_mask_path")))
        source_mask_npz_path    = mask_path.resolve()
        if region_mode == "continuous_wall_npz":
            loaded_continuous_wall_flags = _load_continuous_wall_target_npz(
                mask_path, boundary_geometry
            )
            target_wall_mask = _continuous_wall_flags_to_grid_overlay(
                wall_indices,
                x_um,
                z_um,
                shape,
                boundary_geometry,
                loaded_continuous_wall_flags,
            )
        else:
            x_key = str(
                _config_value(config, ("x_coordinates_key",), default="x_um")
            )
            z_key = str(
                _config_value(config, ("z_coordinates_key",), default="z_um")
            )
            mask_key = str(
                _config_value(config, ("mask_array_key",), default="target_mask")
            )
            region_x, region_z, region_mask = _load_region_mask_npz(
                mask_path,
                x_coordinates_key=x_key,
                z_coordinates_key=z_key,
                mask_array_key=mask_key,
            )

            grid_aligned = (
                region_mask.shape == shape
                and np.array_equal(region_x, x_um)
                and np.array_equal(region_z, z_um)
            )
            if grid_aligned:
                target_wall_mask = np.asarray(region_mask, dtype=bool).copy()
            else:
                wall_site_coordinates = np.column_stack(
                    (x_um[wall_indices[:, 0]], z_um[wall_indices[:, 1]])
                )
                target_flags = _sample_nearest_mask(
                    region_mask,
                    region_x,
                    region_z,
                    wall_site_coordinates,
                )
                target_wall_mask = np.zeros(shape, dtype=bool)
                if wall_indices.shape[0] > 0:
                    target_wall_mask[
                        wall_indices[:, 0], wall_indices[:, 1]
                    ] = target_flags

        # This final intersection guarantees that targets remain on valid solid-wall cells only.
        target_wall_mask &= solid_wall_mask

    # =============================================================================
    # Determine which wall cells are eligible for molecular targets and how many molecules they have.
    # =============================================================================
    # Create the final number-density grid.
    # Non-target cells stay zero, while selected wall cells receive the configured molecule density.
    density_field = np.zeros(shape, dtype=np.float64)
    density_field[target_wall_mask] = density

    # Convert the cell overlay into one immutable flag per authoritative solid
    # face.  Runtime molecular queries then operate on Gamma_w directly instead
    # of sampling a half-grid coordinate where nearest-cell ties can favour one
    # side of an otherwise symmetric vessel.
    if loaded_continuous_wall_flags is not None:
        wall_target_positive = loaded_continuous_wall_flags
    elif region_mask is not None and np.array_equal(
        region_mask, solid_wall_mask
    ):
        wall_target_positive = np.ones(
            wall_coordinates.shape[0], dtype=bool
        )
    else:
        wall_target_positive = (
            np.zeros(wall_coordinates.shape[0], dtype=bool)
            if region_mask is None
            else _sample_nearest_mask(
                region_mask,
                np.asarray(region_x, dtype=np.float64),
                np.asarray(region_z, dtype=np.float64),
                wall_selection_coordinates,
            )
        )
    target_positive_face_indices = np.flatnonzero(wall_target_positive)
    if target_positive_face_indices.size:
        target_positive_wall_coordinates = np.ascontiguousarray(
            wall_coordinates[target_positive_face_indices], dtype=np.float64
        )
        target_positive_wall_tree = cKDTree(
            target_positive_wall_coordinates, copy_data=True
        )
        maximum_target_face_half_length_um = 0.5 * float(
            np.max(wall_length_um[target_positive_face_indices])
        )
    else:
        target_positive_wall_tree = None
        maximum_target_face_half_length_um = float("nan")

    # Pack every prepared value into one object for the particle and bonding code.
    # Contiguous arrays store values next to each other in memory, which makes repeated calculations faster.
    return MolecularTargetField(
        enabled=enabled,
        region_mode=region_mode,
        target_density_molecules_per_m2=float(density),
        spacing_um=spacing_um,
        x_coordinates_um=np.ascontiguousarray(x_um),
        z_coordinates_um=np.ascontiguousarray(z_um),
        solid_wall_mask=solid_wall_mask,
        target_wall_mask=np.ascontiguousarray(target_wall_mask, dtype=bool),
        target_density_field_molecules_per_m2=np.ascontiguousarray(density_field),
        wall_coordinates_xz_um=np.ascontiguousarray(wall_coordinates),
        wall_normal_xz=np.ascontiguousarray(wall_normals),
        wall_axis=np.ascontiguousarray(wall_axis, dtype=np.int8),
        wall_length_um=np.ascontiguousarray(wall_length_um, dtype=np.float64),
        wall_start_xz_um=np.ascontiguousarray(wall_start, dtype=np.float64),
        wall_end_xz_um=np.ascontiguousarray(wall_end, dtype=np.float64),
        wall_tangent_xz=np.ascontiguousarray(wall_tangent, dtype=np.float64),
        wall_ring_index=np.ascontiguousarray(wall_ring_index, dtype=np.int64),
        wall_arclength_start_um=np.ascontiguousarray(
            wall_arclength_start, dtype=np.float64
        ),
        wall_arclength_end_um=np.ascontiguousarray(
            wall_arclength_end, dtype=np.float64
        ),
        boundary_ring_length_um=float(boundary_ring_length),
        wall_target_positive=np.ascontiguousarray(wall_target_positive, dtype=bool),
        open_boundary_mask=np.ascontiguousarray(open_boundary, dtype=bool),
        _maximum_wall_half_length_um=(
            0.5 * float(np.max(wall_length_um))
            if wall_length_um.size
            else 0.0
        ),
        _boundary_geometry=boundary_geometry,
        region_x_coordinates_um=None if region_x is None else np.ascontiguousarray(region_x),
        region_z_coordinates_um=None if region_z is None else np.ascontiguousarray(region_z),
        region_mask=None if region_mask is None else np.ascontiguousarray(region_mask, dtype=bool),
        source_mask_npz_path=source_mask_npz_path,
        _wall_tree=wall_tree,
        _target_positive_wall_tree=target_positive_wall_tree,
        _maximum_target_face_half_length_um=maximum_target_face_half_length_um,
    )


def _config_value(
    config: object | Mapping[str, Any] | None,
    names: tuple[str, ...],
    *,
    default: object = _MISSING,
) -> Any:
    if config is not None:
        for name in names:
            if isinstance(config, Mapping) and name in config:
                return config[name]
            if hasattr(config, name):
                return getattr(config, name)
    if default is not _MISSING:
        return default
    joined = " or ".join(repr(name) for name in names)
    raise ValueError(f"Molecular target configuration is missing {joined}.")


def _project_point_to_wall_segment(
    point_xz_um: np.ndarray,
    start_xz_um: np.ndarray,
    end_xz_um: np.ndarray,
) -> np.ndarray:
    """Project one point onto an arbitrary continuous-wall element."""

    point = np.asarray(point_xz_um, dtype=np.float64)
    start = np.asarray(start_xz_um, dtype=np.float64)
    end = np.asarray(end_xz_um, dtype=np.float64)
    delta = end - start
    fraction = float(np.clip(((point - start) @ delta) / (delta @ delta), 0.0, 1.0))
    return start + fraction * delta


def _load_region_mask_npz(
    path: Path,
    *,
    x_coordinates_key: str,
    z_coordinates_key: str,
    mask_array_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(f"Molecular target mask NPZ does not exist: {path}")
    with np.load(path, allow_pickle=False) as data:
        x_key = _resolve_npz_key(data, x_coordinates_key, ("x_coordinates_um", "x_um"))
        z_key = _resolve_npz_key(data, z_coordinates_key, ("z_coordinates_um", "z_um"))
        mask_key = _resolve_npz_key(data, mask_array_key, ("target_region_mask", "target_mask"))
        x_um = np.asarray(data[x_key], dtype=np.float64)
        z_um = np.asarray(data[z_key], dtype=np.float64)
        raw_mask = np.asarray(data[mask_key])
        if raw_mask.dtype.kind != "b":
            raise ValueError(f"{mask_key} must be stored as a Boolean array.")
        mask = np.asarray(raw_mask, dtype=bool)

    if x_um.ndim != 1 or z_um.ndim != 1 or min(x_um.size, z_um.size) < 2:
        raise ValueError("Mask x_coordinates_um and z_coordinates_um must be one-dimensional with at least two values.")
    if np.any(~np.isfinite(x_um)) or np.any(~np.isfinite(z_um)):
        raise ValueError("Mask physical coordinates must contain only finite values.")
    if np.any(np.diff(x_um) <= 0.0) or np.any(np.diff(z_um) <= 0.0):
        raise ValueError("Mask physical coordinates must be strictly increasing.")
    expected_shape = (x_um.size, z_um.size)
    if mask.shape != expected_shape:
        raise ValueError(f"Target mask shape {mask.shape} does not match coordinate shape {expected_shape}.")
    return np.ascontiguousarray(x_um), np.ascontiguousarray(z_um), np.ascontiguousarray(mask)


def _load_continuous_wall_target_npz(
    path: Path,
    boundary_geometry: ContinuousVesselGeometry,
) -> np.ndarray:
    """Load target flags tied to one immutable v16 boundary arclength."""

    resolved = Path(path).expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(
            f"Continuous molecular target file does not exist: {resolved}"
        )
    required = {
        "target_geometry_schema",
        "continuous_geometry_hash_sha256",
        "wall_ring_index",
        "wall_arclength_start_um",
        "wall_arclength_end_um",
        "wall_target_positive",
    }
    with np.load(resolved, allow_pickle=False) as data:
        missing = sorted(required.difference(data.files))
        if missing:
            raise ValueError(
                "Continuous molecular target NPZ is missing arrays: "
                + ", ".join(missing)
            )
        schema = str(np.asarray(data["target_geometry_schema"]).item())
        if schema != "v16_continuous_wall_arclength_target":
            raise ValueError(
                "Continuous molecular target NPZ has unsupported schema "
                f"{schema!r}."
            )
        expected_hash = str(boundary_geometry.geometry_hash_sha256)
        stored_hash = str(
            np.asarray(data["continuous_geometry_hash_sha256"]).item()
        )
        if stored_hash != expected_hash:
            raise ValueError(
                "Continuous molecular target geometry hash does not match the "
                "active v16 vessel boundary. Regenerate the target artifact."
            )
        ring_index = np.asarray(data["wall_ring_index"], dtype=np.int64)
        arc_start = np.asarray(data["wall_arclength_start_um"], dtype=np.float64)
        arc_end = np.asarray(data["wall_arclength_end_um"], dtype=np.float64)
        flags_raw = np.asarray(data["wall_target_positive"])
    expected_ring = np.asarray(
        boundary_geometry.solid_face_ring_index, dtype=np.int64
    )
    expected_start = np.asarray(
        boundary_geometry.solid_face_arclength_start_um, dtype=np.float64
    )
    expected_end = np.asarray(
        boundary_geometry.solid_face_arclength_end_um, dtype=np.float64
    )
    if flags_raw.dtype.kind != "b":
        raise ValueError("wall_target_positive must be stored as a Boolean array.")
    if (
        ring_index.shape != expected_ring.shape
        or not np.array_equal(ring_index, expected_ring)
        or not np.array_equal(arc_start, expected_start)
        or not np.array_equal(arc_end, expected_end)
    ):
        raise ValueError(
            "Continuous molecular target wall arclength arrays do not exactly "
            "match the active v16 geometry."
        )
    flags = np.asarray(flags_raw, dtype=bool)
    if flags.shape != expected_ring.shape:
        raise ValueError(
            "wall_target_positive must contain one flag per v16 solid-wall element."
        )
    return np.ascontiguousarray(flags)


def _continuous_wall_flags_to_grid_overlay(
    wall_indices: np.ndarray,
    x_um: np.ndarray,
    z_um: np.ndarray,
    shape: tuple[int, int],
    boundary_geometry: ContinuousVesselGeometry,
    wall_target_positive: np.ndarray,
) -> np.ndarray:
    overlay = np.zeros(shape, dtype=bool)
    if wall_indices.shape[0] == 0:
        return overlay
    coordinates = np.column_stack(
        (x_um[wall_indices[:, 0]], z_um[wall_indices[:, 1]])
    )
    state = boundary_geometry.exact_solid_wall_state_xz_um_accelerated(
        coordinates
    )
    owner = np.asarray(state.primary_face_index, dtype=np.int64)
    overlay[wall_indices[:, 0], wall_indices[:, 1]] = wall_target_positive[owner]
    return overlay


def _resolve_npz_key(data: Mapping[str, np.ndarray], preferred: str, aliases: tuple[str, ...]) -> str:
    candidates: list[str] = []
    for key in (preferred, *aliases):
        if key and key not in candidates:
            candidates.append(key)
    for key in candidates:
        if key in data:
            return key
    raise ValueError(f"Molecular target mask NPZ is missing arrays; tried keys: {', '.join(candidates)}.")


def _merge_intervals(
    starts: list[float] | np.ndarray,
    ends: list[float] | np.ndarray,
    lower: float,
    upper: float,
) -> tuple[np.ndarray, np.ndarray]:
    start_array = np.asarray(starts, dtype=np.float64)
    end_array = np.asarray(ends, dtype=np.float64)
    if start_array.size == 0:
        return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
    valid = end_array > start_array
    start_array = np.maximum(start_array[valid], lower)
    end_array = np.minimum(end_array[valid], upper)
    valid = end_array > start_array
    start_array = start_array[valid]
    end_array = end_array[valid]
    if start_array.size == 0:
        return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
    order = np.argsort(start_array, kind="mergesort")
    start_array = start_array[order]
    end_array = end_array[order]
    merged_starts = [float(start_array[0])]
    merged_ends = [float(end_array[0])]
    for start, end in zip(start_array[1:], end_array[1:]):
        if float(start) <= merged_ends[-1]:
            merged_ends[-1] = max(merged_ends[-1], float(end))
        else:
            merged_starts.append(float(start))
            merged_ends.append(float(end))
    return (
        np.asarray(merged_starts, dtype=np.float64),
        np.asarray(merged_ends, dtype=np.float64),
    )


def _sample_nearest_mask(
    mask: np.ndarray | None,
    x_coordinates_um: np.ndarray | None,
    z_coordinates_um: np.ndarray | None,
    points_xz_um: np.ndarray,
) -> np.ndarray:
    points = _as_point_rows(points_xz_um, "points_xz_um")
    if mask is None or x_coordinates_um is None or z_coordinates_um is None or points.shape[0] == 0:
        return np.zeros(points.shape[0], dtype=bool)
    values = np.asarray(mask, dtype=bool)
    x_um = np.asarray(x_coordinates_um, dtype=np.float64)
    z_um = np.asarray(z_coordinates_um, dtype=np.float64)
    ix, valid_x = _nearest_axis_indices(x_um, points[:, 0])
    iz, valid_z = _nearest_axis_indices(z_um, points[:, 1])
    result = np.zeros(points.shape[0], dtype=bool)
    valid = valid_x & valid_z
    result[valid] = values[ix[valid], iz[valid]]
    return result


def _nearest_axis_indices(axis: np.ndarray, coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    insertion = np.searchsorted(axis, coordinates, side="left")
    upper = np.clip(insertion, 0, axis.size - 1)
    lower = np.clip(insertion - 1, 0, axis.size - 1)
    choose_upper = np.abs(axis[upper] - coordinates) < np.abs(coordinates - axis[lower])
    indices = np.where(choose_upper, upper, lower).astype(np.int64)
    lower_edge = axis[0] - 0.5 * (axis[1] - axis[0])
    upper_edge = axis[-1] + 0.5 * (axis[-1] - axis[-2])
    valid = np.isfinite(coordinates) & (coordinates >= lower_edge) & (coordinates <= upper_edge)
    return indices, valid


def _as_point_rows(values: np.ndarray, name: str) -> np.ndarray:
    points = np.asarray(values, dtype=np.float64)
    if points.ndim == 1:
        if points.shape != (2,):
            raise ValueError(f"{name} must have shape (2,) or (N, 2).")
        points = points.reshape(1, 2)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError(f"{name} must have shape (2,) or (N, 2).")
    if np.any(~np.isfinite(points)):
        raise ValueError(f"{name} must contain only finite values.")
    return points


def _broadcast_point_rows(values: np.ndarray, count: int, name: str) -> np.ndarray:
    rows = _as_point_rows(values, name)
    if rows.shape[0] == 1 and count != 1:
        return np.broadcast_to(rows, (count, 2)).copy()
    if rows.shape[0] != count:
        raise ValueError(f"{name} must contain one vector or exactly {count} vectors.")
    return rows


def _broadcast_scalars(values: np.ndarray | float, count: int, name: str) -> np.ndarray:
    scalars = np.asarray(values, dtype=np.float64)
    if scalars.ndim == 0:
        return np.full(count, float(scalars), dtype=np.float64)
    if scalars.ndim != 1 or scalars.size not in (1, count):
        raise ValueError(f"{name} must be scalar or have shape ({count},).")
    if scalars.size == 1 and count != 1:
        return np.full(count, float(scalars[0]), dtype=np.float64)
    return scalars.copy()
