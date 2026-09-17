"""Planar and volumetric domain and vessel-tree spatial queries."""

from __future__ import annotations

import math
from typing import Iterable, Optional

import numpy as np

from ..config.yaml_config import RegionMask
from ..geometry.geometry import _angle_deg
from .models import Vessel
from .tree_index import (
    SegmentIndex,
    build_segment_index,
    point_segment_distances_3d,
    segment_distances_3d,
    segment_intersections_3d,
)


class TreeQueryMixin:
    def _invalidate_tree_cache(self) -> None:
        self._segment_index = None
        self._segment_index_size = -1

    def _segments(self) -> SegmentIndex:
        cached = getattr(self, "_segment_index", None)
        if cached is None or getattr(self, "_segment_index_size", -1) != len(
            self.vessels
        ):
            cached = build_segment_index(self.vessels)
            self._segment_index = cached
            self._segment_index_size = len(self.vessels)
        return cached

    def _terminal_count(self) -> int:
        return sum(1 for vessel in self.vessels if not vessel.children)

    # ----- Domain helpers ---------------------------------------------

    def _domain_lengths(self) -> np.ndarray:
        if self.cfg.geometry_mode == "planar_2d":
            return np.asarray(
                [self.cfg.cube_um, 0.0, self.cfg.cube_um],
                dtype=float,
            )
        return np.asarray(
            [
                self.cfg.domain_size_x_um,
                self.cfg.domain_size_y_um,
                self.cfg.domain_size_z_um,
            ],
            dtype=float,
        )

    def _domain_volume(self) -> float:
        lengths = self._domain_lengths()
        if self.cfg.geometry_mode == "planar_2d":
            return float(lengths[0] * lengths[2])
        return float(np.prod(lengths))

    def _lc(self) -> float:
        if self.cfg.geometry_mode == "planar_2d":
            return math.sqrt(self._domain_volume() / math.pi)
        return (3.0 * self._domain_volume() / (4.0 * math.pi)) ** (1.0 / 3.0)

    def _project(self, p: Iterable[float]) -> np.ndarray:
        arr = np.asarray(p, dtype=float).reshape(-1)
        if arr.size < 3:
            raise ValueError("A vessel point must contain x, y, and z coordinates.")
        if self.cfg.geometry_mode == "planar_2d":
            return np.asarray(
                [float(arr[0]), self.cfg.plane_y_um, float(arr[2])],
                dtype=float,
            )
        return np.asarray(arr[:3], dtype=float)

    def _random_point(self) -> np.ndarray:
        if self.cfg.terminal_distribution == "gaussian":
            x = self.rng.normal(
                self.cfg.gaussian_mean_x_um, self.cfg.gaussian_sigma_x_um
            )
            y = self.rng.normal(
                self.cfg.gaussian_mean_y_um, self.cfg.gaussian_sigma_y_um
            )
            z = self.rng.normal(
                self.cfg.gaussian_mean_z_um, self.cfg.gaussian_sigma_z_um
            )
        elif self.cfg.terminal_distribution == "gaussian_mixture":
            first = self.rng.random() < max(
                0.0, min(1.0, self.cfg.gaussian_mixture_weight)
            )
            if first:
                x = self.rng.normal(
                    self.cfg.gaussian_mean_x_um, self.cfg.gaussian_sigma_x_um
                )
                y = self.rng.normal(
                    self.cfg.gaussian_mean_y_um, self.cfg.gaussian_sigma_y_um
                )
                z = self.rng.normal(
                    self.cfg.gaussian_mean_z_um, self.cfg.gaussian_sigma_z_um
                )
            else:
                x = self.rng.normal(
                    self.cfg.gaussian2_mean_x_um, self.cfg.gaussian2_sigma_x_um
                )
                y = self.rng.normal(
                    self.cfg.gaussian2_mean_y_um, self.cfg.gaussian2_sigma_y_um
                )
                z = self.rng.normal(
                    self.cfg.gaussian2_mean_z_um, self.cfg.gaussian2_sigma_z_um
                )
        else:
            lengths = self._domain_lengths()
            x = self.rng.uniform(0.0, lengths[0])
            y = self.rng.uniform(0.0, lengths[1]) if lengths[1] > 0.0 else 0.0
            z = self.rng.uniform(0.0, lengths[2])

        if self.cfg.geometry_mode == "planar_2d":
            y = self.cfg.plane_y_um
        return np.asarray([x, y, z], dtype=float)

    def _in_bounds(self, p: np.ndarray, margin: float = 0.0) -> bool:
        lengths = self._domain_lengths()
        if self.cfg.geometry_mode == "planar_2d":
            return bool(
                margin <= p[0] <= lengths[0] - margin
                and margin <= p[2] <= lengths[2] - margin
            )
        return bool(
            np.all(np.asarray(p, dtype=float) >= margin)
            and np.all(np.asarray(p, dtype=float) <= lengths - margin)
        )

    def _inside_any(
        self, p: np.ndarray, regions: tuple[RegionMask, ...]
    ) -> bool:
        return any(region.contains(p) for region in regions)

    def _in_extended_domain(self, p: np.ndarray) -> bool:
        return self._in_bounds(p) and not self._inside_any(
            p, self.cfg.avascular_regions
        )

    def _in_perfusion_domain(self, p: np.ndarray) -> bool:
        return self._in_extended_domain(p) and not self._inside_any(
            p, self.cfg.carriage_regions
        )

    def _segment_domain_ok(
        self, a: np.ndarray, b: np.ndarray, role: str
    ) -> bool:
        checker = (
            self._in_extended_domain
            if role in {"transport", "perforator"}
            else self._in_perfusion_domain
        )
        length = float(np.linalg.norm(np.asarray(b) - np.asarray(a)))
        samples = max(9, int(math.ceil(length / 50.0)) + 1)
        for t in np.linspace(0.0, 1.0, samples):
            if not checker((1.0 - t) * a + t * b):
                return False
        return True

    def _bifurcation_domain_ok(self, vessel: Vessel, point: np.ndarray) -> bool:
        if vessel.role == "transport":
            return self._in_extended_domain(point)
        return self._in_perfusion_domain(point)

    def _lmin(self) -> float:
        terminal_count = max(self._terminal_count(), 1)
        if self.cfg.geometry_mode == "planar_2d":
            scale = math.sqrt(self.cfg.nu / (terminal_count + 1))
        else:
            scale = (self.cfg.nu / (terminal_count + 1)) ** (1.0 / 3.0)
        return self._lc() * scale * self._fr_runtime

    # ----- Tree queries ------------------------------------------------

    def _closest_distance(self, p: np.ndarray) -> float:
        if not self.vessels:
            return math.inf
        segments = self._segments()
        distances = point_segment_distances_3d(
            self._project(p), segments.starts, segments.ends
        )
        return float(np.min(distances)) if distances.size else math.inf

    def _branch_source_mask(self, segments: SegmentIndex) -> np.ndarray:
        mask = segments.branchable.copy()
        if self.cfg.main_trunk_only_branching:
            mask &= segments.main_trunk
        return mask

    def _neighborhood(self, p: np.ndarray) -> list[int]:
        radius = float(self.cfg.fn_neighborhood) * self._lc()
        point = self._project(p)
        segments = self._segments()
        branch_source = self._branch_source_mask(segments)
        if segments.ids.size == 0 or not np.any(branch_source):
            return []

        if segments.midpoint_tree is None:
            candidate_idx = np.flatnonzero(branch_source)
        else:
            search_radius = radius + segments.max_half_length
            candidate_idx = np.asarray(
                segments.midpoint_tree.query_ball_point(point, search_radius),
                dtype=int,
            )
            candidate_idx = (
                candidate_idx[branch_source[candidate_idx]]
                if candidate_idx.size
                else np.empty(0, dtype=int)
            )

        ids = []
        if candidate_idx.size:
            distances = point_segment_distances_3d(
                point,
                segments.starts[candidate_idx],
                segments.ends[candidate_idx],
            )
            ids = segments.ids[
                candidate_idx[distances <= radius]
            ].astype(int).tolist()

        if not ids:
            branchable_idx = np.flatnonzero(branch_source)
            k = min(8, branchable_idx.size)
            midpoint_distances = np.linalg.norm(
                segments.midpoints[branchable_idx] - point, axis=1
            )
            nearest_local = np.argpartition(midpoint_distances, k - 1)[:k]
            nearest_idx = branchable_idx[nearest_local]
            order = np.argsort(midpoint_distances[nearest_local])
            ids = segments.ids[nearest_idx[order]].astype(int).tolist()
        return ids

    def _segment_intersects_tree(
        self, a: np.ndarray, b: np.ndarray, exclude: set[int]
    ) -> bool:
        segments = self._segments()
        if segments.ids.size == 0:
            return False
        mask = ~np.isin(segments.ids, list(exclude)) if exclude else np.ones(
            segments.ids.size, dtype=bool
        )
        if not np.any(mask):
            return False
        return bool(
            np.any(
                segment_intersections_3d(
                    a, b, segments.starts[mask], segments.ends[mask]
                )
            )
        )

    def _segment_too_close_to_tree(
        self,
        a: np.ndarray,
        b: np.ndarray,
        exclude: set[int],
        candidate_radius: float | None = None,
    ) -> bool:
        surface_clearance = float(self.cfg.min_vessel_clearance_um)
        if surface_clearance <= 0.0:
            return False
        radius = (
            self.haemodynamics.minimum_size_um
            if candidate_radius is None
            else float(candidate_radius)
        )
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        trimmed_a = a + 0.03 * (b - a)
        trimmed_b = a + 0.97 * (b - a)
        segments = self._segments()
        mask = ~np.isin(segments.ids, list(exclude)) if exclude else np.ones(
            segments.ids.size, dtype=bool
        )
        if not np.any(mask):
            return False
        distances = segment_distances_3d(
            trimmed_a, trimmed_b, segments.starts[mask], segments.ends[mask]
        )
        existing_radii = np.asarray(
            [self.vessels[int(vid)].radius for vid in segments.ids[mask]],
            dtype=float,
        )
        required = radius + existing_radii + surface_clearance
        return bool(np.any(distances < required))

    def _junction_geometry_ok(
        self, vessel_ids: Optional[Iterable[int]] = None
    ) -> bool:
        ids = vessel_ids if vessel_ids is not None else range(len(self.vessels))
        for vid in ids:
            if vid < 0 or vid >= len(self.vessels):
                continue
            vessel = self.vessels[vid]
            if not vessel.children:
                continue
            child_dirs = [
                self.vessels[child].x_d - self.vessels[child].x_p
                for child in vessel.children
            ]
            for i in range(len(child_dirs)):
                for j in range(i + 1, len(child_dirs)):
                    opening = _angle_deg(child_dirs[i], child_dirs[j])
                    if not (
                        self.cfg.min_daughter_angle_deg
                        <= opening
                        <= self.cfg.max_opening_angle_deg
                    ):
                        return False
            parent_dir = vessel.x_d - vessel.x_p
            if any(
                _angle_deg(parent_dir, child_dir)
                > self.cfg.max_daughter_bend_deg
                for child_dir in child_dirs
            ):
                return False
        return True
