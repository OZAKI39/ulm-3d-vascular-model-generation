"""Revised-v20 trajectory-carried categorical vessel ownership.

The regular-grid ``vessel_id`` field is intentionally absent from this module.
An active particle owns one directed vessel segment from birth until its
continuous centre path crosses a parent/child commitment section.  The section
catalogue is constructed from the same grid-independent transition scale used
by the Revised-v16 continuous lumen.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import math
from typing import Iterable

import numpy as np
from shapely.geometry import LineString, Point

try:
    from numba import njit
except ImportError:  # pragma: no cover - production installs include Numba
    njit = None

from ulm_vascular_model_generator.utils.core.models import Vessel


_TOPOLOGY_SCHEMA = "revised_v20_continuous_commitment_sections_v1"


@dataclass(frozen=True, slots=True)
class TopologicalCommitmentCatalog:
    """Immutable directed sections and dense adjacency used by particle steps.

    All vessel IDs in this object are the public trajectory convention:
    ``actual Vessel.vid + 1``.  Zero is therefore never a valid active owner.
    """

    schema: str
    root_vessel_id: np.ndarray
    terminal_vessel_id: np.ndarray
    parent_vessel_id: np.ndarray
    child_vessel_id: np.ndarray
    point_xz_um: np.ndarray
    downstream_normal_xz: np.ndarray
    tangent_xz: np.ndarray
    half_width_um: np.ndarray
    transition_end_distance_um: np.ndarray
    commitment_distance_um: np.ndarray
    child_section_by_vessel_id: np.ndarray
    child_section_offsets_by_parent_id: np.ndarray
    child_section_indices_by_parent_id: np.ndarray
    geometry_hash_sha256: str

    @property
    def section_count(self) -> int:
        return int(self.parent_vessel_id.size)

    def to_metadata(self) -> dict[str, object]:
        """Return compact scalar metadata for a trajectory result."""

        return {
            "particle_vessel_ownership_model": self.schema,
            "particle_vessel_id_convention": "one_based_persistent_topological_state",
            "particle_vessel_id_raster_authority": False,
            "topological_commitment_section_count": self.section_count,
            "topological_commitment_geometry_hash_sha256": self.geometry_hash_sha256,
        }


@dataclass(frozen=True, slots=True)
class TopologicalCrossingBatch:
    """Earliest valid ownership transition on each accepted centre chord."""

    fraction: np.ndarray
    new_vessel_id: np.ndarray
    section_index: np.ndarray
    position_xz_um: np.ndarray

    @property
    def earliest_fraction(self) -> float | None:
        finite = self.fraction[np.isfinite(self.fraction)]
        return float(np.min(finite)) if finite.size else None


def build_topological_commitment_catalog(
    vessels: Iterable[Vessel],
    continuous_geometry: object,
) -> TopologicalCommitmentCatalog:
    """Build one grid-independent directed section for every non-root segment."""

    ordered = sorted(list(vessels), key=lambda vessel: int(vessel.vid))
    if not ordered:
        raise ValueError("Revised-v20 ownership requires at least one vessel.")
    by_id = {int(vessel.vid): vessel for vessel in ordered}
    if len(by_id) != len(ordered) or min(by_id) < 0:
        raise ValueError("Revised-v20 ownership requires unique non-negative vessel IDs.")

    positive_radii = [float(vessel.radius) for vessel in ordered if float(vessel.radius) > 0.0]
    if len(positive_radii) != len(ordered):
        raise ValueError("Every Revised-v20 vessel radius must be positive.")
    fillet_radius = 0.05 * min(positive_radii)
    endpoint_tolerance = _endpoint_tolerance_um(ordered)

    roots = [int(vessel.vid) + 1 for vessel in ordered if int(vessel.parent_id) < 0]
    terminals = [int(vessel.vid) + 1 for vessel in ordered if len(vessel.children) == 0]
    if not roots or not terminals:
        raise ValueError("Revised-v20 ownership requires root and terminal vessels.")

    rows: dict[int, dict[str, object]] = {}
    for child in ordered:
        child_id = int(child.vid)
        parent_id = int(child.parent_id)
        if parent_id < 0:
            continue
        if parent_id not in by_id:
            raise ValueError(f"Vessel {child_id} references missing parent {parent_id}.")
        parent_children = {int(value) for value in by_id[parent_id].children}
        if child_id not in parent_children:
            raise ValueError(
                f"Vessel {child_id} and parent {parent_id} disagree about their directed edge."
            )

        proximal = _xz(child.x_p)
        distal = _xz(child.x_d)
        delta = distal - proximal
        length = float(np.linalg.norm(delta))
        if not math.isfinite(length) or length <= endpoint_tolerance:
            raise ValueError(f"Vessel {child_id} has no positive directed length.")
        normal = delta / length
        tangent = np.asarray([-normal[1], normal[0]], dtype=np.float64)
        incident = [
            vessel
            for vessel in ordered
            if _same_endpoint(_xz(vessel.x_p), proximal, endpoint_tolerance)
            or _same_endpoint(_xz(vessel.x_d), proximal, endpoint_tolerance)
        ]
        junction_radius = max(float(vessel.radius) for vessel in incident)
        transition_end = _transition_length_um(
            length,
            junction_radius,
            float(child.radius),
        )
        base_distance = min(transition_end + fillet_radius, 0.49 * length)
        if base_distance <= 0.0 or base_distance >= length:
            raise ValueError(
                f"Vessel {child_id} has no interior location for a commitment section."
            )
        rows[child_id] = {
            "parent_id": parent_id,
            "proximal": proximal,
            "normal": normal,
            "tangent": tangent,
            "half_width": float(child.radius),
            "length": length,
            "transition_end": transition_end,
            "distance": base_distance,
            "maximum_distance": 0.49 * length,
        }

    # A common-tangent junction may still make sibling cross-sections touch.
    # Move only that sibling group by the smallest shared deterministic fraction
    # needed to separate the sections by the v16 physical fillet scale.
    for parent in ordered:
        child_ids = sorted(int(value) for value in parent.children)
        if len(child_ids) < 2:
            continue
        missing = [value for value in child_ids if value not in rows]
        if missing:
            raise ValueError(
                f"Parent vessel {int(parent.vid)} has unknown children {missing}."
            )
        if _minimum_sibling_distance(child_ids, rows, 0.0) >= fillet_radius:
            continue
        if _minimum_sibling_distance(child_ids, rows, 1.0) < fillet_radius:
            raise ValueError(
                f"Children {child_ids} never form independent commitment sections."
            )
        lower = 0.0
        upper = 1.0
        for _ in range(64):
            middle = 0.5 * (lower + upper)
            if _minimum_sibling_distance(child_ids, rows, middle) >= fillet_radius:
                upper = middle
            else:
                lower = middle
        for child_id in child_ids:
            row = rows[child_id]
            base = float(row["distance"])
            maximum = float(row["maximum_distance"])
            row["distance"] = base + upper * (maximum - base)

    child_ids = sorted(rows)
    parent_vessel_id = np.asarray(
        [int(rows[child_id]["parent_id"]) + 1 for child_id in child_ids],
        dtype=np.int32,
    )
    child_vessel_id = np.asarray([child_id + 1 for child_id in child_ids], dtype=np.int32)
    normals = np.ascontiguousarray(
        np.vstack([np.asarray(rows[child_id]["normal"], dtype=np.float64) for child_id in child_ids])
        if child_ids
        else np.empty((0, 2), dtype=np.float64)
    )
    tangents = np.ascontiguousarray(
        np.vstack([np.asarray(rows[child_id]["tangent"], dtype=np.float64) for child_id in child_ids])
        if child_ids
        else np.empty((0, 2), dtype=np.float64)
    )
    distances = np.asarray([float(rows[child_id]["distance"]) for child_id in child_ids])
    transition_ends = np.asarray(
        [float(rows[child_id]["transition_end"]) for child_id in child_ids],
        dtype=np.float64,
    )
    points = np.ascontiguousarray(
        np.vstack(
            [
                np.asarray(rows[child_id]["proximal"], dtype=np.float64)
                + distances[index] * normals[index]
                for index, child_id in enumerate(child_ids)
            ]
        )
        if child_ids
        else np.empty((0, 2), dtype=np.float64)
    )
    half_widths = np.asarray(
        [float(rows[child_id]["half_width"]) for child_id in child_ids],
        dtype=np.float64,
    )

    _validate_continuous_sections(
        continuous_geometry,
        child_vessel_id,
        points,
        tangents,
        half_widths,
    )
    _validate_nonintersecting_sections(child_vessel_id, points, tangents, half_widths)

    maximum_owner = max(by_id) + 1
    child_section = np.full(maximum_owner + 1, -1, dtype=np.int32)
    for section_index, owner in enumerate(child_vessel_id):
        child_section[int(owner)] = int(section_index)
    grouped: list[list[int]] = [[] for _ in range(maximum_owner + 1)]
    for section_index, owner in enumerate(parent_vessel_id):
        grouped[int(owner)].append(section_index)
    offsets = np.zeros(maximum_owner + 2, dtype=np.int64)
    for owner in range(maximum_owner + 1):
        offsets[owner + 1] = offsets[owner] + len(grouped[owner])
    indices = np.empty(int(offsets[-1]), dtype=np.int32)
    for owner in range(maximum_owner + 1):
        indices[offsets[owner] : offsets[owner + 1]] = grouped[owner]

    digest = sha256()
    for values in (
        parent_vessel_id,
        child_vessel_id,
        points,
        normals,
        tangents,
        half_widths,
        transition_ends,
        distances,
    ):
        array = np.ascontiguousarray(values)
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())

    return TopologicalCommitmentCatalog(
        schema=_TOPOLOGY_SCHEMA,
        root_vessel_id=np.asarray(roots, dtype=np.int32),
        terminal_vessel_id=np.asarray(terminals, dtype=np.int32),
        parent_vessel_id=parent_vessel_id,
        child_vessel_id=child_vessel_id,
        point_xz_um=points,
        downstream_normal_xz=normals,
        tangent_xz=tangents,
        half_width_um=half_widths,
        transition_end_distance_um=transition_ends,
        commitment_distance_um=distances,
        child_section_by_vessel_id=np.ascontiguousarray(child_section),
        child_section_offsets_by_parent_id=np.ascontiguousarray(offsets),
        child_section_indices_by_parent_id=np.ascontiguousarray(indices),
        geometry_hash_sha256=digest.hexdigest(),
    )


def inspect_topological_crossings(
    start_position_xz_um: np.ndarray,
    end_position_xz_um: np.ndarray,
    current_vessel_id: np.ndarray,
    active: np.ndarray,
    catalog: TopologicalCommitmentCatalog,
    *,
    use_numba: bool,
) -> TopologicalCrossingBatch:
    """Find the first directed parent/child transition on each centre chord."""

    start = np.ascontiguousarray(start_position_xz_um, dtype=np.float64)
    end = np.ascontiguousarray(end_position_xz_um, dtype=np.float64)
    owner = np.ascontiguousarray(current_vessel_id, dtype=np.int32)
    live = np.ascontiguousarray(active, dtype=np.bool_)
    count = int(owner.size)
    if start.shape != (count, 2) or end.shape != (count, 2) or live.shape != (count,):
        raise ValueError("Topological crossing arrays have inconsistent shapes.")
    invalid = live & (
        (owner <= 0)
        | (owner >= catalog.child_section_offsets_by_parent_id.size - 1)
    )
    if np.any(invalid):
        raise ValueError("An active particle has no valid persistent topological vessel ID.")

    kernel = (
        _inspect_topological_crossings_numba
        if use_numba and njit is not None
        else _inspect_topological_crossings_kernel
    )
    fraction, new_owner, section = kernel(
        start,
        end,
        owner,
        live,
        catalog.parent_vessel_id,
        catalog.child_vessel_id,
        catalog.point_xz_um,
        catalog.downstream_normal_xz,
        catalog.tangent_xz,
        catalog.half_width_um,
        catalog.child_section_by_vessel_id,
        catalog.child_section_offsets_by_parent_id,
        catalog.child_section_indices_by_parent_id,
    )
    position = np.full((count, 2), np.nan, dtype=np.float64)
    finite = np.isfinite(fraction)
    if np.any(finite):
        position[finite] = start[finite] + fraction[finite, None] * (
            end[finite] - start[finite]
        )
    return TopologicalCrossingBatch(
        fraction=np.asarray(fraction, dtype=np.float64),
        new_vessel_id=np.asarray(new_owner, dtype=np.int32),
        section_index=np.asarray(section, dtype=np.int32),
        position_xz_um=position,
    )


def _inspect_topological_crossings_kernel(
    start,
    end,
    owner,
    active,
    parent_id,
    child_id,
    point,
    normal,
    tangent,
    half_width,
    child_section,
    child_offsets,
    child_indices,
):
    count = owner.size
    fractions = np.full(count, np.nan, dtype=np.float64)
    new_owner = owner.copy()
    sections = np.full(count, -1, dtype=np.int32)
    machine = np.finfo(np.float64).eps
    for lane in range(count):
        if not active[lane]:
            continue
        current = int(owner[lane])
        dx = end[lane, 0] - start[lane, 0]
        dz = end[lane, 1] - start[lane, 1]
        coordinate_scale = max(
            abs(start[lane, 0]),
            abs(start[lane, 1]),
            abs(end[lane, 0]),
            abs(end[lane, 1]),
            1.0,
        )
        epsilon = 2048.0 * machine * coordinate_scale
        best_fraction = 2.0
        best_section = -1
        best_owner = current

        reverse_section = int(child_section[current])
        if reverse_section >= 0:
            candidate = _directed_crossing_fraction(
                start[lane, 0],
                start[lane, 1],
                dx,
                dz,
                reverse_section,
                -1.0,
                point,
                normal,
                tangent,
                half_width,
                epsilon,
            )
            if candidate < best_fraction:
                best_fraction = candidate
                best_section = reverse_section
                best_owner = int(parent_id[reverse_section])

        first = int(child_offsets[current])
        last = int(child_offsets[current + 1])
        for item in range(first, last):
            section_index = int(child_indices[item])
            candidate = _directed_crossing_fraction(
                start[lane, 0],
                start[lane, 1],
                dx,
                dz,
                section_index,
                1.0,
                point,
                normal,
                tangent,
                half_width,
                epsilon,
            )
            if candidate < best_fraction:
                best_fraction = candidate
                best_section = section_index
                best_owner = int(child_id[section_index])

        if best_section >= 0 and best_fraction <= 1.0 + 8.0 * machine:
            fractions[lane] = min(max(best_fraction, 0.0), 1.0)
            new_owner[lane] = best_owner
            sections[lane] = best_section
    return fractions, new_owner, sections


def _directed_crossing_fraction(
    start_x,
    start_z,
    delta_x,
    delta_z,
    section,
    direction_sign,
    point,
    normal,
    tangent,
    half_width,
    epsilon,
):
    relative_x = start_x - point[section, 0]
    relative_z = start_z - point[section, 1]
    phi_start = relative_x * normal[section, 0] + relative_z * normal[section, 1]
    phi_delta = delta_x * normal[section, 0] + delta_z * normal[section, 1]
    directed_start = direction_sign * phi_start
    directed_delta = direction_sign * phi_delta
    directed_end = directed_start + directed_delta
    if directed_delta <= epsilon:
        return 2.0
    if directed_start > epsilon or directed_end < -epsilon:
        return 2.0
    fraction = 0.0 if abs(directed_start) <= epsilon else -directed_start / directed_delta
    if fraction < -epsilon or fraction > 1.0 + epsilon:
        return 2.0
    intersection_x = start_x + fraction * delta_x - point[section, 0]
    intersection_z = start_z + fraction * delta_z - point[section, 1]
    lateral = abs(
        intersection_x * tangent[section, 0]
        + intersection_z * tangent[section, 1]
    )
    if lateral > half_width[section] + epsilon:
        return 2.0
    return min(max(fraction, 0.0), 1.0)


if njit is not None:  # pragma: no branch
    _directed_crossing_fraction = njit(cache=True)(_directed_crossing_fraction)
    _inspect_topological_crossings_numba = njit(cache=True)(
        _inspect_topological_crossings_kernel
    )
else:  # pragma: no cover
    _inspect_topological_crossings_numba = _inspect_topological_crossings_kernel


def _transition_length_um(
    vessel_length_um: float,
    junction_radius_um: float,
    branch_radius_um: float,
) -> float:
    """Return the exact v16 common-tangent transition scale."""

    return min(
        0.45 * float(vessel_length_um),
        max(
            2.0 * float(junction_radius_um),
            float(junction_radius_um) + float(branch_radius_um),
            4.0 * abs(float(junction_radius_um) - float(branch_radius_um)),
        ),
    )


def _minimum_sibling_distance(
    child_ids: list[int], rows: dict[int, dict[str, object]], fraction: float
) -> float:
    lines = [_section_line(rows[child_id], fraction) for child_id in child_ids]
    minimum = float("inf")
    for first in range(len(lines)):
        for second in range(first + 1, len(lines)):
            minimum = min(minimum, float(lines[first].distance(lines[second])))
    return minimum


def _section_line(row: dict[str, object], fraction: float = 0.0) -> LineString:
    base = float(row["distance"])
    maximum = float(row["maximum_distance"])
    distance = base + float(fraction) * (maximum - base)
    point = np.asarray(row["proximal"], dtype=np.float64) + distance * np.asarray(
        row["normal"], dtype=np.float64
    )
    tangent = np.asarray(row["tangent"], dtype=np.float64)
    width = float(row["half_width"])
    return LineString((point - width * tangent, point + width * tangent))


def _validate_continuous_sections(
    continuous_geometry: object,
    child_ids: np.ndarray,
    points: np.ndarray,
    tangents: np.ndarray,
    half_widths: np.ndarray,
) -> None:
    if child_ids.size == 0:
        return
    polygon = getattr(continuous_geometry, "lumen_polygon", None)
    if polygon is None:
        raise ValueError("Revised-v20 ownership requires the continuous lumen polygon.")
    for index, child_id in enumerate(child_ids):
        point = points[index]
        if not bool(polygon.covers(Point(tuple(point)))):
            raise ValueError(
                f"Commitment section for vessel {int(child_id)} is outside the continuous lumen."
            )
        # Trim only the floating-point boundary endpoints; the stored physical
        # aperture remains the exact branch radius used by crossing tests.
        trimmed = half_widths[index] * (1.0 - 1.0e-12)
        line = LineString(
            (
                point - trimmed * tangents[index],
                point + trimmed * tangents[index],
            )
        )
        if not bool(polygon.covers(line)):
            raise ValueError(
                f"Commitment aperture for vessel {int(child_id)} is not contained in the continuous lumen."
            )


def _validate_nonintersecting_sections(
    child_ids: np.ndarray,
    points: np.ndarray,
    tangents: np.ndarray,
    half_widths: np.ndarray,
) -> None:
    lines = [
        LineString(
            (
                points[index] - half_widths[index] * tangents[index],
                points[index] + half_widths[index] * tangents[index],
            )
        )
        for index in range(child_ids.size)
    ]
    for first in range(len(lines)):
        for second in range(first + 1, len(lines)):
            if lines[first].intersects(lines[second]):
                raise ValueError(
                    "Revised-v20 commitment sections intersect for child vessels "
                    f"{int(child_ids[first])} and {int(child_ids[second])}."
                )


def _endpoint_tolerance_um(vessels: list[Vessel]) -> float:
    scale = max(
        max(abs(float(value)) for vessel in vessels for point in (vessel.x_p, vessel.x_d) for value in (point[0], point[2])),
        1.0,
    )
    return 1024.0 * np.finfo(np.float64).eps * scale


def _same_endpoint(first: np.ndarray, second: np.ndarray, tolerance_um: float) -> bool:
    return bool(np.linalg.norm(first - second) <= float(tolerance_um))


def _xz(point: object) -> np.ndarray:
    values = np.asarray(point, dtype=np.float64)
    return np.asarray([values[0], values[2]], dtype=np.float64)
