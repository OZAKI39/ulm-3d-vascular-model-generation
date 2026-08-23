"""v6 C1 port conditioning and continuous-field junction transitions."""

from __future__ import annotations

import time
from dataclasses import replace
from typing import Any

import numpy as np
import trimesh
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .local_implicit_junction import (
    _clean_local_mesh,
    _dense_local_samples,
    define_junction_collars,
    sample_from_junction,
)
from .lumen_builder import balanced_manifold_union, build_variable_radius_tube
from .surface_qc import _section_polygon
from .surface_transition import (
    CFD_DERIVED_EXTENSION,
    CORE_BOUNDARY_EXTENSION_ORIGIN,
    EXPLICIT_BRANCH_FACE,
    JUNCTION_CORE_FACE,
    SOURCE_CENTERLINE,
    TRANSITION_COLLAR_FACE,
    _boundary_loops,
    _collapse_boundary_loops_to_count,
    _loop_at_plane,
    _phase_align,
    _resample_loop,
    _zipper_faces,
    build_open_explicit_tubes,
    clip_implicit_patch_to_transition_inner,
    combine_stitched_surfaces,
    trim_branches_to_transition_outer,
)
from .types import (
    BranchGeometry,
    GeometryValidationError,
    HybridBuildDetails,
    JunctionCollar,
    LocalJunctionPatch,
    PortGeometry,
)


PURE_BRANCH_FACE = np.uint8(4)


def _cleanup_continuous_patch(
    mesh: trimesh.Trimesh, tolerance_um: float
) -> trimesh.Trimesh:
    """Collapse only nanometre-scale MC edges on the local field patch."""

    if tolerance_um <= 0.0:
        return _clean_local_mesh(mesh)
    output = mesh.copy()
    for _ in range(4):
        vertices = np.asarray(output.vertices, dtype=float)
        faces = np.asarray(output.faces, dtype=np.int64)
        edges = np.asarray(output.edges_unique, dtype=np.int64)
        lengths = np.linalg.norm(vertices[edges[:, 1]] - vertices[edges[:, 0]], axis=1)
        short = edges[lengths < float(tolerance_um)]
        if not len(short):
            break
        parent = np.arange(len(vertices), dtype=np.int64)

        def find(value: int) -> int:
            while parent[value] != value:
                parent[value] = parent[parent[value]]
                value = int(parent[value])
            return value

        for first, second in short:
            root_first = find(int(first))
            root_second = find(int(second))
            if root_first != root_second:
                parent[root_second] = root_first
        roots = np.asarray([find(index) for index in range(len(vertices))])
        _, inverse = np.unique(roots, return_inverse=True)
        merged = np.zeros((int(inverse.max()) + 1, 3), dtype=float)
        counts = np.bincount(inverse)
        np.add.at(merged, inverse, vertices)
        merged /= counts[:, None]
        remapped = inverse[faces]
        valid = (
            (remapped[:, 0] != remapped[:, 1])
            & (remapped[:, 1] != remapped[:, 2])
            & (remapped[:, 2] != remapped[:, 0])
        )
        remapped = remapped[valid]
        _, unique_indices = np.unique(
            np.sort(remapped, axis=1), axis=0, return_index=True
        )
        output = trimesh.Trimesh(
            vertices=merged,
            faces=remapped[np.sort(unique_indices)],
            process=False,
        )
        output.remove_unreferenced_vertices()
    return _clean_local_mesh(output)


def _angle_deg(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    cosine = float(
        np.clip(
            np.dot(first, second)
            / max(float(np.linalg.norm(first) * np.linalg.norm(second)), 1.0e-15),
            -1.0,
            1.0,
        )
    )
    return float(np.degrees(np.arccos(cosine)))


def _port_source_samples(
    branch: BranchGeometry,
    endpoint_index: int,
    count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    count = min(max(int(count), 3), len(branch.points_um))
    if endpoint_index == -1:
        points = np.asarray(branch.points_um[-count:], dtype=float)
        radii = np.asarray(branch.radius_um[-count:], dtype=float)
    else:
        points = np.asarray(branch.points_um[:count][::-1], dtype=float)
        radii = np.asarray(branch.radius_um[:count][::-1], dtype=float)
    arc = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))))
    return points, radii, arc


def fit_port_source_state(
    branch: BranchGeometry,
    endpoint_index: int,
    port: PortGeometry,
    fit_points: int,
) -> dict[str, Any]:
    """Fit the source-side tangent and dr/ds from several real centerline points."""

    points, radii, arc = _port_source_samples(branch, endpoint_index, fit_points)
    # Endpoint-biased weighted PCA suppresses single-edge noise without
    # rotating a genuinely curved source branch toward its more distant chord.
    weights = np.geomspace(0.05, 1.0, len(points), dtype=float)
    center = np.average(points, axis=0, weights=weights)
    centered = points - center
    covariance = (centered * weights[:, None]).T @ centered / float(weights.sum())
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    tangent = np.asarray(eigenvectors[:, int(np.argmax(eigenvalues))], dtype=float)
    tangent /= np.linalg.norm(tangent)
    if float(np.dot(tangent, port.outward_tangent)) < 0.0:
        tangent = -tangent
    endpoint_tangent = np.asarray(points[-1] - points[-2], dtype=float)
    endpoint_tangent /= np.linalg.norm(endpoint_tangent)
    # Estimate the endpoint radius derivative from the final three real samples.
    # A quadratic endpoint derivative amplifies the small alternating radius noise
    # present at otherwise smooth ports.  The median of all pairwise secants is the
    # local Theil--Sen slope: it still uses a compliant three-point endpoint fit,
    # preserves genuine monotone taper, and rejects a single radius excursion.
    local_arc = arc[-3:]
    local_radii = radii[-3:]
    pairwise_slopes = np.asarray(
        [
            (local_radii[j] - local_radii[i])
            / max(float(local_arc[j] - local_arc[i]), 1.0e-15)
            for i in range(len(local_arc))
            for j in range(i + 1, len(local_arc))
        ],
        dtype=float,
    )
    slope = float(np.median(pairwise_slopes))
    intercept = float(radii[-1])
    return {
        "fit_point_count": int(len(points)),
        "fit_points_um": points.tolist(),
        "fit_radius_um": radii.tolist(),
        "fit_arc_um": arc.tolist(),
        "source_fitted_tangent": tangent,
        "source_endpoint_segment_tangent": endpoint_tangent,
        "source_radius_slope_um_per_um": float(slope),
        "source_radius_intercept_um": float(intercept),
        "source_radius_fit_method": "last_3_point_theil_sen",
        "source_radius_pairwise_slopes_um_per_um": pairwise_slopes.tolist(),
        "v5_tangent_jump_deg": _angle_deg(tangent, port.outward_tangent),
        "endpoint_segment_jump_deg": _angle_deg(
            endpoint_tangent, port.outward_tangent
        ),
        "v5_radius_slope_jump_um_per_um": abs(float(slope)),
    }


def _hermite_vector(
    fraction: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    start_derivative: np.ndarray,
    end_derivative: np.ndarray,
) -> np.ndarray:
    t = np.asarray(fraction, dtype=float)[:, None]
    h00 = 2.0 * t**3 - 3.0 * t**2 + 1.0
    h10 = t**3 - 2.0 * t**2 + t
    h01 = -2.0 * t**3 + 3.0 * t**2
    h11 = t**3 - t**2
    return h00 * start + h10 * start_derivative + h01 * end + h11 * end_derivative


def _hermite_scalar(
    fraction: np.ndarray,
    start: float,
    end: float,
    start_derivative: float,
    end_derivative: float,
) -> np.ndarray:
    return _hermite_vector(
        fraction,
        np.asarray([start]),
        np.asarray([end]),
        np.asarray([start_derivative]),
        np.asarray([end_derivative]),
    )[:, 0]


def conditioned_port_profile(
    branch: BranchGeometry,
    endpoint_index: int,
    port: PortGeometry,
    config: CFDLumenConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Create an extension-only Hermite tangent/radius conditioning profile."""

    state = fit_port_source_state(
        branch, endpoint_index, port, config.v6.port_source_fit_points
    )
    fitted = np.asarray(state["source_fitted_tangent"], dtype=float)
    target = np.asarray(port.outward_tangent, dtype=float)
    # The fitted tangent is the v6 diagnostic required by the protocol.  Do not
    # bend an already C1 join merely because a curved upstream fit differs from
    # its endpoint segment: conditioning is needed only when both estimates
    # confirm a direction discontinuity at the actual CUT_PORT interface.
    tangent_conditioning_applied = bool(
        state["v5_tangent_jump_deg"]
        >= config.v6.tangent_conditioning_threshold_deg
        and state["endpoint_segment_jump_deg"]
        >= config.v6.tangent_conditioning_threshold_deg
    )
    conditioned_start_tangent = fitted if tangent_conditioning_applied else target
    radius = float(port.radius_um)
    diameter = 2.0 * radius
    base_spacing = float(
        np.clip(
            config.geometry.resample_radius_fraction * radius,
            config.geometry.min_resample_spacing_um,
            config.geometry.max_resample_spacing_um,
        )
    )
    spacing = min(
        base_spacing,
        max(
            config.v6.port_extension_min_spacing_um,
            config.v6.port_extension_spacing_radius_fraction * radius,
        ),
    )
    count = max(2, int(np.ceil(port.extension_length_um / spacing)))
    distances = np.linspace(0.0, port.extension_length_um, count + 1)
    tangent_blend = min(
        config.v6.tangent_blend_diameters * diameter,
        port.extension_length_um,
    )
    points = port.original_position_um[None, :] + distances[:, None] * target[None, :]
    tangent_selected = distances <= tangent_blend + 1.0e-12
    if tangent_blend > 0.0 and np.any(tangent_selected):
        fraction = distances[tangent_selected] / tangent_blend
        blend_end = port.original_position_um + tangent_blend * target
        points[tangent_selected] = _hermite_vector(
            fraction,
            np.asarray(port.original_position_um, dtype=float),
            blend_end,
            tangent_blend * conditioned_start_tangent,
            tangent_blend * target,
        )

    slope = float(state["source_radius_slope_um_per_um"])
    radius_conditioning_applied = bool(
        abs(slope)
        >= config.v6.radius_slope_conditioning_threshold_um_per_um
    )
    radius_blend = min(
        config.v6.radius_blend_diameters * diameter,
        port.extension_length_um,
    )
    if radius_conditioning_applied:
        radius_blend = min(radius_blend, 0.5 * radius / abs(slope))
    else:
        radius_blend = 0.0
    # Return to the immutable CUT_PORT radius before the constant-radius cap
    # segment.  The cubic absorbs the measured source slope entirely inside
    # CFD_EXTENSION while preserving both the CUT_PORT and final port areas.
    radius_end = radius
    radii = np.full(len(distances), radius_end, dtype=float)
    radius_selected = distances <= radius_blend + 1.0e-12
    if radius_blend > 0.0 and np.any(radius_selected):
        fraction = distances[radius_selected] / radius_blend
        radii[radius_selected] = _hermite_scalar(
            fraction,
            radius,
            radius_end,
            radius_blend * slope,
            0.0,
        )
    if np.any(radii <= 0.0):
        raise GeometryValidationError(
            f"V6_PORT_RADIUS_CONDITIONING_FAILED: port {port.cut_port_id} has nonpositive radius"
        )
    initial_tangent = points[1] - points[0]
    initial_tangent /= np.linalg.norm(initial_tangent)
    endpoint_tangent = np.asarray(
        state["source_endpoint_segment_tangent"], dtype=float
    )
    initial_slope = float((radii[1] - radii[0]) / np.linalg.norm(points[1] - points[0]))
    state.update(
        {
            "port_id": port.port_id,
            "cut_port_id": port.cut_port_id,
            "branch_id": branch.branch_id,
            "endpoint_index": endpoint_index,
            "tangent_blend_length_um": tangent_blend,
            "tangent_conditioning_applied": tangent_conditioning_applied,
            "radius_blend_length_um": radius_blend,
            "radius_conditioning_applied": radius_conditioning_applied,
            "radius_slope_conditioning_threshold_um_per_um": (
                config.v6.radius_slope_conditioning_threshold_um_per_um
            ),
            "radius_after_conditioning_um": radius_end,
            "v6_tangent_jump_deg": _angle_deg(fitted, initial_tangent),
            "v5_interface_tangent_jump_deg": state[
                "endpoint_segment_jump_deg"
            ],
            "v6_interface_tangent_jump_deg": _angle_deg(
                endpoint_tangent, initial_tangent
            ),
            "v6_radius_slope_jump_um_per_um": abs(slope - initial_slope),
            "radius_profile_before": np.full(len(distances), radius).tolist(),
            "radius_profile_after": radii.tolist(),
            "distance_profile_um": distances.tolist(),
            "v5_extension_tangent": target.tolist(),
            "v6_initial_extension_tangent": initial_tangent.tolist(),
        }
    )
    return points, radii, distances, state


def extend_branches_to_cfd_ports_v6(
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
) -> tuple[list[BranchGeometry], list[dict[str, Any]], list[dict[str, Any]]]:
    """Condition only CFD-derived extension points; source samples remain byte-identical."""

    port_by_node = {int(port.local_node_id): port for port in ports}
    output: list[BranchGeometry] = []
    point_rows: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    for branch in branches:
        extended = replace(branch)
        points = np.asarray(branch.points_um, dtype=float).copy()
        radii = np.asarray(branch.radius_um, dtype=float).copy()
        source_points = points.copy()
        source_radii = radii.copy()
        point_type = np.full(len(points), SOURCE_CENTERLINE, dtype=np.uint8)
        cut_ids = [""] * len(points)
        edge_ids = np.full(len(points), -1, dtype=np.int64)
        core_distance = np.full(len(points), np.nan, dtype=float)
        core_positions = np.full((len(points), 3), np.nan, dtype=float)
        for endpoint in (-1, 0):
            port = port_by_node.get(int(branch.local_node_ids[endpoint]))
            if port is None:
                continue
            extension_points, extension_radii, distances, state = conditioned_port_profile(
                branch, endpoint, port, config
            )
            diagnostics.append(state)
            original_core = np.asarray(
                port.original_core_cut_position_um
                if port.original_core_cut_position_um is not None
                else port.original_position_um,
                dtype=float,
            )
            source_cut_id = port.source_core_cut_port_id or port.cut_port_id
            for index, (distance, point, point_radius) in enumerate(
                zip(distances, extension_points, extension_radii)
            ):
                point_rows.append(
                    {
                        "branch_id": branch.branch_id,
                        "extension_point_index": index,
                        "point_type": (
                            "CFD_EXTENSION_ORIGIN" if index == 0 else "CFD_EXTENSION"
                        ),
                        "source_cut_port_id": source_cut_id,
                        "active_cfd_port_id": port.cut_port_id,
                        "source_global_edge_id": port.global_edge_id,
                        "distance_from_cfd_extension_origin_um": float(distance),
                        "distance_from_core_boundary_um": float(
                            port.source_core_to_cfd_cut_length_um + distance
                        ),
                        "x_um": float(point[0]),
                        "y_um": float(point[1]),
                        "z_um": float(point[2]),
                        "radius_um": float(point_radius),
                        "original_core_cut_x_um": float(original_core[0]),
                        "original_core_cut_y_um": float(original_core[1]),
                        "original_core_cut_z_um": float(original_core[2]),
                        "original_source_boundary_type": "CORE_ROI_BOUNDARY",
                        "extension_origin_type": "CFD_DOMAIN_SOURCE_CUT",
                        "boundary_type_at_end": "CFD_BOUNDARY_PORT",
                        "point_boundary_type": (
                            "CFD_BOUNDARY_PORT"
                            if index == len(distances) - 1
                            else "CFD_DOMAIN_SOURCE_CUT"
                            if index == 0
                            else "NONE"
                        ),
                    }
                )
            derived_type = np.full(len(distances), CFD_DERIVED_EXTENSION, dtype=np.uint8)
            derived_type[0] = CORE_BOUNDARY_EXTENSION_ORIGIN
            derived_ids = [source_cut_id] * len(distances)
            derived_edges = np.full(len(distances), port.global_edge_id, dtype=np.int64)
            derived_distances = port.source_core_to_cfd_cut_length_um + distances
            derived_core = np.repeat(original_core[None, :], len(distances), axis=0)
            if endpoint == -1:
                points = np.vstack((points[:-1], extension_points))
                radii = np.concatenate((radii[:-1], extension_radii))
                point_type = np.concatenate((point_type[:-1], derived_type))
                cut_ids = [*cut_ids[:-1], *derived_ids]
                edge_ids = np.concatenate((edge_ids[:-1], derived_edges))
                core_distance = np.concatenate((core_distance[:-1], derived_distances))
                core_positions = np.vstack((core_positions[:-1], derived_core))
            else:
                points = np.vstack((extension_points[::-1], points[1:]))
                radii = np.concatenate((extension_radii[::-1], radii[1:]))
                point_type = np.concatenate((derived_type[::-1], point_type[1:]))
                cut_ids = [*derived_ids[::-1], *cut_ids[1:]]
                edge_ids = np.concatenate((derived_edges[::-1], edge_ids[1:]))
                core_distance = np.concatenate((derived_distances[::-1], core_distance[1:]))
                core_positions = np.vstack((derived_core[::-1], core_positions[1:]))
        if not np.array_equal(np.asarray(branch.points_um), source_points):
            raise GeometryValidationError("V6 modified source branch centerline samples")
        if not np.array_equal(np.asarray(branch.radius_um), source_radii):
            raise GeometryValidationError("V6 modified source branch radii")
        extended.points_um = points
        extended.radius_um = radii
        extended.arc_length_um = np.concatenate(
            ([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
        )
        extended.construction_point_type = point_type
        extended.construction_source_cut_port_id = tuple(cut_ids)
        extended.construction_source_global_edge_id = edge_ids
        extended.construction_distance_from_core_boundary_um = core_distance
        extended.construction_original_core_cut_position_um = core_positions
        output.append(extended)
    return output, point_rows, diagnostics


def _project_to_polyline(
    query: np.ndarray,
    points: np.ndarray,
    radii: np.ndarray,
    arc: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    starts = points[:-1]
    vectors = points[1:] - starts
    squared = np.sum(vectors**2, axis=1)
    valid = squared > 1.0e-20
    starts = starts[valid]
    vectors = vectors[valid]
    squared = squared[valid]
    start_radii = radii[:-1][valid]
    radius_delta = (radii[1:] - radii[:-1])[valid]
    start_arc = arc[:-1][valid]
    segment_length = np.sqrt(squared)
    tree = cKDTree(starts + 0.5 * vectors)
    candidate_count = min(8, len(starts))
    _, candidates = tree.query(query, k=candidate_count, workers=-1)
    if candidate_count == 1:
        candidates = candidates[:, None]
    candidate_starts = starts[candidates]
    candidate_vectors = vectors[candidates]
    fraction = np.sum(
        (query[:, None, :] - candidate_starts) * candidate_vectors,
        axis=2,
    ) / squared[candidates]
    fraction = np.clip(fraction, 0.0, 1.0)
    projected = candidate_starts + fraction[:, :, None] * candidate_vectors
    projected_radius = start_radii[candidates] + fraction * radius_delta[candidates]
    phi = np.linalg.norm(query[:, None, :] - projected, axis=2) - projected_radius
    selected = np.argmin(phi, axis=1)
    rows = np.arange(len(query))
    distance = start_arc[candidates[rows, selected]] + (
        fraction[rows, selected] * segment_length[candidates[rows, selected]]
    )
    return phi[rows, selected], distance


def _quintic_smoothstep(value: np.ndarray) -> np.ndarray:
    value = np.clip(np.asarray(value, dtype=float), 0.0, 1.0)
    return 6.0 * value**5 - 15.0 * value**4 + 10.0 * value**3


def build_continuous_junction_patch(
    roi: ROIRecord,
    branches_by_id: dict[int, BranchGeometry],
    junction_node_id: int,
    collars: list[JunctionCollar],
    config: CFDLumenConfig,
    *,
    cells_across_min_diameter: int | None = None,
) -> LocalJunctionPatch:
    """Extract one zero-set from a quintic core/branch scalar-field blend."""

    started = time.perf_counter()
    cells = int(
        cells_across_min_diameter
        or config.junction.implicit.cells_across_min_diameter
    )
    minimum_diameter = 2.0 * min(collar.collar_radius_um for collar in collars)
    grid_spacing = minimum_diameter / cells
    sample_spacing = 0.5 * grid_spacing
    point_blocks: list[np.ndarray] = []
    radius_blocks: list[np.ndarray] = []
    arc_blocks: list[np.ndarray] = []
    for collar in collars:
        points, radii = _dense_local_samples(
            branches_by_id[collar.branch_id], collar, sample_spacing
        )
        point_blocks.append(points)
        radius_blocks.append(radii)
        arc_blocks.append(
            np.concatenate(
                ([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
            )
        )
    sample_points = np.vstack(point_blocks)
    sample_radii = np.concatenate(radius_blocks)
    max_radius = float(sample_radii.max())
    padding = config.junction.bbox_padding_radius * max_radius
    minimum = np.min(sample_points - sample_radii[:, None], axis=0) - padding
    maximum = np.max(sample_points + sample_radii[:, None], axis=0) + padding
    junction_center = np.asarray(
        roi.local_node_positions_um[junction_node_id], dtype=float
    )
    junction_radius = float(roi.local_node_radius_um[junction_node_id])
    minimum = np.minimum(minimum, junction_center - junction_radius - padding)
    maximum = np.maximum(maximum, junction_center + junction_radius + padding)
    axes = [
        np.arange(minimum[axis], maximum[axis] + 0.5 * grid_spacing, grid_spacing)
        for axis in range(3)
    ]
    nx, ny, nz = map(len, axes)
    grid_cells = int(nx * ny * nz)
    if grid_cells > config.junction.implicit.max_grid_cells:
        raise GeometryValidationError(
            f"V6 local junction {junction_node_id} grid has {grid_cells:,} cells"
        )
    field = np.empty((nz, ny, nx), dtype=np.float32)
    core_tree = cKDTree(sample_points)
    core_k = min(config.junction.implicit.k_nearest, len(sample_points))
    plane_size = nx * ny
    slab = max(1, config.junction.implicit.chunk_size // max(1, plane_size))
    for z_start in range(0, nz, slab):
        z_end = min(nz, z_start + slab)
        zz, yy, xx = np.meshgrid(
            axes[2][z_start:z_end], axes[1], axes[0], indexing="ij"
        )
        query = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        distances, indices = core_tree.query(query, k=core_k, workers=-1)
        if core_k == 1:
            phi_junction = distances - sample_radii[indices]
        else:
            phi_junction = np.min(distances - sample_radii[indices], axis=1)
        branch_phi: list[np.ndarray] = []
        branch_arc: list[np.ndarray] = []
        for points, radii, arc in zip(point_blocks, radius_blocks, arc_blocks):
            phi, projected_arc = _project_to_polyline(query, points, radii, arc)
            branch_phi.append(phi)
            branch_arc.append(projected_arc)
        phi_stack = np.vstack(branch_phi)
        arc_stack = np.vstack(branch_arc)
        owner = np.argmin(phi_stack, axis=0)
        columns = np.arange(len(query))
        phi_branch = phi_stack[owner, columns]
        projected = arc_stack[owner, columns]
        transition_start = np.asarray(
            [collar.explicit_cap_distance_um for collar in collars], dtype=float
        )[owner]
        transition_end = np.asarray(
            [collar.collar_distance_um for collar in collars], dtype=float
        )[owner]
        coordinate = (projected - transition_start) / np.maximum(
            transition_end - transition_start, 1.0e-12
        )
        weight = _quintic_smoothstep(coordinate)
        phi = (1.0 - weight) * phi_junction + weight * phi_branch
        field[z_start:z_end] = phi.reshape((z_end - z_start, ny, nx)).astype(
            np.float32
        )
    if not (float(field.min()) < 0.0 < float(field.max())):
        raise GeometryValidationError(
            f"V6 local junction {junction_node_id} field does not bracket zero"
        )
    vertices_zyx, faces, _, _ = marching_cubes(
        field,
        level=0.0,
        spacing=(grid_spacing, grid_spacing, grid_spacing),
    )
    vertices = vertices_zyx[:, ::-1] + minimum
    raw = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    clean = _cleanup_continuous_patch(
        raw, config.v6.local_patch_cleanup_tolerance_um
    )
    if not clean.is_watertight:
        raise GeometryValidationError(
            f"V6 local junction {junction_node_id} patch is not watertight"
        )
    metadata: dict[str, Any] = {
        "junction_node_id": junction_node_id,
        "method": "continuous_implicit_field",
        "equation": "(1-w(s))*phi_J + w(s)*phi_B",
        "transition_coordinate": "nearest branch centerline segment projection",
        "smoothstep": "6t^5 - 15t^4 + 10t^3",
        "one_marching_cubes_extraction": True,
        "surface_loop_stitching": False,
        "global_smoothing": False,
        "local_sliver_cleanup_tolerance_um": (
            config.v6.local_patch_cleanup_tolerance_um
        ),
        "dtype": "float32",
        "cells_across_min_diameter": cells,
        "grid_spacing_um": grid_spacing,
        "grid_dimensions_xyz": [nx, ny, nz],
        "grid_cell_count": grid_cells,
        "bbox_min_um": minimum.tolist(),
        "bbox_max_um": maximum.tolist(),
        "centerline_sample_count": int(len(sample_points)),
        "raw_triangle_count": int(len(raw.faces)),
        "clean_triangle_count": int(len(clean.faces)),
        "degenerate_triangles_removed": int(len(raw.faces) - len(clean.faces)),
        "local_volume_um3": float(abs(clean.volume)),
        "runtime_s": time.perf_counter() - started,
        "collars": [
            {
                "branch_id": collar.branch_id,
                "endpoint_index": collar.endpoint_index,
                "transition_start_um": collar.explicit_cap_distance_um,
                "transition_end_um": collar.collar_distance_um,
                "pure_branch_overlap_end_um": collar.implicit_extent_um,
                "pure_branch_overlap_length_um": collar.overlap_length_um,
                "collar_radius_um": collar.collar_radius_um,
            }
            for collar in collars
        ],
    }
    return LocalJunctionPatch(
        junction_node_id=junction_node_id,
        raw_mesh=raw,
        clean_mesh=clean,
        collars=collars,
        centerline_points_um=sample_points,
        centerline_radius_um=sample_radii,
        metadata=metadata,
    )


def _merge_distance(collar: JunctionCollar) -> float:
    return float(collar.collar_distance_um + 0.5 * collar.overlap_length_um)


def _trim_branches_to_merge_planes(
    branches: list[BranchGeometry],
    collars_by_junction: dict[int, list[JunctionCollar]],
) -> tuple[list[BranchGeometry], dict[int, list[JunctionCollar]]]:
    by_branch: dict[int, list[JunctionCollar]] = {}
    for collars in collars_by_junction.values():
        for collar in collars:
            by_branch.setdefault(collar.branch_id, []).append(collar)
    output: list[BranchGeometry] = []
    for branch in branches:
        lower = 0.0
        upper = float(branch.arc_length_um[-1])
        for collar in by_branch.get(branch.branch_id, []):
            if collar.endpoint_index == 0:
                lower = max(lower, _merge_distance(collar))
            else:
                upper = min(
                    upper,
                    float(branch.arc_length_um[-1]) - _merge_distance(collar),
                )
        if upper - lower <= 1.0e-8:
            raise GeometryValidationError(
                f"V6_PURE_BRANCH_OVERLAP_FAILED: branch {branch.branch_id} has no explicit region"
            )
        inside = branch.arc_length_um[
            (branch.arc_length_um > lower) & (branch.arc_length_um < upper)
        ]
        targets = np.concatenate(([lower], inside, [upper]))
        trimmed = replace(branch)
        trimmed.points_um = np.column_stack(
            [
                np.interp(targets, branch.arc_length_um, branch.points_um[:, axis])
                for axis in range(3)
            ]
        )
        trimmed.radius_um = np.interp(targets, branch.arc_length_um, branch.radius_um)
        trimmed.arc_length_um = targets - lower
        output.append(trimmed)
    return output, by_branch


def _clip_patch_at_pure_branch_planes(
    patch: LocalJunctionPatch,
    branches_by_id: dict[int, BranchGeometry],
    tube_sides: int,
) -> tuple[
    trimesh.Trimesh,
    dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]],
]:
    clipped = patch.clean_mesh.copy()
    planes: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for collar in patch.collars:
        center, _, tangent = sample_from_junction(
            branches_by_id[collar.branch_id],
            collar.endpoint_index,
            _merge_distance(collar),
        )
        candidate = trimesh.intersections.slice_mesh_plane(
            clipped,
            plane_normal=-tangent,
            plane_origin=center,
            cap=False,
        )
        if candidate is None or not len(candidate.faces):
            raise GeometryValidationError(
                f"V6_PURE_BRANCH_CLIP_FAILED: J{patch.junction_node_id}/B{collar.branch_id}"
            )
        clipped = candidate
        clipped.merge_vertices(digits_vertex=8)
        clipped.remove_unreferenced_vertices()
        planes[collar.branch_id] = (center, tangent)
    loops = _boundary_loops(clipped)
    tolerance = max(float(patch.metadata["grid_spacing_um"]) * 1.0e-3, 1.0e-6)
    for center, tangent in planes.values():
        _loop_at_plane(clipped, loops, center, tangent, tolerance)
    clipped = _collapse_boundary_loops_to_count(clipped, loops, tube_sides)
    loops = _boundary_loops(clipped)
    output: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for branch_id, (center, tangent) in planes.items():
        loop = _loop_at_plane(clipped, loops, center, tangent, tolerance)
        if len(loop) != tube_sides:
            raise GeometryValidationError(
                f"V6_LOOP_WELD_FAILED: local loop has {len(loop)} vertices, expected {tube_sides}"
            )
        output[branch_id] = (loop, center, tangent)
    clipped = _condition_pure_branch_surface(
        clipped, patch.collars, branches_by_id
    )
    return clipped, output


def _condition_pure_branch_surface(
    mesh: trimesh.Trimesh,
    collars: list[JunctionCollar],
    branches_by_id: dict[int, BranchGeometry],
) -> trimesh.Trimesh:
    """Project only the pure-field overlap onto its analytic branch zero-set.

    Marching cubes samples the correct ``phi_B=0`` surface but its last grid
    row does not generally share the analytic tube normal.  A quintic,
    extension-side-only projection removes that discretization error before
    the same-section weld; the core and transition vertices remain untouched.
    """

    vertices = np.asarray(mesh.vertices, dtype=float).copy()
    for collar in collars:
        branch = branches_by_id[collar.branch_id]
        points, radii, arc = (
            (branch.points_um, branch.radius_um, branch.arc_length_um)
            if collar.endpoint_index == 0
            else (
                branch.points_um[::-1],
                branch.radius_um[::-1],
                float(branch.arc_length_um[-1]) - branch.arc_length_um[::-1],
            )
        )
        _, distance = _project_to_polyline(vertices, points, radii, arc)
        start = float(collar.collar_distance_um)
        merge = _merge_distance(collar)
        conditioning_end = merge
        selected = (distance >= start) & (distance <= merge + 1.0e-6)
        for vertex_id in np.flatnonzero(selected):
            center, radius, tangent = sample_from_junction(
                branch, collar.endpoint_index, float(distance[vertex_id])
            )
            relative = vertices[vertex_id] - center
            radial = relative - float(np.dot(relative, tangent)) * tangent
            norm = float(np.linalg.norm(radial))
            if norm <= 1.0e-12:
                continue
            target = center + radius * radial / norm
            weight = float(
                _quintic_smoothstep(
                    np.asarray(
                        [
                            (distance[vertex_id] - start)
                            / max(conditioning_end - start, 1.0e-12)
                        ]
                    )
                )[0]
            )
            vertices[vertex_id] = (1.0 - weight) * vertices[vertex_id] + weight * target
    conditioned = trimesh.Trimesh(
        vertices=vertices, faces=np.asarray(mesh.faces), process=False
    )
    trimesh.repair.fix_normals(conditioned, multibody=True)
    return conditioned


def _loop_metrics(
    first: np.ndarray,
    second: np.ndarray,
    center: np.ndarray,
) -> dict[str, float]:
    first_center = np.mean(first, axis=0)
    second_center = np.mean(second, axis=0)
    first_radius = float(np.mean(np.linalg.norm(first - center[None, :], axis=1)))
    second_radius = float(np.mean(np.linalg.norm(second - center[None, :], axis=1)))
    tree_first = cKDTree(first)
    tree_second = cKDTree(second)
    hausdorff = max(
        float(tree_first.query(second, k=1)[0].max()),
        float(tree_second.query(first, k=1)[0].max()),
    )
    return {
        "loop_center_error_um": float(np.linalg.norm(first_center - second_center)),
        "loop_radius_error_um": abs(first_radius - second_radius),
        "loop_hausdorff_um": hausdorff,
    }


def _combine_welded(
    patches: dict[int, trimesh.Trimesh],
    patch_loops: dict[tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]],
    explicit_meshes: dict[int, trimesh.Trimesh],
    explicit_loops: dict[tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]],
    collars: dict[int, list[JunctionCollar]],
    config: CFDLumenConfig,
) -> tuple[trimesh.Trimesh, list[dict[str, Any]]]:
    loop_rows: list[dict[str, Any]] = []
    moved_patches = {node_id: mesh.copy() for node_id, mesh in patches.items()}
    for node_id, node_collars in collars.items():
        for collar in node_collars:
            key = (node_id, collar.branch_id)
            local_loop, center, _ = patch_loops[key]
            explicit_loop, _, _ = explicit_loops[key]
            local_points = np.asarray(moved_patches[node_id].vertices[local_loop], dtype=float)
            explicit_points = np.asarray(
                explicit_meshes[collar.branch_id].vertices[explicit_loop], dtype=float
            )
            aligned, shift, reversed_order, phase_error = _phase_align(
                local_points, explicit_points
            )
            metrics = _loop_metrics(local_points, aligned, center)
            radius = float(collar.collar_radius_um)
            eligible = (
                metrics["loop_center_error_um"]
                <= config.v6.loop_center_tolerance_radius_fraction * radius
                and metrics["loop_radius_error_um"]
                <= config.v6.loop_radius_tolerance_fraction * radius
                and metrics["loop_hausdorff_um"]
                <= config.v6.loop_hausdorff_tolerance_radius_fraction * radius
            )
            loop_rows.append(
                {
                    "junction_node_id": node_id,
                    "branch_id": collar.branch_id,
                    "merge_region": "PURE_BRANCH",
                    "merge_distance_um": _merge_distance(collar),
                    "source_radius_um": radius,
                    **metrics,
                    "phase_shift": int(shift),
                    "explicit_order_reversed": bool(reversed_order),
                    "phase_alignment_squared_error_um2": float(phase_error),
                    "clip_weld_eligible": bool(eligible),
                }
            )
            if not eligible:
                raise GeometryValidationError(
                    f"V6_LOOP_WELD_TOLERANCE_FAILED: J{node_id}/B{collar.branch_id}"
                )
            vertices = np.asarray(moved_patches[node_id].vertices, dtype=float).copy()
            vertices[local_loop] = aligned
            moved_patches[node_id] = trimesh.Trimesh(
                vertices=vertices,
                faces=np.asarray(moved_patches[node_id].faces),
                process=False,
            )
    combined = trimesh.util.concatenate(
        [*moved_patches.values(), *explicit_meshes.values()]
    )
    combined.merge_vertices(digits_vertex=8)
    keep = combined.nondegenerate_faces(height=1.0e-12) & combined.unique_faces()
    combined.update_faces(keep)
    combined.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(combined, multibody=True)
    if combined.is_watertight and combined.volume < 0.0:
        combined.invert()
    if not combined.is_watertight:
        raise GeometryValidationError("V6_LOOP_WELD_TOPOLOGY_FAILED: surface is open")
    return combined, loop_rows


def _analytic_ring_from_previous(
    previous: np.ndarray,
    previous_center: np.ndarray,
    center: np.ndarray,
    tangent: np.ndarray,
    radius: float,
) -> np.ndarray:
    radial = np.asarray(previous, dtype=float) - previous_center[None, :]
    radial -= (radial @ tangent)[:, None] * tangent[None, :]
    norms = np.linalg.norm(radial, axis=1)
    if np.any(norms <= 1.0e-12):
        raise GeometryValidationError("V6_PURE_FIELD_RING_FAILED: zero radial vector")
    return center[None, :] + radius * radial / norms[:, None]


def _pure_field_adapter(
    patch_mesh: trimesh.Trimesh,
    patch_loop: np.ndarray,
    explicit_mesh: trimesh.Trimesh,
    explicit_loop: np.ndarray,
    branch: BranchGeometry,
    collar: JunctionCollar,
    tube_sides: int,
) -> tuple[trimesh.Trimesh, np.ndarray, dict[str, Any]]:
    """Align both merge loops through the same analytic branch field.

    The adapter lies wholly after ``w=1``.  It is not a junction transition:
    every station is projected to the unchanged centerline/radius tube field.
    """

    start_distance = float(collar.collar_distance_um)
    end_distance = float(collar.implicit_extent_um)
    start_center, start_radius, start_tangent = sample_from_junction(
        branch, collar.endpoint_index, start_distance
    )
    end_center, _, _ = sample_from_junction(
        branch, collar.endpoint_index, end_distance
    )
    patch_vertices = np.asarray(patch_mesh.vertices, dtype=float).copy()
    patch_points = patch_vertices[patch_loop]
    # The first common ring is the exact MC zero-set boundary.  Moving it to a
    # perfect circle creates a normal kink in its incident patch triangles.
    # Analytic projection therefore begins at the next pure-branch station.
    snapped_start = patch_points.copy()
    patch_mesh = trimesh.Trimesh(
        vertices=patch_vertices,
        faces=np.asarray(patch_mesh.faces),
        process=False,
    )

    distance_one = start_distance + (end_distance - start_distance) / 3.0
    center_one, radius_one, tangent_one = sample_from_junction(
        branch, collar.endpoint_index, distance_one
    )
    ring_one = _analytic_ring_from_previous(
        snapped_start, start_center, center_one, tangent_one, radius_one
    )
    distance_two = start_distance + 2.0 * (end_distance - start_distance) / 3.0
    center_two, radius_two, tangent_two = sample_from_junction(
        branch, collar.endpoint_index, distance_two
    )
    ring_two_dense = _analytic_ring_from_previous(
        ring_one, center_one, center_two, tangent_two, radius_two
    )
    ring_two = _resample_loop(ring_two_dense, tube_sides)
    explicit_points = np.asarray(explicit_mesh.vertices[explicit_loop], dtype=float)
    explicit_aligned, shift, reversed_order, phase_error = _phase_align(
        ring_two, explicit_points
    )

    blocks = (snapped_start, ring_one, ring_two, explicit_aligned)
    offsets = np.cumsum((0, *[len(block) for block in blocks[:-1]]))
    faces: list[tuple[int, int, int]] = []
    for block_index in range(len(blocks) - 1):
        first_count = len(blocks[block_index])
        second_count = len(blocks[block_index + 1])
        zipped = _zipper_faces(first_count, second_count)
        zipped = np.where(
            zipped < first_count,
            zipped + offsets[block_index],
            zipped - first_count + offsets[block_index + 1],
        )
        faces.extend(map(tuple, zipped.tolist()))
    adapter = trimesh.Trimesh(
        vertices=np.vstack(blocks),
        faces=np.asarray(faces, dtype=np.int64),
        process=False,
    )
    loop_metrics = _loop_metrics(snapped_start, ring_one, start_center)
    return patch_mesh, adapter, {
        "junction_node_id": collar.junction_node_id,
        "branch_id": collar.branch_id,
        "merge_region": "PURE_BRANCH",
        "merge_backend": "common_target_loop_vertex_weld",
        "merge_start_um": start_distance,
        "merge_end_um": end_distance,
        "pure_branch_field_station_count": 4,
        "patch_loop_vertex_count": int(len(patch_loop)),
        "target_loop_vertex_count": int(tube_sides),
        "phase_shift": int(shift),
        "explicit_order_reversed": bool(reversed_order),
        "phase_alignment_squared_error_um2": float(phase_error),
        **loop_metrics,
    }


def _build_common_target_weld(
    patches: dict[int, LocalJunctionPatch],
    construction_branches: list[BranchGeometry],
    collars: dict[int, list[JunctionCollar]],
    tube_sides: int,
) -> tuple[
    trimesh.Trimesh,
    dict[int, trimesh.Trimesh],
    list[BranchGeometry],
    list[dict[str, Any]],
    np.ndarray,
]:
    branch_by_id = {branch.branch_id: branch for branch in construction_branches}
    pure_collars = {
        node_id: [
            replace(
                collar,
                explicit_cap_distance_um=collar.collar_distance_um,
                implicit_extent_um=collar.implicit_extent_um,
            )
            for collar in node_collars
        ]
        for node_id, node_collars in collars.items()
    }
    clipped_patches: dict[int, trimesh.Trimesh] = {}
    implicit_loops: dict[
        tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]
    ] = {}
    for node_id, patch in patches.items():
        tolerance = max(float(patch.metadata["grid_spacing_um"]) * 1.0e-3, 1.0e-6)
        clipped, loops = clip_implicit_patch_to_transition_inner(
            patch.clean_mesh,
            pure_collars[node_id],
            branch_by_id,
            tolerance,
            tube_sides,
        )
        clipped_patches[node_id] = clipped
        for branch_id, payload in loops.items():
            implicit_loops[(node_id, branch_id)] = payload
    trimmed, collars_by_branch = trim_branches_to_transition_outer(
        construction_branches, pure_collars
    )
    explicit_meshes, explicit_loops = build_open_explicit_tubes(
        trimmed, collars_by_branch, tube_sides
    )
    adapters: list[trimesh.Trimesh] = []
    rows: list[dict[str, Any]] = []
    for node_id, node_collars in pure_collars.items():
        for collar in node_collars:
            key = (node_id, collar.branch_id)
            patch_loop, _, _ = implicit_loops[key]
            explicit_loop, _, _ = explicit_loops[key]
            moved_patch, adapter, row = _pure_field_adapter(
                clipped_patches[node_id],
                patch_loop,
                explicit_meshes[collar.branch_id],
                explicit_loop,
                branch_by_id[collar.branch_id],
                collar,
                tube_sides,
            )
            clipped_patches[node_id] = moved_patch
            adapters.append(adapter)
            rows.append(row)
    mesh, provenance_labels = combine_stitched_surfaces(
        list(clipped_patches.values()), adapters, list(explicit_meshes.values())
    )
    if not mesh.is_watertight:
        raise GeometryValidationError("V6_COMMON_TARGET_WELD_FAILED: surface is open")
    provenance_labels[provenance_labels == TRANSITION_COLLAR_FACE] = PURE_BRANCH_FACE
    geometric_labels = classify_v6_face_regions(mesh, patches, construction_branches)
    patch_faces = provenance_labels == JUNCTION_CORE_FACE
    provenance_labels[patch_faces] = geometric_labels[patch_faces]
    return mesh, clipped_patches, trimmed, rows, provenance_labels


def _boolean_fallback(
    patches: dict[int, LocalJunctionPatch],
    branches: list[BranchGeometry],
    collars: dict[int, list[JunctionCollar]],
    tube_sides: int,
) -> trimesh.Trimesh:
    trimmed, _ = _trim_branches_to_merge_planes(branches, collars)
    tubes = [build_variable_radius_tube(branch, tube_sides) for branch in trimmed]
    return balanced_manifold_union(
        [*[patch.clean_mesh for patch in patches.values()], *tubes]
    )


def classify_v6_face_regions(
    mesh: trimesh.Trimesh,
    patches: dict[int, LocalJunctionPatch],
    branches: list[BranchGeometry],
) -> np.ndarray:
    centers = np.asarray(mesh.triangles_center, dtype=float)
    labels = np.full(len(mesh.faces), EXPLICIT_BRANCH_FACE, dtype=np.uint8)
    for patch in patches.values():
        for collar in patch.collars:
            branch = next(item for item in branches if item.branch_id == collar.branch_id)
            points, _, arc = (
                (branch.points_um, branch.radius_um, branch.arc_length_um)
                if collar.endpoint_index == 0
                else (
                    branch.points_um[::-1],
                    branch.radius_um[::-1],
                    float(branch.arc_length_um[-1]) - branch.arc_length_um[::-1],
                )
            )
            _, distance = _project_to_polyline(centers, points, np.zeros(len(points)), arc)
            radial = np.linalg.norm(
                centers
                - points[np.argmin(
                    np.linalg.norm(centers[:, None, :] - points[None, :, :], axis=2),
                    axis=1,
                )],
                axis=1,
            )
            local = (distance <= _merge_distance(collar) + collar.collar_radius_um) & (
                radial <= 2.5 * collar.collar_radius_um
            )
            labels[local & (distance <= collar.explicit_cap_distance_um)] = (
                JUNCTION_CORE_FACE
            )
            labels[
                local
                & (distance > collar.explicit_cap_distance_um)
                & (distance < collar.collar_distance_um)
            ] = TRANSITION_COLLAR_FACE
            labels[
                local
                & (distance >= collar.collar_distance_um)
                & (distance <= _merge_distance(collar) + 1.0e-9)
            ] = PURE_BRANCH_FACE
    return labels


def build_continuous_field_hybrid(
    branches: list[BranchGeometry],
    roi: ROIRecord,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    tube_sides: int | None = None,
    cells_across_min_diameter: int | None = None,
) -> tuple[trimesh.Trimesh, HybridBuildDetails, list[dict[str, Any]]]:
    """Build the v6 surface with one field per junction and pure-branch merging."""

    started = time.perf_counter()
    sides = int(tube_sides or config.geometry.tube_sides)
    collars = define_junction_collars(roi, branches, config)
    construction_branches, port_rows, port_diagnostics = (
        extend_branches_to_cfd_ports_v6(branches, ports, config)
    )
    branch_by_id = {branch.branch_id: branch for branch in construction_branches}
    field_started = time.perf_counter()
    patches = {
        node_id: build_continuous_junction_patch(
            roi,
            branch_by_id,
            node_id,
            node_collars,
            config,
            cells_across_min_diameter=cells_across_min_diameter,
        )
        for node_id, node_collars in collars.items()
    }
    field_runtime = time.perf_counter() - field_started
    merge_started = time.perf_counter()
    merge_backend = "common_target_loop_vertex_weld"
    fallback_reason: str | None = None
    try:
        mesh, clipped_patches, trimmed, loop_rows, face_region = (
            _build_common_target_weld(
            patches,
            construction_branches,
            collars,
            sides,
            )
        )
    except Exception as exc:
        fallback_reason = f"{type(exc).__name__}: {exc}"
        merge_backend = "pure_branch_local_manifold_boolean"
        mesh = _boolean_fallback(
            patches, construction_branches, collars, sides
        )
        clipped_patches = {
            node_id: patch.clean_mesh for node_id, patch in patches.items()
        }
        trimmed, _ = _trim_branches_to_merge_planes(
            construction_branches, collars
        )
        face_region = classify_v6_face_regions(
            mesh, patches, construction_branches
        )
        loop_rows = [
            {
                "junction_node_id": node_id,
                "branch_id": collar.branch_id,
                "merge_region": "PURE_BRANCH",
                "merge_distance_um": _merge_distance(collar),
                "clip_weld_eligible": False,
                "fallback_reason": fallback_reason,
            }
            for node_id, node_collars in collars.items()
            for collar in node_collars
        ]
    merge_runtime = time.perf_counter() - merge_started
    for patch in patches.values():
        patch.metadata["merge_backend"] = merge_backend
        patch.metadata["merge_fallback_reason"] = fallback_reason
    details = HybridBuildDetails(
        patches=patches,
        merged_junction_meshes={
            node_id: clipped_patches.get(node_id, patch.clean_mesh)
            for node_id, patch in patches.items()
        },
        trimmed_branches=trimmed,
        merge_steps=[
            {
                "junction_node_id": node_id,
                "merge_backend": merge_backend,
                "merge_region": "PURE_BRANCH",
                "incident_branch_ids": [
                    collar.branch_id for collar in node_collars
                ],
            }
            for node_id, node_collars in collars.items()
        ],
        explicit_cap_planes=[],
        runtime_s={
            "continuous_local_field": field_runtime,
            "pure_branch_merge": merge_runtime,
            "reconstruction_total": time.perf_counter() - started,
            "wall_total_with_step_qc": time.perf_counter() - started,
        },
        transition_backend="continuous_implicit_field",
        face_region=face_region,
        transition_rows=loop_rows,
        constructed_branches=construction_branches,
        port_extension_rows=port_rows,
        transition_fallback_reason=fallback_reason,
    )
    return mesh, details, port_diagnostics


def pure_branch_loop_section_metrics(
    patch: trimesh.Trimesh,
    explicit: trimesh.Trimesh,
    branch: BranchGeometry,
    collar: JunctionCollar,
) -> dict[str, Any]:
    center, source_radius, tangent = sample_from_junction(
        branch, collar.endpoint_index, _merge_distance(collar)
    )
    first = _section_polygon(patch, center, tangent)
    second = _section_polygon(explicit, center, tangent)
    return {
        "junction_node_id": collar.junction_node_id,
        "branch_id": collar.branch_id,
        "source_radius_um": source_radius,
        "patch_area_um2": first[0] if first else None,
        "explicit_area_um2": second[0] if second else None,
    }
