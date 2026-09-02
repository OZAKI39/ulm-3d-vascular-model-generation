"""
Two-dimensional open-edge boundary construction for the X-Z lumen solver.

The finite-volume projection solves incompressibility on face fluxes, so inlet
and outlet conditions must also be expressed on faces.  This module therefore
does not paint a thick velocity band inside the lumen.  Instead, it finds the
actual mask-boundary faces near root and terminal vessel endpoints, assigns
each selected face a normal vector, a length, a label, and a prescribed signed
flux, then passes those faces to the projection code unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from ulm_vascular_model_generator.utils.core.models import Vessel

from ..core.types import GridDomain, RasterizedVessels


@dataclass(frozen=True)
class BoundaryFluxFields:
    """
    Discrete inlet and outlet open edges used by the finite-volume projection.

    There are two parallel representations in this container:

    * cell-centered helper arrays, such as inlet_label and boundary_velocity,
      which are convenient for diagnostics, plotting, and saved NPZ fields;
    * explicit face arrays, such as open_face_index_ij and open_face_flux,
      which are the authoritative data used by the finite-volume solver.

    The face arrays are the physical boundary condition.  The cell arrays are
    derived annotations and must not be used as a separate flux definition.
    """

    inlet_label: np.ndarray
    outlet_label: np.ndarray
    boundary_velocity_xz_um_s: np.ndarray
    boundary_normal_xz: np.ndarray
    boundary_weight: np.ndarray
    boundary_edge_length_um: np.ndarray
    open_boundary_flux_um2_s: np.ndarray
    inlet_mask: np.ndarray
    outlet_mask: np.ndarray
    target_inlet_q2d_um2_s: float
    target_outlet_q2d_um2_s: float
    inlet_target_by_label_um2_s: np.ndarray
    outlet_target_by_label_um2_s: np.ndarray
    inlet_ids: tuple[int, ...]
    outlet_ids: tuple[int, ...]
    open_face_cell_ij: np.ndarray
    open_face_index_ij: np.ndarray
    open_face_axis: np.ndarray
    open_face_normal_xz: np.ndarray
    open_face_center_xz_um: np.ndarray
    open_face_length_um: np.ndarray
    open_face_label: np.ndarray
    open_face_kind: np.ndarray
    open_face_flux_um2_s: np.ndarray
    open_section_point_xz_um: np.ndarray
    open_section_outward_normal_xz: np.ndarray
    open_section_tangent_xz: np.ndarray
    open_section_half_width_um: np.ndarray
    open_section_label: np.ndarray
    open_section_kind: np.ndarray
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class BoundaryFaceCatalog:
    """
    Flat list of all lumen-to-non-lumen faces in the current mask.

    Building this once avoids repeatedly scanning the entire mask for every
    root and terminal vessel.  Each row represents one candidate open edge with
    a cell location, a face-array index, a normal, a physical center position,
    and a face length.
    """

    cell_ij: np.ndarray
    face_index_ij: np.ndarray
    axis: np.ndarray
    normal_xz: np.ndarray
    center_xz_um: np.ndarray
    length_um: np.ndarray


@dataclass(frozen=True)
class _BoundaryProfile:
    """
    Selected faces and flux distribution for one vessel endpoint.

    The profile is built independently for each root inlet or terminal outlet.
    It stores only faces selected for that endpoint, plus the local parabolic
    weighting that distributes the vessel's 2D target flux over those faces.
    """

    cell_ij: np.ndarray
    face_index_ij: np.ndarray
    axis: np.ndarray
    normal_xz: np.ndarray
    center_xz_um: np.ndarray
    length_um: np.ndarray
    weight: np.ndarray
    signed_flux_um2_s: np.ndarray
    target_q2d_um2_s: float


def build_flux_boundaries(
    domain: GridDomain,
    raster: RasterizedVessels,
    vessels: list[Vessel] | tuple[Vessel, ...],
    *,
    effective_thickness_um: float,
    depth_cells: float,
) -> BoundaryFluxFields:
    """
    Build root-inlet and terminal-outlet open edges from 2D flux targets.

    For a format-v3 planar graph, the loader first maps exported true 2-D flux
    to the pipeline's explicit extrusion depth. Legacy graphs already provide
    3-D-equivalent flow. Thus the internal vessel value is divided by the same
    effective thickness here to recover the authoritative 2-D target flux.
    That target is spread over actual mask-boundary faces near the endpoint.
    """
    # These grid arrays are not the primary boundary condition.  
    # They are compact, human-readable fields used for saved outputs, visualization, and diagnostics.  
    # If multiple open faces touch the same cell, their weights, normals, edge lengths, and signed fluxes are accumulated here.
    shape           = domain.shape
    inlet_label     = np.zeros(shape, dtype=np.int32)
    outlet_label    = np.zeros(shape, dtype=np.int32)
    velocity        = np.zeros((*shape, 2), dtype=np.float32)
    normal          = np.zeros((*shape, 2), dtype=np.float32)
    normal_sum      = np.zeros((*shape, 2), dtype=np.float64)
    weight          = np.zeros(shape, dtype=np.float32)
    edge_length     = np.zeros(shape, dtype=np.float32)
    open_flux       = np.zeros(shape, dtype=np.float64)

    target_inlet        = 0.0
    target_outlet       = 0.0
    inlet_ids: list[int] = []
    outlet_ids: list[int] = []
    inlet_targets       = [0.0]
    outlet_targets      = [0.0]
    skipped_inlets: list[int] = []
    skipped_outlets: list[int] = []

    # Convert 3D flow Q [um^3/s] to the planar 2D flux q [um^2/s] used by the X-Z finite-volume system.  
    # The epsilon guard avoids a division-by-zero failure if a malformed config passes zero thickness.
    h0 = max(float(effective_thickness_um), np.finfo(float).eps)

    # ====================================================================================================
    # ======= Catalog every lumen boundary face once. 
    # ====================================================================================================
    # Search and record all outer edge surfaces (faces) of the entire fluid channel (cavity)
    # Enumerate every lumen face adjacent to non-lumen space. These faces are possible open boundaries, but most of them are solid walls.
    # Identify all exposed surfaces and compile them into a catalog
    face_catalog = enumerate_lumen_boundary_faces(domain, raster.lumen_mask)

    # These lists hold the authoritative face-based boundary condition.  
    # They are concatenated at the end into one row per selected open face.
    open_cells: list[np.ndarray] = []
    open_face_indices: list[np.ndarray] = []
    open_axes: list[np.ndarray] = []
    open_normals: list[np.ndarray] = []
    open_centers: list[np.ndarray] = []
    open_lengths: list[np.ndarray] = []
    open_labels: list[np.ndarray] = []
    open_kinds: list[np.ndarray] = []
    open_fluxes: list[np.ndarray] = []
    section_points: list[np.ndarray] = []
    section_normals: list[np.ndarray] = []
    section_tangents: list[np.ndarray] = []
    section_half_widths: list[float] = []
    section_labels: list[int] = []
    section_kinds: list[int] = []

    label = 1
    for vessel in _root_vessels(vessels):
        # Rationally distribute the total water flow to every pixel on the cut surface.
        profile = _boundary_profile(
            domain,
            raster,
            face_catalog,
            vessel,
            q2d_um2_s=max(float(vessel.flow_rate), 0.0) / h0,
            depth_cells=depth_cells,
            at_distal_end=False,
        )

        if profile.cell_ij.size == 0:
            # Missing inlet faces usually means the mask or endpoint geometry
            # no longer exposes a clear boundary at this root.  Keep the ID in
            # metadata so the failure can be audited.
            skipped_inlets.append(int(vessel.vid))
            continue
        _accumulate_profile(
            profile=profile,
            label=label,
            kind=-1,
            label_grid=inlet_label,
            velocity=velocity,
            normal_sum=normal_sum,
            weight=weight,
            edge_length=edge_length,
            open_flux=open_flux,
        )
        _append_open_faces(
            profile,
            label=label,
            kind=-1,
            open_cells=open_cells,
            open_face_indices=open_face_indices,
            open_axes=open_axes,
            open_normals=open_normals,
            open_centers=open_centers,
            open_lengths=open_lengths,
            open_labels=open_labels,
            open_kinds=open_kinds,
            open_fluxes=open_fluxes,
        )
        _append_open_section(
            profile,
            label=label,
            kind=-1,
            section_points=section_points,
            section_normals=section_normals,
            section_tangents=section_tangents,
            section_half_widths=section_half_widths,
            section_labels=section_labels,
            section_kinds=section_kinds,
        )
        target_inlet += profile.target_q2d_um2_s
        inlet_ids.append(int(vessel.vid))
        inlet_targets.append(float(profile.target_q2d_um2_s))
        label += 1

    label = 1
    for vessel in _terminal_vessels(vessels):
        # Terminal outlets use the distal endpoint.  
        # Their outward normal follows the parent-to-child vessel direction.
        profile = _boundary_profile(
            domain,
            raster,
            face_catalog,
            vessel,
            q2d_um2_s=max(float(vessel.flow_rate), 0.0) / h0,
            depth_cells=depth_cells,
            at_distal_end=True,
        )
        if profile.cell_ij.size == 0:
            # A skipped outlet is serious because per-outlet flux controls the
            # branch ratios that particles experience downstream.  Store it in
            # metadata rather than silently dropping the problem.
            skipped_outlets.append(int(vessel.vid))
            continue
        _accumulate_profile(
            profile=profile,
            label=label,
            kind=1,
            label_grid=outlet_label,
            velocity=velocity,
            normal_sum=normal_sum,
            weight=weight,
            edge_length=edge_length,
            open_flux=open_flux,
        )
        _append_open_faces(
            profile,
            label=label,
            kind=1,
            open_cells=open_cells,
            open_face_indices=open_face_indices,
            open_axes=open_axes,
            open_normals=open_normals,
            open_centers=open_centers,
            open_lengths=open_lengths,
            open_labels=open_labels,
            open_kinds=open_kinds,
            open_fluxes=open_fluxes,
        )
        _append_open_section(
            profile,
            label=label,
            kind=1,
            section_points=section_points,
            section_normals=section_normals,
            section_tangents=section_tangents,
            section_half_widths=section_half_widths,
            section_labels=section_labels,
            section_kinds=section_kinds,
        )
        target_outlet += profile.target_q2d_um2_s
        outlet_ids.append(int(vessel.vid))
        outlet_targets.append(float(profile.target_q2d_um2_s))
        label += 1

    # A cell can receive several selected open faces, for example around a
    # corner-shaped raster endpoint.  Normalize the accumulated normal so the
    # helper field stores direction rather than flux magnitude.
    norm    = np.linalg.norm(normal_sum, axis=-1, keepdims=True)
    normal  = np.divide(normal_sum, np.maximum(norm, np.finfo(float).eps), out=np.zeros_like(normal_sum), where=norm > 0.0).astype(
        np.float32
    )
    inlet_mask  = inlet_label > 0
    outlet_mask = outlet_label > 0

    # This is only the requested target-balance error.  Final physical
    # acceptance later uses actual face fluxes after projection.
    target_balance_error = abs(target_inlet - target_outlet) / max(target_inlet, np.finfo(float).eps)

    return BoundaryFluxFields(
        inlet_label=inlet_label,
        outlet_label=outlet_label,
        boundary_velocity_xz_um_s=velocity,
        boundary_normal_xz=normal,
        boundary_weight=weight,
        boundary_edge_length_um=edge_length,
        open_boundary_flux_um2_s=open_flux.astype(np.float32),
        inlet_mask=inlet_mask,
        outlet_mask=outlet_mask,
        target_inlet_q2d_um2_s=float(target_inlet),
        target_outlet_q2d_um2_s=float(target_outlet),
        inlet_target_by_label_um2_s=np.asarray(inlet_targets, dtype=np.float64),
        outlet_target_by_label_um2_s=np.asarray(outlet_targets, dtype=np.float64),
        inlet_ids=tuple(inlet_ids),
        outlet_ids=tuple(outlet_ids),
        open_face_cell_ij=_concat_2d(open_cells),
        open_face_index_ij=_concat_2d(open_face_indices),
        open_face_axis=_concat_1d(open_axes, dtype=np.int8),
        open_face_normal_xz=_concat_2d(open_normals).astype(np.float32, copy=False),
        open_face_center_xz_um=_concat_2d_float(open_centers),
        open_face_length_um=_concat_1d(open_lengths, dtype=np.float64),
        open_face_label=_concat_1d(open_labels, dtype=np.int32),
        open_face_kind=_concat_1d(open_kinds, dtype=np.int8),
        open_face_flux_um2_s=_concat_1d(open_fluxes, dtype=np.float64),
        open_section_point_xz_um=_stack_section_vectors(section_points),
        open_section_outward_normal_xz=_stack_section_vectors(section_normals),
        open_section_tangent_xz=_stack_section_vectors(section_tangents),
        open_section_half_width_um=np.asarray(section_half_widths, dtype=np.float64),
        open_section_label=np.asarray(section_labels, dtype=np.int32),
        open_section_kind=np.asarray(section_kinds, dtype=np.int8),
        metadata={
            "n_inlet_boundaries": int(len(inlet_ids)),
            "n_outlet_boundaries": int(len(outlet_ids)),
            "n_open_boundary_edges": int(sum(x.shape[0] for x in open_cells)),
            "inlet_boundary_cells": int(np.count_nonzero(inlet_mask)),
            "outlet_boundary_cells": int(np.count_nonzero(outlet_mask)),
            "target_flux_balance_relative_error": float(target_balance_error),
            "skipped_inlet_vessel_ids": tuple(skipped_inlets),
            "skipped_outlet_vessel_ids": tuple(skipped_outlets),
            "boundary_representation": "mask_boundary_open_edges",
        },
    )


def apply_open_boundary_fluxes(
    flux_x_um2_s: np.ndarray,
    flux_z_um2_s: np.ndarray,
    boundaries: BoundaryFluxFields,
    *,
    outlet_actual_by_label_um2_s: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return face fluxes with fixed inlet/outlet open-edge fluxes applied.

    The projection is allowed to correct only internal fluid-fluid faces.  Open
    boundary faces must keep their prescribed fluxes, so this function overwrites
    those exact face entries immediately before pressure projection.
    """

    # Work on copies so the caller can still inspect the unmodified internal
    # fluxes if needed for staged diagnostics.
    flux_x = np.asarray(flux_x_um2_s, dtype=np.float64).copy()
    flux_z = np.asarray(flux_z_um2_s, dtype=np.float64).copy()
    for row in range(boundaries.open_face_flux_um2_s.size):
        axis = int(boundaries.open_face_axis[row])
        i, j = (int(x) for x in boundaries.open_face_index_ij[row])
        signed_outflow = _realized_open_face_flux(
            boundaries, row, outlet_actual_by_label_um2_s
        )
        normal = boundaries.open_face_normal_xz[row]
        # Face arrays store flux in the positive grid-axis direction.  The
        # boundary profile stores signed flux relative to the outward normal.
        # Convert between these conventions using the face normal sign.
        if axis == 0:
            global_flux = signed_outflow if float(normal[0]) > 0.0 else -signed_outflow
            flux_x[i, j] = global_flux
        else:
            global_flux = signed_outflow if float(normal[1]) > 0.0 else -signed_outflow
            flux_z[i, j] = global_flux
    return flux_x, flux_z


def realized_open_boundary_flux_field(
    boundaries: BoundaryFluxFields,
    outlet_actual_by_label_um2_s: np.ndarray | None = None,
) -> np.ndarray:
    """Return the cell helper field with pressure-solved outlet totals."""

    field = np.zeros_like(boundaries.open_boundary_flux_um2_s, dtype=np.float64)
    for row, (i, j) in enumerate(boundaries.open_face_cell_ij):
        field[int(i), int(j)] += _realized_open_face_flux(
            boundaries, row, outlet_actual_by_label_um2_s
        )
    return field.astype(np.float32)


def _realized_open_face_flux(
    boundaries: BoundaryFluxFields,
    row: int,
    outlet_actual_by_label_um2_s: np.ndarray | None,
) -> float:
    target_face_flux = float(boundaries.open_face_flux_um2_s[row])
    if int(boundaries.open_face_kind[row]) < 0 or outlet_actual_by_label_um2_s is None:
        return target_face_flux
    label = int(boundaries.open_face_label[row])
    actual = np.asarray(outlet_actual_by_label_um2_s, dtype=np.float64)
    target = np.asarray(boundaries.outlet_target_by_label_um2_s, dtype=np.float64)
    if label < 1 or label >= actual.size or label >= target.size:
        raise ValueError(f"Missing realized pressure-outlet flux for label {label}.")
    denominator = float(target[label])
    if denominator <= np.finfo(np.float64).eps:
        raise ValueError(f"Outlet label {label} has a non-positive reference flux.")
    return target_face_flux * float(actual[label]) / denominator


def _boundary_profile(
    domain: GridDomain,
    raster: RasterizedVessels,
    face_catalog: BoundaryFaceCatalog,
    vessel: Vessel,
    *,
    q2d_um2_s: float,
    depth_cells: float,
    at_distal_end: bool,
) -> _BoundaryProfile:
    """
    From thousands of walls, accurately find that specific inlet or outlet 
    and arrange the water flow in a way that best conforms to the laws of real physics
    """
    p0 = np.asarray([vessel.x_p[0], vessel.x_p[2]], dtype=float)
    p1 = np.asarray([vessel.x_d[0], vessel.x_d[2]], dtype=float)
    segment = p1 - p0
    length = float(np.linalg.norm(segment))
    if length <= np.finfo(float).eps:
        return _empty_profile(float(q2d_um2_s))
    unit = segment / length

    # The outward normal points out of the computational lumen.  
    # At a distal outlet it follows flow direction; at a root inlet it points upstream, opposite to the parent-to-child direction.
    outward = unit if at_distal_end else -unit
    radius = max(float(vessel.radius), np.finfo(float).eps)
    spacing = float(domain.spacing_um)
    max_search_um = max(float(depth_cells) * spacing, 0.75 * spacing)
    endpoint = p1 if at_distal_end else p0
    tangent = np.asarray([-unit[1], unit[0]], dtype=float)

    center = face_catalog.center_xz_um
    if center.size == 0:
        return _empty_profile(float(q2d_um2_s))

    # Use a cheap bounding box first to keep the later vector filters focused on faces that could plausibly belong to this endpoint.
    bbox = radius + max_search_um + 2.0 * spacing
    near = (np.abs(center[:, 0] - endpoint[0]) <= bbox) & (np.abs(center[:, 1] - endpoint[1]) <= bbox)

    # alignment > 0 means the mask face points generally out of this endpoint.
    # A positive threshold rejects side-wall faces that are close in space but
    # should remain solid walls.
    alignment = face_catalog.normal_xz @ outward

    # Axial and lateral coordinates are measured in the local vessel frame.  The
    # endpoint plane is axial=0 for root inlets and axial=length for terminal
    # outlets.
    axial = (center - p0) @ unit
    lateral = np.abs((center - p0) @ tangent)
    target_axial = length if at_distal_end else 0.0

    # First pass: strict alignment and endpoint distance.  This usually selects
    # the flat cap produced by the geometry rasterizer.
    selected = near & (alignment > 0.15) & (np.abs(axial - target_axial) <= max_search_um) & (lateral <= radius + 0.75 * spacing)
    if not np.any(selected):
        # Fallback: allow a wider axial and lateral window.  This handles
        # slanted, staircase-like endpoints after rasterization without
        # selecting the whole terminal branch as a thick velocity band.
        selected = near & (alignment > 0.0) & (np.abs(axial - target_axial) <= 2.0 * max_search_um) & (lateral <= radius + spacing)
    if not np.any(selected):
        return _empty_profile(float(q2d_um2_s))

    selected_lateral = lateral[selected]

    # A parabolic lateral weight approximates a planar Poiseuille profile across
    # the open edge.  Multiplying by normal alignment reduces the contribution
    # of stair-step faces that point only partly along the endpoint normal.
    profile = np.maximum(0.0, 1.0 - (selected_lateral / max(radius, np.finfo(float).eps)) ** 2)
    profile *= np.maximum(alignment[selected], 0.0)
    if not np.any(profile > 0.0):
        # Degenerate raster geometry can make all parabolic weights vanish.
        # Uniform weights preserve the target flux instead of dropping the
        # boundary condition.
        profile = np.ones_like(profile, dtype=float)
    lengths = face_catalog.length_um[selected].astype(float, copy=False)

    selected_centers = face_catalog.center_xz_um[selected].astype(np.float64, copy=False)

    # Normalize by profile * length so integrating signed_flux over selected
    # faces gives exactly the requested q2d target before projection.
    denominator = float(np.sum(profile * lengths))
    if denominator <= np.finfo(float).eps:
        return _empty_profile(float(q2d_um2_s))

    # Positive signed flux means outflow through the outward normal.  Root
    # inlets are therefore negative because flow enters the computational domain.
    signed_target = float(q2d_um2_s) if at_distal_end else -float(q2d_um2_s)
    signed_flux = signed_target * profile * lengths / denominator

    return _BoundaryProfile(
        cell_ij=face_catalog.cell_ij[selected].astype(np.int32, copy=True),
        face_index_ij=face_catalog.face_index_ij[selected].astype(np.int32, copy=True),
        axis=face_catalog.axis[selected].astype(np.int8, copy=True),
        normal_xz=face_catalog.normal_xz[selected].astype(np.float32, copy=True),
        center_xz_um=selected_centers.astype(np.float64, copy=True),
        length_um=lengths.astype(np.float64, copy=True),
        weight=profile.astype(np.float64, copy=True),
        signed_flux_um2_s=signed_flux.astype(np.float64, copy=True),
        target_q2d_um2_s=float(q2d_um2_s),
    )


def enumerate_lumen_boundary_faces(domain: GridDomain, lumen_mask: np.ndarray) -> BoundaryFaceCatalog:
    """
    Enumerate every lumen face adjacent to non-lumen space.

    These faces are possible open boundaries, but most of them are solid walls.
    Later endpoint-specific filtering chooses the small subset that should be
    treated as root inlets or terminal outlets.
    """

    lumen = np.asarray(lumen_mask, dtype=bool)
    spacing = float(domain.spacing_um)
    cells: list[tuple[int, int]] = []
    face_indices: list[tuple[int, int]] = []
    axes: list[int] = []
    normals: list[tuple[float, float]] = []
    centers: list[tuple[float, float]] = []
    lengths: list[float] = []
    nx, nz = lumen.shape
    for i, j in np.argwhere(lumen):
        x = float(domain.x_coordinates_um[i])
        z = float(domain.z_coordinates_um[j])

        # For x-faces, face_index_ij stores the index in the flux_x array, which
        # has shape (nx + 1, nz).  A face at i separates cell i-1 and i.
        if i == 0 or not lumen[i - 1, j]:
            _append_face(cells, face_indices, axes, normals, centers, lengths, i, j, i, j, 0, (-1.0, 0.0), (x - 0.5 * spacing, z), spacing)
        if i == nx - 1 or not lumen[i + 1, j]:
            _append_face(cells, face_indices, axes, normals, centers, lengths, i, j, i + 1, j, 0, (1.0, 0.0), (x + 0.5 * spacing, z), spacing)

        # For z-faces, face_index_ij stores the index in the flux_z array, which
        # has shape (nx, nz + 1).  A face at j separates cell j-1 and j.
        if j == 0 or not lumen[i, j - 1]:
            _append_face(cells, face_indices, axes, normals, centers, lengths, i, j, i, j, 1, (0.0, -1.0), (x, z - 0.5 * spacing), spacing)
        if j == nz - 1 or not lumen[i, j + 1]:
            _append_face(cells, face_indices, axes, normals, centers, lengths, i, j, i, j + 1, 1, (0.0, 1.0), (x, z + 0.5 * spacing), spacing)

    if not cells:
        return BoundaryFaceCatalog(
            cell_ij=np.zeros((0, 2), dtype=np.int32),
            face_index_ij=np.zeros((0, 2), dtype=np.int32),
            axis=np.zeros(0, dtype=np.int8),
            normal_xz=np.zeros((0, 2), dtype=np.float32),
            center_xz_um=np.zeros((0, 2), dtype=np.float64),
            length_um=np.zeros(0, dtype=np.float64),
        )
    return BoundaryFaceCatalog(
        cell_ij=np.asarray(cells, dtype=np.int32),
        face_index_ij=np.asarray(face_indices, dtype=np.int32),
        axis=np.asarray(axes, dtype=np.int8),
        normal_xz=np.asarray(normals, dtype=np.float32),
        center_xz_um=np.asarray(centers, dtype=np.float64),
        length_um=np.asarray(lengths, dtype=np.float64),
    )


def _append_face(
    cells: list[tuple[int, int]],
    face_indices: list[tuple[int, int]],
    axes: list[int],
    normals: list[tuple[float, float]],
    centers: list[tuple[float, float]],
    lengths: list[float],
    cell_i: int,
    cell_j: int,
    face_i: int,
    face_j: int,
    axis: int,
    normal: tuple[float, float],
    center: tuple[float, float],
    length: float,
) -> None:
    """Append one candidate boundary face to the catalog lists."""

    cells.append((int(cell_i), int(cell_j)))
    face_indices.append((int(face_i), int(face_j)))
    axes.append(int(axis))
    normals.append(normal)
    centers.append(center)
    lengths.append(float(length))


def _accumulate_profile(
    *,
    profile: _BoundaryProfile,
    label: int,
    kind: int,
    label_grid: np.ndarray,
    velocity: np.ndarray,
    normal_sum: np.ndarray,
    weight: np.ndarray,
    edge_length: np.ndarray,
    open_flux: np.ndarray,
) -> None:
    """
    Accumulate one endpoint profile into grid-shaped helper arrays.

    Multiple selected faces can touch the same cell.  The explicit face arrays
    remain one row per face, while these helper grids combine values by cell for
    convenient visualization and NPZ output.
    """

    for row, (i, j) in enumerate(profile.cell_ij):
        signed_flux = float(profile.signed_flux_um2_s[row])
        length = max(float(profile.length_um[row]), np.finfo(float).eps)
        normal = profile.normal_xz[row].astype(float, copy=False)
        label_grid[int(i), int(j)] = int(label)

        # Convert signed flux per face length into a velocity-like helper
        # vector.  This is not used as a hard velocity overwrite in the solver.
        velocity[int(i), int(j)] += (signed_flux / length * normal).astype(np.float32)

        # Weight normals by absolute flux so the displayed normal reflects the
        # dominant open face if a cell has multiple selected faces.
        normal_sum[int(i), int(j)] += normal * abs(signed_flux)
        weight[int(i), int(j)] += float(profile.weight[row])
        edge_length[int(i), int(j)] += float(profile.length_um[row])
        open_flux[int(i), int(j)] += signed_flux


def _append_open_faces(
    profile: _BoundaryProfile,
    *,
    label: int,
    kind: int,
    open_cells: list[np.ndarray],
    open_face_indices: list[np.ndarray],
    open_axes: list[np.ndarray],
    open_normals: list[np.ndarray],
    open_centers: list[np.ndarray],
    open_lengths: list[np.ndarray],
    open_labels: list[np.ndarray],
    open_kinds: list[np.ndarray],
    open_fluxes: list[np.ndarray],
) -> None:
    """Append one endpoint profile to the authoritative row-wise face arrays."""

    n = int(profile.cell_ij.shape[0])
    open_cells.append(profile.cell_ij)
    open_face_indices.append(profile.face_index_ij)
    open_axes.append(profile.axis)
    open_normals.append(profile.normal_xz)
    open_centers.append(profile.center_xz_um)
    open_lengths.append(profile.length_um)
    open_labels.append(np.full(n, int(label), dtype=np.int32))
    open_kinds.append(np.full(n, int(kind), dtype=np.int8))
    open_fluxes.append(profile.signed_flux_um2_s)


def _append_open_section(
    profile: _BoundaryProfile,
    *,
    label: int,
    kind: int,
    section_points: list[np.ndarray],
    section_normals: list[np.ndarray],
    section_tangents: list[np.ndarray],
    section_half_widths: list[float],
    section_labels: list[int],
    section_kinds: list[int],
) -> None:
    """Append directed sections that exactly cover the selected open faces.

    A slanted raster endpoint is a staircase rather than one straight line.
    Treating its anatomical endpoint as a single plane can place the lifecycle
    plane upstream or downstream of the actual CFD opening.  The canonical
    representation therefore merges only contiguous, collinear faces.  A
    staircase may produce several small sections, but every one is exactly on
    the authoritative open face that it represents.
    """

    sections = derive_coplanar_open_sections(
        profile.face_index_ij,
        profile.axis,
        profile.normal_xz,
        profile.center_xz_um,
        profile.length_um,
        np.full(profile.axis.size, int(label), dtype=np.int32),
        np.full(profile.axis.size, int(kind), dtype=np.int8),
    )
    section_points.extend(sections["point_xz_um"])
    section_normals.extend(sections["normal_xz"])
    section_tangents.extend(sections["tangent_xz"])
    section_half_widths.extend(sections["half_width_um"].tolist())
    section_labels.extend(sections["label"].tolist())
    section_kinds.extend(sections["kind"].tolist())


def derive_coplanar_open_sections(
    face_index_ij: np.ndarray,
    face_axis: np.ndarray,
    face_normal_xz: np.ndarray,
    face_center_xz_um: np.ndarray,
    face_length_um: np.ndarray,
    face_label: np.ndarray,
    face_kind: np.ndarray,
) -> dict[str, np.ndarray]:
    """Merge contiguous open faces only when they lie on the same plane.

    The raster face index supplies an exact grouping key: X-normal faces share
    a plane when their first face index is equal, while Z-normal faces share a
    plane when their second face index is equal.  The returned finite sections
    are therefore an exact union of the supplied CFD open faces, not a nearby
    anatomical surrogate.
    """

    indices = np.asarray(face_index_ij, dtype=np.int32)
    axes = np.asarray(face_axis, dtype=np.int8)
    normals = np.asarray(face_normal_xz, dtype=np.float64)
    centers = np.asarray(face_center_xz_um, dtype=np.float64)
    lengths = np.asarray(face_length_um, dtype=np.float64)
    labels = np.asarray(face_label, dtype=np.int32)
    kinds = np.asarray(face_kind, dtype=np.int8)
    count = int(axes.size)
    expected = {
        "face_index_ij": (count, 2),
        "face_normal_xz": (count, 2),
        "face_center_xz_um": (count, 2),
        "face_length_um": (count,),
        "face_label": (count,),
        "face_kind": (count,),
    }
    actual = {
        "face_index_ij": indices.shape,
        "face_normal_xz": normals.shape,
        "face_center_xz_um": centers.shape,
        "face_length_um": lengths.shape,
        "face_label": labels.shape,
        "face_kind": kinds.shape,
    }
    for name, shape in expected.items():
        if actual[name] != shape:
            raise ValueError(f"{name} must have shape {shape}, got {actual[name]}.")
    if count == 0:
        return {
            "point_xz_um": np.empty((0, 2), dtype=np.float64),
            "normal_xz": np.empty((0, 2), dtype=np.float64),
            "tangent_xz": np.empty((0, 2), dtype=np.float64),
            "half_width_um": np.empty(0, dtype=np.float64),
            "label": np.empty(0, dtype=np.int32),
            "kind": np.empty(0, dtype=np.int8),
        }

    points_out: list[np.ndarray] = []
    normals_out: list[np.ndarray] = []
    tangents_out: list[np.ndarray] = []
    widths_out: list[float] = []
    labels_out: list[int] = []
    kinds_out: list[int] = []
    keys: dict[tuple[int, int, int, int, int], list[int]] = {}
    for row in range(count):
        axis = int(axes[row])
        if axis not in (0, 1):
            raise ValueError("Open-face axes must use 0 for X faces or 1 for Z faces.")
        normal_length = float(np.linalg.norm(normals[row]))
        if not np.isfinite(normal_length) or normal_length <= 0.0:
            raise ValueError("Every open face must have a finite non-zero normal.")
        unit_normal = normals[row] / normal_length
        normal_sign = int(np.sign(unit_normal[axis]))
        if normal_sign == 0 or abs(float(unit_normal[1 - axis])) > 1.0e-8:
            raise ValueError("An authoritative raster face normal must align with its face axis.")
        plane_index = int(indices[row, axis])
        key = (
            int(kinds[row]),
            int(labels[row]),
            axis,
            plane_index,
            normal_sign,
        )
        keys.setdefault(key, []).append(row)

    for key in sorted(keys):
        kind, label, axis, _, _ = key
        rows = np.asarray(keys[key], dtype=np.int64)
        normal = np.mean(normals[rows], axis=0)
        normal /= float(np.linalg.norm(normal))
        tangent = np.asarray([-normal[1], normal[0]], dtype=np.float64)
        projected_centers = centers[rows] @ tangent
        projected_half_lengths = 0.5 * lengths[rows]
        starts = projected_centers - projected_half_lengths
        ends = projected_centers + projected_half_lengths
        order = np.argsort(starts, kind="stable")
        scale = max(
            1.0,
            float(np.max(np.abs(centers[rows]), initial=0.0)),
            float(np.max(lengths[rows], initial=0.0)),
        )
        tolerance = 128.0 * np.finfo(np.float64).eps * scale
        run_start = float(starts[order[0]])
        run_end = float(ends[order[0]])
        plane_offset = float(np.mean(centers[rows] @ normal))
        for local_index in order[1:]:
            next_start = float(starts[local_index])
            next_end = float(ends[local_index])
            if next_start <= run_end + tolerance:
                run_end = max(run_end, next_end)
                continue
            midpoint = 0.5 * (run_start + run_end)
            points_out.append(normal * plane_offset + tangent * midpoint)
            normals_out.append(normal.copy())
            tangents_out.append(tangent.copy())
            widths_out.append(0.5 * (run_end - run_start))
            labels_out.append(label)
            kinds_out.append(kind)
            run_start = next_start
            run_end = next_end
        midpoint = 0.5 * (run_start + run_end)
        points_out.append(normal * plane_offset + tangent * midpoint)
        normals_out.append(normal.copy())
        tangents_out.append(tangent.copy())
        widths_out.append(0.5 * (run_end - run_start))
        labels_out.append(label)
        kinds_out.append(kind)

    return {
        "point_xz_um": np.ascontiguousarray(points_out, dtype=np.float64),
        "normal_xz": np.ascontiguousarray(normals_out, dtype=np.float64),
        "tangent_xz": np.ascontiguousarray(tangents_out, dtype=np.float64),
        "half_width_um": np.ascontiguousarray(widths_out, dtype=np.float64),
        "label": np.ascontiguousarray(labels_out, dtype=np.int32),
        "kind": np.ascontiguousarray(kinds_out, dtype=np.int8),
    }


def _empty_profile(target_q2d_um2_s: float) -> _BoundaryProfile:
    """Return a valid zero-face profile for skipped or degenerate endpoints."""

    return _BoundaryProfile(
        cell_ij=np.zeros((0, 2), dtype=np.int32),
        face_index_ij=np.zeros((0, 2), dtype=np.int32),
        axis=np.zeros(0, dtype=np.int8),
        normal_xz=np.zeros((0, 2), dtype=np.float32),
        center_xz_um=np.zeros((0, 2), dtype=np.float64),
        length_um=np.zeros(0, dtype=np.float64),
        weight=np.zeros(0, dtype=np.float64),
        signed_flux_um2_s=np.zeros(0, dtype=np.float64),
        target_q2d_um2_s=float(target_q2d_um2_s),
    )


def _concat_1d(parts: list[np.ndarray], *, dtype: object) -> np.ndarray:
    """Concatenate optional 1D profile arrays without special-casing callers."""

    if not parts:
        return np.zeros(0, dtype=dtype)
    return np.concatenate([np.asarray(part, dtype=dtype).reshape(-1) for part in parts]).astype(dtype, copy=False)


def _concat_2d(parts: list[np.ndarray]) -> np.ndarray:
    """Concatenate optional 2D profile arrays without special-casing callers."""

    if not parts:
        return np.zeros((0, 2), dtype=np.int32)
    return np.vstack(parts)


def _concat_2d_float(parts: list[np.ndarray]) -> np.ndarray:
    """Concatenate physical X-Z rows as float64, including the empty case."""

    if not parts:
        return np.zeros((0, 2), dtype=np.float64)
    return np.vstack([np.asarray(part, dtype=np.float64) for part in parts])


def _stack_section_vectors(parts: list[np.ndarray]) -> np.ndarray:
    """Stack one X-Z vector per physical open section."""

    if not parts:
        return np.zeros((0, 2), dtype=np.float64)
    return np.vstack([np.asarray(part, dtype=np.float64).reshape(1, 2) for part in parts])


def _root_vessels(vessels: list[Vessel] | tuple[Vessel, ...]) -> list[Vessel]:
    """Root vessels define inlet endpoints because they have no parent."""

    return [v for v in vessels if int(v.parent_id) < 0]


def _terminal_vessels(vessels: list[Vessel] | tuple[Vessel, ...]) -> list[Vessel]:
    """Terminal vessels define outlet endpoints because they have no children."""

    return [v for v in vessels if len(getattr(v, "children", [])) == 0]
