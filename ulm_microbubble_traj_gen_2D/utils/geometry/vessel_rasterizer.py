"""
Rasterize vessel corridors into a 2D X-Z lumen mask and physical metadata.

The rasterizer deliberately separates geometry construction from physical
attribute assignment. The continuous Shapely lumen is intersected exactly
with every candidate Cartesian cell, producing a physical area fraction before
the Boolean solver mask is derived. Vessel IDs, radii, flows, directions, and
seed velocities are assigned afterward without changing the accepted mask.
"""

import math
import os
import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import shapely
from scipy import ndimage
from shapely import make_valid
from shapely.geometry import Point
from shapely.ops import unary_union

from ulm_vascular_model_generator.utils.core.models import Vessel
from ..core.config import DomainConfig
from .continuous_vessel_geometry import ContinuousVesselGeometry
from ..runtime.console_output import print_warning
from ..core.types import GridDomain, RasterizedVessels


# Exact constructive intersections are substantially more expensive than
# prepared spatial predicates.  Keep tiles small enough to follow thin vessel
# corridors closely while amortizing the scalar tile-query overhead.
_EXACT_COVERAGE_TILE_EDGE_CELLS = 32
_EXACT_COVERAGE_MAX_WORKERS = 4
_coverage_thread_state = threading.local()


def rasterize_vessels(
    vessels: list[Vessel],
    domain: GridDomain,
    cfg: DomainConfig,
    *,
    effective_thickness_um: float,
    continuous_geometry: ContinuousVesselGeometry,
    dynamic_viscosity_mpas: float = 3.0,
) -> RasterizedVessels:
    """
    Convert vessel centerlines and radii to grid arrays used by the solver.

    The returned object provides the Cartesian sampling and diagnostic support
    for the boundary-fitted DOLFINx solution. A valid result contains exact
    cell-area fractions, a centre-valid Boolean mask, cell-wise vessel
    metadata, and continuous-wall distance/normal samples.
    """
    # ===================================================================================================
    # ======= Allocate all output arrays up front so every later step writes into one shared grid.  
    # ===================================================================================================
    # The mask is produced by continuous geometry first; 
    # the remaining arrays are filled only for cells that survive that mask step.
    nx, nz       = domain.shape
    if continuous_geometry.shape != tuple(int(value) for value in domain.shape):
        raise ValueError(
            "continuous_geometry must be built for the supplied domain."
        )
    vessel_id    = np.full((nx, nz), -1, dtype=np.int32)
    radius_um    = np.zeros((nx, nz), dtype=np.float32)
    flow_rate    = np.zeros((nx, nz), dtype=np.float32)
    q2d_flow     = np.zeros((nx, nz), dtype=np.float32)
    viscosity    = np.zeros((nx, nz), dtype=np.float32)
    direction    = np.zeros((nx, nz, 2), dtype=np.float32)
    
    distance_to_centerline = np.full((nx, nz), np.inf, dtype=np.float32)    # Distance from each cell to the nearest vessel centerline

    min_lumen_radius_cells = float(cfg.min_lumen_radius_cells)
    if not math.isfinite(min_lumen_radius_cells) or min_lumen_radius_cells < 0.0:
        raise ValueError(
            "domain.min_lumen_radius_cells must be finite and non-negative."
        )
    minimum_assignment_radius = min_lumen_radius_cells * float(domain.spacing_um)

    effective_thickness = float(effective_thickness_um)
    if not math.isfinite(effective_thickness) or effective_thickness <= 0.0:
        raise ValueError("effective_thickness_um must be finite and positive.")
    viscosity_value = float(dynamic_viscosity_mpas)
    if not math.isfinite(viscosity_value) or viscosity_value <= 0.0:
        raise ValueError("dynamic_viscosity_mpas must be finite and positive.")

    for vessel in vessels:
        flow = float(vessel.flow_rate)
        if not math.isfinite(flow) or flow < 0.0:
            raise ValueError(f"Vessel {vessel.vid} flow_rate must be finite and non-negative.")

    # ===================================================================================================
    # ======= Convert all blood vessel shapes into a pixel mask and assign metadata to each cell.
    # ===================================================================================================
    lumen_fraction = _exact_polygon_cell_coverage(
        continuous_geometry.lumen_polygon,
        domain,
    )
    lumen_mask = _center_validated_majority_mask(
        lumen_fraction,
        domain,
        continuous_geometry,
    )
    lumen_mask = _cleanup_lumen_mask(lumen_mask)
    junction_core_mask = _continuous_junction_core_mask(
        vessels, domain, lumen_mask
    )

    # ===================================================================================================

    for vessel in vessels:
        # Assign metadata to existing lumen cells.  This pass must not add new
        # lumen cells; otherwise the exact-coverage Boolean mask
        # would be replaced by a simpler centerline-distance rasterization.
        _rasterize_single_vessel(
            vessel=vessel,
            domain=domain,
            minimum_assignment_radius_um=minimum_assignment_radius,
            effective_thickness_um=effective_thickness,
            dynamic_viscosity_mpas=viscosity_value,
            lumen_mask=lumen_mask,
            vessel_id=vessel_id,
            radius_um=radius_um,
            flow_rate=flow_rate,
            q2d_flow=q2d_flow,
            viscosity=viscosity,
            direction=direction,
            distance_to_centerline=distance_to_centerline,
        )

    if not np.any(lumen_mask):
        raise ValueError(
            "Exact Shapely cell coverage produced an empty lumen mask; "
            "decrease grid_spacing_um."
        )

    # Junction polygons can create cells that belong to the fluid domain but do
    # not lie close enough to any single segment centerline to receive metadata
    # in the nearest-segment pass.  They still need flow, radius, direction, and
    # viscosity values for initialization and particle output, so they inherit
    # the nearest already-assigned lumen cell.
    _fill_added_lumen_cells(
        lumen_mask=lumen_mask,
        vessel_id=vessel_id,
        radius_um=radius_um,
        flow_rate=flow_rate,
        q2d_flow=q2d_flow,
        viscosity=viscosity,
        direction=direction,
        distance_to_centerline=distance_to_centerline,
    )
    _warn_if_under_resolved(vessels, domain, cfg)

    # Sample distance and inward normal from the same continuous wall used by
    # DOLFINx and particle contact. No distance transform or pixel-gradient
    # reconstruction is used.
    distance_to_wall, wall_normal, wall_mask = _continuous_wall_fields_on_grid(
        domain,
        lumen_mask,
        continuous_geometry,
    )
    distance_to_centerline[~lumen_mask] = np.nan

    return RasterizedVessels(
        lumen_mask=lumen_mask,
        wall_mask=wall_mask,
        vessel_id=vessel_id,
        radius_um=radius_um,
        flow_rate_um3_s=flow_rate,
        q2d_flow_um2_s=q2d_flow,
        viscosity_mpas=viscosity,
        direction_xz=direction,
        distance_to_centerline_um=distance_to_centerline,
        distance_to_wall_um=distance_to_wall.astype(np.float32),
        wall_normal_xz=wall_normal.astype(np.float32),
        lumen_fraction=lumen_fraction,
        junction_core_mask=junction_core_mask,
    )


def _continuous_junction_core_mask(
    vessels: list[Vessel],
    domain: GridDomain,
    lumen_mask: np.ndarray,
) -> np.ndarray:
    """Sample grid-independent physical junction disks for diagnostics only."""

    open_keys = set()
    for vessel in vessels:
        if int(vessel.parent_id) < 0:
            open_keys.add(
                (round(float(vessel.x_p[0]), 6), round(float(vessel.x_p[2]), 6))
            )
        if not getattr(vessel, "children", []):
            open_keys.add(
                (round(float(vessel.x_d[0]), 6), round(float(vessel.x_d[2]), 6))
            )

    radii_by_key: dict[tuple[float, float], list[float]] = {}
    for vessel in vessels:
        for point in (vessel.x_p, vessel.x_d):
            key = (round(float(point[0]), 6), round(float(point[2]), 6))
            radii_by_key.setdefault(key, []).append(float(vessel.radius))
    disks = [
        Point(key).buffer(max(radii), quad_segs=16)
        for key, radii in radii_by_key.items()
        if key not in open_keys and len(radii) >= 2
    ]
    if not disks:
        return np.zeros(domain.shape, dtype=bool)
    polygon = unary_union(disks)
    if not polygon.is_valid:
        polygon = make_valid(polygon)
    coverage = _exact_polygon_cell_coverage(polygon, domain)
    return np.asarray(
        (coverage >= 0.5) & np.asarray(lumen_mask, dtype=bool),
        dtype=bool,
    )


def _exact_polygon_cell_coverage(
    polygon,
    domain: GridDomain,
) -> np.ndarray:
    """Return exact Shapely intersection area divided by Cartesian cell area.

    A polygon's global bounding box is a poor candidate region for a sparse
    vascular tree: it includes most of the solid tissue between branches.
    Prepared tile and cell predicates conservatively reject those zero-area
    cells before the expensive constructive intersections are evaluated.
    Every cell that touches the polygon still goes through the original exact
    GEOS intersection, preserving the physical area-fraction contract.
    """

    shape = tuple(int(value) for value in domain.shape)
    spacing = float(domain.spacing_um)
    if len(shape) != 2 or min(shape) < 1:
        raise ValueError("Exact lumen coverage requires a non-empty 2D grid.")
    if not math.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("domain.spacing_um must be finite and positive.")
    if polygon is None or polygon.is_empty:
        return np.zeros(shape, dtype=np.float32)
    if not polygon.is_valid:
        polygon = make_valid(polygon)
    if polygon.is_empty:
        return np.zeros(shape, dtype=np.float32)
    shapely.prepare(polygon)

    coverage = np.zeros(shape, dtype=np.float64)
    half = 0.5 * spacing
    minimum_x, minimum_z, maximum_x, maximum_z = polygon.bounds
    x = np.asarray(domain.x_coordinates_um, dtype=np.float64)
    z = np.asarray(domain.z_coordinates_um, dtype=np.float64)
    ix = np.flatnonzero(
        (x + half >= minimum_x) & (x - half <= maximum_x)
    )
    iz = np.flatnonzero(
        (z + half >= minimum_z) & (z - half <= maximum_z)
    )
    if ix.size == 0 or iz.size == 0:
        return coverage.astype(np.float32)

    cell_area = spacing * spacing
    tile_edge = _EXACT_COVERAGE_TILE_EDGE_CELLS
    tile_tasks = []
    for ix_start in range(0, int(ix.size), tile_edge):
        tile_ix = ix[ix_start : ix_start + tile_edge]
        tile_minimum_x = float(x[tile_ix[0]] - half)
        tile_maximum_x = float(x[tile_ix[-1]] + half)
        for iz_start in range(0, int(iz.size), tile_edge):
            tile_iz = iz[iz_start : iz_start + tile_edge]
            tile = shapely.box(
                tile_minimum_x,
                float(z[tile_iz[0]] - half),
                tile_maximum_x,
                float(z[tile_iz[-1]] + half),
            )
            if not bool(shapely.intersects(polygon, tile)):
                continue
            tile_tasks.append(
                (
                    x,
                    z,
                    tile_ix,
                    tile_iz,
                    half,
                    cell_area,
                )
            )

    worker_count = min(
        _EXACT_COVERAGE_MAX_WORKERS,
        os.cpu_count() or 1,
        len(tile_tasks),
    )
    if worker_count <= 1:
        tile_results = (
            _exact_coverage_for_tile((polygon, *task))
            for task in tile_tasks
        )
        for tile_ix, tile_iz, local_coverage in tile_results:
            coverage[np.ix_(tile_ix, tile_iz)] = local_coverage
    else:
        with ThreadPoolExecutor(
            max_workers=worker_count,
            thread_name_prefix="vessel-coverage",
            initializer=_initialize_exact_coverage_worker,
            initargs=(shapely.to_wkb(polygon),),
        ) as executor:
            for tile_ix, tile_iz, local_coverage in executor.map(
                _exact_coverage_for_worker_tile,
                tile_tasks,
            ):
                coverage[np.ix_(tile_ix, tile_iz)] = local_coverage
    return np.ascontiguousarray(coverage, dtype=np.float32)


def _initialize_exact_coverage_worker(polygon_wkb: bytes) -> None:
    """Give every GEOS worker its own prepared polygon instance."""

    polygon = shapely.from_wkb(polygon_wkb)
    shapely.prepare(polygon)
    _coverage_thread_state.polygon = polygon


def _exact_coverage_for_worker_tile(
    task,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate one tile with the calling worker's private GEOS geometry."""

    return _exact_coverage_for_tile(
        (_coverage_thread_state.polygon, *task)
    )


def _exact_coverage_for_tile(
    task,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate exact fractions for polygon-intersecting cells in one tile."""

    polygon, x, z, tile_ix, tile_iz, half, cell_area = task
    center_x, center_z = np.meshgrid(
        x[tile_ix],
        z[tile_iz],
        indexing="ij",
    )
    cells = shapely.box(
        center_x - half,
        center_z - half,
        center_x + half,
        center_z + half,
    )
    intersects = np.asarray(
        shapely.intersects(polygon, cells),
        dtype=bool,
    )
    if not np.any(intersects):
        return tile_ix, tile_iz, np.zeros(cells.shape, dtype=np.float64)

    intersection_area = np.asarray(
        shapely.area(shapely.intersection(cells[intersects], polygon)),
        dtype=np.float64,
    )
    local_coverage = np.zeros(cells.shape, dtype=np.float64)
    local_coverage[intersects] = np.clip(
        intersection_area / cell_area,
        0.0,
        1.0,
    )
    return tile_ix, tile_iz, local_coverage


def _continuous_wall_fields_on_grid(
    domain: GridDomain,
    lumen_mask: np.ndarray,
    continuous_geometry: ContinuousVesselGeometry,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample exact solid-wall distance and inward normal at lumen centres."""

    lumen = np.asarray(lumen_mask, dtype=bool)
    if lumen.shape != tuple(domain.shape):
        raise ValueError("lumen_mask must match domain.shape.")
    distance = np.zeros(domain.shape, dtype=np.float64)
    inward_normal = np.zeros((*domain.shape, 2), dtype=np.float64)
    indices = np.argwhere(lumen)
    chunk_size = 100_000
    for start in range(0, int(indices.shape[0]), chunk_size):
        rows = indices[start : start + chunk_size]
        points = np.column_stack(
            (
                domain.x_coordinates_um[rows[:, 0]],
                domain.z_coordinates_um[rows[:, 1]],
            )
        )
        state = continuous_geometry.exact_solid_wall_state_xz_um_accelerated(
            points
        )
        sampled_distance = np.asarray(
            state.distance_um, dtype=np.float64
        ).reshape(-1)
        sampled_normal = np.asarray(
            state.inward_normal_xz, dtype=np.float64
        ).reshape(-1, 2)
        distance[rows[:, 0], rows[:, 1]] = sampled_distance
        inward_normal[rows[:, 0], rows[:, 1]] = sampled_normal

    # A Cartesian cell belongs to the solid-wall overlay when its centre is no
    # farther than the cell half-diagonal from the continuous solid wall.
    half_diagonal = math.sqrt(0.5) * float(domain.spacing_um)
    tolerance = 64.0 * np.finfo(np.float64).eps * max(
        1.0,
        float(np.max(np.abs(domain.origin_um))),
        float(domain.spacing_um),
    )
    wall_mask = lumen & (distance <= half_diagonal + tolerance)
    return distance, inward_normal, wall_mask


def _center_validated_majority_mask(
    lumen_fraction: np.ndarray,
    domain: GridDomain,
    continuous_geometry: ContinuousVesselGeometry,
) -> np.ndarray:
    """Keep majority-covered cells whose sampling centre is inside the lumen."""

    fraction = np.asarray(lumen_fraction, dtype=np.float64)
    if fraction.shape != tuple(domain.shape):
        raise ValueError("lumen_fraction must match domain.shape.")
    candidates = np.argwhere(fraction >= 0.5)
    mask = np.zeros(domain.shape, dtype=bool)
    chunk_size = 100_000
    for start in range(0, int(candidates.shape[0]), chunk_size):
        rows = candidates[start : start + chunk_size]
        points = np.column_stack(
            (
                domain.x_coordinates_um[rows[:, 0]],
                domain.z_coordinates_um[rows[:, 1]],
            )
        )
        inside = np.asarray(
            continuous_geometry.contains_xz_um(points),
            dtype=bool,
        ).reshape(-1)
        accepted = rows[inside]
        mask[accepted[:, 0], accepted[:, 1]] = True
    return mask


def _cleanup_lumen_mask(lumen_mask: np.ndarray) -> np.ndarray:
    """Remove tiny fragments without adding cells rejected by exact geometry."""

    lumen = np.asarray(lumen_mask, dtype=bool).copy()
    if not np.any(lumen):
        return lumen

    # Do not call binary_fill_holes here.  A one-cell exterior pocket beside a
    # concave vessel junction can be enclosed in the Cartesian mask even though
    # its centre lies outside the authoritative continuous lumen.  Filling such
    # a pocket breaks the centre-valid mask contract and later makes DOLFINx
    # reject the added sampling point.
    labels, count = ndimage.label(lumen, structure=np.ones((3, 3), dtype=bool))
    if count <= 1:
        return lumen

    # Remove only tiny isolated fragments.  Larger disconnected components are
    # intentionally left for the connectivity validator to catch and stop.
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    tiny_labels = np.flatnonzero((sizes <= 4) & (np.arange(sizes.size) > 0))
    if tiny_labels.size:
        lumen[np.isin(labels, tiny_labels)] = False
    return lumen


def _rasterize_single_vessel(
    vessel: Vessel,
    domain: GridDomain,
    minimum_assignment_radius_um: float,
    effective_thickness_um: float,
    dynamic_viscosity_mpas: float,
    lumen_mask: np.ndarray,
    vessel_id: np.ndarray,
    radius_um: np.ndarray,
    flow_rate: np.ndarray,
    q2d_flow: np.ndarray,
    viscosity: np.ndarray,
    direction: np.ndarray,
    distance_to_centerline: np.ndarray,
) -> None:
    """
    Rasterize one vessel segment and keep the nearest segment per cell.

    Despite the name, this function no longer decides whether a cell is lumen.
    The lumen decision has already been made by the continuous geometry and
    exact-coverage stage. This function only assigns the nearest vessel's
    physical attributes to cells that are already inside the accepted lumen.
    """
    # Compute the properties of the vessel segment in the X-Z plane.
    p0      = np.asarray([vessel.x_p[0], vessel.x_p[2]], dtype=float)
    p1      = np.asarray([vessel.x_d[0], vessel.x_d[2]], dtype=float)
    seg     = p1 - p0
    length  = float(np.linalg.norm(seg))

    unit = seg / length
    physical_radius = float(vessel.radius)
    assignment_radius = max(physical_radius, minimum_assignment_radius_um)

    # Determine the bounding box of the vessel segment in grid coordinates, expanded by the mask radius.
    # Restricting work to this box keeps the nearest-centerline assignment fast
    # even on million-cell grids.
    ix0 = max(0, int(math.floor((min(p0[0], p1[0]) - assignment_radius - domain.origin_um[0]) / domain.spacing_um)))
    ix1 = min(domain.shape[0] - 1, int(math.ceil((max(p0[0], p1[0]) + assignment_radius - domain.origin_um[0]) / domain.spacing_um)))
    iz0 = max(0, int(math.floor((min(p0[1], p1[1]) - assignment_radius - domain.origin_um[2]) / domain.spacing_um)))
    iz1 = min(domain.shape[1] - 1, int(math.ceil((max(p0[1], p1[1]) + assignment_radius - domain.origin_um[2]) / domain.spacing_um)))
    if ix1 < ix0 or iz1 < iz0:
        return

    # Compute the distance from each grid cell in the bounding box to the vessel centerline, 
    # and update the raster arrays for cells that are closer than any previous vessel segment.
    xs = domain.x_coordinates_um[ix0 : ix1 + 1]
    zs = domain.z_coordinates_um[iz0 : iz1 + 1]
    xx, zz = np.meshgrid(xs, zs, indexing="ij")
    points = np.stack([xx, zz], axis=-1)

    rel = points - p0
    t = np.clip((rel[..., 0] * seg[0] + rel[..., 1] * seg[1]) / (length * length), 0.0, 1.0)
    closest = p0 + t[..., None] * seg
    dist = np.linalg.norm(points - closest, axis=-1)
    view = np.s_[ix0 : ix1 + 1, iz0 : iz1 + 1]
    local_lumen = lumen_mask[view]

    # Do not use this candidate mask to create new lumen cells.  It only limits
    # which accepted lumen cells may be attributed to this vessel segment.
    candidate = (dist <= assignment_radius) & local_lumen
    if not np.any(candidate):
        return

    # Overlapping corridors are expected at junctions and nearby branches.  The
    # nearest centerline wins so every ordinary lumen cell has one stable vessel
    # ID and one set of physical parameters.
    nearest = candidate & (dist < distance_to_centerline[view])
    if not np.any(nearest):
        return

    local_vid = vessel_id[view]
    local_radius = radius_um[view]
    local_flow = flow_rate[view]
    local_q2d_flow = q2d_flow[view]
    local_viscosity = viscosity[view]
    local_direction = direction[view]
    local_distance = distance_to_centerline[view]

    local_vid[nearest] = int(vessel.vid)
    local_radius[nearest] = physical_radius
    local_flow[nearest] = float(vessel.flow_rate)

    # Convert the exported 3D flow rate to the 2D flux used by the planar flow
    # solver.  The effective thickness is a modeling parameter, not a geometric
    # change to the lumen mask.
    local_q2d_flow[nearest] = float(vessel.flow_rate) / effective_thickness_um
    local_viscosity[nearest] = float(dynamic_viscosity_mpas)
    local_direction[nearest] = unit

    local_distance[nearest] = dist[nearest]


def _fill_added_lumen_cells(
    *,
    lumen_mask: np.ndarray,
    vessel_id: np.ndarray,
    radius_um: np.ndarray,
    flow_rate: np.ndarray,
    q2d_flow: np.ndarray,
    viscosity: np.ndarray,
    direction: np.ndarray,
    distance_to_centerline: np.ndarray,
) -> None:
    """Assign inherited metadata to lumen cells created only by junction filling."""

    added = lumen_mask & (vessel_id < 0)
    if not np.any(added):
        return
    source = vessel_id >= 0
    if not np.any(source):
        return

    # Distance-transform indices provide the nearest already-attributed source
    # cell for every unassigned lumen cell.  This keeps junction cells usable by
    # the flow and particle code without changing the mask or inventing a new
    # vessel object.
    _, nearest = ndimage.distance_transform_edt(~source, return_indices=True)
    ni = nearest[0][added]
    nj = nearest[1][added]
    vessel_id[added] = vessel_id[ni, nj]
    radius_um[added] = radius_um[ni, nj]
    flow_rate[added] = flow_rate[ni, nj]
    q2d_flow[added] = q2d_flow[ni, nj]
    viscosity[added] = viscosity[ni, nj]
    direction[added] = direction[ni, nj]
    distance_to_centerline[added] = distance_to_centerline[ni, nj]


def _warn_if_under_resolved(vessels: list[Vessel], domain: GridDomain, cfg: DomainConfig) -> None:
    """Warn when the physical graph is too thin for the chosen grid spacing."""

    min_diameter_cells = (
        2.0 * min(float(v.radius) for v in vessels) / float(domain.spacing_um)
    )
    required = float(cfg.min_resolved_diameter_cells)
    if not math.isfinite(required) or required < 0.0:
        raise ValueError(
            "domain.min_resolved_diameter_cells must be finite and non-negative."
        )
    if required > 0.0 and min_diameter_cells < required:
        # This is intentionally a warning, not a mask-radius inflation.  The
        # physically correct response for wall-shear-sensitive runs is a finer
        # grid or a smaller ROI, not pretending the DCCO radius is larger.
        print_warning(
            "The smallest physical vessel diameter spans "
            f"{min_diameter_cells:.2f} grid cells, below domain.min_resolved_diameter_cells={required:.2f}. "
            "Decrease domain.grid_spacing_um for wall-shear-sensitive runs."
        )
