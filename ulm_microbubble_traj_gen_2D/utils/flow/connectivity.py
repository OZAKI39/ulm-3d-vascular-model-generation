"""
Connectivity checks for rasterized fluid domains.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy import ndimage

from ..runtime.console_output import print_error, print_key_values, print_success, print_warning
from ..core.types import GridDomain, RasterizedVessels

if TYPE_CHECKING:
    from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry


@dataclass(frozen=True)
class FluidConnectivityReport:
    """
    Summary of connected components inside the rasterized lumen mask.
    """
    n_components: int
    total_lumen_cells: int
    largest_component_cells: int
    disconnected_cells: int
    disconnected_fraction: float
    component_sizes: tuple[int, ...]
    disconnected_vessel_ids: tuple[int, ...]
    hole_cells: int


def validate_fluid_connectivity(
    raster: RasterizedVessels,
    domain: GridDomain,
    *,
    continuous_geometry: "ContinuousVesselGeometry | None" = None,
    max_reported_components: int = 8,
) -> FluidConnectivityReport:
    """
    Stop the run if the rasterized lumen is not one connected fluid domain.
    """
    # Use 4-connectivity to label connected components in the lumen mask.
    labels, n_components    = _label_lumen_components(raster.lumen_mask)
    component_sizes         = _component_sizes(labels, n_components)
    total_lumen_cells       = int(np.count_nonzero(raster.lumen_mask))
    if total_lumen_cells <= 0:
        raise ValueError("Connectivity check failed: the rasterized lumen mask is empty.")

    largest_component_cells = int(component_sizes[0]) if component_sizes else 0
    disconnected_cells      = int(total_lumen_cells - largest_component_cells)
    disconnected_fraction   = float(disconnected_cells / total_lumen_cells)
    disconnected_vessel_ids = _disconnected_vessel_ids(raster, labels, component_sizes)     # find the original vessel IDs of the disconnected components
    hole_cells = _hole_cell_count(
        raster.lumen_mask,
        domain=domain,
        continuous_geometry=continuous_geometry,
    )

    report = FluidConnectivityReport(
        n_components=int(n_components),
        total_lumen_cells=total_lumen_cells,
        largest_component_cells=largest_component_cells,
        disconnected_cells=disconnected_cells,
        disconnected_fraction=disconnected_fraction,
        component_sizes=tuple(int(x) for x in component_sizes[: int(max_reported_components)]),
        disconnected_vessel_ids=tuple(int(x) for x in disconnected_vessel_ids[: int(max_reported_components)]),
        hole_cells=int(hole_cells),
    )

    if int(n_components) != 1 or int(hole_cells) != 0:
        _print_connectivity_failure(report, domain)
        raise ValueError(_connectivity_error_message(report, domain))

    print_success(
        "Fluid connectivity passed: "
        f"one lumen domain, no holes, {total_lumen_cells} cells at {domain.spacing_um:.3f} um spacing."
    )
    return report


def _label_lumen_components(lumen_mask: np.ndarray) -> tuple[np.ndarray, int]:
    """
    two pixels must be "edge-to-edge" to be considered connected; 
    if they are only "corner-to-corner" (diagonally adjacent), the fluid cannot pass through and is considered disconnected.
    """
    # define the connectivity structure for 2D labeling (4-connectivity)
    structure = np.asarray([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
    
    # label the connected components in the lumen mask
    labels, n_components = ndimage.label(np.asarray(lumen_mask, dtype=bool), structure=structure)
    return labels.astype(np.int32, copy=False), int(n_components)


def _component_sizes(labels: np.ndarray, n_components: int) -> list[int]:
    if int(n_components) <= 0:
        return []
    counts = np.bincount(labels.ravel(), minlength=int(n_components) + 1)[1:]
    return sorted((int(x) for x in counts), reverse=True)


def _disconnected_vessel_ids(raster: RasterizedVessels, labels: np.ndarray, component_sizes: list[int]) -> list[int]:
    """
    Find all connected blocks except the largest one (main road), 
    and find out which blood vessels these disconnected blocks originally belonged to by comparing with the original data.
    """
    if len(component_sizes) <= 1:
        return []
    counts          = np.bincount(labels.ravel())
    largest_label   = int(np.argmax(counts[1:]) + 1)
    disconnected    = (labels > 0) & (labels != largest_label)
    vessel_ids      = np.unique(raster.vessel_id[disconnected])
    vessel_ids      = vessel_ids[vessel_ids >= 0]
    return sorted(int(x) for x in vessel_ids)


def _hole_cell_count(
    lumen_mask: np.ndarray,
    *,
    domain: GridDomain | None = None,
    continuous_geometry: "ContinuousVesselGeometry | None" = None,
) -> int:
    """
    Count enclosed missing cells whose centres belong to the continuous lumen.

    A concave exterior pocket can look enclosed after Cartesian sampling even
    though the authoritative continuous polygon places it outside the lumen.
    Such a cell is not an internal fluid hole and must not be added to the mask.
    """
    lumen = np.asarray(lumen_mask, dtype=bool)
    filled = ndimage.binary_fill_holes(lumen, structure=ndimage.generate_binary_structure(2, 1))
    hole_rows = np.argwhere(filled & ~lumen)
    if hole_rows.size == 0:
        return 0
    if continuous_geometry is None:
        return int(hole_rows.shape[0])
    if domain is None:
        raise ValueError(
            "domain is required when continuous_geometry is supplied."
        )
    points = np.column_stack(
        (
            np.asarray(domain.x_coordinates_um)[hole_rows[:, 0]],
            np.asarray(domain.z_coordinates_um)[hole_rows[:, 1]],
        )
    )
    inside = np.asarray(
        continuous_geometry.contains_xz_um(points), dtype=bool
    ).reshape(-1)
    return int(np.count_nonzero(inside))


def _print_connectivity_failure(report: FluidConnectivityReport, domain: GridDomain) -> None:
    print_error("Fluid connectivity check failed.")
    print_key_values(
        [
            ("Connected lumen domains", report.n_components),
            ("Total lumen cells", report.total_lumen_cells),
            ("Largest component cells", report.largest_component_cells),
            (
                "Disconnected lumen cells",
                f"{report.disconnected_cells} ({100.0 * report.disconnected_fraction:.2f}%)",
            ),
            ("Internal hole cells", report.hole_cells),
            ("Largest component sizes", list(report.component_sizes)),
            ("Disconnected vessel IDs", list(report.disconnected_vessel_ids)),
            ("Current grid spacing", f"{domain.spacing_um:.3f} um"),
        ]
    )
    print_warning(
        "Suggested fix: decrease domain.grid_spacing_um or increase "
        "domain.min_lumen_radius_cells."
    )


def _connectivity_error_message(report: FluidConnectivityReport, domain: GridDomain) -> str:
    return (
        "Rasterized lumen failed mask-quality checks. "
        f"Found {report.n_components} components and {report.hole_cells} hole cells at grid_spacing_um={domain.spacing_um:.3f}. "
        "The flow solve was stopped before velocity-field generation."
    )
