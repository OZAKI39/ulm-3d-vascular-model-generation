"""Interactive viewer for the meshes used by the hybrid velocity scheme.

The viewer deliberately concentrates on discretisation:

* Cartesian cells and how the continuous lumen boundary cuts them;
* DOLFINx triangles, their hybrid-region role, size, and shape quality;
* the continuous solid wall as the geometric reference.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Mapping

import numpy as np


PACKAGE_DIR = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PACKAGE_DIR / "results"
FIELD_FILENAME = "velocity_and_wall_shear_field.npz"

_VELOCITY_REGION_NAMES = {
    0: "outside",
    1: "finite element",
    2: "transition",
    3: "regular grid",
}
_CARTESIAN_CELL_NAMES = {
    0: "outside",
    1: "cut, centre outside",
    2: "cut, centre inside",
    3: "full lumen",
}
_TRIANGLE_SIZE_NAMES = {
    1: "finer than grid",
    2: "similar to grid",
    3: "coarser than grid",
}


@dataclass(frozen=True)
class HybridVelocityQualityData:
    """Validated Cartesian mesh, DOLFINx mesh, and continuous wall."""

    source_path: Path
    x_coordinates_um: np.ndarray
    z_coordinates_um: np.ndarray
    spacing_um: float
    cartesian_fields: Mapping[str, np.ndarray]
    finite_element_fields: Mapping[str, np.ndarray]
    continuous_wall_start_xz_um: np.ndarray
    continuous_wall_end_xz_um: np.ndarray
    finite_element_cell_vertices_xz_um: np.ndarray
    finite_element_distance_um: float
    regular_grid_distance_um: float

    @property
    def shape(self) -> tuple[int, int]:
        return (
            int(self.x_coordinates_um.size),
            int(self.z_coordinates_um.size),
        )

    @property
    def extent(self) -> tuple[float, float, float, float]:
        half = 0.5 * self.spacing_um
        return (
            float(self.x_coordinates_um[0] - half),
            float(self.x_coordinates_um[-1] + half),
            float(self.z_coordinates_um[0] - half),
            float(self.z_coordinates_um[-1] + half),
        )


@dataclass(frozen=True)
class _Layer:
    key: str
    label: str
    unit: str
    values: np.ndarray
    cmap: str
    mesh: str
    fixed_limits: tuple[float, float] | None = None
    category_names: Mapping[int, str] | None = None


def _resolve_field_archive(path: str | Path | None) -> Path:
    if path is not None:
        candidate = Path(path).expanduser().resolve()
        if candidate.is_dir():
            candidate = candidate / FIELD_FILENAME
        if not candidate.is_file():
            raise FileNotFoundError(f"Field archive does not exist: {candidate}")
        return candidate

    candidates = sorted(
        DEFAULT_RESULTS_DIR.glob(f"*/{FIELD_FILENAME}"),
        key=lambda item: item.stat().st_mtime_ns,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            f"No {FIELD_FILENAME} was found below {DEFAULT_RESULTS_DIR}."
        )
    return candidates[0].resolve()


def _required_array(
    archive: np.lib.npyio.NpzFile,
    key: str,
    *,
    ndim: int | None = None,
) -> np.ndarray:
    if key not in archive.files:
        raise ValueError(f"Field archive is missing required array {key!r}.")
    value = np.asarray(archive[key])
    if ndim is not None and value.ndim != ndim:
        raise ValueError(
            f"Field array {key!r} must have {ndim} dimensions, found {value.shape}."
        )
    return value


def _grid_array(
    archive: np.lib.npyio.NpzFile,
    key: str,
    shape: tuple[int, int],
) -> np.ndarray:
    value = _required_array(archive, key)
    if value.shape != shape:
        raise ValueError(
            f"Field array {key!r} has shape {value.shape}; expected {shape}."
        )
    return np.ascontiguousarray(value)


def _point_to_segments_distance(
    points_xz_um: np.ndarray,
    segment_start_xz_um: np.ndarray,
    segment_end_xz_um: np.ndarray,
) -> np.ndarray:
    """Return exact distances to the nearest continuous-wall segment."""

    import shapely
    from shapely.strtree import STRtree

    points = np.asarray(points_xz_um, dtype=np.float64)
    starts = np.asarray(segment_start_xz_um, dtype=np.float64)
    ends = np.asarray(segment_end_xz_um, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("Distance points must have shape (N, 2).")
    if starts.ndim != 2 or starts.shape[1] != 2 or starts.shape != ends.shape:
        raise ValueError("Wall start/end arrays must both have shape (N, 2).")
    if starts.shape[0] == 0:
        raise ValueError("At least one continuous solid-wall segment is required.")

    segment_coordinates = np.stack((starts, ends), axis=1)
    wall_segments = shapely.linestrings(segment_coordinates)
    point_geometries = shapely.points(points)
    _, distances = STRtree(wall_segments).query_nearest(
        point_geometries,
        all_matches=False,
        return_distance=True,
    )
    return np.ascontiguousarray(distances, dtype=np.float64)


def _cartesian_mesh_fields(
    lumen_mask: np.ndarray,
    lumen_fraction: np.ndarray,
    hybrid_region: np.ndarray,
    finite_element_weight: np.ndarray,
    distance_to_wall_um: np.ndarray,
) -> dict[str, np.ndarray]:
    lumen = np.asarray(lumen_mask, dtype=bool)
    fraction = np.asarray(lumen_fraction, dtype=np.float64)
    tolerance = 1.0e-6
    covered = fraction > tolerance
    cut = covered & (fraction < 1.0 - tolerance)

    cell_type = np.zeros(lumen.shape, dtype=np.uint8)
    cell_type[cut & ~lumen] = 1
    cell_type[cut & lumen] = 2
    cell_type[covered & ~cut] = 3

    return {
        "cartesian_cell_type": cell_type,
        "cartesian_hybrid_region": np.ascontiguousarray(
            hybrid_region, dtype=np.uint8
        ),
        "cartesian_finite_element_weight": np.ascontiguousarray(
            finite_element_weight, dtype=np.float64
        ),
        "cartesian_lumen_fraction": np.ascontiguousarray(
            fraction, dtype=np.float64
        ),
        "cartesian_wall_distance_um": np.ascontiguousarray(
            distance_to_wall_um, dtype=np.float64
        ),
        "lumen_mask": np.ascontiguousarray(lumen),
    }


def _finite_element_mesh_fields(
    triangles_xz_um: np.ndarray,
    wall_start_xz_um: np.ndarray,
    wall_end_xz_um: np.ndarray,
    finite_element_distance_um: float,
    regular_grid_distance_um: float,
    cartesian_spacing_um: float,
) -> dict[str, np.ndarray]:
    triangles = np.asarray(triangles_xz_um, dtype=np.float64)
    if triangles.ndim != 3 or triangles.shape[1:] != (3, 2):
        raise ValueError(
            "fem_cell_vertices_xz_um must have shape (cell_count, 3, 2)."
        )
    if triangles.shape[0] == 0:
        raise ValueError("The DOLFINx mesh must contain at least one triangle.")

    edge_vectors = np.roll(triangles, -1, axis=1) - triangles
    edge_lengths = np.linalg.norm(edge_vectors, axis=2)
    twice_area = np.abs(
        edge_vectors[:, 0, 0] * edge_vectors[:, 1, 1]
        - edge_vectors[:, 0, 1] * edge_vectors[:, 1, 0]
    )
    area = 0.5 * twice_area
    edge_square_sum = np.sum(edge_lengths * edge_lengths, axis=1)
    quality = np.divide(
        4.0 * math.sqrt(3.0) * area,
        edge_square_sum,
        out=np.zeros_like(area),
        where=edge_square_sum > 0.0,
    )

    vertex_distances = _point_to_segments_distance(
        triangles.reshape(-1, 2),
        wall_start_xz_um,
        wall_end_xz_um,
    ).reshape(-1, 3)
    centroid = np.mean(triangles, axis=1)
    centroid_distance = _point_to_segments_distance(
        centroid,
        wall_start_xz_um,
        wall_end_xz_um,
    )

    region = np.full(triangles.shape[0], 2, dtype=np.uint8)
    region[np.all(vertex_distances <= finite_element_distance_um, axis=1)] = 1
    region[np.all(vertex_distances >= regular_grid_distance_um, axis=1)] = 3

    maximum_edge = np.max(edge_lengths, axis=1)
    size_class = np.full(triangles.shape[0], 2, dtype=np.uint8)
    size_class[maximum_edge < 0.8 * cartesian_spacing_um] = 1
    size_class[maximum_edge > 1.25 * cartesian_spacing_um] = 3

    return {
        "fem_triangle_region": region,
        "fem_triangle_size_class": size_class,
        "fem_triangle_area_um2": np.ascontiguousarray(area),
        "fem_triangle_minimum_edge_um": np.ascontiguousarray(
            np.min(edge_lengths, axis=1)
        ),
        "fem_triangle_maximum_edge_um": np.ascontiguousarray(maximum_edge),
        "fem_triangle_shape_quality": np.ascontiguousarray(quality),
        "fem_triangle_wall_distance_um": np.ascontiguousarray(centroid_distance),
        "fem_triangle_centroid_xz_um": np.ascontiguousarray(centroid),
    }


def _assemble_data(
    *,
    source_path: Path,
    x_coordinates_um: np.ndarray,
    z_coordinates_um: np.ndarray,
    spacing_um: float,
    lumen_mask: np.ndarray,
    lumen_fraction: np.ndarray,
    hybrid_region: np.ndarray,
    finite_element_weight: np.ndarray,
    distance_to_wall_um: np.ndarray,
    wall_start_xz_um: np.ndarray,
    wall_end_xz_um: np.ndarray,
    triangles_xz_um: np.ndarray,
    finite_element_distance_um: float,
    regular_grid_distance_um: float,
) -> HybridVelocityQualityData:
    x = np.asarray(x_coordinates_um, dtype=np.float64)
    z = np.asarray(z_coordinates_um, dtype=np.float64)
    spacing = float(spacing_um)
    if x.ndim != 1 or z.ndim != 1 or x.size < 2 or z.size < 2:
        raise ValueError("Mesh viewing requires one-dimensional coordinates and a 2 x 2 grid.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(z)):
        raise ValueError("Grid coordinates must be finite.")
    if np.any(np.diff(x) <= 0.0) or np.any(np.diff(z) <= 0.0):
        raise ValueError("Grid coordinates must be strictly increasing.")
    if not math.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("spacing_um must be finite and positive.")
    tolerance = 1.0e-7 * max(spacing, 1.0)
    if not np.allclose(np.diff(x), spacing, rtol=0.0, atol=tolerance):
        raise ValueError("X coordinates are inconsistent with spacing_um.")
    if not np.allclose(np.diff(z), spacing, rtol=0.0, atol=tolerance):
        raise ValueError("Z coordinates are inconsistent with spacing_um.")

    shape = (int(x.size), int(z.size))
    grid_inputs = {
        "lumen_mask": lumen_mask,
        "lumen_fraction": lumen_fraction,
        "hybrid_velocity_region": hybrid_region,
        "hybrid_finite_element_weight": finite_element_weight,
        "distance_to_wall_um": distance_to_wall_um,
    }
    for key, value in grid_inputs.items():
        if np.asarray(value).shape != shape:
            raise ValueError(
                f"Field array {key!r} has shape {np.asarray(value).shape}; "
                f"expected {shape}."
            )

    starts = np.asarray(wall_start_xz_um, dtype=np.float64)
    ends = np.asarray(wall_end_xz_um, dtype=np.float64)
    if starts.ndim != 2 or starts.shape[1] != 2 or starts.shape != ends.shape:
        raise ValueError(
            "Continuous-wall start/end arrays must both have shape (N, 2)."
        )
    near = float(finite_element_distance_um)
    far = float(regular_grid_distance_um)
    if (
        not math.isfinite(near)
        or not math.isfinite(far)
        or near < 0.0
        or far <= near
    ):
        raise ValueError(
            "Hybrid distances must satisfy 0 <= finite_element_distance "
            "< regular_grid_distance."
        )

    triangles = np.ascontiguousarray(triangles_xz_um, dtype=np.float64)
    cartesian_fields = _cartesian_mesh_fields(
        lumen_mask,
        lumen_fraction,
        hybrid_region,
        finite_element_weight,
        distance_to_wall_um,
    )
    finite_element_fields = _finite_element_mesh_fields(
        triangles,
        starts,
        ends,
        near,
        far,
        spacing,
    )
    return HybridVelocityQualityData(
        source_path=source_path,
        x_coordinates_um=np.ascontiguousarray(x),
        z_coordinates_um=np.ascontiguousarray(z),
        spacing_um=spacing,
        cartesian_fields=cartesian_fields,
        finite_element_fields=finite_element_fields,
        continuous_wall_start_xz_um=np.ascontiguousarray(starts),
        continuous_wall_end_xz_um=np.ascontiguousarray(ends),
        finite_element_cell_vertices_xz_um=triangles,
        finite_element_distance_um=near,
        regular_grid_distance_um=far,
    )


def load_hybrid_velocity_quality_data(
    path: str | Path | None = None,
) -> HybridVelocityQualityData:
    """Load mesh-partition data from one saved hybrid result."""

    field_path = _resolve_field_archive(path)
    with np.load(field_path, allow_pickle=False) as archive:
        x = np.asarray(_required_array(archive, "x_coordinates_um", ndim=1))
        z = np.asarray(_required_array(archive, "z_coordinates_um", ndim=1))
        shape = (int(x.size), int(z.size))
        spacing_values = np.asarray(_required_array(archive, "spacing_um")).reshape(-1)
        if spacing_values.size != 1:
            raise ValueError("spacing_um must contain exactly one value.")

        lumen = _grid_array(archive, "lumen_mask", shape)
        fraction = _grid_array(archive, "lumen_fraction", shape)
        region = _grid_array(archive, "hybrid_velocity_region", shape)
        weight = _grid_array(archive, "hybrid_finite_element_weight", shape)
        wall_distance = _grid_array(archive, "distance_to_wall_um", shape)
        triangles = _required_array(archive, "fem_cell_vertices_xz_um")
        wall_start = _required_array(
            archive, "continuous_wall_start_xz_um", ndim=2
        )
        wall_end = _required_array(archive, "continuous_wall_end_xz_um", ndim=2)
        near = np.asarray(
            _required_array(archive, "hybrid_finite_element_distance_um"),
            dtype=np.float64,
        ).reshape(-1)
        far = np.asarray(
            _required_array(archive, "hybrid_regular_grid_distance_um"),
            dtype=np.float64,
        ).reshape(-1)
        if near.size != 1 or far.size != 1:
            raise ValueError("Each hybrid distance must contain exactly one value.")

    return _assemble_data(
        source_path=field_path,
        x_coordinates_um=x,
        z_coordinates_um=z,
        spacing_um=float(spacing_values[0]),
        lumen_mask=lumen,
        lumen_fraction=fraction,
        hybrid_region=region,
        finite_element_weight=weight,
        distance_to_wall_um=wall_distance,
        wall_start_xz_um=wall_start,
        wall_end_xz_um=wall_end,
        triangles_xz_um=triangles,
        finite_element_distance_um=float(near[0]),
        regular_grid_distance_um=float(far[0]),
    )


def hybrid_velocity_quality_data_from_objects(
    domain: object,
    raster: object,
    flow: object,
    *,
    continuous_geometry: object,
    source_path: str | Path = "in_memory_hybrid_velocity",
) -> HybridVelocityQualityData:
    """Build mesh-partition data from one in-memory hybrid solution."""

    from ulm_microbubble_traj_gen_2D.utils.flow.hybrid_velocity import (
        hybrid_region_map,
    )

    x = np.asarray(getattr(domain, "x_coordinates_um"), dtype=np.float64)
    z = np.asarray(getattr(domain, "z_coordinates_um"), dtype=np.float64)
    shape = (int(x.size), int(z.size))
    if shape != tuple(int(value) for value in getattr(domain, "shape")):
        raise ValueError("Domain coordinates disagree with domain.shape.")

    lumen = np.asarray(getattr(raster, "lumen_mask"))
    fraction = np.asarray(getattr(raster, "lumen_fraction"))
    wall_distance = np.asarray(getattr(raster, "distance_to_wall_um"))
    hybrid = getattr(flow, "hybrid_velocity")
    region, weight = hybrid_region_map(wall_distance, lumen, hybrid)
    finite_element = hybrid.finite_element

    return _assemble_data(
        source_path=Path(source_path).expanduser().resolve(),
        x_coordinates_um=x,
        z_coordinates_um=z,
        spacing_um=float(getattr(domain, "spacing_um")),
        lumen_mask=lumen,
        lumen_fraction=fraction,
        hybrid_region=region,
        finite_element_weight=weight,
        distance_to_wall_um=wall_distance,
        wall_start_xz_um=np.asarray(
            getattr(continuous_geometry, "solid_face_start_xz_um")
        ),
        wall_end_xz_um=np.asarray(
            getattr(continuous_geometry, "solid_face_end_xz_um")
        ),
        triangles_xz_um=np.asarray(finite_element.cell_vertices_xz_um),
        finite_element_distance_um=float(hybrid.finite_element_distance_um),
        regular_grid_distance_um=float(hybrid.regular_grid_distance_um),
    )


def build_hybrid_velocity_quality_data_from_config(
    config_path: str | Path,
    *,
    quick_test: bool = False,
) -> HybridVelocityQualityData:
    """Run the configured solve and build its mesh-partition view."""

    from ulm_microbubble_traj_gen_2D.utils.core.config import load_config
    from ulm_microbubble_traj_gen_2D.utils.flow.dolfinx_gmsh_solver import (
        solve_dolfinx_stokes_gmsh_2d,
    )
    from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
        build_continuous_vessel_geometry,
    )
    from ulm_microbubble_traj_gen_2D.utils.geometry.grid_domain import (
        build_domain_from_vessels,
    )
    from ulm_microbubble_traj_gen_2D.utils.geometry.vessel_rasterizer import (
        rasterize_vessels,
    )
    from ulm_microbubble_traj_gen_2D.utils.io.vascular_io import load_physics_input

    cfg = load_config(config_path, quick_test=quick_test)
    physics_input = load_physics_input(
        cfg.model_dir,
        planar_extrusion_depth_um=cfg.field.effective_thickness_um,
    )
    domain = build_domain_from_vessels(physics_input.vessels, cfg.domain)
    continuous_geometry = build_continuous_vessel_geometry(
        physics_input.vessels,
        domain,
        maximum_boundary_element_length_um=(
            cfg.domain.continuous_boundary_maximum_element_length_um
        ),
    )
    raster = rasterize_vessels(
        physics_input.vessels,
        domain,
        cfg.domain,
        effective_thickness_um=cfg.field.effective_thickness_um,
        continuous_geometry=continuous_geometry,
        dynamic_viscosity_mpas=(
            cfg.field.kinematic_viscosity_um2_s
            * cfg.field.blood_density_kg_m3
            * 1.0e-9
        ),
    )
    flow = solve_dolfinx_stokes_gmsh_2d(
        domain,
        raster,
        cfg.field,
        physics_input.vessels,
        continuous_geometry,
        vessel_metadata=physics_input.vessel_metadata,
    )
    return hybrid_velocity_quality_data_from_objects(
        domain,
        raster,
        flow,
        continuous_geometry=continuous_geometry,
        source_path=cfg.source_path,
    )


def _build_layers(data: HybridVelocityQualityData) -> list[_Layer]:
    cartesian = data.cartesian_fields
    finite_element = data.finite_element_fields
    return [
        _Layer(
            "cartesian_cell_type",
            "[Cartesian] Cell type",
            "",
            cartesian["cartesian_cell_type"],
            "tab10",
            "cartesian",
            (-0.5, 3.5),
            _CARTESIAN_CELL_NAMES,
        ),
        _Layer(
            "cartesian_hybrid_region",
            "[Cartesian] Velocity region",
            "",
            cartesian["cartesian_hybrid_region"],
            "Set1",
            "cartesian",
            (-0.5, 3.5),
            _VELOCITY_REGION_NAMES,
        ),
        _Layer(
            "cartesian_finite_element_weight",
            "[Cartesian] FEM blend weight",
            "fraction",
            cartesian["cartesian_finite_element_weight"],
            "magma",
            "cartesian",
            (0.0, 1.0),
        ),
        _Layer(
            "cartesian_lumen_fraction",
            "[Cartesian] Lumen coverage",
            "fraction",
            cartesian["cartesian_lumen_fraction"],
            "viridis",
            "cartesian",
            (0.0, 1.0),
        ),
        _Layer(
            "cartesian_wall_distance_um",
            "[Cartesian] Wall distance",
            "µm",
            cartesian["cartesian_wall_distance_um"],
            "magma",
            "cartesian",
        ),
        _Layer(
            "fem_triangle_region",
            "[DOLFINx] Triangle region",
            "",
            finite_element["fem_triangle_region"],
            "Set1",
            "finite_element",
            (0.5, 3.5),
            {key: value for key, value in _VELOCITY_REGION_NAMES.items() if key},
        ),
        _Layer(
            "fem_triangle_size_class",
            "[DOLFINx] Relative size",
            "",
            finite_element["fem_triangle_size_class"],
            "Accent",
            "finite_element",
            (0.5, 3.5),
            _TRIANGLE_SIZE_NAMES,
        ),
        _Layer(
            "fem_triangle_area_um2",
            "[DOLFINx] Triangle area",
            "µm²",
            finite_element["fem_triangle_area_um2"],
            "viridis",
            "finite_element",
        ),
        _Layer(
            "fem_triangle_minimum_edge_um",
            "[DOLFINx] Minimum edge",
            "µm",
            finite_element["fem_triangle_minimum_edge_um"],
            "viridis",
            "finite_element",
        ),
        _Layer(
            "fem_triangle_maximum_edge_um",
            "[DOLFINx] Maximum edge",
            "µm",
            finite_element["fem_triangle_maximum_edge_um"],
            "viridis",
            "finite_element",
        ),
        _Layer(
            "fem_triangle_shape_quality",
            "[DOLFINx] Shape quality",
            "0 to 1",
            finite_element["fem_triangle_shape_quality"],
            "RdYlGn",
            "finite_element",
            (0.0, 1.0),
        ),
        _Layer(
            "fem_triangle_wall_distance_um",
            "[DOLFINx] Centroid-wall distance",
            "µm",
            finite_element["fem_triangle_wall_distance_um"],
            "magma",
            "finite_element",
        ),
    ]


class HybridVelocityQualityViewer:
    """Matplotlib viewer for Cartesian and finite-element mesh partitions."""

    def __init__(self, data: HybridVelocityQualityData):
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection, PolyCollection
        from matplotlib.widgets import CheckButtons, RadioButtons
        from scipy.spatial import cKDTree

        self.data = data
        self.layers = _build_layers(data)
        self.layer_by_label = {layer.label: layer for layer in self.layers}
        self.current_layer = self.layers[0]
        self._show_cell_grid = True
        self._show_continuous_wall = True
        self._show_finite_element_mesh = True
        self._updating_limits = False
        centroids = data.finite_element_fields["fem_triangle_centroid_xz_um"]
        self._triangle_centroid_tree = cKDTree(centroids)

        self.figure = plt.figure(figsize=(16.2, 9.2))
        self.axes = self.figure.add_axes((0.055, 0.09, 0.59, 0.84))
        self.radio_axes = self.figure.add_axes((0.79, 0.39, 0.205, 0.54))
        self.check_axes = self.figure.add_axes((0.79, 0.29, 0.205, 0.075))
        self.colorbar_axes = self.figure.add_axes((0.665, 0.09, 0.012, 0.84))

        self.image = self.axes.imshow(
            self._display_values(self.current_layer).T,
            origin="lower",
            extent=data.extent,
            interpolation="nearest",
            aspect="equal",
            cmap=self._layer_cmap(self.current_layer),
            zorder=1,
        )
        self.finite_element_surface = PolyCollection(
            data.finite_element_cell_vertices_xz_um,
            closed=True,
            edgecolors="none",
            linewidths=0.0,
            zorder=2,
        )
        self.finite_element_surface.set_visible(False)
        self.axes.add_collection(self.finite_element_surface)
        self._active_mappable = self.image
        self._apply_layer_limits(self.image, self.current_layer)
        self.colorbar = self.figure.colorbar(self.image, cax=self.colorbar_axes)
        self._configure_colorbar()

        self.grid_collection = LineCollection(
            [], colors="#555555", linewidths=0.35, alpha=0.55, zorder=5
        )
        self.axes.add_collection(self.grid_collection)
        self.wall_collection = LineCollection(
            self._wall_segments(),
            colors="#080808",
            linewidths=0.75,
            alpha=0.95,
            zorder=7,
        )
        self.axes.add_collection(self.wall_collection)
        self.finite_element_collection = LineCollection(
            [],
            colors="#00c8d2",
            linewidths=0.38,
            alpha=0.72,
            zorder=6,
        )
        self.axes.add_collection(self.finite_element_collection)

        self.radio = RadioButtons(
            self.radio_axes,
            [layer.label for layer in self.layers],
            active=0,
            activecolor="#1f77b4",
        )
        for label in self.radio.labels:
            label.set_fontsize(8.2)
        self.radio.on_clicked(self._on_layer_selected)
        self.radio_axes.set_title("Mesh partition view", fontsize=10, loc="left")

        self.checks = CheckButtons(
            self.check_axes,
            (
                "Continuous wall",
                "DOLFINx triangle edges",
                "Cartesian cell edges",
            ),
            (True, True, True),
        )
        for label in self.checks.labels:
            label.set_fontsize(8.5)
        self.checks.on_clicked(self._on_overlay_toggled)

        self.summary_text = self.figure.text(
            0.79,
            0.265,
            "",
            ha="left",
            va="top",
            fontsize=8.1,
            family="monospace",
        )
        self.hover_text = self.figure.text(
            0.06,
            0.025,
            "Move over a cell or triangle for exact mesh values.",
            ha="left",
            va="center",
            fontsize=9,
        )

        self.axes.set_xlabel("X (µm)")
        self.axes.set_ylabel("Z (µm)")
        self.axes.set_xlim(data.extent[:2])
        self.axes.set_ylim(data.extent[2:])
        self.axes.callbacks.connect("xlim_changed", self._on_limits_changed)
        self.axes.callbacks.connect("ylim_changed", self._on_limits_changed)
        self.figure.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.figure.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.figure.canvas.mpl_connect("key_press_event", self._on_key)
        self._update_title_and_summary()
        self._refresh_dynamic_overlays()

    @staticmethod
    def _display_values(layer: _Layer) -> np.ma.MaskedArray:
        values = np.asarray(layer.values, dtype=np.float64)
        return np.ma.array(values, mask=~np.isfinite(values))

    @staticmethod
    def _layer_cmap(layer: _Layer) -> object:
        from matplotlib.colors import ListedColormap

        if layer.key == "cartesian_cell_type":
            return ListedColormap(("#eeeeee", "#d62728", "#d66fba", "#b8bd22"))
        if layer.key in {"cartesian_hybrid_region", "fem_triangle_region"}:
            colors = ("#eeeeee", "#377eb8", "#ff9f1c", "#2ca02c")
            return ListedColormap(colors if layer.mesh == "cartesian" else colors[1:])
        if layer.key == "fem_triangle_size_class":
            return ListedColormap(("#377eb8", "#f2c14e", "#d62728"))
        return layer.cmap

    def _apply_layer_limits(self, mappable: object, layer: _Layer) -> None:
        if layer.fixed_limits is not None:
            mappable.set_clim(*layer.fixed_limits)
            return
        values = self._display_values(layer).compressed()
        if values.size == 0:
            mappable.set_clim(0.0, 1.0)
            return
        minimum = float(np.min(values))
        maximum = float(np.max(values))
        if minimum == maximum:
            padding = max(abs(minimum) * 0.05, 1.0e-12)
            minimum -= padding
            maximum += padding
        mappable.set_clim(minimum, maximum)

    def _configure_colorbar(self) -> None:
        layer = self.current_layer
        self.colorbar.ax.yaxis.set_label_position("left")
        self.colorbar.set_label(
            layer.label if not layer.unit else f"{layer.label} ({layer.unit})"
        )
        if layer.category_names:
            ticks = sorted(layer.category_names)
            self.colorbar.set_ticks(ticks)
            self.colorbar.set_ticklabels(
                [layer.category_names[value] for value in ticks]
            )
        else:
            from matplotlib.ticker import AutoLocator

            self.colorbar.locator = AutoLocator()
            self.colorbar.update_ticks()

    def _wall_segments(self) -> np.ndarray:
        return np.stack(
            (
                self.data.continuous_wall_start_xz_um,
                self.data.continuous_wall_end_xz_um,
            ),
            axis=1,
        )

    def _on_layer_selected(self, label: str) -> None:
        self.current_layer = self.layer_by_label[str(label)]
        if self.current_layer.mesh == "cartesian":
            self.image.set_data(self._display_values(self.current_layer).T)
            self.image.set_cmap(self._layer_cmap(self.current_layer))
            self.image.set_visible(True)
            self.finite_element_surface.set_visible(False)
            self._active_mappable = self.image
        else:
            self.finite_element_surface.set_array(
                self._display_values(self.current_layer)
            )
            self.finite_element_surface.set_cmap(
                self._layer_cmap(self.current_layer)
            )
            self.finite_element_surface.set_visible(True)
            self.image.set_visible(False)
            self._active_mappable = self.finite_element_surface
        self._apply_layer_limits(self._active_mappable, self.current_layer)
        self.colorbar.update_normal(self._active_mappable)
        self._configure_colorbar()
        self._update_title_and_summary()
        self._refresh_dynamic_overlays()
        self.figure.canvas.draw_idle()

    def select_layer(self, key: str) -> None:
        """Select a mesh layer by its stable key."""

        for index, layer in enumerate(self.layers):
            if layer.key == key:
                self.radio.set_active(index)
                return
        available = ", ".join(layer.key for layer in self.layers)
        raise KeyError(f"Unknown mesh layer {key!r}; available: {available}")

    def _on_overlay_toggled(self, label: str) -> None:
        if label == "Continuous wall":
            self._show_continuous_wall = not self._show_continuous_wall
        elif label == "DOLFINx triangle edges":
            self._show_finite_element_mesh = not self._show_finite_element_mesh
        elif label == "Cartesian cell edges":
            self._show_cell_grid = not self._show_cell_grid
        self._refresh_dynamic_overlays()
        self.figure.canvas.draw_idle()

    def _on_limits_changed(self, _axes: object) -> None:
        if not self._updating_limits:
            self._refresh_dynamic_overlays()

    def _visible_index_bounds(self) -> tuple[int, int, int, int]:
        xmin, xmax = sorted(self.axes.get_xlim())
        zmin, zmax = sorted(self.axes.get_ylim())
        x = self.data.x_coordinates_um
        z = self.data.z_coordinates_um
        half = 0.5 * self.data.spacing_um
        ix0 = max(0, int(np.searchsorted(x + half, xmin, side="left")))
        ix1 = min(x.size, int(np.searchsorted(x - half, xmax, side="right")))
        iz0 = max(0, int(np.searchsorted(z + half, zmin, side="left")))
        iz1 = min(z.size, int(np.searchsorted(z - half, zmax, side="right")))
        return ix0, max(ix0, ix1), iz0, max(iz0, iz1)

    def _refresh_dynamic_overlays(self) -> None:
        ix0, ix1, iz0, iz1 = self._visible_index_bounds()
        nx = ix1 - ix0
        nz = iz1 - iz0
        h = self.data.spacing_um

        segments = []
        if self._show_cell_grid and nx > 0 and nz > 0 and nx * nz <= 12_000:
            x_edges = self.data.x_coordinates_um[ix0:ix1] - 0.5 * h
            z_edges = self.data.z_coordinates_um[iz0:iz1] - 0.5 * h
            x_edges = np.append(
                x_edges, self.data.x_coordinates_um[ix1 - 1] + 0.5 * h
            )
            z_edges = np.append(
                z_edges, self.data.z_coordinates_um[iz1 - 1] + 0.5 * h
            )
            segments.extend(
                ((float(value), float(z_edges[0])), (float(value), float(z_edges[-1])))
                for value in x_edges
            )
            segments.extend(
                ((float(x_edges[0]), float(value)), (float(x_edges[-1]), float(value)))
                for value in z_edges
            )
        self.grid_collection.set_segments(segments)
        self.wall_collection.set_visible(self._show_continuous_wall)
        self.finite_element_collection.set_segments(
            self._visible_finite_element_edges()
        )
        self.finite_element_collection.set_visible(
            self._show_finite_element_mesh
        )

    def _visible_finite_element_edges(self) -> np.ndarray:
        if not self._show_finite_element_mesh:
            return np.empty((0, 2, 2), dtype=np.float64)
        triangles = self.data.finite_element_cell_vertices_xz_um
        xmin, xmax = sorted(self.axes.get_xlim())
        zmin, zmax = sorted(self.axes.get_ylim())
        minimum = np.min(triangles, axis=1)
        maximum = np.max(triangles, axis=1)
        visible = (
            (maximum[:, 0] >= xmin)
            & (minimum[:, 0] <= xmax)
            & (maximum[:, 1] >= zmin)
            & (minimum[:, 1] <= zmax)
        )
        selected = triangles[visible]
        if selected.shape[0] > 20_000:
            return np.empty((0, 2, 2), dtype=np.float64)
        return np.concatenate(
            (
                selected[:, (0, 1)],
                selected[:, (1, 2)],
                selected[:, (2, 0)],
            ),
            axis=0,
        )

    def _on_scroll(self, event: object) -> None:
        if getattr(event, "inaxes", None) is not self.axes:
            return
        xdata = getattr(event, "xdata", None)
        zdata = getattr(event, "ydata", None)
        if xdata is None or zdata is None:
            return
        factor = 0.80 if getattr(event, "button", None) == "up" else 1.25
        xmin, xmax = self.axes.get_xlim()
        zmin, zmax = self.axes.get_ylim()
        width = (xmax - xmin) * factor
        height = (zmax - zmin) * factor
        x_fraction = (float(xdata) - xmin) / max(xmax - xmin, np.finfo(float).eps)
        z_fraction = (float(zdata) - zmin) / max(zmax - zmin, np.finfo(float).eps)
        self._updating_limits = True
        self.axes.set_xlim(
            float(xdata) - x_fraction * width,
            float(xdata) + (1.0 - x_fraction) * width,
        )
        self.axes.set_ylim(
            float(zdata) - z_fraction * height,
            float(zdata) + (1.0 - z_fraction) * height,
        )
        self._updating_limits = False
        self._refresh_dynamic_overlays()
        self.figure.canvas.draw_idle()

    def _on_motion(self, event: object) -> None:
        if getattr(event, "inaxes", None) is not self.axes:
            return
        x_value = getattr(event, "xdata", None)
        z_value = getattr(event, "ydata", None)
        if x_value is None or z_value is None:
            return
        if self.current_layer.mesh == "finite_element":
            self._show_triangle_hover(float(x_value), float(z_value))
        else:
            self._show_cartesian_hover(float(x_value), float(z_value))
        self.figure.canvas.draw_idle()

    def _show_cartesian_hover(self, x_value: float, z_value: float) -> None:
        i = int(np.floor((x_value - self.data.extent[0]) / self.data.spacing_um))
        j = int(np.floor((z_value - self.data.extent[2]) / self.data.spacing_um))
        if not (0 <= i < self.data.shape[0] and 0 <= j < self.data.shape[1]):
            return
        fields = self.data.cartesian_fields
        cell_type = int(fields["cartesian_cell_type"][i, j])
        region = int(fields["cartesian_hybrid_region"][i, j])
        details = (
            f"cell=({i}, {j})  |  "
            f"centre=({self.data.x_coordinates_um[i]:.6g}, "
            f"{self.data.z_coordinates_um[j]:.6g}) µm  |  "
            f"type={_CARTESIAN_CELL_NAMES[cell_type]}  |  "
            f"region={_VELOCITY_REGION_NAMES[region]}  |  "
            f"coverage={float(fields['cartesian_lumen_fraction'][i, j]):.6g}  |  "
            f"FEM weight={float(fields['cartesian_finite_element_weight'][i, j]):.6g}"
        )
        self.hover_text.set_text(details)

    def _show_triangle_hover(self, x_value: float, z_value: float) -> None:
        _, index = self._triangle_centroid_tree.query((x_value, z_value))
        i = int(index)
        fields = self.data.finite_element_fields
        centroid = fields["fem_triangle_centroid_xz_um"][i]
        region = int(fields["fem_triangle_region"][i])
        size_class = int(fields["fem_triangle_size_class"][i])
        self.hover_text.set_text(
            f"nearest triangle={i}  |  "
            f"centroid=({centroid[0]:.6g}, {centroid[1]:.6g}) µm  |  "
            f"region={_VELOCITY_REGION_NAMES[region]}  |  "
            f"size={_TRIANGLE_SIZE_NAMES[size_class]}  |  "
            f"area={fields['fem_triangle_area_um2'][i]:.6g} µm²  |  "
            f"edges={fields['fem_triangle_minimum_edge_um'][i]:.6g}.."
            f"{fields['fem_triangle_maximum_edge_um'][i]:.6g} µm  |  "
            f"quality={fields['fem_triangle_shape_quality'][i]:.6g}"
        )

    def _on_key(self, event: object) -> None:
        key = str(getattr(event, "key", "")).lower()
        if key in {"r", "home"}:
            self.axes.set_xlim(self.data.extent[:2])
            self.axes.set_ylim(self.data.extent[2:])
            self.figure.canvas.draw_idle()
        elif key == "g":
            self.checks.set_active(2)
        elif key == "m":
            self.checks.set_active(1)
        elif key == "w":
            self.checks.set_active(0)

    def _update_title_and_summary(self) -> None:
        self.axes.set_title(
            f"{self.current_layer.label} — {self._source_label()}",
            loc="left",
        )
        cartesian = self.data.cartesian_fields
        finite_element = self.data.finite_element_fields
        cell_types = cartesian["cartesian_cell_type"]
        regions = cartesian["cartesian_hybrid_region"]
        triangle_regions = finite_element["fem_triangle_region"]
        sizes = finite_element["fem_triangle_size_class"]
        quality = finite_element["fem_triangle_shape_quality"]
        area = finite_element["fem_triangle_area_um2"]
        values = self._display_values(self.current_layer).compressed()
        value_summary = "no finite values"
        if values.size:
            value_summary = (
                f"layer min/med/max:\n"
                f"{np.min(values):.5g} / {np.median(values):.5g} / "
                f"{np.max(values):.5g}"
            )
        self.summary_text.set_text(
            "CARTESIAN MESH\n"
            f"shape: {self.data.shape[0]} x {self.data.shape[1]}\n"
            f"spacing: {self.data.spacing_um:.6g} µm\n"
            f"outside/cut-out/cut-in/full:\n"
            f"{self._counts(cell_types, (0, 1, 2, 3))}\n"
            f"FEM/blend/grid cells:\n"
            f"{self._counts(regions, (1, 2, 3))}\n\n"
            "DOLFINX TRIANGLES\n"
            f"total: {area.size}\n"
            f"FEM/blend/grid triangles:\n"
            f"{self._counts(triangle_regions, (1, 2, 3))}\n"
            f"fine/similar/coarse:\n"
            f"{self._counts(sizes, (1, 2, 3))}\n"
            f"area min/med/max (µm²):\n"
            f"{np.min(area):.5g} / {np.median(area):.5g} / {np.max(area):.5g}\n"
            f"quality min/med: {np.min(quality):.5g} / {np.median(quality):.5g}\n"
            f"distance limits: {self.data.finite_element_distance_um:.5g}, "
            f"{self.data.regular_grid_distance_um:.5g} µm\n\n"
            f"{value_summary}"
        )

    @staticmethod
    def _counts(values: np.ndarray, keys: tuple[int, ...]) -> str:
        return " / ".join(str(int(np.count_nonzero(values == key))) for key in keys)

    def _source_label(self) -> str:
        if self.data.source_path.name == FIELD_FILENAME:
            return self.data.source_path.parent.name
        return self.data.source_path.stem

    def save_snapshot(self, path: str | Path, *, dpi: int = 180) -> Path:
        """Save the current viewport as PNG, PDF, or SVG."""

        output_path = Path(path).expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self.figure.savefig(output_path, dpi=int(dpi), bbox_inches="tight")
        return output_path

    def show(self) -> None:
        import matplotlib.pyplot as plt

        plt.show()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect Cartesian cell types and DOLFINx triangle partitions "
            "against the continuous vessel wall."
        )
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--field",
        type=Path,
        default=None,
        help=(
            "A velocity_and_wall_shear_field.npz or its result directory. "
            "The newest result is used when omitted."
        ),
    )
    source.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Run the configured DOLFINx solve and inspect its meshes.",
    )
    parser.add_argument(
        "--quick-test",
        action="store_true",
        help="Apply YAML quick-test overrides when --config is used.",
    )
    parser.add_argument(
        "--layer",
        default="cartesian_cell_type",
        help="Initial stable mesh-layer key.",
    )
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=None,
        help="Optionally save the initial/current view to PNG, PDF, or SVG.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Build the viewer without opening a GUI.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    data = (
        build_hybrid_velocity_quality_data_from_config(
            args.config, quick_test=bool(args.quick_test)
        )
        if args.config is not None
        else load_hybrid_velocity_quality_data(args.field)
    )
    viewer = HybridVelocityQualityViewer(data)
    viewer.select_layer(str(args.layer))
    if args.snapshot is not None:
        saved = viewer.save_snapshot(args.snapshot)
        print(f"Saved mesh-partition snapshot: {saved}")
    print(f"Hybrid field archive: {data.source_path}")
    print(
        "Controls: mouse wheel=zoom, toolbar=pan/home, "
        "M=DOLFINx edges, G=Cartesian edges, W=continuous wall."
    )
    if not args.no_show:
        viewer.show()


if __name__ == "__main__":
    main()
