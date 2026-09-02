"""Academic interactive desktop viewer for validated production CFD fields.

This module is deliberately a viewer, not a solver or scientific post-processor.
It validates and reads the accepted production VTU without modifying it on disk.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pyvista as pv
import yaml
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from scipy.spatial import cKDTree

from utils.cfd_flow.io import read_json, sha256_file, write_json
from utils.cfd_flow.visualization import _nearest_cell_streamlines


PROJECT_ROOT = Path(__file__).resolve().parent
RUN_PATTERN = "production_tau1_base_promotion_anchor003274_*"
EXPECTED_RUN_STATUS = "CFD_FLOW_PRODUCTION_TAU1_INTEGRATION_AND_VISUAL_REGRESSION_PASS"
REQUIRED_ARRAYS = (
    "velocity_phy",
    "velocity_magnitude_m_s",
    "velocity_magnitude_mm_s",
    "pressure_gauge_pa",
    "pressure_absolute_solver_pa",
    "rho_lattice",
)
FIELD_ORDER = ("velocity", "pressure", "rho")
PORT_ORDER = ("inlet", "outlet_01", "outlet_02", "outlet_03")
PORT_LABELS = {
    "inlet": "INLET",
    "outlet_01": "OUTLET 01",
    "outlet_02": "OUTLET 02",
    "outlet_03": "OUTLET 03",
}
PORT_COLORS = {
    "inlet": "#176B68",
    "outlet_01": "#4C78A8",
    "outlet_02": "#D18F2F",
    "outlet_03": "#7A5195",
}
BACKGROUND = "#FAF8F3"
CONTEXT_COLOR = "#AEB5BA"


class VisualizerInputError(RuntimeError):
    """Raised when production artifacts do not satisfy the viewing contract."""


@dataclass(frozen=True, slots=True)
class VisualConfig:
    """User-selected display options."""

    width: int = 1800
    height: int = 1100
    initial_scalar: str = "velocity"
    streamline_seeds: int = 24
    build_streamlines: bool = True
    show_ports: bool = True
    full_range: bool = False
    debug_cells: bool = False
    numerical_pressure_debug: bool = False


@dataclass(frozen=True, slots=True)
class FieldRange:
    """Raw and robust display range for one scalar."""

    raw_min: float
    raw_max: float
    percentile_min: float
    percentile_max: float

    def selected(self, full: bool) -> tuple[float, float]:
        return (
            (self.raw_min, self.raw_max)
            if full
            else (self.percentile_min, self.percentile_max)
        )


@dataclass(frozen=True, slots=True)
class FieldSpec:
    """Presentation metadata for a validated VTU cell array."""

    array: str
    title: str
    units: str
    cmap: Any


@dataclass(slots=True)
class VisualData:
    """Validated production evidence and cached display geometry."""

    run_dir: Path
    vtu_path: Path
    vtu_sha256: str
    manifest: dict[str, Any]
    metrics: dict[str, Any]
    steady_qc: dict[str, Any]
    run_summary: dict[str, Any]
    physical_flux: dict[str, Any]
    plane_contract: dict[str, Any]
    grid_um: pv.UnstructuredGrid
    surface_um: pv.PolyData
    centers_um: np.ndarray
    center_tree_um: cKDTree
    fields: dict[str, FieldSpec]
    ranges: dict[str, FieldRange]
    streamlines_um: pv.PolyData | None
    valid_streamline_count: int


def build_parser() -> argparse.ArgumentParser:
    """Build the small, desktop-oriented command-line interface."""

    parser = argparse.ArgumentParser(
        description="Academic interactive viewer for accepted production CFD fields."
    )
    parser.add_argument("--run-dir", type=Path, help="Explicit production run directory")
    parser.add_argument("--vtu", type=Path, help="Explicit production VTU (highest priority)")
    parser.add_argument(
        "--scalar",
        choices=("velocity", "pressure", "rho", "numerical-pressure"),
        default="velocity",
    )
    parser.add_argument("--window-width", type=int, default=1800)
    parser.add_argument("--window-height", type=int, default=1100)
    parser.add_argument("--streamline-seeds", type=int, default=24)
    parser.add_argument("--no-streamlines", action="store_true")
    parser.add_argument("--no-ports", action="store_true")
    parser.add_argument("--full-range", action="store_true")
    parser.add_argument("--publication-screenshot", type=Path)
    parser.add_argument("--debug-cells", action="store_true")
    parser.add_argument("--show-numerical-pressure-debug", action="store_true")
    parser.add_argument("--off-screen", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser


def _load_json_object(path: Path, role: str) -> dict[str, Any]:
    if not path.is_file():
        raise VisualizerInputError(f"{role} missing: {path}")
    value = read_json(path)
    if not isinstance(value, dict):
        raise VisualizerInputError(f"{role} is not a JSON object: {path}")
    return value


def _candidate_timestamp(run_dir: Path) -> float:
    summary_path = run_dir / "qc" / "run_summary.json"
    try:
        completed = str(read_json(summary_path).get("completed_at", ""))
        return datetime.fromisoformat(completed).timestamp()
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return run_dir.stat().st_mtime


def _candidate_is_complete(run_dir: Path) -> bool:
    """Cheap but cryptographically strict discovery check."""

    try:
        flow = run_dir / "flow"
        vtu = flow / "production_steady_flow_field.vtu"
        manifest = _load_json_object(
            flow / "production_steady_flow_field_manifest.json", "VTU manifest"
        )
        summary = _load_json_object(run_dir / "qc" / "run_summary.json", "run summary")
        steady = _load_json_object(
            run_dir / "qc" / "production_steady_qc.json", "steady QC"
        )
        metrics = run_dir / "qc" / "production_primary_metrics.json"
        flux = run_dir / "steady_replay" / "physical_port_flux.json"
        return bool(
            manifest.get("status") == "PASS"
            and summary.get("status") == EXPECTED_RUN_STATUS
            and steady.get("status") == "PASS"
            and metrics.is_file()
            and flux.is_file()
            and vtu.is_file()
            and int(manifest.get("cell_count", 0)) > 0
            and sha256_file(vtu) == manifest.get("sha256")
        )
    except (OSError, ValueError, TypeError, VisualizerInputError):
        return False


def locate_latest_run(project_root: Path = PROJECT_ROOT) -> Path:
    """Return the newest complete, PASS production promotion run."""

    output = Path(project_root).resolve() / "outputs" / "cfd_flow"
    candidates = sorted(
        (path for path in output.glob(RUN_PATTERN) if path.is_dir()),
        key=_candidate_timestamp,
        reverse=True,
    )
    for candidate in candidates:
        if _candidate_is_complete(candidate):
            return candidate.resolve()
    raise VisualizerInputError(f"No complete PASS production run found under {output}")


def resolve_run_and_vtu(
    *,
    project_root: Path = PROJECT_ROOT,
    explicit_run_dir: Path | None = None,
    explicit_vtu: Path | None = None,
) -> tuple[Path, Path]:
    """Resolve inputs with VTU > run directory > latest-run priority."""

    if explicit_vtu is not None:
        vtu = explicit_vtu.expanduser().resolve()
        return vtu.parent.parent, vtu
    if explicit_run_dir is not None:
        run_dir = explicit_run_dir.expanduser().resolve()
    else:
        run_dir = locate_latest_run(project_root)
    return run_dir, run_dir / "flow" / "production_steady_flow_field.vtu"


def calculate_field_range(values: np.ndarray) -> FieldRange:
    """Calculate honest raw and p1-p99 visualization limits."""

    array = np.asarray(values, dtype=np.float64).reshape(-1)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise VisualizerInputError("Scalar field is empty or contains non-finite values")
    raw_min, raw_max = (
        float(value) for value in (np.min(array), np.max(array))
    )
    low, high = (float(value) for value in np.percentile(array, (1.0, 99.0)))
    if not raw_min < raw_max:
        delta = max(abs(raw_min), 1.0) * 1.0e-12
        raw_min, raw_max = raw_min - delta, raw_max + delta
    if not low < high:
        low, high = raw_min, raw_max
    return FieldRange(raw_min, raw_max, low, high)


def _pressure_colormap(low: float, high: float) -> Any:
    """Place neutral grey at zero without forcing a symmetric range."""

    if not low < 0.0 < high:
        return "cividis"
    zero_fraction = float(np.clip(-low / (high - low), 0.05, 0.95))
    base = colormaps["coolwarm"]
    return LinearSegmentedColormap.from_list(
        "gauge_pressure_zero",
        [(0.0, base(0.02)), (zero_fraction, base(0.5)), (1.0, base(0.98))],
    )


def validate_required_fields(grid: pv.DataSet, manifest: dict[str, Any]) -> None:
    """Validate array names, shapes, finiteness, and redundant units."""

    if manifest.get("status") != "PASS":
        raise VisualizerInputError("VTU manifest status is not PASS")
    expected_cells = int(manifest.get("cell_count", -1))
    if expected_cells <= 0 or grid.n_cells != expected_cells:
        raise VisualizerInputError(
            f"VTU cell count mismatch: grid={grid.n_cells}, manifest={expected_cells}"
        )
    missing = [name for name in REQUIRED_ARRAYS if name not in grid.cell_data]
    if missing:
        raise VisualizerInputError(f"VTU missing required cell arrays: {missing}")
    for name in REQUIRED_ARRAYS:
        values = np.asarray(grid.cell_data[name])
        if values.shape[0] != grid.n_cells or not np.all(np.isfinite(values)):
            raise VisualizerInputError(f"Invalid or non-finite VTU array: {name}")
    velocity = np.asarray(grid.cell_data["velocity_phy"], dtype=np.float64)
    if velocity.ndim != 2 or velocity.shape[1] != 3:
        raise VisualizerInputError("velocity_phy must have exactly three components")
    speed = np.linalg.norm(velocity, axis=1)
    if not np.allclose(
        speed, grid.cell_data["velocity_magnitude_m_s"], rtol=2.0e-12, atol=1.0e-18
    ):
        raise VisualizerInputError("velocity magnitude is inconsistent with velocity_phy")
    if not np.allclose(
        speed * 1.0e3,
        grid.cell_data["velocity_magnitude_mm_s"],
        rtol=2.0e-12,
        atol=1.0e-15,
    ):
        raise VisualizerInputError("velocity mm/s field is inconsistent with velocity_phy")


def _plane_contract_candidates(run_dir: Path, project_root: Path) -> list[Path]:
    candidates: list[Path] = []
    snapshots = sorted((run_dir / "input").glob("*cfd_flow*.yaml"))
    snapshots += [project_root / "configs" / "cfd_flow_promotion_regression.yaml"]
    for snapshot in snapshots:
        try:
            value = yaml.safe_load(snapshot.read_text(encoding="utf-8"))
            configured = Path(value["paths"]["physical_plane_contract"])
            candidates.append(
                configured if configured.is_absolute() else project_root / configured
            )
        except (OSError, TypeError, KeyError, yaml.YAMLError):
            continue
    candidates.extend(
        (project_root / "outputs" / "cfd_flow").glob(
            "*/qc/physical_port_flux_plane_contract_v3.json"
        )
    )
    return list(dict.fromkeys(path.resolve() for path in candidates))


def _validate_plane_contract(
    contract: dict[str, Any], expected_contract_sha256: str
) -> None:
    if contract.get("status") != "PASS":
        raise VisualizerInputError("Physical plane contract status is not PASS")
    if contract.get("contract_sha256") != expected_contract_sha256:
        raise VisualizerInputError("Physical plane contract provenance mismatch")
    ports = contract.get("ports", {})
    if tuple(ports) != PORT_ORDER and set(ports) != set(PORT_ORDER):
        raise VisualizerInputError("Physical plane contract must define exactly four ports")
    for label in PORT_ORDER:
        try:
            plane = ports[label]["planes"]["central"]
            origin = np.asarray(plane["origin_m"], dtype=float)
            normal = np.asarray(plane["unit_normal"], dtype=float)
            basis_u = np.asarray(plane["basis_u"], dtype=float)
            basis_v = np.asarray(plane["basis_v"], dtype=float)
            contour = np.asarray(plane["physical_aperture_contour_uv_m"], dtype=float)
        except (KeyError, TypeError, ValueError) as error:
            raise VisualizerInputError(f"Malformed physical plane for {label}") from error
        if (
            origin.shape != (3,)
            or normal.shape != (3,)
            or basis_u.shape != (3,)
            or basis_v.shape != (3,)
            or contour.ndim != 2
            or contour.shape[1] != 2
            or len(contour) < 8
            or not np.all(np.isfinite(contour))
            or not math.isclose(float(np.linalg.norm(normal)), 1.0, rel_tol=0.0, abs_tol=1e-9)
        ):
            raise VisualizerInputError(f"Invalid physical plane geometry for {label}")


def load_plane_contract(
    run_dir: Path,
    physical_flux: dict[str, Any],
    *,
    project_root: Path = PROJECT_ROOT,
) -> dict[str, Any]:
    """Load the exact validated port-plane contract, never estimated coordinates."""

    expected = str(physical_flux.get("plane_contract_sha256", ""))
    for path in _plane_contract_candidates(run_dir, Path(project_root).resolve()):
        if not path.is_file():
            continue
        contract = _load_json_object(path, "physical plane contract")
        if contract.get("contract_sha256") == expected:
            _validate_plane_contract(contract, expected)
            contract["_loaded_from"] = str(path)
            return contract
    raise VisualizerInputError(
        f"Validated physical plane contract {expected or '<missing hash>'} not found"
    )


def _load_production_records(run_dir: Path) -> tuple[dict[str, Any], ...]:
    metrics = _load_json_object(
        run_dir / "qc" / "production_primary_metrics.json", "primary metrics"
    )
    steady = _load_json_object(
        run_dir / "qc" / "production_steady_qc.json", "steady QC"
    )
    summary = _load_json_object(run_dir / "qc" / "run_summary.json", "run summary")
    flux = _load_json_object(
        run_dir / "steady_replay" / "physical_port_flux.json", "physical port flux"
    )
    if steady.get("status") != "PASS" or flux.get("status") != "PASS":
        raise VisualizerInputError("Production steady or physical-flux QC is not PASS")
    if summary.get("status") != EXPECTED_RUN_STATUS:
        raise VisualizerInputError("Production run summary is not an accepted PASS run")
    if summary.get("steady_solution_source") != "VALIDATED_RESEARCH_BASE_ACCEPTED_RESTART":
        raise VisualizerInputError("Unexpected production steady source classification")
    required_metrics = (
        "iteration",
        "rho_mean",
        "Qin_m3_s",
        "Qout_m3_s",
        "physical_volume_closure",
        "inlet_gauge_pressure_pa",
    )
    if any(name not in metrics for name in required_metrics):
        raise VisualizerInputError("Production primary metrics are incomplete")
    if not np.all(np.isfinite([float(metrics[name]) for name in required_metrics])):
        raise VisualizerInputError("Production primary metrics contain non-finite values")
    if not math.isclose(
        float(metrics["Qin_m3_s"]),
        float(flux["Qin_m3_s"]),
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        raise VisualizerInputError("Primary metrics and physical-flux Qin disagree")
    return metrics, steady, summary, flux


def _inlet_seeds(plane: dict[str, Any], count: int) -> np.ndarray:
    """Create deterministic interior seeds on the validated inlet aperture."""

    origin = np.asarray(plane["origin_m"], dtype=float)
    basis_u = np.asarray(plane["basis_u"], dtype=float)
    basis_v = np.asarray(plane["basis_v"], dtype=float)
    radius = math.sqrt(float(plane["aperture_physical_area_m2"]) / math.pi)
    indices = np.arange(count, dtype=float)
    radial = 0.68 * np.sqrt((indices + 0.5) / count)
    angles = indices * math.pi * (3.0 - math.sqrt(5.0))
    return (
        origin
        + (radius * radial * np.cos(angles))[:, None] * basis_u
        + (radius * radial * np.sin(angles))[:, None] * basis_v
    )


def build_streamlines(
    centers_m: np.ndarray,
    velocity_m_s: np.ndarray,
    inlet_plane: dict[str, Any],
    *,
    dx_m: float,
    seed_count: int,
) -> tuple[pv.PolyData, int]:
    """Render the existing nearest-cell/wall-stop integration as VTK polylines."""

    seeds = _inlet_seeds(inlet_plane, seed_count)
    lines = _nearest_cell_streamlines(
        centers_m, velocity_m_s, seeds, dx_m=dx_m
    )
    if not lines:
        return pv.PolyData(), 0
    nearest = cKDTree(centers_m)
    points: list[np.ndarray] = []
    connectivity: list[int] = []
    speed_values: list[np.ndarray] = []
    offset = 0
    speed_mm_s = np.linalg.norm(velocity_m_s, axis=1) * 1.0e3
    for line in lines:
        count = len(line)
        points.append(line * 1.0e6)
        connectivity.extend((count, *range(offset, offset + count)))
        _, indices = nearest.query(line, k=1)
        speed_values.append(speed_mm_s[np.asarray(indices, dtype=int)])
        offset += count
    poly = pv.PolyData(np.vstack(points))
    poly.lines = np.asarray(connectivity, dtype=np.int64)
    poly.point_data["velocity_magnitude_mm_s"] = np.concatenate(speed_values)
    return poly, len(lines)


def load_and_validate_data(
    run_dir: Path,
    vtu_path: Path,
    config: VisualConfig,
    *,
    project_root: Path = PROJECT_ROOT,
) -> VisualData:
    """Read, validate, and cache all geometry required by the viewer."""

    run_dir = run_dir.resolve()
    vtu_path = vtu_path.resolve()
    manifest = _load_json_object(
        vtu_path.parent / "production_steady_flow_field_manifest.json", "VTU manifest"
    )
    if not vtu_path.is_file():
        raise VisualizerInputError(f"VTU missing: {vtu_path}")
    actual_sha = sha256_file(vtu_path)
    if manifest.get("sha256") and actual_sha != manifest["sha256"]:
        raise VisualizerInputError(
            f"VTU SHA256 mismatch: expected {manifest['sha256']}, got {actual_sha}"
        )
    metrics, steady, summary, flux = _load_production_records(run_dir)
    plane_contract = load_plane_contract(
        run_dir, flux, project_root=Path(project_root).resolve()
    )
    grid = pv.read(vtu_path)
    if not isinstance(grid, pv.UnstructuredGrid):
        grid = grid.cast_to_unstructured_grid()
    validate_required_fields(grid, manifest)
    if not math.isclose(
        float(np.mean(grid.cell_data["rho_lattice"])),
        float(metrics["rho_mean"]),
        rel_tol=2.0e-13,
        abs_tol=0.0,
    ):
        raise VisualizerInputError("VTU rho mean is inconsistent with production QC")
    pressure_offset = np.asarray(grid.cell_data["pressure_absolute_solver_pa"]) - np.asarray(
        grid.cell_data["pressure_gauge_pa"]
    )
    expected_offset = float(summary["numerical_contract"]["pressure_reference_pa"])
    if not np.allclose(pressure_offset, expected_offset, rtol=2.0e-13, atol=1.0e-6):
        raise VisualizerInputError("Gauge and numerical pressure arrays are inconsistent")

    centers_m = np.asarray(grid.cell_centers().points, dtype=np.float64)
    velocity = np.asarray(grid.cell_data["velocity_phy"], dtype=np.float64)
    dx_m = float(summary["numerical_contract"]["dx_m"])
    streamlines: pv.PolyData | None = None
    valid_streamlines = 0
    if config.build_streamlines:
        streamlines, valid_streamlines = build_streamlines(
            centers_m,
            velocity,
            plane_contract["ports"]["inlet"]["planes"]["central"],
            dx_m=dx_m,
            seed_count=config.streamline_seeds,
        )

    # Centralized display-only coordinate conversion; the source VTU is never saved.
    grid.points = np.asarray(grid.points, dtype=np.float64) * 1.0e6
    centers_um = centers_m * 1.0e6
    surface = grid.extract_surface(algorithm="dataset_surface")
    pressure_range = calculate_field_range(grid.cell_data["pressure_gauge_pa"])
    pressure_clim = pressure_range.selected(config.full_range)
    fields: dict[str, FieldSpec] = {
        "velocity": FieldSpec(
            "velocity_magnitude_mm_s", "Velocity magnitude", "mm/s", "viridis"
        ),
        "pressure": FieldSpec(
            "pressure_gauge_pa",
            "Gauge pressure",
            "Pa",
            _pressure_colormap(*pressure_clim),
        ),
        "rho": FieldSpec("rho_lattice", "Lattice density", "dimensionless", "cividis"),
    }
    if config.numerical_pressure_debug:
        fields["numerical-pressure"] = FieldSpec(
            "pressure_absolute_solver_pa",
            "Numerical pressure (OFFSET INCLUDED — NOT PHYSIOLOGICAL)",
            "Pa",
            "cividis",
        )
    ranges = {
        key: calculate_field_range(grid.cell_data[spec.array])
        for key, spec in fields.items()
    }
    return VisualData(
        run_dir=run_dir,
        vtu_path=vtu_path,
        vtu_sha256=actual_sha,
        manifest=manifest,
        metrics=metrics,
        steady_qc=steady,
        run_summary=summary,
        physical_flux=flux,
        plane_contract=plane_contract,
        grid_um=grid,
        surface_um=surface,
        centers_um=centers_um,
        center_tree_um=cKDTree(centers_um),
        fields=fields,
        ranges=ranges,
        streamlines_um=streamlines,
        valid_streamline_count=valid_streamlines,
    )


def _camera_record(plotter: pv.Plotter) -> dict[str, list[float]]:
    position, focal_point, view_up = plotter.camera_position
    return {
        "position_um": [float(value) for value in position],
        "focal_point_um": [float(value) for value in focal_point],
        "view_up": [float(value) for value in view_up],
    }


def _image_qc(path: Path) -> dict[str, Any]:
    with Image.open(path) as image:
        pixels = np.asarray(image.convert("RGB"), dtype=np.float64)
        width, height = image.size
    standard_deviation = float(np.std(pixels))
    checks = {
        "minimum_width": width >= 1600,
        "minimum_height": height >= 1000,
        "not_all_white": bool(np.any(pixels < 250.0)),
        "not_all_black": bool(np.any(pixels > 5.0)),
        "nonblank": standard_deviation > 1.0,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "path": str(path),
        "sha256": sha256_file(path),
        "width_px": width,
        "height_px": height,
        "pixel_standard_deviation": standard_deviation,
        "checks": checks,
    }


class AcademicCFDViewer:
    """PyVista/VTK desktop viewer with deliberately compact academic controls."""

    def __init__(
        self,
        data: VisualData,
        config: VisualConfig,
        *,
        off_screen: bool = False,
        publication: bool = False,
    ) -> None:
        self.data = data
        self.config = config
        self.publication = publication
        self.current_field = config.initial_scalar
        self.full_range = config.full_range
        self.plane_mode = "clip"
        self.plane_origin = np.asarray(data.grid_um.center, dtype=float)
        self.plane_normal = np.array((1.0, 0.0, 0.0))
        self.show_vectors = False
        self.show_streamlines = bool(
            config.build_streamlines and data.valid_streamline_count > 0
        )
        self.show_outlet_split = False
        self.show_pressure_drop = False
        self.show_help = not publication
        self.show_edges = False
        size = (
            max(config.width, 2400) if publication else config.width,
            max(config.height, 1600) if publication else config.height,
        )
        self.plotter = pv.Plotter(off_screen=off_screen or publication, window_size=size)
        self.plotter.set_background(BACKGROUND)
        self.field_actor: Any = None
        self.context_actor: Any = None
        self.glyph_mesh: pv.PolyData | None = None
        self._build_scene()

    def _current_limits(self) -> tuple[float, float]:
        return self.data.ranges[self.current_field].selected(self.full_range)

    def _field_subset(self) -> pv.DataSet:
        if self.plane_mode == "slice":
            return self.data.grid_um.slice(
                normal=self.plane_normal, origin=self.plane_origin
            )
        return self.data.grid_um.clip(
            normal=self.plane_normal, origin=self.plane_origin, invert=False
        )

    def _build_scene(self) -> None:
        self.context_actor = self.plotter.add_mesh(
            self.data.surface_um,
            color=CONTEXT_COLOR,
            opacity=0.18,
            show_edges=False,
            smooth_shading=False,
            pickable=False,
            name="vessel_context",
        )
        self._replace_field_actor(render=False)
        if self.config.show_ports:
            self._add_ports()
        if self.show_streamlines:
            self._add_streamlines()
        if self.config.debug_cells:
            self.plotter.add_points(
                self.data.centers_um[:: max(1, self.data.grid_um.n_cells // 12_000)],
                color="#555555",
                opacity=0.16,
                point_size=2,
                name="debug_cell_centres",
                pickable=False,
            )
        self.plotter.add_axes(line_width=2, color="#333333")
        self._update_overlays()
        self._set_isometric()
        if not self.publication:
            self._add_interactions()

    def _replace_field_actor(self, *, render: bool = True) -> None:
        subset = self._field_subset()
        if subset.n_cells == 0 and subset.n_points == 0:
            return
        try:
            self.plotter.remove_scalar_bar(render=False)
        except (IndexError, StopIteration):
            pass
        field = self.data.fields[self.current_field]
        clim = self._current_limits()
        cmap = (
            _pressure_colormap(*clim)
            if self.current_field == "pressure"
            else field.cmap
        )
        self.field_actor = self.plotter.add_mesh(
            subset,
            scalars=field.array,
            preference="cell",
            cmap=cmap,
            clim=clim,
            show_edges=False,
            smooth_shading=False,
            nan_color="#BBBBBB",
            scalar_bar_args={
                "title": f"{field.title} ({field.units})",
                "vertical": True,
                "position_x": 0.86,
                "position_y": 0.18,
                "width": 0.09,
                "height": 0.58,
                "title_font_size": 16,
                "label_font_size": 13,
                "color": "#222222",
                "fmt": "%.4g",
            },
            name="cfd_field",
            reset_camera=False,
        )
        self._update_field_overlay(render=render)

    def _plane_callback(self, normal: Sequence[float], origin: Sequence[float]) -> None:
        self.plane_normal = np.asarray(normal, dtype=float)
        self.plane_origin = np.asarray(origin, dtype=float)
        subset = self._field_subset()
        if self.field_actor is not None and (subset.n_cells or subset.n_points):
            self.field_actor.mapper.dataset = subset
            self.plotter.render()

    def _add_interactions(self) -> None:
        self.plotter.add_plane_widget(
            self._plane_callback,
            normal=self.plane_normal,
            origin=self.plane_origin,
            bounds=self.data.grid_um.bounds,
            factor=1.05,
            color="#666666",
            tubing=False,
            outline_translation=False,
            implicit=True,
            interaction_event="always",
        )
        callbacks = {
            "1": lambda: self.set_field("velocity"),
            "2": lambda: self.set_field("pressure"),
            "3": lambda: self.set_field("rho"),
            "c": lambda: self.set_plane_mode("clip"),
            "l": lambda: self.set_plane_mode("slice"),
            "v": self.toggle_vectors,
            "t": self.toggle_streamlines,
            "o": self.toggle_outlet_split,
            "d": self.toggle_pressure_drop,
            "i": self._set_isometric,
            "x": lambda: self._view_axis((1.0, 0.0, 0.0)),
            "y": lambda: self._view_axis((0.0, 1.0, 0.0)),
            "z": lambda: self._view_axis((0.0, 0.0, 1.0)),
            "r": self.reset_camera,
            "f": self.focus_vessel,
            "a": lambda: self.set_full_range(False),
            "m": lambda: self.set_full_range(True),
            "s": self.save_interactive_screenshot,
            "h": self.toggle_help,
            "e": self.toggle_edges,
            "q": self.plotter.close,
        }
        if self.config.numerical_pressure_debug:
            callbacks["4"] = lambda: self.set_field("numerical-pressure")
        for key, callback in callbacks.items():
            self.plotter.add_key_event(key, callback)
        self.plotter.enable_surface_point_picking(
            callback=self._pick_callback,
            show_message=False,
            show_point=True,
            point_size=10,
            color="#222222",
            left_clicking=False,
            picker="cell",
        )

    def _add_ports(self) -> None:
        label_points: list[np.ndarray] = []
        labels: list[str] = []
        for label in PORT_ORDER:
            plane = self.data.plane_contract["ports"][label]["planes"]["central"]
            origin = np.asarray(plane["origin_m"], dtype=float) * 1.0e6
            normal = np.asarray(plane["unit_normal"], dtype=float)
            basis_u = np.asarray(plane["basis_u"], dtype=float)
            basis_v = np.asarray(plane["basis_v"], dtype=float)
            contour = np.asarray(plane["physical_aperture_contour_uv_m"], dtype=float)
            contour_points = (
                origin
                + contour[:, :1] * 1.0e6 * basis_u
                + contour[:, 1:] * 1.0e6 * basis_v
            )
            disk = pv.PolyData(
                contour_points,
                faces=np.r_[len(contour_points), np.arange(len(contour_points))],
            )
            outline = pv.lines_from_points(
                np.vstack((contour_points, contour_points[0])), close=False
            )
            self.plotter.add_mesh(
                disk,
                color=PORT_COLORS[label],
                opacity=0.32,
                lighting=False,
                name=f"port_disk_{label}",
                pickable=False,
            )
            self.plotter.add_mesh(
                outline,
                color=PORT_COLORS[label],
                line_width=4,
                name=f"port_outline_{label}",
                pickable=False,
            )
            flow_direction = -normal if label == "inlet" else normal
            diameter_um = float(plane["local_hydraulic_diameter_m"]) * 1.0e6
            arrow_length = max(2.2 * diameter_um, 2.0)
            arrow = pv.Arrow(
                start=origin,
                direction=flow_direction,
                scale=arrow_length,
                tip_length=0.28,
                tip_radius=0.09,
                shaft_radius=0.025,
            )
            self.plotter.add_mesh(
                arrow,
                color=PORT_COLORS[label],
                name=f"port_arrow_{label}",
                pickable=False,
            )
            label_points.append(origin + flow_direction * arrow_length * 1.2)
            labels.append(PORT_LABELS[label])
        self.plotter.add_point_labels(
            np.vstack(label_points),
            labels,
            font_size=13,
            text_color="#222222",
            font_family="arial",
            show_points=False,
            shape="rounded_rect",
            shape_color="#F7F5EF",
            shape_opacity=0.85,
            always_visible=True,
            name="port_labels",
            pickable=False,
        )

    def _add_streamlines(self) -> None:
        if self.data.streamlines_um is None or self.data.streamlines_um.n_points == 0:
            return
        clim = self.data.ranges["velocity"].selected(self.full_range)
        self.plotter.add_mesh(
            self.data.streamlines_um,
            scalars="velocity_magnitude_mm_s",
            cmap="viridis",
            clim=clim,
            line_width=4,
            render_lines_as_tubes=True,
            show_scalar_bar=False,
            name="streamlines",
            pickable=False,
        )

    def _build_glyphs(self) -> pv.PolyData:
        target = min(2200, self.data.grid_um.n_cells)
        indices = np.linspace(0, self.data.grid_um.n_cells - 1, target, dtype=int)
        velocity = np.asarray(self.data.grid_um.cell_data["velocity_phy"])[indices]
        magnitude = np.linalg.norm(velocity, axis=1)
        nonzero = magnitude > np.finfo(float).tiny
        points = self.data.centers_um[indices][nonzero]
        directions = velocity[nonzero] / magnitude[nonzero, None]
        cloud = pv.PolyData(points)
        cloud.point_data["velocity_direction"] = directions
        cloud.point_data["velocity_magnitude_mm_s"] = magnitude[nonzero] * 1.0e3
        dx_um = float(self.data.run_summary["numerical_contract"]["dx_m"]) * 1.0e6
        arrow_size = max(3.0 * dx_um, float(self.data.grid_um.length) / 150.0)
        return cloud.glyph(
            orient="velocity_direction", scale=False, factor=arrow_size, tolerance=None
        )

    def _update_field_overlay(self, *, render: bool = True) -> None:
        field = self.data.fields[self.current_field]
        limits = self.data.ranges[self.current_field]
        display_min, display_max = limits.selected(self.full_range)
        warning = ""
        if self.current_field == "pressure":
            warning = (
                "\nGauge pressure shown. The ~3.39 MPa reference is an LBM numerical offset."
            )
        elif self.current_field == "numerical-pressure":
            warning = "\nNUMERICAL OFFSET INCLUDED — NOT PHYSIOLOGICAL PRESSURE."
        text = (
            f"CURRENT FIELD: {field.title} [{field.units}]\n"
            f"RAW RANGE: {limits.raw_min:.7g} to {limits.raw_max:.7g}\n"
            f"DISPLAY RANGE: {display_min:.7g} to {display_max:.7g} "
            f"({'FULL' if self.full_range else 'p1–p99'})"
            f"{warning}"
        )
        self.plotter.add_text(
            text,
            position=(15, self.plotter.window_size[1] - 115),
            font_size=11,
            color="#222222",
            font="arial",
            name="field_information",
            render=render,
        )

    def _scientific_summary(self) -> str:
        metrics = self.data.metrics
        contract = self.data.run_summary["numerical_contract"]
        target = float(contract["target_volume_flow_m3_s"]) * 60.0e12
        measured = float(metrics["Qin_m3_s"]) * 60.0e12
        return (
            "PRODUCTION CFD — VALIDATED BASE\n"
            f"iteration: {int(metrics['iteration'])}\n"
            f"dx: {float(contract['dx_m']) * 1e6:.2f} µm   tau: {float(contract['tau']):.1f}\n"
            f"Q target: {target:.10f} nL/min\n"
            f"Q measured: {measured:.10f} nL/min\n"
            f"closure: {float(metrics['physical_volume_closure']):.7g}\n"
            f"rho mean: {float(metrics['rho_mean']):.10f}\n"
            "STEADY SOURCE: VALIDATED RESEARCH BASE ACCEPTED RESTART\n"
            "fresh full production steady solve: NO\n"
            "Resolution: Coarse→Base sensitivity PASS; formal C/B/F GCI not completed\n"
            "WSS: not part of this validated visualization."
        )

    def _help_text(self) -> str:
        fourth = "\n4  Numerical pressure DEBUG" if self.config.numerical_pressure_debug else ""
        return (
            "SHORTCUTS\n"
            "1 Velocity   2 Gauge pressure   3 Density"
            f"{fourth}\n"
            "C Clip   L Slice   V Vectors   T Streamlines\n"
            "O Outlet split   D Pressure drops   E Context edges\n"
            "I/X/Y/Z Camera   R Reset   F Focus\n"
            "A Auto p1–p99   M Full range   S Screenshot\n"
            "Right-click Pick cell   H Help   Q Quit"
        )

    def _update_overlays(self) -> None:
        self.plotter.add_text(
            "ACADEMIC PRODUCTION CFD VIEWER V3",
            position="upper_edge",
            font_size=16,
            color="#1F2933",
            font="arial",
            name="main_title",
        )
        if not self.publication:
            self.plotter.add_text(
                self._scientific_summary(),
                position=(15, self.plotter.window_size[1] - 390),
                font_size=9,
                color="#30343B",
                font="arial",
                name="scientific_summary",
            )
        else:
            metrics = self.data.metrics
            contract = self.data.run_summary["numerical_contract"]
            footer = (
                f"Validated Tau1 Base CFD | dx={float(contract['dx_m']) * 1e6:.2f} µm, "
                f"tau={float(contract['tau']):.1f} | Gauge pressure / physical velocity | "
                f"Accepted iteration {int(metrics['iteration'])}"
            )
            self.plotter.add_text(
                footer,
                position="lower_edge",
                font_size=10,
                color="#30343B",
                font="arial",
                name="publication_footer",
            )
        self._update_field_overlay(render=False)
        self._update_help()
        self._update_optional_summary()

    def _update_help(self) -> None:
        self.plotter.remove_actor("help_overlay", render=False)
        if self.show_help:
            self.plotter.add_text(
                self._help_text(),
                position=(15, 20),
                font_size=9,
                color="#222222",
                font="arial",
                name="help_overlay",
                render=False,
            )

    def _update_optional_summary(self) -> None:
        self.plotter.remove_actor("optional_summary", render=False)
        blocks: list[str] = []
        metrics = self.data.metrics
        if self.show_outlet_split:
            fractions = metrics["flow_fractions"]
            blocks.append(
                "OUTLET FLOW FRACTIONS\n"
                + "\n".join(
                    f"Outlet {index:02d}: {100.0 * float(fractions[f'outlet_{index:02d}']):.3f}%"
                    for index in range(1, 4)
                )
            )
        if self.show_pressure_drop:
            drops = metrics["pressure_drops_pa"]
            blocks.append(
                "GAUGE PRESSURE SUMMARY\n"
                f"Inlet: {float(metrics['inlet_gauge_pressure_pa']):.4f} Pa\n"
                + "\n".join(
                    f"ΔP{index:02d}: {float(drops[f'outlet_{index:02d}']):.4f} Pa"
                    for index in range(1, 4)
                )
            )
        if blocks:
            self.plotter.add_text(
                "\n\n".join(blocks),
                position=(self.plotter.window_size[0] - 390, self.plotter.window_size[1] - 300),
                font_size=10,
                color="#222222",
                font="arial",
                name="optional_summary",
                render=False,
            )

    def _pick_callback(self, point: Sequence[float] | None) -> None:
        if point is None:
            return
        coordinates = np.asarray(point, dtype=float)
        if coordinates.shape != (3,) or not np.all(np.isfinite(coordinates)):
            return
        _, index = self.data.center_tree_um.query(coordinates, k=1)
        index = int(index)
        velocity = np.asarray(self.data.grid_um.cell_data["velocity_phy"])[index] * 1.0e3
        magnitude = float(self.data.grid_um.cell_data["velocity_magnitude_mm_s"][index])
        pressure = float(self.data.grid_um.cell_data["pressure_gauge_pa"][index])
        rho = float(self.data.grid_um.cell_data["rho_lattice"][index])
        text = (
            f"PICKED CELL {index}\n"
            f"x/y/z: {coordinates[0]:.4f}, {coordinates[1]:.4f}, {coordinates[2]:.4f} µm\n"
            f"velocity: [{velocity[0]:.6g}, {velocity[1]:.6g}, {velocity[2]:.6g}] mm/s\n"
            f"|velocity|: {magnitude:.7g} mm/s\n"
            f"gauge pressure: {pressure:.7g} Pa\n"
            f"rho_lattice: {rho:.10g}"
        )
        self.plotter.add_text(
            text,
            position=(self.plotter.window_size[0] - 560, 25),
            font_size=10,
            color="#222222",
            font="arial",
            name="pick_information",
        )

    def set_field(self, field: str) -> None:
        if field not in self.data.fields:
            raise VisualizerInputError(f"Field is unavailable or debug-disabled: {field}")
        self.current_field = field
        self._replace_field_actor()

    def set_plane_mode(self, mode: str) -> None:
        if mode not in {"clip", "slice"}:
            raise ValueError(mode)
        self.plane_mode = mode
        self._replace_field_actor()

    def set_full_range(self, enabled: bool) -> None:
        self.full_range = enabled
        self._replace_field_actor()
        if self.show_streamlines:
            self._add_streamlines()

    def toggle_vectors(self) -> None:
        self.show_vectors = not self.show_vectors
        if self.show_vectors:
            if self.glyph_mesh is None:
                self.glyph_mesh = self._build_glyphs()
            self.plotter.add_mesh(
                self.glyph_mesh,
                scalars="velocity_magnitude_mm_s",
                cmap="viridis",
                clim=self.data.ranges["velocity"].selected(self.full_range),
                show_scalar_bar=False,
                name="velocity_glyphs",
                pickable=False,
            )
        else:
            self.plotter.remove_actor("velocity_glyphs")

    def toggle_streamlines(self) -> None:
        if self.data.streamlines_um is None or self.data.valid_streamline_count == 0:
            return
        self.show_streamlines = not self.show_streamlines
        if self.show_streamlines:
            self._add_streamlines()
        else:
            self.plotter.remove_actor("streamlines")

    def toggle_outlet_split(self) -> None:
        self.show_outlet_split = not self.show_outlet_split
        self._update_optional_summary()
        self.plotter.render()

    def toggle_pressure_drop(self) -> None:
        self.show_pressure_drop = not self.show_pressure_drop
        self._update_optional_summary()
        self.plotter.render()

    def toggle_help(self) -> None:
        self.show_help = not self.show_help
        self._update_help()
        self.plotter.render()

    def toggle_edges(self) -> None:
        self.show_edges = not self.show_edges
        self.context_actor = self.plotter.add_mesh(
            self.data.surface_um,
            color=CONTEXT_COLOR,
            opacity=0.20,
            show_edges=self.show_edges,
            edge_color="#777777",
            line_width=1,
            name="vessel_context",
            pickable=False,
        )

    def _set_isometric(self) -> None:
        self.plotter.view_isometric()
        self.plotter.camera.zoom(1.12)
        self.plotter.render()

    def _view_axis(self, vector: tuple[float, float, float]) -> None:
        view_up = (0.0, 0.0, 1.0) if vector != (0.0, 0.0, 1.0) else (0.0, 1.0, 0.0)
        self.plotter.view_vector(vector, viewup=view_up)
        self.plotter.render()

    def reset_camera(self) -> None:
        self.plotter.reset_camera()
        self.plotter.render()

    def focus_vessel(self) -> None:
        self.plotter.reset_camera(bounds=self.data.surface_um.bounds)
        self.plotter.camera.zoom(1.08)
        self.plotter.render()

    def _screenshot_metadata(self, path: Path, purpose: str) -> dict[str, Any]:
        field = self.data.fields[self.current_field]
        return {
            "status": "PASS",
            "purpose": purpose,
            "created_at": datetime.now().isoformat(),
            "screenshot": str(path),
            "screenshot_sha256": sha256_file(path),
            "scalar": self.current_field,
            "array": field.array,
            "units": field.units,
            "raw_range": [
                self.data.ranges[self.current_field].raw_min,
                self.data.ranges[self.current_field].raw_max,
            ],
            "display_range": list(self._current_limits()),
            "range_mode": "full" if self.full_range else "p1-p99",
            "camera": _camera_record(self.plotter),
            "source_vtu": str(self.data.vtu_path),
            "source_vtu_sha256": self.data.vtu_sha256,
            "source_classification": self.data.run_summary["steady_solution_source"],
            "source_iteration": int(self.data.metrics["iteration"]),
        }

    def save_interactive_screenshot(self) -> Path:
        output = self.data.run_dir / "visualization" / "interactive_v3"
        output.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = output / f"screenshot_{stamp}.png"
        self.plotter.screenshot(path, scale=2)
        write_json(
            path.with_suffix(".json"),
            self._screenshot_metadata(path, "interactive_v3_user_screenshot"),
        )
        print(f"Screenshot saved: {path}")
        return path

    def save_publication_screenshot(self, path: Path) -> Path:
        path = path.expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        self.show_help = False
        self._update_help()
        self._set_isometric()
        self.plotter.show(auto_close=False, interactive=False)
        self.plotter.screenshot(path)
        write_json(
            path.with_suffix(".json"),
            self._screenshot_metadata(path, "academic_publication_screenshot"),
        )
        self.plotter.close()
        return path

    def show(self) -> None:
        self.plotter.show(
            title="Academic Production CFD Viewer V3",
            interactive=True,
            auto_close=True,
        )


def run_self_test(data: VisualData, config: VisualConfig) -> tuple[dict[str, Any], Path]:
    """Exercise all three field updates and verify real off-screen renders."""

    output = data.run_dir / "visualization" / "interactive_v3"
    output.mkdir(parents=True, exist_ok=True)
    test_config = VisualConfig(
        width=max(config.width, 1600),
        height=max(config.height, 1000),
        initial_scalar="velocity",
        streamline_seeds=config.streamline_seeds,
        build_streamlines=True,
        show_ports=True,
        full_range=config.full_range,
        debug_cells=False,
        numerical_pressure_debug=config.numerical_pressure_debug,
    )
    viewer = AcademicCFDViewer(data, test_config, off_screen=True, publication=True)
    viewer.plotter.show(auto_close=False, interactive=False)
    images: dict[str, Any] = {}
    for field in FIELD_ORDER:
        viewer.set_field(field)
        viewer.plotter.render()
        path = output / f"self_test_{field}.png"
        viewer.plotter.screenshot(path)
        images[field] = _image_qc(path)
    viewer.plotter.close()
    field_ranges = {
        key: {
            "units": data.fields[key].units,
            "raw": [value.raw_min, value.raw_max],
            "p1_p99": [value.percentile_min, value.percentile_max],
        }
        for key, value in data.ranges.items()
        if key in FIELD_ORDER
    }
    checks = {
        "manifest_pass": data.manifest.get("status") == "PASS",
        "cell_count": data.grid_um.n_cells == int(data.manifest["cell_count"]),
        "surface_nonempty": data.surface_um.n_cells > 0,
        "four_ports": len(data.plane_contract["ports"]) == 4,
        "several_valid_streamlines": data.valid_streamline_count >= 4,
        "three_field_renders": all(item["status"] == "PASS" for item in images.values()),
        "finite_ranges": all(
            np.all(np.isfinite(record["raw"] + record["p1_p99"]))
            for record in field_ranges.values()
        ),
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "solver_calls": 0,
        "run_dir": str(data.run_dir),
        "vtu": str(data.vtu_path),
        "vtu_sha256": data.vtu_sha256,
        "arrays": list(REQUIRED_ARRAYS),
        "cells": data.grid_um.n_cells,
        "surface_cells": data.surface_um.n_cells,
        "valid_streamlines": data.valid_streamline_count,
        "port_count": len(data.plane_contract["ports"]),
        "field_ranges": field_ranges,
        "renders": images,
        "checks": checks,
    }
    report_path = output / "interactive_v3_self_test.json"
    write_json(report_path, report)
    if report["status"] != "PASS":
        raise RuntimeError(f"Interactive V3 self-test failed: {checks}")
    return report, report_path


def _print_startup(data: VisualData) -> None:
    print("CFD Visualizer V3")
    print(f"Run: {data.run_dir}")
    print(f"VTU: {data.vtu_path}")
    print(f"VTU SHA256: {data.vtu_sha256}")
    print(f"Cells: {data.grid_um.n_cells}")
    print(f"Source: {data.run_summary['steady_solution_source']}")
    print(f"Iteration: {int(data.metrics['iteration'])}")
    print("Fields: velocity / gauge pressure / rho")


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not 1200 <= args.window_width <= 7680 or not 800 <= args.window_height <= 4320:
        print("CFD Visualizer V3 input invalid: unreasonable window dimensions")
        print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_INPUT_INVALID")
        return 2
    if not 8 <= args.streamline_seeds <= 128:
        print("CFD Visualizer V3 input invalid: --streamline-seeds must be 8..128")
        print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_INPUT_INVALID")
        return 2
    if args.scalar == "numerical-pressure" and not args.show_numerical_pressure_debug:
        print("CFD Visualizer V3 input invalid: numerical pressure requires debug opt-in")
        print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_INPUT_INVALID")
        return 2
    try:
        run_dir, vtu = resolve_run_and_vtu(
            project_root=PROJECT_ROOT,
            explicit_run_dir=args.run_dir,
            explicit_vtu=args.vtu,
        )
        config = VisualConfig(
            width=args.window_width,
            height=args.window_height,
            initial_scalar=args.scalar,
            streamline_seeds=args.streamline_seeds,
            build_streamlines=not args.no_streamlines or args.self_test,
            show_ports=not args.no_ports or args.self_test,
            full_range=args.full_range,
            debug_cells=args.debug_cells,
            numerical_pressure_debug=args.show_numerical_pressure_debug,
        )
        data = load_and_validate_data(run_dir, vtu, config, project_root=PROJECT_ROOT)
        _print_startup(data)
        try:
            if args.self_test:
                report, path = run_self_test(data, config)
                print(f"Self-test: {report['status']}")
                print(f"Valid streamlines: {report['valid_streamlines']}")
                print(f"Report: {path}")
                print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_READY")
                return 0
            if args.publication_screenshot:
                viewer = AcademicCFDViewer(data, config, off_screen=True, publication=True)
                path = viewer.save_publication_screenshot(args.publication_screenshot)
                print(f"Publication screenshot: {path}")
                print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_READY")
                return 0
            if args.off_screen:
                viewer = AcademicCFDViewer(data, config, off_screen=True, publication=True)
                output = (
                    data.run_dir
                    / "visualization"
                    / "interactive_v3"
                    / "off_screen_preview.png"
                )
                viewer.save_publication_screenshot(output)
                print(f"Off-screen preview: {output}")
                print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_READY")
                return 0
            viewer = AcademicCFDViewer(data, config)
            print("Interactive window opening...")
            viewer.show()
        except Exception as error:  # VTK backend exceptions differ by platform
            print(f"Interactive PyVista window could not open: {error}")
            print(
                "Use a local desktop Python environment or "
                "--publication-screenshot for off-screen rendering."
            )
            print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_GUI_ENVIRONMENT_BLOCKED")
            return 3
        print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_READY")
        return 0
    except (VisualizerInputError, OSError, ValueError, KeyError) as error:
        print(f"CFD Visualizer V3 input invalid: {error}")
        print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_INPUT_INVALID")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
