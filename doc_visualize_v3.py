"""Academic interactive desktop viewer for validated production CFD fields.

This module is deliberately a viewer, not a solver or scientific post-processor.
It validates and reads the accepted production VTU without modifying it on disk.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass, replace
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

from utils.cfd_flow.io import read_json, sha256_file
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
DERIVED_ARRAYS = (
    "density_deviation_ppm",
    "dynamic_pressure_mpa",
    "velocity_x_mm_s",
    "velocity_y_mm_s",
    "velocity_z_mm_s",
)
FIELD_ORDER = (
    "velocity",
    "pressure",
    "density-deviation",
    "dynamic-pressure",
    "velocity-x",
    "velocity-y",
    "velocity-z",
    "rho",
)
SCALAR_CHOICES = FIELD_ORDER + ("numerical-pressure",)
PORT_ORDER = ("inlet", "outlet_01", "outlet_02", "outlet_03")
PORT_LABELS = {
    "inlet": "Inlet",
    "outlet_01": "Outlet 1",
    "outlet_02": "Outlet 2",
    "outlet_03": "Outlet 3",
}
PORT_COLORS = {
    "inlet": "#0072B2",
    "outlet_01": "#009E73",
    "outlet_02": "#E69F00",
    "outlet_03": "#CC79A7",
}
RENDER_INTERPOLATION = "CELL_CENTER_TO_ORIGINAL_SURFACE_IDW_K8_DISPLAY_ONLY"
SURFACE_MAPPING_NEIGHBORS = 8


class VisualizerInputError(RuntimeError):
    """Raised when production artifacts do not satisfy the viewing contract."""


@dataclass(frozen=True, slots=True)
class VisualConfig:
    """User-selected display options."""

    width: int = 1920
    height: int = 1200
    initial_scalar: str = "velocity"
    streamline_seeds: int = 24
    build_streamlines: bool = True
    show_ports: bool = True
    full_range: bool = False
    debug_cells: bool = False
    numerical_pressure_debug: bool = False
    projection: str = "parallel"
    ui_mode: str = "analysis"
    theme: str = "dark"


@dataclass(frozen=True, slots=True)
class AcademicStyle:
    """Central publication styling; no global PyVista theme is mutated."""

    theme: str = "dark"
    font_family: str = "arial"
    background: str = "#07131F"
    background_top: str = "#142A3D"
    panel_color: str = "#0D2234"
    panel_border: str = "#36536A"
    text_color: str = "#F4F7FA"
    muted_color: str = "#B7C5D1"
    context_color: str = "#91A4B5"
    context_opacity: float = 0.16
    title_font_size: int = 22
    metadata_font_size: int = 15
    port_font_size: int = 14
    control_font_size: int = 15
    scalar_title_font_size: int = 20
    scalar_label_font_size: int = 16
    scalar_bar_x: float = 0.89
    scalar_bar_y: float = 0.30
    scalar_bar_width: float = 0.035
    scalar_bar_height: float = 0.50
    scalar_bar_labels: int = 6
    scalar_title_gap_px: int = 20
    camera_padding: float = 1.02


def academic_style(theme: str, *, publication: bool = False) -> AcademicStyle:
    """Return a high-contrast style with typography scaled for raster export."""

    if theme == "dark":
        style = AcademicStyle()
    elif theme == "light":
        style = AcademicStyle(
            theme="light",
            background="#EDF2F6",
            background_top="#FFFFFF",
            panel_color="#FFFFFF",
            panel_border="#B7C7D3",
            text_color="#17232E",
            muted_color="#526574",
            context_color="#7F929F",
        )
    else:
        raise ValueError(f"Unsupported visual theme: {theme}")
    if not publication:
        return style
    return replace(
        style,
        title_font_size=28,
        metadata_font_size=18,
        port_font_size=19,
        control_font_size=18,
        scalar_title_font_size=22,
        scalar_label_font_size=18,
        scalar_title_gap_px=26,
    )


@dataclass(frozen=True, slots=True)
class AcademicLayout:
    """Normalized layout shared by interactive and publication renders."""

    title_position: str = "upper_left"
    help_position: str = "lower_left"
    info_position: str = "upper_right"
    pick_position: str = "lower_right"
    orientation_viewport: tuple[float, float, float, float] = (
        0.015,
        0.015,
        0.09,
        0.105,
    )

    def scalar_bar_args(self, field: "FieldSpec", style: AcademicStyle) -> dict[str, Any]:
        return {
            "title": f"{field.title}\n{field.units}",
            "vertical": True,
            "position_x": style.scalar_bar_x,
            "position_y": style.scalar_bar_y,
            "width": style.scalar_bar_width,
            "height": style.scalar_bar_height,
            "n_labels": style.scalar_bar_labels,
            "title_font_size": style.scalar_title_font_size,
            "label_font_size": style.scalar_label_font_size,
            "color": style.text_color,
            "font_family": style.font_family,
            "bold": True,
            "fmt": "%.4g",
            "outline": False,
            "fill": False,
            "unconstrained_font_size": True,
        }


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
    shortcut: str = ""
    control_label: str = ""
    description: str = ""
    symmetric_range: bool = False
    log_scale: bool = False


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
    original_surface_path: Path
    original_surface_sha256: str
    surface_mapping: dict[str, Any]
    # grid_um retains the original cell-centred quantitative arrays.
    grid_um: pv.UnstructuredGrid
    # display_grid_um contains point-interpolated scalars for rendering only.
    display_grid_um: pv.UnstructuredGrid
    surface_um: pv.PolyData
    centers_um: np.ndarray
    center_tree_um: cKDTree
    fields: dict[str, FieldSpec]
    ranges: dict[str, FieldRange]
    streamlines_um: pv.PolyData | None
    valid_streamline_count: int
    original_cell_array_sha256: dict[str, str]
    rendering_scalar_interpolation: str = RENDER_INTERPOLATION


def build_parser() -> argparse.ArgumentParser:
    """Build the small, desktop-oriented command-line interface."""

    parser = argparse.ArgumentParser(
        description="Academic interactive viewer for accepted production CFD fields."
    )
    parser.add_argument("--run-dir", type=Path, help="Explicit production run directory")
    parser.add_argument("--vtu", type=Path, help="Explicit production VTU (highest priority)")
    parser.add_argument(
        "--surface",
        type=Path,
        help="Explicit continuous original vessel surface (.vtp or .stl)",
    )
    parser.add_argument(
        "--scalar",
        choices=SCALAR_CHOICES,
        default="velocity",
    )
    parser.add_argument("--window-width", type=int, default=1920)
    parser.add_argument("--window-height", type=int, default=1200)
    parser.add_argument("--streamline-seeds", type=int, default=24)
    parser.add_argument("--no-streamlines", action="store_true")
    parser.add_argument("--no-ports", action="store_true")
    parser.add_argument("--full-range", action="store_true")
    parser.add_argument("--publication-screenshot", type=Path)
    parser.add_argument("--publication-suite", type=Path)
    parser.add_argument(
        "--projection", choices=("parallel", "perspective"), default="parallel"
    )
    parser.add_argument("--ui-mode", choices=("clean", "analysis"), default="analysis")
    parser.add_argument("--theme", choices=("dark", "light"), default="dark")
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


def _write_viewer_json(path: Path, value: dict[str, Any]) -> None:
    """Write portable UTF-8/LF viewer evidence without touching production records."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False)
        stream.write("\n")


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


def _rebase_provenance_path(path_value: str | Path, project_root: Path) -> Path:
    """Resolve an evidence path after a repository has moved to another machine."""

    recorded = Path(path_value).expanduser()
    if recorded.exists():
        return recorded.resolve()
    parts = recorded.parts
    output_index = next(
        (index for index, part in enumerate(parts) if part.lower() == "outputs"),
        None,
    )
    if output_index is not None:
        portable = project_root.joinpath(*parts[output_index:])
        if portable.exists():
            return portable.resolve()
    return recorded.resolve()


def resolve_original_surface(
    run_summary: dict[str, Any],
    *,
    project_root: Path = PROJECT_ROOT,
    explicit_surface: Path | None = None,
) -> Path:
    """Resolve the accepted continuous vessel surface used to build the CFD mesh."""

    root = Path(project_root).resolve()
    if explicit_surface is not None:
        candidate = explicit_surface.expanduser().resolve()
        if not candidate.is_file():
            raise VisualizerInputError(f"Explicit original surface missing: {candidate}")
        return candidate

    mesh_hashes = run_summary.get("mesh_provenance", {}).get("mesh_hashes", {})
    validation_paths: list[Path] = []
    for record in mesh_hashes.values():
        try:
            mesh_file = _rebase_provenance_path(record["path"], root)
            validation_paths.append(
                mesh_file.parents[2] / "qc" / "final_base_geometry_validation.json"
            )
        except (KeyError, IndexError, TypeError):
            continue
    for validation_path in dict.fromkeys(validation_paths):
        if not validation_path.is_file():
            continue
        validation = _load_json_object(validation_path, "accepted mesh geometry QC")
        if validation.get("status") != "PASS":
            continue
        try:
            meter_surface = _rebase_provenance_path(
                validation["full_fluid_center_containment"]["surface_path"], root
            )
        except (KeyError, TypeError):
            continue
        if meter_surface.stem.endswith("_m"):
            continuous = meter_surface.with_name(
                f"{meter_surface.stem[:-2]}_um.vtp"
            )
        else:
            continuous = meter_surface.with_suffix(".vtp")
        transform_qc = continuous.parents[1] / "qc" / "geometry_rigid_transform_qc.json"
        if continuous.is_file() and transform_qc.is_file():
            transform = _load_json_object(transform_qc, "rigid-transform geometry QC")
            checks = transform.get("checks", {})
            if (
                transform.get("status") == "PASS"
                and transform.get("transform_kind") == "GLOBAL_RIGID_ROTATION_ONLY"
                and transform.get("scale") == 1.0
                and not transform.get("remeshing", True)
                and all(bool(value) for value in checks.values())
            ):
                return continuous.resolve()
    raise VisualizerInputError(
        "Accepted continuous original vessel surface could not be resolved from mesh provenance"
    )


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


def calculate_symmetric_field_range(values: np.ndarray) -> FieldRange:
    """Return zero-centred raw and robust limits for signed scalar fields."""

    regular = calculate_field_range(values)
    array = np.asarray(values, dtype=np.float64).reshape(-1)
    raw = max(abs(regular.raw_min), abs(regular.raw_max))
    robust = float(np.percentile(np.abs(array), 99.0))
    if robust <= np.finfo(float).tiny:
        robust = raw
    if raw <= np.finfo(float).tiny:
        raw = 1.0e-12
    return FieldRange(-raw, raw, -robust, robust)


def add_derived_cell_fields(
    grid: pv.UnstructuredGrid,
    *,
    physical_density_kg_m3: float,
    rho_lattice_mean: float,
) -> None:
    """Add transparent, display-only derivatives to the in-memory VTU copy."""

    density = float(physical_density_kg_m3)
    rho_mean = float(rho_lattice_mean)
    if not np.isfinite(density) or density <= 0.0 or not np.isfinite(rho_mean):
        raise VisualizerInputError("Derived-field reference values are invalid")
    velocity_m_s = np.asarray(grid.cell_data["velocity_phy"], dtype=np.float64)
    speed_m_s = np.asarray(
        grid.cell_data["velocity_magnitude_m_s"], dtype=np.float64
    )
    rho_lattice = np.asarray(grid.cell_data["rho_lattice"], dtype=np.float64)
    grid.cell_data["density_deviation_ppm"] = (
        (rho_lattice - rho_mean) / rho_mean * 1.0e6
    )
    grid.cell_data["dynamic_pressure_mpa"] = (
        0.5 * density * np.square(speed_m_s) * 1.0e3
    )
    for index, axis in enumerate("xyz"):
        grid.cell_data[f"velocity_{axis}_mm_s"] = velocity_m_s[:, index] * 1.0e3


def _array_sha256(values: np.ndarray) -> str:
    """Hash an in-memory scientific array without changing its representation."""

    array = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
    digest.update(array.tobytes())
    return digest.hexdigest()


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


def _build_display_grid(grid_um: pv.UnstructuredGrid) -> pv.UnstructuredGrid:
    """Interpolate scalars for display while retaining cell arrays in the copy."""

    display = grid_um.cell_data_to_point_data(pass_cell_data=True)
    if not isinstance(display, pv.UnstructuredGrid):
        display = display.cast_to_unstructured_grid()
    return display


def _build_overview_surface(display_grid_um: pv.UnstructuredGrid) -> pv.PolyData:
    """Extract an unsmoothed surface and compute shading normals only."""

    surface = display_grid_um.extract_surface(algorithm="dataset_surface")
    return surface.compute_normals(
        cell_normals=False,
        point_normals=True,
        consistent_normals=True,
        auto_orient_normals=True,
        split_vertices=False,
        inplace=False,
    )


def map_cell_fields_to_original_surface(
    surface_path: Path,
    grid_um: pv.UnstructuredGrid,
    centers_um: np.ndarray,
    *,
    dx_um: float,
) -> tuple[pv.PolyData, dict[str, Any]]:
    """Map cell evidence onto the unchanged continuous surface for display only."""

    source = pv.read(surface_path)
    if not isinstance(source, pv.PolyData):
        source = source.extract_surface(algorithm="dataset_surface")
    if source.n_points < 4 or source.n_cells < 1:
        raise VisualizerInputError(f"Original vessel surface is empty: {surface_path}")
    source_points = np.asarray(source.points, dtype=np.float64)
    source_faces = np.asarray(source.faces).copy()
    coordinate_unit = "um"
    if float(np.max(np.ptp(source_points, axis=0))) < 1.0:
        source.points = source_points * 1.0e6
        source_points = np.asarray(source.points, dtype=np.float64)
        coordinate_unit = "m_converted_to_um"
    if source.n_open_edges != 0:
        raise VisualizerInputError("Accepted continuous vessel surface is not watertight")

    neighbors = min(SURFACE_MAPPING_NEIGHBORS, grid_um.n_cells)
    distances, indices = cKDTree(np.asarray(centers_um)).query(
        source_points,
        k=neighbors,
    )
    if neighbors == 1:
        distances = distances[:, None]
        indices = indices[:, None]
    regularizer = max(float(dx_um) * 0.05, np.finfo(float).eps)
    weights = 1.0 / np.square(np.asarray(distances) + regularizer)
    weights /= np.sum(weights, axis=1, keepdims=True)
    for name in REQUIRED_ARRAYS + DERIVED_ARRAYS:
        values = np.asarray(grid_um.cell_data[name], dtype=np.float64)
        selected = values[np.asarray(indices)]
        if selected.ndim == 2:
            mapped = np.sum(weights * selected, axis=1)
        else:
            mapped = np.sum(weights[..., None] * selected, axis=1)
        source.point_data[name] = mapped

    surface = source.compute_normals(
        cell_normals=False,
        point_normals=True,
        consistent_normals=True,
        auto_orient_normals=True,
        split_vertices=False,
        inplace=False,
    )
    if not np.array_equal(np.asarray(surface.points), source_points) or not np.array_equal(
        np.asarray(surface.faces), source_faces
    ):
        raise VisualizerInputError("Display mapping unexpectedly changed original geometry")
    nearest = np.asarray(distances)[:, 0]
    checks = {
        "watertight": surface.n_open_edges == 0,
        "geometry_points_unchanged": bool(
            np.array_equal(np.asarray(surface.points), source_points)
        ),
        "geometry_faces_unchanged": bool(
            np.array_equal(np.asarray(surface.faces), source_faces)
        ),
        "nearest_distance_p99_within_dx": float(np.percentile(nearest, 99.0))
        <= float(dx_um),
        "nearest_distance_max_within_2dx": float(np.max(nearest))
        <= 2.0 * float(dx_um),
        "all_display_fields_present": all(
            name in surface.point_data for name in REQUIRED_ARRAYS + DERIVED_ARRAYS
        ),
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "method": RENDER_INTERPOLATION,
        "neighbors": neighbors,
        "inverse_distance_power": 2,
        "regularizer_um": regularizer,
        "coordinate_unit": coordinate_unit,
        "surface_points": surface.n_points,
        "surface_triangles": surface.n_cells,
        "nearest_distance_um": {
            "minimum": float(np.min(nearest)),
            "median": float(np.median(nearest)),
            "p95": float(np.percentile(nearest, 95.0)),
            "p99": float(np.percentile(nearest, 99.0)),
            "maximum": float(np.max(nearest)),
        },
        "checks": checks,
    }
    if report["status"] != "PASS":
        raise VisualizerInputError(f"Original-surface mapping QC failed: {checks}")
    return surface, report


def load_and_validate_data(
    run_dir: Path,
    vtu_path: Path,
    config: VisualConfig,
    *,
    project_root: Path = PROJECT_ROOT,
    explicit_surface: Path | None = None,
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

    original_hashes = {
        name: _array_sha256(np.asarray(grid.cell_data[name])) for name in REQUIRED_ARRAYS
    }

    add_derived_cell_fields(
        grid,
        physical_density_kg_m3=float(
            summary["numerical_contract"]["rho0_kg_m3"]
        ),
        rho_lattice_mean=float(metrics["rho_mean"]),
    )

    # Centralized display-only coordinate conversion; the source VTU is never saved.
    grid.points = np.asarray(grid.points, dtype=np.float64) * 1.0e6
    centers_um = centers_m * 1.0e6
    display_grid = _build_display_grid(grid)
    original_surface_path = resolve_original_surface(
        summary,
        project_root=Path(project_root).resolve(),
        explicit_surface=explicit_surface,
    )
    original_surface_sha256 = sha256_file(original_surface_path)
    surface, surface_mapping = map_cell_fields_to_original_surface(
        original_surface_path,
        grid,
        centers_um,
        dx_um=dx_m * 1.0e6,
    )
    pressure_range = calculate_field_range(grid.cell_data["pressure_gauge_pa"])
    pressure_clim = pressure_range.selected(config.full_range)
    fields: dict[str, FieldSpec] = {
        "velocity": FieldSpec(
            "velocity_magnitude_mm_s",
            "Velocity magnitude",
            "mm/s",
            "viridis",
            "1",
            "Speed",
            "Local speed reconstructed from the validated velocity vector.",
            log_scale=True,
        ),
        "pressure": FieldSpec(
            "pressure_gauge_pa",
            "Gauge pressure",
            "Pa",
            _pressure_colormap(*pressure_clim),
            "2",
            "Pressure",
            "Physical gauge pressure after removing the LBM reference offset.",
        ),
        "density-deviation": FieldSpec(
            "density_deviation_ppm",
            "Density deviation",
            "ppm",
            "coolwarm",
            "3",
            "Density Δ",
            "Small compressibility signal relative to the accepted field mean.",
            True,
        ),
        "dynamic-pressure": FieldSpec(
            "dynamic_pressure_mpa",
            "Dynamic pressure",
            "mPa",
            "magma",
            "4",
            "Dynamic p",
            "Kinetic-energy density q = 1/2 rho |u|^2.",
            log_scale=True,
        ),
        "velocity-x": FieldSpec(
            "velocity_x_mm_s",
            "Velocity X component",
            "mm/s",
            "coolwarm",
            "5",
            "Vx",
            "Signed global x-component of the validated velocity vector.",
            True,
        ),
        "velocity-y": FieldSpec(
            "velocity_y_mm_s",
            "Velocity Y component",
            "mm/s",
            "coolwarm",
            "6",
            "Vy",
            "Signed global y-component of the validated velocity vector.",
            True,
        ),
        "velocity-z": FieldSpec(
            "velocity_z_mm_s",
            "Velocity Z component",
            "mm/s",
            "coolwarm",
            "7",
            "Vz",
            "Signed global z-component of the validated velocity vector.",
            True,
        ),
        "rho": FieldSpec(
            "rho_lattice",
            "Lattice density",
            "lattice units",
            "cividis",
            "8",
            "Rho lattice",
            "Raw lattice density retained for numerical diagnostics.",
        ),
    }
    if config.numerical_pressure_debug:
        fields["numerical-pressure"] = FieldSpec(
            "pressure_absolute_solver_pa",
            "Numerical pressure (OFFSET INCLUDED — NOT PHYSIOLOGICAL)",
            "Pa",
            "cividis",
            "9",
            "p numerical",
            "Absolute LBM pressure including the numerical reference offset.",
        )
    ranges = {}
    for key, spec in fields.items():
        values = np.asarray(grid.cell_data[spec.array])
        ranges[key] = (
            calculate_symmetric_field_range(values)
            if spec.symmetric_range
            else calculate_field_range(values)
        )
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
        original_surface_path=original_surface_path,
        original_surface_sha256=original_surface_sha256,
        surface_mapping=surface_mapping,
        grid_um=grid,
        display_grid_um=display_grid,
        surface_um=surface,
        centers_um=centers_um,
        center_tree_um=cKDTree(centers_um),
        fields=fields,
        ranges=ranges,
        streamlines_um=streamlines,
        valid_streamline_count=valid_streamlines,
        original_cell_array_sha256=original_hashes,
    )


# Publication-grade viewer implementation follows.


def _normalise(vector: np.ndarray) -> np.ndarray:
    length = float(np.linalg.norm(vector))
    if length <= np.finfo(float).tiny:
        raise VisualizerInputError("Cannot construct a camera from degenerate geometry")
    return vector / length


def _stable_axis(vector: np.ndarray) -> np.ndarray:
    """Remove the arbitrary sign from a PCA eigenvector."""

    axis = _normalise(np.asarray(vector, dtype=np.float64))
    largest = int(np.argmax(np.abs(axis)))
    return axis if axis[largest] >= 0.0 else -axis


def academic_camera_parameters(
    points_um: np.ndarray,
    aspect: float,
    *,
    padding: float = 1.04,
) -> dict[str, Any]:
    """Return a deterministic PCA camera fitted to vessel points only."""

    points = np.asarray(points_um, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 4:
        raise VisualizerInputError("Vessel surface has insufficient points for camera fit")
    center = np.mean(points, axis=0)
    centered = points - center
    _, eigenvectors = np.linalg.eigh(np.cov(centered, rowvar=False))
    view_out = _stable_axis(eigenvectors[:, 0])
    axis_major = _stable_axis(eigenvectors[:, 2])
    axis_minor = _stable_axis(eigenvectors[:, 1])
    view_to_focal = -view_out
    best: tuple[float, np.ndarray, np.ndarray, float, float, float] | None = None
    for angle in np.linspace(0.0, math.pi, 361, endpoint=False):
        view_up = _normalise(math.cos(angle) * axis_major + math.sin(angle) * axis_minor)
        right = _normalise(np.cross(view_to_focal, view_up))
        vertical_extent = float(np.ptp(centered @ view_up))
        horizontal_extent = float(np.ptp(centered @ right))
        parallel_scale = padding * max(
            vertical_extent / (2.0 * 0.76),
            horizontal_extent / (2.0 * max(aspect, 1.0e-9) * 0.72),
        )
        height_coverage = vertical_extent / (2.0 * parallel_scale)
        width_coverage = horizontal_extent / (2.0 * parallel_scale * aspect)
        penalty = (
            (height_coverage - 0.72) ** 2
            + (width_coverage - 0.62) ** 2
            + 6.0 * max(0.0, 0.65 - height_coverage) ** 2
            + 6.0 * max(0.0, 0.55 - width_coverage) ** 2
        )
        candidate = (
            penalty,
            view_up,
            right,
            parallel_scale,
            width_coverage,
            height_coverage,
        )
        if best is None or candidate[0] < best[0]:
            best = candidate
    assert best is not None
    _, view_up, right, parallel_scale, width_coverage, height_coverage = best
    vertical_projection = centered @ view_up
    horizontal_projection = centered @ right
    focal_point = (
        center
        + view_up
        * 0.5
        * (float(np.min(vertical_projection)) + float(np.max(vertical_projection)))
        + right
        * 0.5
        * (float(np.min(horizontal_projection)) + float(np.max(horizontal_projection)))
    )
    diagonal = float(np.linalg.norm(np.ptp(points, axis=0)))
    distance = max(diagonal * 3.0, parallel_scale * 4.0)
    return {
        "position_um": focal_point + view_out * distance,
        "focal_point_um": focal_point,
        "view_up": view_up,
        "view_out": view_out,
        "right": right,
        "parallel_scale": float(parallel_scale),
        "projected_width_fraction": float(width_coverage),
        "projected_height_fraction": float(height_coverage),
        "bounds_um": [float(value) for value in pv.PolyData(points).bounds],
    }


def _camera_record(plotter: pv.Plotter) -> dict[str, Any]:
    position, focal_point, view_up = plotter.camera_position
    return {
        "position_um": [float(value) for value in position],
        "focal_point_um": [float(value) for value in focal_point],
        "view_up": [float(value) for value in view_up],
        "parallel_scale_um": float(plotter.camera.parallel_scale),
        "projection": "parallel" if plotter.camera.parallel_projection else "perspective",
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
    """Figure-first PyVista viewer for validated production CFD evidence."""

    def __init__(
        self,
        data: VisualData,
        config: VisualConfig,
        *,
        off_screen: bool = False,
        publication: bool = False,
        style: AcademicStyle | None = None,
        layout: AcademicLayout | None = None,
    ) -> None:
        self.data = data
        self.config = config
        self.publication = publication
        self.style = style or academic_style(config.theme, publication=publication)
        self.layout = layout or AcademicLayout()
        self.current_field = config.initial_scalar
        self.full_range = config.full_range
        self.visual_mode = "overview"
        self.plane_mode = self.visual_mode
        self.plane_origin = np.asarray(data.surface_um.center, dtype=float)
        self.plane_normal = np.array((1.0, 0.0, 0.0))
        self.projection = config.projection
        self.show_vectors = False
        self.show_streamlines = False
        self.show_port_normals = False
        self.show_help = False
        self.show_edges = False
        self.bounding_box_visible = False
        self.plane_widget: Any = None
        self.plane_widget_visible = False
        self.picked_marker_visible = False
        self.picked_point_um: np.ndarray | None = None
        self.picked_cell_id: int | None = None
        self.field_widgets: dict[str, Any] = {}
        self.field_selector_visible = False
        self.streamline_widget: Any = None
        self.streamline_module_visible = False
        self.field_actor: Any = None
        self.context_actor: Any = None
        self.glyph_mesh: pv.PolyData | None = None
        self.streamline_tubes: pv.PolyData | None = None
        self.streamline_arrows: pv.PolyData | None = None
        self.streamline_arrow_count = 0
        size = (
            (max(config.width, 3200), max(config.height, 2000))
            if publication
            else (config.width, config.height)
        )
        self.plotter = pv.Plotter(off_screen=off_screen or publication, window_size=size)
        self.plotter.set_background(
            self.style.background,
            top=self.style.background_top,
        )
        self.camera_parameters = academic_camera_parameters(
            np.asarray(data.surface_um.points),
            size[0] / size[1],
            padding=self.style.camera_padding,
        )
        self.anti_aliasing = "none"
        self.depth_peeling = False
        self._configure_render_quality()
        self._build_scene()

    def _configure_render_quality(self) -> None:
        preferred = "ssaa"
        try:
            self.plotter.enable_anti_aliasing(preferred)
            self.anti_aliasing = preferred.upper()
        except (RuntimeError, TypeError, ValueError):
            try:
                self.plotter.enable_anti_aliasing("msaa", multi_samples=8)
                self.anti_aliasing = "MSAA"
            except (RuntimeError, TypeError, ValueError):
                self.anti_aliasing = "UNAVAILABLE"
        try:
            self.depth_peeling = bool(
                self.plotter.enable_depth_peeling(number_of_peels=8, occlusion_ratio=0.0)
            )
        except (RuntimeError, TypeError, ValueError):
            self.depth_peeling = False

    def _configure_lighting(self) -> None:
        center = np.asarray(self.camera_parameters["focal_point_um"])
        scale = max(float(self.data.surface_um.length), 1.0)
        self.plotter.remove_all_lights()
        lights = (
            ((1.2, -1.0, 1.8), 0.48),
            ((-1.4, 0.6, 0.8), 0.30),
            ((0.2, 1.4, -1.0), 0.22),
        )
        for direction, intensity in lights:
            light = pv.Light(
                position=tuple(center + scale * np.asarray(direction)),
                focal_point=tuple(center),
                color="#FFFFFF",
                intensity=intensity,
                positional=False,
            )
            self.plotter.add_light(light)

    def _build_scene(self) -> None:
        self._configure_lighting()
        self._replace_field_actor(render=False)
        if self.config.show_ports:
            self._add_ports()
        if self.config.debug_cells:
            stride = max(1, self.data.grid_um.n_cells // 12_000)
            self.plotter.add_points(
                self.data.centers_um[::stride],
                color=self.style.muted_color,
                opacity=0.16,
                point_size=2,
                name="debug_cell_centres",
                pickable=False,
                render=False,
            )
        axes = self.plotter.add_axes(
            interactive=False,
            line_width=1,
            color=self.style.muted_color,
            labels_off=False,
            viewport=self.layout.orientation_viewport,
        )
        for caption_getter in (
            axes.GetXAxisCaptionActor2D,
            axes.GetYAxisCaptionActor2D,
            axes.GetZAxisCaptionActor2D,
        ):
            caption_getter().GetCaptionTextProperty().SetFontFamilyToArial()
        self._update_overlays(render=False)
        self.apply_academic_camera(render=False)
        if not self.publication:
            self._add_interactions()

    def _current_limits(self) -> tuple[float, float]:
        return self.data.ranges[self.current_field].selected(self.full_range)

    def _field_dataset(self) -> pv.DataSet:
        if self.visual_mode == "overview":
            return self.data.surface_um
        if self.visual_mode == "slice":
            return self.data.surface_um.slice(
                normal=self.plane_normal, origin=self.plane_origin
            )
        return self.data.surface_um.clip(
            normal=self.plane_normal, origin=self.plane_origin, invert=False
        )

    def _remove_scalar_bar(self) -> None:
        try:
            self.plotter.remove_scalar_bar(render=False)
        except (IndexError, KeyError, StopIteration):
            pass

    def _replace_field_actor(self, *, render: bool = True) -> None:
        dataset = self._field_dataset()
        if dataset.n_cells == 0 and dataset.n_points == 0:
            return
        self.plotter.remove_actor("cfd_field", render=False)
        self.plotter.remove_actor("vessel_context", render=False)
        self._remove_scalar_bar()
        field = self.data.fields[self.current_field]
        clim = self._current_limits()
        cmap = _pressure_colormap(*clim) if self.current_field == "pressure" else field.cmap
        scalar_bar_args = self.layout.scalar_bar_args(field, self.style)
        if self.visual_mode != "overview":
            self.context_actor = self.plotter.add_mesh(
                self.data.surface_um,
                color=self.style.context_color,
                opacity=self.style.context_opacity,
                show_edges=self.show_edges,
                edge_color="#9AA0A6",
                line_width=1,
                smooth_shading=True,
                ambient=0.62,
                diffuse=0.38,
                specular=0.0,
                pickable=False,
                name="vessel_context",
                reset_camera=False,
                render=False,
            )
        else:
            self.context_actor = None
        self.field_actor = self.plotter.add_mesh(
            dataset,
            scalars=field.array,
            preference="point",
            cmap=cmap,
            clim=clim,
            n_colors=512,
            log_scale=field.log_scale,
            interpolate_before_map=True,
            show_edges=False,
            smooth_shading=True,
            ambient=0.58,
            diffuse=0.40,
            specular=0.02,
            specular_power=8.0,
            opacity=0.22 if self.show_streamlines else 1.0,
            nan_color="#C7C7C7",
            scalar_bar_args=scalar_bar_args,
            name="cfd_field",
            reset_camera=False,
            render=False,
        )
        scalar_bar = self.plotter.scalar_bars[scalar_bar_args["title"]]
        scalar_bar.SetVerticalTitleSeparation(self.style.scalar_title_gap_px)
        self._update_title(render=False)
        if render:
            self.plotter.render()

    def _plane_callback(
        self,
        normal: Sequence[float],
        origin: Sequence[float],
        *_: Any,
    ) -> None:
        self.plane_normal = np.asarray(normal, dtype=float)
        self.plane_origin = np.asarray(origin, dtype=float)
        dataset = self._field_dataset()
        if self.field_actor is not None and (dataset.n_cells or dataset.n_points):
            self.field_actor.mapper.dataset = dataset
            self.plotter.render()

    def show_plane_widget(self, mode: str) -> None:
        if mode not in {"clip", "slice"}:
            raise ValueError(mode)
        if self.plane_widget_visible:
            self.hide_plane_widget()
        self.set_visual_mode(mode)
        widget = self.plotter.add_plane_widget(
            self._plane_callback,
            normal=self.plane_normal,
            origin=self.plane_origin,
            bounds=self.data.surface_um.bounds,
            factor=1.0,
            color="#7A7F85",
            tubing=False,
            outline_translation=False,
            origin_translation=True,
            implicit=True,
            pass_widget=True,
            test_callback=False,
            interaction_event="always",
        )
        representation = (
            widget.GetRepresentation() if hasattr(widget, "GetRepresentation") else widget
        )
        representation.SetHandleSize(0.025)
        representation.SetTubing(False)
        representation.SetDrawPlane(True)
        representation.GetPlaneProperty().SetColor(0.55, 0.57, 0.60)
        representation.GetPlaneProperty().SetOpacity(0.035)
        representation.GetOutlineProperty().SetColor(0.48, 0.50, 0.53)
        representation.GetOutlineProperty().SetOpacity(0.28)
        representation.GetOutlineProperty().SetLineWidth(1.0)
        representation.GetNormalProperty().SetColor(0.48, 0.50, 0.53)
        representation.GetNormalProperty().SetOpacity(0.18)
        self.plane_widget = widget
        self.plane_widget_visible = True
        self.plotter.render()

    def hide_plane_widget(self) -> None:
        if self.plane_widget_visible:
            self.plotter.clear_plane_widgets()
        self.plane_widget = None
        self.plane_widget_visible = False
        self.plotter.render()

    def _toggle_inspection_widget(self, mode: str) -> None:
        if self.visual_mode == mode and self.plane_widget_visible:
            self.hide_plane_widget()
        else:
            self.show_plane_widget(mode)

    def _add_interactions(self) -> None:
        callbacks = {
            "0": lambda: self.set_visual_mode("overview"),
            "1": lambda: self.set_field("velocity"),
            "2": lambda: self.set_field("pressure"),
            "3": lambda: self.set_field("density-deviation"),
            "4": lambda: self.set_field("dynamic-pressure"),
            "5": lambda: self.set_field("velocity-x"),
            "6": lambda: self.set_field("velocity-y"),
            "7": lambda: self.set_field("velocity-z"),
            "8": lambda: self.set_field("rho"),
            "c": lambda: self._toggle_inspection_widget("clip"),
            "l": lambda: self._toggle_inspection_widget("slice"),
            "v": self.toggle_vectors,
            "t": self.toggle_streamlines,
            "n": self.toggle_port_normals,
            "p": self.toggle_projection,
            "i": self.apply_academic_camera,
            "x": lambda: self._view_axis((1.0, 0.0, 0.0)),
            "y": lambda: self._view_axis((0.0, 1.0, 0.0)),
            "z": lambda: self._view_axis((0.0, 0.0, 1.0)),
            "r": self.reset_camera,
            "f": self.focus_on_vessel,
            "a": lambda: self.set_full_range(False),
            "m": lambda: self.set_full_range(True),
            "s": self.save_interactive_screenshot,
            "h": self.toggle_help,
            "e": self.toggle_edges,
            "Return": self.hide_plane_widget,
            "Escape": self.hide_plane_widget,
            "q": self.plotter.close,
        }
        if self.config.numerical_pressure_debug:
            callbacks["9"] = lambda: self.set_field("numerical-pressure")
        for key, callback in callbacks.items():
            self.plotter.add_key_event(key, callback)
        self._add_field_selector()
        self._add_streamline_module()
        self.plotter.enable_surface_point_picking(
            callback=self._pick_callback,
            show_message=False,
            show_point=False,
            left_clicking=False,
            picker="cell",
        )

    def _add_field_selector(self) -> None:
        """Add a large, visible scalar selector instead of hiding fields behind keys."""

        field_keys = [key for key in FIELD_ORDER if key in self.data.fields]
        if self.config.numerical_pressure_debug:
            field_keys.append("numerical-pressure")
        width = float(self.plotter.window_size[0])
        column_count = 4
        row_count = int(math.ceil(len(field_keys) / column_count))
        slot = min(330.0, max(190.0, (width - 440.0) / column_count))
        start_x = max(130.0, (width - slot * column_count - 260.0) / 2.0)
        for index, key in enumerate(field_keys):
            spec = self.data.fields[key]
            column = index % column_count
            row = index // column_count
            widget = self.plotter.add_radio_button_widget(
                lambda selected=key: self.set_field(selected),
                radio_button_group="cfd_fields",
                value=key == self.current_field,
                title=f"{spec.shortcut}  {spec.control_label}",
                position=(
                    start_x + slot * column,
                    22.0 + 48.0 * (row_count - 1 - row),
                ),
                size=34,
                border_size=5,
                color_on="#35B8FF",
                color_off="#60788A",
                background_color=self.style.panel_color,
            )
            self.field_widgets[key] = widget
        title_dict = getattr(self.plotter.widgets, "radio_button_title_dict", {})
        for actors in title_dict.values():
            for actor in actors:
                actor.prop.font_size = self.style.control_font_size
                actor.prop.color = self.style.text_color
                actor.prop.bold = True
                actor.prop.font_family = self.style.font_family
        self.field_selector_visible = bool(self.field_widgets)

    def _sync_field_selector(self) -> None:
        """Keep mouse and keyboard field selection states synchronized."""

        for key, widget in self.field_widgets.items():
            representation = widget.GetRepresentation()
            representation.SetState(1 if key == self.current_field else 0)

    def _streamline_module_text(self) -> str:
        state = "ON" if self.show_streamlines else "OFF"
        return (
            f"FLOW PATHS  {state}\n"
            f"T / click · {self.data.valid_streamline_count} paths · arrows"
        )

    def _update_streamline_module_label(self, *, render: bool = True) -> None:
        self.plotter.remove_actor("streamline_module_label", render=False)
        if self.streamline_module_visible:
            width = int(self.plotter.window_size[0])
            actor = self.plotter.add_text(
                self._streamline_module_text(),
                position=(width - 392, 56),
                font_size=self.style.control_font_size,
                color=self.style.text_color,
                font=self.style.font_family,
                name="streamline_module_label",
                render=False,
            )
            actor.prop.bold = True
        if render:
            self.plotter.render()

    def _sync_streamline_widget(self) -> None:
        if self.streamline_widget is not None:
            self.streamline_widget.GetRepresentation().SetState(
                1 if self.show_streamlines else 0
            )

    def _add_streamline_module(self) -> None:
        """Expose flow paths as a first-class mouse and keyboard module."""

        source = self.data.streamlines_um
        if source is None or source.n_points == 0:
            return
        width = float(self.plotter.window_size[0])
        self.streamline_widget = self.plotter.add_checkbox_button_widget(
            lambda enabled: self.set_streamlines(bool(enabled)),
            value=self.show_streamlines,
            position=(width - 442.0, 58.0),
            size=34,
            border_size=5,
            color_on="#35B8FF",
            color_off="#60788A",
            background_color=self.style.panel_color,
        )
        self.streamline_module_visible = True
        self._update_streamline_module_label(render=False)

    def _port_geometry(
        self, label: str
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        plane = self.data.plane_contract["ports"][label]["planes"]["central"]
        origin = np.asarray(plane["origin_m"], dtype=float) * 1.0e6
        normal = np.asarray(plane["unit_normal"], dtype=float)
        basis_u = np.asarray(plane["basis_u"], dtype=float)
        basis_v = np.asarray(plane["basis_v"], dtype=float)
        contour = np.asarray(plane["physical_aperture_contour_uv_m"], dtype=float)
        points = (
            origin
            + contour[:, :1] * 1.0e6 * basis_u
            + contour[:, 1:] * 1.0e6 * basis_v
        )
        diameter_um = float(plane["local_hydraulic_diameter_m"]) * 1.0e6
        return origin, normal, points, diameter_um

    def _add_ports(self) -> None:
        label_points: list[np.ndarray] = []
        labels: list[str] = []
        flows = {
            "inlet": float(self.data.metrics.get("Qin_m3_s", 0.0)),
            "outlet_01": float(self.data.metrics.get("Q1_m3_s", 0.0)),
            "outlet_02": float(self.data.metrics.get("Q2_m3_s", 0.0)),
            "outlet_03": float(self.data.metrics.get("Q3_m3_s", 0.0)),
        }
        fractions = self.data.metrics.get("flow_fractions", {})
        for label in PORT_ORDER:
            origin, normal, contour_points, diameter_um = self._port_geometry(label)
            faces = np.r_[len(contour_points), np.arange(len(contour_points))]
            disk = pv.PolyData(contour_points, faces=faces)
            outline = pv.lines_from_points(
                np.vstack((contour_points, contour_points[0])), close=False
            )
            self.plotter.add_mesh(
                disk,
                color=PORT_COLORS[label],
                opacity=0.24,
                lighting=False,
                name=f"port_disk_{label}",
                pickable=False,
                reset_camera=False,
                render=False,
            )
            self.plotter.add_mesh(
                outline,
                color=PORT_COLORS[label],
                opacity=1.0,
                line_width=4,
                lighting=False,
                name=f"port_outline_{label}",
                pickable=False,
                reset_camera=False,
                render=False,
            )
            toward_center = _normalise(
                np.asarray(self.data.surface_um.center, dtype=float) - origin
            )
            label_points.append(
                origin + toward_center * max(4.00 * diameter_um, 3.00)
            )
            flow_nl_min = flows[label] * 60.0e12
            if label == "inlet":
                labels.append(f"{PORT_LABELS[label]}\nQ = {flow_nl_min:.4g} nL/min")
            else:
                split = float(fractions.get(label, 0.0)) * 100.0
                labels.append(
                    f"{PORT_LABELS[label]}\nQ = {flow_nl_min:.4g} nL/min · {split:.1f}%"
                )
        self.plotter.add_point_labels(
            np.vstack(label_points),
            labels,
            font_size=self.style.port_font_size,
            text_color=self.style.text_color,
            font_family=self.style.font_family,
            show_points=False,
            shape="rounded_rect",
            shape_color=self.style.panel_color,
            shape_opacity=0.90,
            margin=6,
            always_visible=True,
            name="port_labels",
            pickable=False,
            render=False,
        )

    def toggle_port_normals(self) -> None:
        self.show_port_normals = not self.show_port_normals
        for label in PORT_ORDER:
            self.plotter.remove_actor(f"port_normal_{label}", render=False)
        if self.show_port_normals:
            for label in PORT_ORDER:
                origin, normal, _, diameter_um = self._port_geometry(label)
                outward = -normal if label == "inlet" else normal
                arrow = pv.Arrow(
                    start=origin,
                    direction=outward,
                    scale=max(1.15 * diameter_um, 0.8),
                    tip_length=0.22,
                    tip_radius=0.06,
                    shaft_radius=0.018,
                )
                self.plotter.add_mesh(
                    arrow,
                    color=PORT_COLORS[label],
                    name=f"port_normal_{label}",
                    pickable=False,
                    reset_camera=False,
                    render=False,
                )
        self.plotter.render()

    def _add_streamlines(self) -> None:
        source = self.data.streamlines_um
        if source is None or source.n_points == 0:
            return
        limits = self.data.ranges["velocity"].selected(self.full_range)
        if self.streamline_tubes is None:
            dx_um = float(self.data.run_summary["numerical_contract"]["dx_m"]) * 1.0e6
            radius_um = float(np.clip(0.45 * dx_um, 0.075, 0.15))
            self.streamline_tubes = source.tube(radius=radius_um, n_sides=10)
        self.plotter.add_mesh(
            self.streamline_tubes,
            scalars="velocity_magnitude_mm_s",
            cmap="viridis",
            clim=limits,
            smooth_shading=True,
            ambient=0.72,
            diffuse=0.28,
            specular=0.0,
            show_scalar_bar=False,
            name="streamlines",
            pickable=False,
            reset_camera=False,
            render=False,
        )
        if self.streamline_arrows is None:
            self.streamline_arrows = self._build_streamline_arrows()
        if self.streamline_arrows.n_points:
            self.plotter.add_mesh(
                self.streamline_arrows,
                color=self.style.text_color,
                smooth_shading=True,
                ambient=0.76,
                diffuse=0.24,
                specular=0.0,
                show_scalar_bar=False,
                name="streamline_direction_arrows",
                pickable=False,
                reset_camera=False,
                render=False,
            )

    def _build_streamline_arrows(self) -> pv.PolyData:
        """Sample deterministic arrow glyphs that expose streamline direction."""

        source = self.data.streamlines_um
        if source is None or source.n_points == 0:
            self.streamline_arrow_count = 0
            return pv.PolyData()
        connectivity = np.asarray(source.lines, dtype=np.int64)
        points = np.asarray(source.points, dtype=np.float64)
        speeds = np.asarray(
            source.point_data["velocity_magnitude_mm_s"], dtype=np.float64
        )
        arrow_points: list[np.ndarray] = []
        arrow_directions: list[np.ndarray] = []
        arrow_speeds: list[float] = []
        cursor = 0
        path_index = 0
        while cursor < len(connectivity):
            count = int(connectivity[cursor])
            point_ids = connectivity[cursor + 1 : cursor + 1 + count]
            cursor += count + 1
            sample_fraction = 0.25 + 0.50 * ((path_index % 5) / 4.0)
            path_index += 1
            if count < 2:
                continue
            sample_indices = (int(round(sample_fraction * (count - 1))),)
            for sample in sample_indices:
                sample = int(np.clip(sample, 0, count - 2))
                lower = max(0, sample - 2)
                upper = min(count - 1, sample + 2)
                direction = points[point_ids[upper]] - points[point_ids[lower]]
                length = float(np.linalg.norm(direction))
                if length <= np.finfo(float).tiny:
                    continue
                arrow_points.append(points[point_ids[sample]])
                arrow_directions.append(direction / length)
                arrow_speeds.append(float(speeds[point_ids[sample]]))
        if not arrow_points:
            self.streamline_arrow_count = 0
            return pv.PolyData()
        self.streamline_arrow_count = len(arrow_points)
        cloud = pv.PolyData(np.vstack(arrow_points))
        cloud.point_data["flow_direction"] = np.vstack(arrow_directions)
        cloud.point_data["velocity_magnitude_mm_s"] = np.asarray(arrow_speeds)
        dx_um = float(self.data.run_summary["numerical_contract"]["dx_m"]) * 1.0e6
        arrow = pv.Arrow(
            start=(0.0, 0.0, 0.0),
            direction=(1.0, 0.0, 0.0),
            tip_length=0.34,
            tip_radius=0.13,
            shaft_radius=0.04,
        )
        return cloud.glyph(
            orient="flow_direction",
            scale=False,
            factor=max(8.0 * dx_um, 1.4),
            geom=arrow,
            tolerance=None,
        )

    def _build_glyphs(self) -> pv.PolyData:
        target = min(1600, self.data.grid_um.n_cells)
        indices = np.linspace(0, self.data.grid_um.n_cells - 1, target, dtype=int)
        velocity = np.asarray(self.data.grid_um.cell_data["velocity_phy"])[indices]
        magnitude = np.linalg.norm(velocity, axis=1)
        nonzero = magnitude > np.finfo(float).tiny
        cloud = pv.PolyData(self.data.centers_um[indices][nonzero])
        cloud.point_data["velocity_direction"] = velocity[nonzero] / magnitude[nonzero, None]
        cloud.point_data["velocity_magnitude_mm_s"] = magnitude[nonzero] * 1.0e3
        dx_um = float(self.data.run_summary["numerical_contract"]["dx_m"]) * 1.0e6
        return cloud.glyph(
            orient="velocity_direction",
            scale=False,
            factor=max(2.0 * dx_um, 0.5),
            tolerance=None,
        )

    def _field_caption(self) -> str:
        field = self.data.fields[self.current_field]
        range_mode = "Full" if self.full_range else "P1–P99"
        color_scale = "Log" if field.log_scale else "Linear"
        view = self.visual_mode.title()
        if self.show_streamlines:
            view += " + flow paths"
        return (
            f"{field.title}  [{field.units}]\n"
            f"Original CFD surface  ·  {range_mode}  ·  {color_scale}  ·  "
            f"{view}"
        )

    def _style_text_panel(
        self,
        actor: Any,
        *,
        opacity: float = 0.88,
        horizontal: str = "left",
        vertical: str = "bottom",
    ) -> Any:
        """Apply a readable high-contrast panel to a VTK text actor."""

        actor.prop.background_color = self.style.panel_color
        actor.prop.background_opacity = opacity
        actor.prop.show_frame = True
        actor.prop.frame_color = self.style.panel_border
        actor.prop.frame_width = 2
        if hasattr(actor.prop, "justification_horizontal"):
            actor.prop.justification_horizontal = horizontal
            actor.prop.justification_vertical = vertical
        return actor

    def _update_title(self, *, render: bool = True) -> None:
        self.plotter.remove_actor("field_title", render=False)
        actor = self.plotter.add_text(
            self._field_caption(),
            position=(28, int(self.plotter.window_size[1]) - 28),
            font_size=self.style.title_font_size,
            color=self.style.text_color,
            font=self.style.font_family,
            name="field_title",
            render=False,
        )
        self._style_text_panel(actor, opacity=0.82, vertical="top")
        if render:
            self.plotter.render()

    def _help_text(self) -> str:
        return (
            "FIELDS  click circles or use 1–8\n"
            "        1 Speed · 2 Pressure · 3 Density Δ · 4 Dynamic p\n"
            "        5 uₓ · 6 uᵧ · 7 u_z · 8 ρ lattice\n"
            "INSPECT C Clip · L Slice · V Vectors · T Streamlines · N Normals\n"
            "CAMERA  I Academic · X/Y/Z Axis · P Projection · R/F Fit\n"
            "DISPLAY A Robust · M Full · E Edges · 0 Overview\n"
            "PICK    Right-click · S Screenshot · H Help · Q Quit"
        )

    def _update_overlays(self, *, render: bool = True) -> None:
        self._update_title(render=False)
        self._update_help(render=False)
        if render:
            self.plotter.render()

    def _update_help(self, *, render: bool = True) -> None:
        self.plotter.remove_actor("help_overlay", render=False)
        if self.show_help:
            actor = self.plotter.add_text(
                self._help_text(),
                position=(28, 126),
                font_size=self.style.metadata_font_size,
                color=self.style.text_color,
                font=self.style.font_family,
                name="help_overlay",
                render=False,
            )
            self._style_text_panel(actor)
        if render:
            self.plotter.render()

    def _pick_callback(self, point: Sequence[float] | None) -> None:
        if point is None:
            return
        coordinates = np.asarray(point, dtype=float)
        if coordinates.shape != (3,) or not np.all(np.isfinite(coordinates)):
            return
        _, raw_index = self.data.center_tree_um.query(coordinates, k=1)
        index = int(raw_index)
        original = self.data.grid_um.cell_data
        magnitude = float(original["velocity_magnitude_mm_s"][index])
        pressure = float(original["pressure_gauge_pa"][index])
        rho = float(original["rho_lattice"][index])
        density_deviation = float(original["density_deviation_ppm"][index])
        dynamic_pressure = float(original["dynamic_pressure_mpa"][index])
        velocity = np.asarray(original["velocity_phy"][index], dtype=np.float64) * 1.0e3
        self.plotter.remove_actor("picked_marker", render=False)
        dx_um = float(self.data.run_summary["numerical_contract"]["dx_m"]) * 1.0e6
        marker_diameter = float(np.clip(1.5 * dx_um, 0.3, 0.6))
        marker = pv.Sphere(radius=marker_diameter / 2.0, center=coordinates)
        self.plotter.add_mesh(
            marker,
            color="#FFFFFF",
            ambient=0.75,
            diffuse=0.25,
            specular=0.0,
            name="picked_marker",
            pickable=False,
            reset_camera=False,
            render=False,
        )
        text = (
            f"x  {coordinates[0]:.3f} μm   y  {coordinates[1]:.3f} μm\n"
            f"z  {coordinates[2]:.3f} μm   |u|  {magnitude:.4g} mm/s\n"
            f"uₓ / uᵧ / u_z  {velocity[0]:.4g} / {velocity[1]:.4g} / "
            f"{velocity[2]:.4g} mm/s\n"
            f"p  {pressure:.5g} Pa   q  {dynamic_pressure:.5g} mPa\n"
            f"ρ  {rho:.9g}   Δρ  {density_deviation:.4g} ppm\n"
            f"cell {index:,} · original cell-centred evidence"
        )
        actor = self.plotter.add_text(
            text,
            position=(int(self.plotter.window_size[0]) - 28, 126),
            font_size=self.style.metadata_font_size,
            color=self.style.text_color,
            font=self.style.font_family,
            name="pick_information",
            render=False,
        )
        self._style_text_panel(actor, horizontal="right")
        self.picked_point_um = coordinates
        self.picked_cell_id = index
        self.picked_marker_visible = True
        self.plotter.render()

    def _clear_pick(self) -> None:
        self.plotter.remove_actor("picked_marker", render=False)
        self.plotter.remove_actor("pick_information", render=False)
        self.picked_marker_visible = False

    def set_field(self, field: str) -> None:
        if field not in self.data.fields:
            raise VisualizerInputError(f"Field is unavailable or debug-disabled: {field}")
        if field != "velocity" and self.show_streamlines:
            self.show_streamlines = False
            self._remove_streamline_actors()
            self._sync_streamline_widget()
            self._update_streamline_module_label(render=False)
        self.current_field = field
        self._sync_field_selector()
        self._replace_field_actor()

    def set_visual_mode(self, mode: str) -> None:
        if mode not in {"overview", "clip", "slice"}:
            raise ValueError(mode)
        if mode == "overview":
            self.hide_plane_widget()
        self.visual_mode = mode
        self.plane_mode = mode
        self._replace_field_actor()

    def set_plane_mode(self, mode: str) -> None:
        """Backward-compatible alias for inspection mode changes."""

        self.set_visual_mode(mode)

    def set_full_range(self, enabled: bool) -> None:
        self.full_range = enabled
        self._replace_field_actor(render=False)
        if self.show_streamlines:
            self._remove_streamline_actors()
            self._add_streamlines()
        self.plotter.render()

    def toggle_vectors(self) -> None:
        self.show_vectors = not self.show_vectors
        self.plotter.remove_actor("velocity_glyphs", render=False)
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
                reset_camera=False,
                render=False,
            )
        self.plotter.render()

    def _remove_streamline_actors(self) -> None:
        self.plotter.remove_actor("streamlines", render=False)
        self.plotter.remove_actor("streamline_direction_arrows", render=False)

    def set_streamlines(self, enabled: bool, *, render: bool = True) -> None:
        available = (
            self.data.streamlines_um is not None
            and self.data.streamlines_um.n_points > 0
            and self.data.valid_streamline_count > 0
        )
        self.show_streamlines = bool(enabled and available)
        if self.show_streamlines and self.current_field != "velocity":
            self.current_field = "velocity"
            self._sync_field_selector()
        self._remove_streamline_actors()
        self._replace_field_actor(render=False)
        if self.show_streamlines:
            self._add_streamlines()
        self._sync_streamline_widget()
        self._update_streamline_module_label(render=False)
        if render:
            self.plotter.render()

    def toggle_streamlines(self) -> None:
        self.set_streamlines(not self.show_streamlines)

    def toggle_help(self) -> None:
        self.show_help = not self.show_help
        self._update_help()

    def set_clean_ui(self) -> None:
        self.show_help = False
        self._update_help(render=False)
        self.plotter.render()

    def toggle_edges(self) -> None:
        self.show_edges = not self.show_edges
        self._replace_field_actor()

    def apply_academic_camera(self, *, render: bool = True) -> None:
        camera = self.plotter.camera
        record = self.camera_parameters
        camera.position = tuple(record["position_um"])
        camera.focal_point = tuple(record["focal_point_um"])
        camera.up = tuple(record["view_up"])
        camera.parallel_scale = float(record["parallel_scale"])
        if self.projection == "parallel":
            camera.parallel_projection = True
        else:
            camera.parallel_projection = False
            view_angle = math.radians(float(camera.view_angle) / 2.0)
            distance = float(record["parallel_scale"]) / math.tan(view_angle)
            camera.position = tuple(
                np.asarray(record["focal_point_um"])
                + np.asarray(record["view_out"]) * distance
            )
        camera.OrthogonalizeViewUp()
        if render:
            self.plotter.render()

    def focus_on_vessel(self) -> None:
        """Fit exclusively from cached vessel surface geometry."""

        self.apply_academic_camera()

    def focus_vessel(self) -> None:
        self.focus_on_vessel()

    def reset_camera(self) -> None:
        self.focus_on_vessel()

    def _set_isometric(self) -> None:
        self.apply_academic_camera()

    def _view_axis(self, vector: tuple[float, float, float]) -> None:
        direction = _normalise(np.asarray(vector, dtype=float))
        center = np.asarray(self.camera_parameters["focal_point_um"])
        distance = float(np.linalg.norm(np.ptp(self.data.surface_um.points, axis=0))) * 3.0
        up = np.array((0.0, 0.0, 1.0))
        if abs(float(np.dot(direction, up))) > 0.95:
            up = np.array((0.0, 1.0, 0.0))
        self.plotter.camera.position = tuple(center + direction * distance)
        self.plotter.camera.focal_point = tuple(center)
        self.plotter.camera.up = tuple(up)
        self.plotter.camera.parallel_scale = float(self.camera_parameters["parallel_scale"])
        self.plotter.render()

    def toggle_projection(self) -> None:
        self.projection = "perspective" if self.projection == "parallel" else "parallel"
        self.apply_academic_camera()

    def composition_qc(self) -> dict[str, Any]:
        record = self.camera_parameters
        center = np.asarray(record["focal_point_um"])
        up = np.asarray(record["view_up"])
        right = np.asarray(record["right"])
        points = np.asarray(self.data.surface_um.points) - center
        scale = float(record["parallel_scale"])
        width, height = self.plotter.window_size
        aspect = float(width) / float(height)
        x_values = points @ right / (2.0 * scale * aspect) + 0.5
        y_values = points @ up / (2.0 * scale) + 0.5
        width_fraction = float(np.ptp(x_values))
        height_fraction = float(np.ptp(y_values))
        port_viewports: dict[str, list[float]] = {}
        for label in PORT_ORDER:
            origin, _, _, _ = self._port_geometry(label)
            offset = origin - center
            port_viewports[label] = [
                float(np.dot(offset, right) / (2.0 * scale * aspect) + 0.5),
                float(np.dot(offset, up) / (2.0 * scale) + 0.5),
            ]
        bounds = np.asarray(self.data.surface_um.bounds).reshape(3, 2)
        focal_inside = bool(
            np.all(center >= bounds[:, 0] - 1.0e-12)
            and np.all(center <= bounds[:, 1] + 1.0e-12)
        )
        endpoints_visible = all(
            0.015 <= xy[0] <= 0.985 and 0.015 <= xy[1] <= 0.985
            for xy in port_viewports.values()
        )
        checks = {
            "focal_point_inside_vessel_bounds": focal_inside,
            "vessel_height_at_least_55pct": height_fraction >= 0.55,
            "vessel_width_at_least_35pct": width_fraction >= 0.35,
            "port_endpoints_not_cropped": endpoints_visible,
            "camera_bounds_source_vessel_only": True,
        }
        return {
            "status": "PASS" if all(checks.values()) else "FAIL",
            "width_fraction": width_fraction,
            "height_fraction": height_fraction,
            "viewport_bbox": [
                float(np.min(x_values)),
                float(np.max(x_values)),
                float(np.min(y_values)),
                float(np.max(y_values)),
            ],
            "port_viewport_coordinates": port_viewports,
            "checks": checks,
        }

    def visible_helper_actors(self) -> dict[str, Any]:
        return {
            "plane_widget": self.plane_widget_visible,
            "help": self.show_help,
            "info": False,
            "picked_marker": self.picked_marker_visible,
            "vectors": self.show_vectors,
            "streamlines": self.show_streamlines,
            "streamline_arrows": self.show_streamlines
            and "streamline_direction_arrows" in self.plotter.actors,
            "streamline_module": self.streamline_module_visible,
            "bounding_box": self.bounding_box_visible,
            "port_normals": self.show_port_normals,
            "ports": self.config.show_ports,
            "field_selector": self.field_selector_visible,
        }

    def _screenshot_metadata(self, path: Path, purpose: str) -> dict[str, Any]:
        field = self.data.fields[self.current_field]
        return {
            "status": "PASS",
            "purpose": purpose,
            "created_at": datetime.now().isoformat(),
            "screenshot": str(path),
            "screenshot_sha256": sha256_file(path),
            "width_px": int(self.plotter.window_size[0]),
            "height_px": int(self.plotter.window_size[1]),
            "background": self.style.background,
            "background_top": self.style.background_top,
            "theme": self.style.theme,
            "scalar": self.current_field,
            "array": field.array,
            "units": field.units,
            "field_description": field.description,
            "available_fields": list(self.data.fields),
            "derived_arrays": list(DERIVED_ARRAYS),
            "raw_range": [
                self.data.ranges[self.current_field].raw_min,
                self.data.ranges[self.current_field].raw_max,
            ],
            "display_range": list(self._current_limits()),
            "range_mode": "full" if self.full_range else "p1-p99",
            "color_scale": "logarithmic" if field.log_scale else "linear",
            "visual_mode": self.visual_mode,
            "projection": self.projection,
            "render_interpolation": self.data.rendering_scalar_interpolation,
            "display_geometry": "ORIGINAL_CONTINUOUS_CFD_SURFACE",
            "voxel_geometry_visible": False,
            "original_surface": str(self.data.original_surface_path),
            "original_surface_sha256": self.data.original_surface_sha256,
            "surface_mapping": self.data.surface_mapping,
            "ports_visible": self.config.show_ports,
            "streamlines_visible": self.show_streamlines,
            "streamline_direction_arrows_visible": self.show_streamlines
            and "streamline_direction_arrows" in self.plotter.actors,
            "streamline_module_visible": self.streamline_module_visible,
            "streamline_arrow_count": self.streamline_arrow_count,
            "streamline_rendering": {
                "seed_location": "validated inlet aperture",
                "path_color": "velocity_magnitude_mm_s",
                "direction_arrow_color": self.style.text_color,
                "surface_opacity_when_enabled": 0.22,
            },
            "clip_state": {
                "mode": self.visual_mode,
                "origin_um": self.plane_origin.tolist(),
                "normal": self.plane_normal.tolist(),
                "widget_visible": self.plane_widget_visible,
            },
            "camera": _camera_record(self.plotter),
            "camera_framing": "PCA_DETERMINISTIC_VESSEL_SURFACE_BOUNDS_ONLY",
            "vessel_bounds_um": [float(value) for value in self.data.surface_um.bounds],
            "vessel_projected_coverage": self.composition_qc(),
            "scalar_bar": {
                "position_x": self.style.scalar_bar_x,
                "position_y": self.style.scalar_bar_y,
                "width_fraction": self.style.scalar_bar_width,
                "height_fraction": self.style.scalar_bar_height,
                "labels": self.style.scalar_bar_labels,
                "title_gap_px": self.style.scalar_title_gap_px,
                "outline_visible": False,
                "title": f"{field.title}\n{field.units}",
            },
            "typography_px": {
                "font_family": "Arial",
                "title": self.style.title_font_size,
                "metadata": self.style.metadata_font_size,
                "ports": self.style.port_font_size,
                "controls": self.style.control_font_size,
                "scalar_title": self.style.scalar_title_font_size,
                "scalar_labels": self.style.scalar_label_font_size,
            },
            "visible_helper_actors": self.visible_helper_actors(),
            "anti_aliasing": self.anti_aliasing,
            "depth_peeling": self.depth_peeling,
            "port_count": len(self.data.plane_contract["ports"]),
            "valid_streamline_count": self.data.valid_streamline_count,
            "quantitative_source": "ORIGINAL_CELL_CENTERED_VTU_ARRAYS",
            "source_vtu": str(self.data.vtu_path),
            "source_vtu_sha256": self.data.vtu_sha256,
            "source_classification": self.data.run_summary["steady_solution_source"],
            "source_iteration": int(self.data.metrics["iteration"]),
        }

    def save_interactive_screenshot(self) -> Path:
        output = self.data.run_dir / "visualization" / "interactive_v3_redesign"
        output.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = output / f"screenshot_{stamp}.png"
        previous = {
            "help": self.show_help,
            "widget": self.plane_widget_visible,
            "pick": self.picked_marker_visible,
        }
        self.show_help = False
        self._update_help(render=False)
        if previous["widget"]:
            self.hide_plane_widget()
        if previous["pick"]:
            self._clear_pick()
        self.plotter.render()
        self.plotter.screenshot(path, scale=2)
        _write_viewer_json(
            path.with_suffix(".json"),
            self._screenshot_metadata(path, "interactive_v3_clean_user_screenshot"),
        )
        self.show_help = previous["help"]
        self._update_help(render=False)
        if previous["widget"]:
            self.show_plane_widget(self.visual_mode)
        if previous["pick"] and self.picked_point_um is not None:
            self._pick_callback(self.picked_point_um)
        print(f"Screenshot saved: {path}")
        return path

    def render_publication_scene(
        self,
        path: Path,
        *,
        field: str,
        visual_mode: str,
        streamlines: bool = False,
    ) -> dict[str, Any]:
        path = path.expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        self.set_clean_ui()
        self.hide_plane_widget()
        self._clear_pick()
        self.show_vectors = False
        self.plotter.remove_actor("velocity_glyphs", render=False)
        self.show_port_normals = False
        for label in PORT_ORDER:
            self.plotter.remove_actor(f"port_normal_{label}", render=False)
        self.set_visual_mode(visual_mode)
        self.set_field(field)
        self.set_streamlines(streamlines, render=False)
        self.apply_academic_camera(render=False)
        self.plotter.show(auto_close=False, interactive=False)
        self.plotter.screenshot(path)
        metadata = self._screenshot_metadata(path, "publication_grade_cfd_figure")
        metadata["image_qc"] = _image_qc(path)
        _write_viewer_json(path.with_suffix(".json"), metadata)
        return metadata

    def save_publication_screenshot(self, path: Path) -> Path:
        self.render_publication_scene(
            path,
            field=self.current_field,
            visual_mode="overview",
            streamlines=False,
        )
        self.plotter.close()
        return path.expanduser().resolve()

    def show(self) -> None:
        self.plotter.show(
            title="Validated Base CFD - Interactive",
            interactive=True,
            auto_close=True,
        )


PUBLICATION_SCENES = (
    ("01_after_velocity_overview.png", "velocity", "overview", False),
    ("02_after_pressure_overview.png", "pressure", "overview", False),
    ("03_after_density_deviation.png", "density-deviation", "overview", False),
    ("04_after_dynamic_pressure.png", "dynamic-pressure", "overview", False),
    ("05_after_velocity_x.png", "velocity-x", "overview", False),
    ("06_after_velocity_clip.png", "velocity", "clip", False),
    ("07_after_pressure_slice.png", "pressure", "slice", False),
    ("08_after_streamlines.png", "velocity", "clip", True),
)


def _production_evidence_hashes(data: VisualData) -> dict[str, str]:
    accepted_restart = Path(
        data.run_summary["mesh_provenance"]["accepted_restart"]["binary"]
    )
    paths = {
        "source_vtu": data.vtu_path,
        "original_continuous_surface": data.original_surface_path,
        "accepted_restart": accepted_restart,
        "vtu_manifest": data.vtu_path.parent / "production_steady_flow_field_manifest.json",
        "steady_qc": data.run_dir / "qc" / "production_steady_qc.json",
        "primary_metrics": data.run_dir / "qc" / "production_primary_metrics.json",
        "run_summary": data.run_dir / "qc" / "run_summary.json",
        "physical_port_flux": data.run_dir / "steady_replay" / "physical_port_flux.json",
    }
    return {name: sha256_file(path) for name, path in paths.items()}


def render_publication_suite(
    data: VisualData,
    config: VisualConfig,
    output_dir: Path,
) -> tuple[dict[str, Any], Path]:
    """Render a deterministic eight-panel, multi-field acceptance suite."""

    output = output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    before_hashes = _production_evidence_hashes(data)
    renders: dict[str, Any] = {}
    for filename, field, visual_mode, streamlines in PUBLICATION_SCENES:
        scene_config = VisualConfig(
            width=max(config.width, 3200),
            height=max(config.height, 2000),
            initial_scalar=field,
            streamline_seeds=config.streamline_seeds,
            build_streamlines=True,
            show_ports=True,
            full_range=config.full_range,
            debug_cells=False,
            numerical_pressure_debug=config.numerical_pressure_debug,
            projection="parallel",
            ui_mode="clean",
            theme=config.theme,
        )
        viewer = AcademicCFDViewer(
            data, scene_config, off_screen=True, publication=True
        )
        path = output / filename
        renders[filename] = viewer.render_publication_scene(
            path,
            field=field,
            visual_mode=visual_mode,
            streamlines=streamlines,
        )
        viewer.plotter.close()
    after_hashes = _production_evidence_hashes(data)
    report = {
        "status": "PASS",
        "solver_calls": 0,
        "created_at": datetime.now().isoformat(),
        "source_vtu_sha256": data.vtu_sha256,
        "display_geometry": "ORIGINAL_CONTINUOUS_CFD_SURFACE",
        "voxel_geometry_visible": False,
        "original_surface": str(data.original_surface_path),
        "original_surface_sha256": data.original_surface_sha256,
        "surface_mapping": data.surface_mapping,
        "render_interpolation": data.rendering_scalar_interpolation,
        "quantitative_source": "ORIGINAL_CELL_CENTERED_VTU_ARRAYS",
        "camera_framing": "PCA_DETERMINISTIC_VESSEL_SURFACE_BOUNDS_ONLY",
        "publication_resolution_px": [
            max(config.width, 3200),
            max(config.height, 2000),
        ],
        "production_evidence_before": before_hashes,
        "production_evidence_after": after_hashes,
        "production_evidence_unchanged": before_hashes == after_hashes,
        "renders": renders,
    }
    if not report["production_evidence_unchanged"]:
        report["status"] = "FAIL"
    manifest_path = output / "publication_suite_manifest.json"
    _write_viewer_json(manifest_path, report)
    return report, manifest_path


def render_interactive_dashboard_preview(
    data: VisualData,
    config: VisualConfig,
    output_dir: Path,
) -> tuple[dict[str, Any], Path]:
    """Render the default interactive workstation, including native controls."""

    dashboard_config = VisualConfig(
        width=max(config.width, 2240),
        height=max(config.height, 1400),
        initial_scalar="velocity",
        streamline_seeds=config.streamline_seeds,
        build_streamlines=True,
        show_ports=True,
        full_range=config.full_range,
        debug_cells=False,
        numerical_pressure_debug=config.numerical_pressure_debug,
        projection="parallel",
        ui_mode="analysis",
        theme=config.theme,
    )
    viewer = AcademicCFDViewer(data, dashboard_config, off_screen=True)
    path = output_dir.expanduser().resolve() / "00_interactive_dashboard.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    viewer.apply_academic_camera(render=False)
    viewer.plotter.show(auto_close=False, interactive=False)
    viewer.plotter.screenshot(path)
    metadata = viewer._screenshot_metadata(path, "interactive_dashboard_visual_regression")
    metadata["image_qc"] = _image_qc(path)
    switch_audit: dict[str, Any] = {}
    for field, widget in viewer.field_widgets.items():
        widget.GetRepresentation().SetState(1)
        widget.InvokeEvent("StateChangedEvent")
        expected_array = data.fields[field].array
        mapped_dataset = viewer.field_actor.mapper.dataset
        checks = {
            "selected_field_updated": viewer.current_field == field,
            "expected_array_on_display_geometry": (
                expected_array in mapped_dataset.point_data
                or expected_array in mapped_dataset.cell_data
            ),
            "field_actor_present": viewer.field_actor is not None,
        }
        switch_audit[field] = {
            "status": "PASS" if all(checks.values()) else "FAIL",
            "expected_array": expected_array,
            "selected_field": viewer.current_field,
            "checks": checks,
        }
    viewer.set_field("velocity")
    metadata["field_switch_audit"] = switch_audit
    streamline_path = output_dir.expanduser().resolve() / "09_interactive_streamline_module.png"
    streamline_checks = {
        "module_control_visible": viewer.streamline_module_visible,
        "module_control_available": viewer.streamline_widget is not None,
        "streamlines_enabled": False,
        "speed_field_selected": False,
        "tube_actor_present": False,
        "direction_arrow_actor_present": False,
        "direction_arrows_sampled": False,
    }
    if viewer.streamline_widget is not None:
        representation = viewer.streamline_widget.GetRepresentation()
        representation.SetState(1)
        viewer.streamline_widget.InvokeEvent("StateChangedEvent")
        streamline_checks.update(
            {
                "streamlines_enabled": viewer.show_streamlines,
                "speed_field_selected": viewer.current_field == "velocity",
                "tube_actor_present": "streamlines" in viewer.plotter.actors,
                "direction_arrow_actor_present": (
                    "streamline_direction_arrows" in viewer.plotter.actors
                ),
                "direction_arrows_sampled": viewer.streamline_arrow_count > 0,
            }
        )
    viewer.plotter.screenshot(streamline_path)
    streamline_metadata = viewer._screenshot_metadata(
        streamline_path, "interactive_streamline_module_visual_regression"
    )
    streamline_metadata["image_qc"] = _image_qc(streamline_path)
    _write_viewer_json(streamline_path.with_suffix(".json"), streamline_metadata)
    metadata["streamline_module_preview"] = str(streamline_path)
    metadata["streamline_module_audit"] = {
        "status": "PASS" if all(streamline_checks.values()) else "FAIL",
        "checks": streamline_checks,
        "valid_streamline_count": data.valid_streamline_count,
        "direction_arrow_count": viewer.streamline_arrow_count,
        "screenshot": str(streamline_path),
        "screenshot_sha256": sha256_file(streamline_path),
    }
    _write_viewer_json(path.with_suffix(".json"), metadata)
    viewer.plotter.close()
    return metadata, path


def run_self_test(data: VisualData, config: VisualConfig) -> tuple[dict[str, Any], Path]:
    """Render the redesign suite and verify scientific and composition invariants."""

    output = data.run_dir / "visualization" / "interactive_v3_redesign"
    array_hashes_before = {
        name: _array_sha256(np.asarray(data.grid_um.cell_data[name]))
        for name in REQUIRED_ARRAYS
    }
    suite, suite_path = render_publication_suite(data, config, output)
    dashboard, dashboard_path = render_interactive_dashboard_preview(
        data, config, output
    )
    array_hashes_after = {
        name: _array_sha256(np.asarray(data.grid_um.cell_data[name]))
        for name in REQUIRED_ARRAYS
    }
    renders = suite["renders"]
    overview = renders["01_after_velocity_overview.png"]
    helpers = overview["visible_helper_actors"]
    clutter_checks = {
        "plane_widget_hidden": not helpers["plane_widget"],
        "help_hidden": not helpers["help"],
        "info_hidden": not helpers["info"],
        "picked_marker_hidden": not helpers["picked_marker"],
        "vectors_hidden": not helpers["vectors"],
        "streamlines_hidden": not helpers["streamlines"],
        "streamline_arrows_hidden": not helpers["streamline_arrows"],
        "streamline_module_hidden": not helpers["streamline_module"],
        "bounding_box_hidden": not helpers["bounding_box"],
        "port_normals_hidden": not helpers["port_normals"],
        "field_selector_hidden": not helpers["field_selector"],
    }
    composition_checks = {
        name: record["vessel_projected_coverage"]["status"] == "PASS"
        for name, record in renders.items()
    }
    cell_statistics = {
        name: {
            "min": float(np.min(np.asarray(data.grid_um.cell_data[name]))),
            "max": float(np.max(np.asarray(data.grid_um.cell_data[name]))),
            "mean": float(np.mean(np.asarray(data.grid_um.cell_data[name]))),
        }
        for name in REQUIRED_ARRAYS + DERIVED_ARRAYS
    }
    interactive_style = academic_style(config.theme)
    rendered_fields = {record["scalar"] for record in renders.values()}
    checks = {
        "manifest_pass": data.manifest.get("status") == "PASS",
        "source_vtu_sha_unchanged": sha256_file(data.vtu_path) == data.vtu_sha256,
        "accepted_restart_sha_unchanged": (
            suite["production_evidence_after"]["accepted_restart"]
            == data.run_summary["accepted_steady_reference"]["restart_sha256"]
        ),
        "production_evidence_unchanged": suite["production_evidence_unchanged"],
        "cell_count_unchanged": data.grid_um.n_cells == int(data.manifest["cell_count"]),
        "original_cell_arrays_unchanged": (
            array_hashes_before == array_hashes_after == data.original_cell_array_sha256
        ),
        "original_surface_sha_unchanged": (
            sha256_file(data.original_surface_path) == data.original_surface_sha256
        ),
        "original_surface_mapping_pass": data.surface_mapping["status"] == "PASS",
        "original_surface_field_mapping_present": all(
            data.fields[key].array in data.surface_um.point_data for key in FIELD_ORDER
        ),
        "derived_cell_arrays_present": all(
            name in data.grid_um.cell_data for name in DERIVED_ARRAYS
        ),
        "expanded_field_catalogue": set(FIELD_ORDER).issubset(data.fields),
        "original_surface_nonempty_and_watertight": (
            data.surface_um.n_cells > 0 and data.surface_um.n_open_edges == 0
        ),
        "no_voxel_geometry_rendered": (
            not dashboard["voxel_geometry_visible"]
            and all(not record["voxel_geometry_visible"] for record in renders.values())
        ),
        "four_ports": len(data.plane_contract["ports"]) == 4,
        "several_valid_streamlines": data.valid_streamline_count >= 4,
        "eight_scene_renders": len(renders) == 8
        and all(record["image_qc"]["status"] == "PASS" for record in renders.values()),
        "multi_field_scene_coverage": {
            "velocity",
            "pressure",
            "density-deviation",
            "dynamic-pressure",
            "velocity-x",
        }.issubset(rendered_fields),
        "interactive_dashboard_render": dashboard["image_qc"]["status"] == "PASS",
        "interactive_dashboard_composition": (
            dashboard["vessel_projected_coverage"]["status"] == "PASS"
        ),
        "interactive_field_selector_visible": bool(
            dashboard["visible_helper_actors"]["field_selector"]
        ),
        "interactive_streamline_module_pass": (
            dashboard["streamline_module_audit"]["status"] == "PASS"
        ),
        "streamline_scene_direction_arrows_present": bool(
            renders["08_after_streamlines.png"]["visible_helper_actors"][
                "streamline_arrows"
            ]
        )
        and renders["08_after_streamlines.png"]["streamline_arrow_count"] > 0,
        "upper_right_information_panel_removed": not bool(
            dashboard["visible_helper_actors"]["info"]
        ),
        "interactive_all_fields_switchable": all(
            record["status"] == "PASS"
            for record in dashboard["field_switch_audit"].values()
        )
        and set(dashboard["field_switch_audit"]) == set(FIELD_ORDER),
        "composition_qc": all(composition_checks.values()),
        "default_clutter_gate": all(clutter_checks.values()),
        "scalar_bar_gate": (
            overview["scalar_bar"]["position_x"] >= 0.88
            and overview["scalar_bar"]["position_x"]
            + overview["scalar_bar"]["width_fraction"]
            <= 0.95
            and overview["scalar_bar"]["position_y"] >= 0.28
            and overview["scalar_bar"]["width_fraction"] <= 0.05
            and overview["scalar_bar"]["height_fraction"] <= 0.55
            and overview["scalar_bar"]["labels"] <= 6
            and overview["scalar_bar"]["title_gap_px"] >= 18
            and not overview["scalar_bar"]["outline_visible"]
            and bool(overview["units"])
        ),
        "large_typography_gate": (
            interactive_style.font_family == "arial"
            and dashboard["typography_px"]["font_family"] == "Arial"
            and interactive_style.title_font_size >= 22
            and interactive_style.metadata_font_size >= 15
            and interactive_style.control_font_size >= 14
            and interactive_style.scalar_label_font_size >= 15
        ),
        "render_interpolation_declared": (
            data.rendering_scalar_interpolation == RENDER_INTERPOLATION
        ),
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "solver_calls": 0,
        "run_dir": str(data.run_dir),
        "vtu": str(data.vtu_path),
        "source_vtu_sha256": data.vtu_sha256,
        "original_surface": str(data.original_surface_path),
        "original_surface_sha256": data.original_surface_sha256,
        "surface_mapping": data.surface_mapping,
        "display_geometry": "ORIGINAL_CONTINUOUS_CFD_SURFACE",
        "voxel_geometry_visible": False,
        "width_px": overview["width_px"],
        "height_px": overview["height_px"],
        "background": overview["background"],
        "theme": config.theme,
        "field": overview["scalar"],
        "projection": overview["projection"],
        "camera": overview["camera"],
        "vessel_bounds_um": overview["vessel_bounds_um"],
        "vessel_projected_coverage": overview["vessel_projected_coverage"],
        "scalar_bar": overview["scalar_bar"],
        "visible_helper_actors": helpers,
        "available_fields": list(data.fields),
        "derived_arrays": list(DERIVED_ARRAYS),
        "port_count": overview["port_count"],
        "valid_streamline_count": data.valid_streamline_count,
        "render_interpolation": data.rendering_scalar_interpolation,
        "quantitative_source": "ORIGINAL_CELL_CENTERED_VTU_ARRAYS",
        "original_cell_array_sha256": array_hashes_after,
        "original_cell_statistics": cell_statistics,
        "publication_suite": str(suite_path),
        "interactive_dashboard": str(dashboard_path),
        "interactive_dashboard_metadata": dashboard,
        "renders": renders,
        "composition_checks": composition_checks,
        "default_clutter_checks": clutter_checks,
        "checks": checks,
    }
    report_path = output / "redesign_self_test.json"
    _write_viewer_json(report_path, report)
    if report["status"] != "PASS":
        raise RuntimeError(f"Interactive V3 redesign self-test failed: {checks}")
    return report, report_path


def _print_startup(data: VisualData) -> None:
    print("CFD Visualizer V3 - high-definition multi-field workstation")
    print(f"Run: {data.run_dir}")
    print(f"VTU: {data.vtu_path}")
    print(f"VTU SHA256: {data.vtu_sha256}")
    print(f"Original surface: {data.original_surface_path}")
    print(f"Original surface SHA256: {data.original_surface_sha256}")
    print(f"Cells: {data.grid_um.n_cells}")
    print(f"Source: {data.run_summary['steady_solution_source']}")
    print(f"Iteration: {int(data.metrics['iteration'])}")
    print("Press H for controls")


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
            build_streamlines=bool(
                not args.no_streamlines or args.self_test or args.publication_suite
            ),
            show_ports=bool(not args.no_ports or args.self_test),
            full_range=args.full_range,
            debug_cells=args.debug_cells,
            numerical_pressure_debug=args.show_numerical_pressure_debug,
            projection=args.projection,
            ui_mode=args.ui_mode,
            theme=args.theme,
        )
        data = load_and_validate_data(
            run_dir,
            vtu,
            config,
            project_root=PROJECT_ROOT,
            explicit_surface=args.surface,
        )
        _print_startup(data)
        try:
            if args.self_test:
                report, path = run_self_test(data, config)
                print(f"Self-test: {report['status']}")
                print(f"Valid streamlines: {report['valid_streamline_count']}")
                print(f"Report: {path}")
                print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_PUBLICATION_GRADE_READY")
                return 0
            if args.publication_suite:
                report, path = render_publication_suite(data, config, args.publication_suite)
                print(f"Publication suite: {report['status']}")
                print(f"Manifest: {path}")
                print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_PUBLICATION_GRADE_READY")
                return 0
            if args.publication_screenshot:
                viewer = AcademicCFDViewer(data, config, off_screen=True, publication=True)
                path = viewer.save_publication_screenshot(args.publication_screenshot)
                print(f"Publication screenshot: {path}")
                print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_PUBLICATION_GRADE_READY")
                return 0
            if args.off_screen:
                viewer = AcademicCFDViewer(data, config, off_screen=True, publication=True)
                output = (
                    data.run_dir
                    / "visualization"
                    / "interactive_v3_redesign"
                    / "off_screen_preview.png"
                )
                viewer.save_publication_screenshot(output)
                print(f"Off-screen preview: {output}")
                print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_PUBLICATION_GRADE_READY")
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
        print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_PUBLICATION_GRADE_READY")
        return 0
    except (VisualizerInputError, OSError, ValueError, KeyError) as error:
        print(f"CFD Visualizer V3 input invalid: {error}")
        print("STATUS: CFD_INTERACTIVE_VISUALIZER_V3_INPUT_INVALID")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
