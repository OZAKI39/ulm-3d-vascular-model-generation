"""Diagnostic outputs for 2D lumen flow failures and accepted fields."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import yaml
from scipy import ndimage
from ulm_vascular_model_generator.utils.core.models import Vessel

from ..core.types import FlowField, GridDomain, RasterizedVessels
from ..geometry.lumen_width import NarrowLumenResult, classify_narrow_lumen_cells
from ..runtime.console_output import print_key_values, print_section


def write_flow_diagnostics(
    output_dir: Path,
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    vessels: list[Vessel] | tuple[Vessel, ...] = (),
) -> dict[str, Any]:
    """Save divergence overlays and scalar diagnostics for a computed flow field."""

    diagnostic_dir = Path(output_dir) / "flow_diagnostics"
    diagnostic_dir.mkdir(parents=True, exist_ok=True)

    masks = _diagnostic_masks(raster, flow)
    narrow_result = classify_narrow_lumen_cells(
        raster.lumen_mask,
        float(domain.spacing_um),
        raster.radius_um,
        min_diameter_px=8.0,
        junction_mask=raster.junction_core_mask,
        open_boundary_mask=masks["open_boundary"],
    )
    divergence = np.asarray(flow.divergence_s_inv, dtype=float)
    summary: dict[str, Any] = {
        "stage_divergence": _stage_divergence_from_metadata(flow),
        "momentum_balance": _momentum_balance_from_metadata(flow),
        "divergence_stats": {
            "core": _array_stats(divergence[masks["core"]]),
            "boundary": _array_stats(divergence[masks["boundary"]]),
            "lumen": _array_stats(divergence[raster.lumen_mask]),
        },
        "mask_quality": _mask_quality(raster, domain, narrow_result, masks),
        "final_velocity_stage": str(
            flow.solver_metadata.get("saved_velocity_stage", "unknown")
        ),
        "post_projection_velocity_modified": bool(
            flow.solver_metadata.get("post_projection_velocity_modified", True)
        ),
        "solver_metadata": _scalar_metadata(flow.solver_metadata),
    }

    outlet_rows = _write_outlet_flux_csv(
        diagnostic_dir / "outlet_flux_errors.csv", flow
    )
    summary["outlet_flux"] = {
        "n_outlets": int(max(0, len(outlet_rows) - 1)),
        "max_relative_error": float(
            max((abs(row["relative_error"]) for row in outlet_rows[1:]), default=0.0)
        ),
    }
    narrow_segment_rows = _write_narrow_segments_csv(
        diagnostic_dir / "narrow_segments.csv",
        raster,
        domain,
        narrow_result,
        vessels,
    )
    summary["narrow_segments"] = _narrow_segments_summary(narrow_segment_rows)

    heatmap_path = diagnostic_dir / "divergence_heatmap_overlay.png"
    _save_divergence_heatmap(heatmap_path, domain, raster, flow, masks)
    summary["files"] = {
        "divergence_heatmap_overlay": str(heatmap_path),
        "outlet_flux_errors": str(diagnostic_dir / "outlet_flux_errors.csv"),
        "narrow_segments": str(diagnostic_dir / "narrow_segments.csv"),
    }

    yaml_path = diagnostic_dir / "flow_diagnostics.yaml"
    yaml_path.write_text(yaml.safe_dump(summary, sort_keys=False), encoding="utf-8")
    _print_diagnostics(summary)
    return summary


def _diagnostic_masks(
    raster: RasterizedVessels, flow: FlowField
) -> dict[str, np.ndarray]:
    lumen = np.asarray(raster.lumen_mask, dtype=bool)
    open_boundary = np.zeros(lumen.shape, dtype=bool)
    if flow.inlet_label is not None:
        open_boundary |= np.asarray(flow.inlet_label > 0, dtype=bool)
    if flow.outlet_label is not None:
        open_boundary |= np.asarray(flow.outlet_label > 0, dtype=bool)
    boundary = (
        ndimage.binary_dilation(
            open_boundary, structure=np.ones((3, 3), dtype=bool), iterations=3
        )
        & lumen
    )
    wall = np.asarray(raster.wall_mask, dtype=bool) & ~open_boundary
    core = lumen & ~boundary
    if not np.any(core):
        core = lumen.copy()
    return {
        "core": core,
        "boundary": boundary,
        "wall": wall,
        "open_boundary": open_boundary,
    }


def _array_stats(values: np.ndarray) -> dict[str, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {
            "count": 0,
            "min": 0.0,
            "mean": 0.0,
            "rms": 0.0,
            "max_abs": 0.0,
            "p95_abs": 0.0,
            "p99_abs": 0.0,
        }
    return {
        "count": int(finite.size),
        "min": float(np.min(finite)),
        "mean": float(np.mean(finite)),
        "rms": float(np.sqrt(np.mean(finite * finite))),
        "max_abs": float(np.max(np.abs(finite))),
        "p95_abs": float(np.percentile(np.abs(finite), 95.0)),
        "p99_abs": float(np.percentile(np.abs(finite), 99.0)),
    }


def _stage_divergence_from_metadata(flow: FlowField) -> dict[str, float]:
    keys = (
        "div_initial",
        "div_after_boundary_apply",
        "div_after_projection",
        "div_final_saved",
    )
    return {key: float(flow.solver_metadata.get(key, float("nan"))) for key in keys}


def _momentum_balance_from_metadata(flow: FlowField) -> dict[str, Any]:
    """Return the DOLFINx linear-system convergence record."""
    metadata = flow.solver_metadata
    acceptance_schema = str(metadata.get("physical_acceptance_schema", "missing"))
    return {
        "acceptance_schema": acceptance_schema,
        "definition": "boundary_fitted_fem_linear_system_convergence",
        "converged": bool(metadata.get("physical_converged", False)),
        "normalized_residual": float(
            metadata.get("linear_solver_relative_residual", float("nan"))
        ),
        "tolerance": float(
            metadata.get(
                "linear_solver_residual_tolerance",
                metadata.get("dolfinx_ksp_rtol", float("nan")),
            )
        ),
        "residual_rms_um_s2": float("nan"),
        "residual_max_um_s2": float("nan"),
        "viscous_acceleration_rms_um_s2": float("nan"),
        "pressure_acceleration_rms_um_s2": float("nan"),
        "boundary_acceleration_rms_um_s2": float("nan"),
        "cell_count": int(metadata.get("mesh_cell_count", 0)),
        "linear_solver_backend": str(metadata.get("linear_solver_backend", "petsc")),
        "ksp_converged_reason": int(
            metadata.get(
                "linear_solver_converged_reason",
                metadata.get("dolfinx_ksp_converged_reason", 0),
            )
        ),
        "ksp_iterations": int(
            metadata.get(
                "linear_solver_iterations",
                metadata.get("dolfinx_ksp_iterations", 0),
            )
        ),
    }


def _mask_quality(
    raster: RasterizedVessels,
    domain: GridDomain,
    narrow_result: NarrowLumenResult,
    diagnostic_masks: dict[str, np.ndarray],
) -> dict[str, Any]:
    lumen = np.asarray(raster.lumen_mask, dtype=bool)
    labels, n_components = ndimage.label(
        lumen, structure=np.asarray([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
    )
    counts = (
        np.bincount(labels.ravel())[1:]
        if n_components > 0
        else np.asarray([], dtype=int)
    )
    filled = ndimage.binary_fill_holes(lumen)
    holes = filled & ~lumen
    narrow = np.asarray(narrow_result.narrow_mask, dtype=bool)
    junction = (
        np.asarray(raster.junction_core_mask, dtype=bool)
        if raster.junction_core_mask is not None
        else np.zeros(lumen.shape, dtype=bool)
    )
    open_boundary = np.asarray(diagnostic_masks["open_boundary"], dtype=bool)
    finite_effective = np.asarray(
        narrow_result.effective_diameter_px[lumen], dtype=float
    )
    finite_effective = finite_effective[np.isfinite(finite_effective)]
    finite_graph = np.asarray(narrow_result.graph_diameter_px[lumen], dtype=float)
    finite_graph = finite_graph[np.isfinite(finite_graph)]
    finite_mask = np.asarray(narrow_result.mask_diameter_px[lumen], dtype=float)
    finite_mask = finite_mask[np.isfinite(finite_mask)]
    narrow_fraction = float(np.count_nonzero(narrow) / max(np.count_nonzero(lumen), 1))
    return {
        "n_components": int(n_components),
        "component_sizes_desc": sorted((int(x) for x in counts), reverse=True)[:10],
        "hole_cells": int(np.count_nonzero(holes)),
        "narrow_cells": int(np.count_nonzero(narrow)),
        "narrow_fraction": narrow_fraction,
        "narrow_lumen_fraction": narrow_fraction,
        "narrow_lumen_definition": "min(D_graph_px, D_mask_medial_axis_px) < 8",
        "narrow_lumen_diameter_threshold_px": float(narrow_result.min_diameter_px),
        "narrow_in_inlet_or_outlet": int(np.count_nonzero(narrow & open_boundary)),
        "narrow_in_junction_core": int(np.count_nonzero(narrow & junction)),
        "narrow_in_main_junction_core": int(np.count_nonzero(narrow & junction)),
        "junction_core_cells": int(np.count_nonzero(junction)),
        "min_resolved_diameter_px": float(np.min(finite_effective))
        if finite_effective.size
        else 0.0,
        "p1_resolved_diameter_px": float(np.percentile(finite_effective, 1.0))
        if finite_effective.size
        else 0.0,
        "p5_resolved_diameter_px": float(np.percentile(finite_effective, 5.0))
        if finite_effective.size
        else 0.0,
        "median_resolved_diameter_px": float(np.median(finite_effective))
        if finite_effective.size
        else 0.0,
        "min_graph_diameter_px": float(np.min(finite_graph))
        if finite_graph.size
        else 0.0,
        "p5_graph_diameter_px": float(np.percentile(finite_graph, 5.0))
        if finite_graph.size
        else 0.0,
        "min_mask_diameter_px": float(np.min(finite_mask)) if finite_mask.size else 0.0,
        "p5_mask_diameter_px": float(np.percentile(finite_mask, 5.0))
        if finite_mask.size
        else 0.0,
        "minimum_continuous_wall_distance_um": (
            float(np.nanmin(raster.distance_to_wall_um[lumen]))
            if np.any(lumen)
            else 0.0
        ),
        "grid_spacing_um": float(domain.spacing_um),
    }


def _write_outlet_flux_csv(path: Path, flow: FlowField) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    header = "label,target_um2_s,actual_um2_s,relative_error\n"
    if (
        flow.outlet_target_by_label_um2_s is None
        or flow.outlet_actual_by_label_um2_s is None
    ):
        path.write_text(header, encoding="utf-8")
        return rows

    targets = np.asarray(flow.outlet_target_by_label_um2_s, dtype=float)
    actuals = np.asarray(flow.outlet_actual_by_label_um2_s, dtype=float)
    lines = [header]
    rows.append(
        {"label": 0.0, "target_um2_s": 0.0, "actual_um2_s": 0.0, "relative_error": 0.0}
    )
    for label in range(1, min(targets.size, actuals.size)):
        actual = float(actuals[label])
        target = float(targets[label])
        rel = (actual - target) / max(abs(target), np.finfo(float).eps)
        rows.append(
            {
                "label": float(label),
                "target_um2_s": target,
                "actual_um2_s": actual,
                "relative_error": float(rel),
            }
        )
        lines.append(f"{label},{target:.12g},{actual:.12g},{rel:.12g}\n")
    path.write_text("".join(lines), encoding="utf-8")
    return rows


def _write_narrow_segments_csv(
    path: Path,
    raster: RasterizedVessels,
    domain: GridDomain,
    narrow_result: NarrowLumenResult,
    vessels: list[Vessel] | tuple[Vessel, ...],
) -> list[dict[str, float | int | bool]]:
    header = (
        "segment_id,radius_um,diameter_px_graph,diameter_px_mask_median,diameter_px_mask_min,"
        "diameter_px_effective_median,diameter_px_effective_min,lumen_cells,narrow_cells,"
        "length_um,flow,is_inlet_or_outlet,is_terminal\n"
    )
    lumen = np.asarray(raster.lumen_mask, dtype=bool)
    vessel_id = np.asarray(raster.vessel_id, dtype=int)
    vessel_by_id = {int(v.vid): v for v in vessels}
    rows: list[dict[str, float | int | bool]] = []
    lines = [header]
    for segment_id in sorted(
        int(v) for v in np.unique(vessel_id[lumen]) if int(v) >= 0
    ):
        mask = lumen & (vessel_id == segment_id)
        if not np.any(mask):
            continue
        vessel = vessel_by_id.get(segment_id)
        radius_um = (
            float(vessel.radius)
            if vessel is not None
            else _finite_median(raster.radius_um[mask])
        )
        graph_diameter = (
            2.0 * radius_um / max(float(domain.spacing_um), np.finfo(float).eps)
        )
        mask_values = _finite_values(narrow_result.mask_diameter_px[mask])
        effective_values = _finite_values(narrow_result.effective_diameter_px[mask])
        narrow_cells = int(np.count_nonzero(narrow_result.narrow_mask & mask))
        length_um = float(vessel.length()) if vessel is not None else 0.0
        flow = (
            float(vessel.flow_rate)
            if vessel is not None
            else _finite_median(raster.flow_rate_um3_s[mask])
        )
        is_terminal = bool(
            vessel is not None and len(getattr(vessel, "children", [])) == 0
        )
        is_inlet_or_outlet = bool(
            vessel is not None and (int(vessel.parent_id) < 0 or is_terminal)
        )
        row: dict[str, float | int | bool] = {
            "segment_id": int(segment_id),
            "radius_um": float(radius_um),
            "diameter_px_graph": float(graph_diameter),
            "diameter_px_mask_median": float(np.median(mask_values))
            if mask_values.size
            else 0.0,
            "diameter_px_mask_min": float(np.min(mask_values))
            if mask_values.size
            else 0.0,
            "diameter_px_effective_median": float(np.median(effective_values))
            if effective_values.size
            else 0.0,
            "diameter_px_effective_min": float(np.min(effective_values))
            if effective_values.size
            else 0.0,
            "lumen_cells": int(np.count_nonzero(mask)),
            "narrow_cells": int(narrow_cells),
            "length_um": float(length_um),
            "flow": float(flow),
            "is_inlet_or_outlet": bool(is_inlet_or_outlet),
            "is_terminal": bool(is_terminal),
        }
        rows.append(row)
        lines.append(
            f"{row['segment_id']},{row['radius_um']:.12g},{row['diameter_px_graph']:.12g},"
            f"{row['diameter_px_mask_median']:.12g},{row['diameter_px_mask_min']:.12g},"
            f"{row['diameter_px_effective_median']:.12g},{row['diameter_px_effective_min']:.12g},"
            f"{row['lumen_cells']},{row['narrow_cells']},{row['length_um']:.12g},{row['flow']:.12g},"
            f"{int(row['is_inlet_or_outlet'])},{int(row['is_terminal'])}\n"
        )
    path.write_text("".join(lines), encoding="utf-8")
    return rows


def _narrow_segments_summary(
    rows: list[dict[str, float | int | bool]],
) -> dict[str, float | int]:
    if not rows:
        return {
            "narrow_segments_count": 0,
            "narrow_segments_fraction": 0.0,
            "segments_count": 0,
            "p5_resolved_diameter_px": 0.0,
            "median_resolved_diameter_px": 0.0,
        }
    narrow_count = sum(1 for row in rows if int(row["narrow_cells"]) > 0)
    effective_min = np.asarray(
        [float(row["diameter_px_effective_min"]) for row in rows], dtype=float
    )
    return {
        "narrow_segments_count": int(narrow_count),
        "narrow_segments_fraction": float(narrow_count / max(len(rows), 1)),
        "segments_count": int(len(rows)),
        "p5_resolved_diameter_px": float(np.percentile(effective_min, 5.0)),
        "median_resolved_diameter_px": float(np.median(effective_min)),
    }


def _finite_values(values: np.ndarray) -> np.ndarray:
    finite = np.asarray(values, dtype=float)
    return finite[np.isfinite(finite)]


def _finite_median(values: np.ndarray) -> float:
    finite = _finite_values(values)
    return float(np.median(finite)) if finite.size else 0.0


def _save_divergence_heatmap(
    path: Path,
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    masks: dict[str, np.ndarray],
) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    divergence = np.ma.array(flow.divergence_s_inv, mask=~raster.lumen_mask)
    finite_abs = np.abs(divergence.compressed())
    vmax = float(np.percentile(finite_abs, 99.0)) if finite_abs.size else 1.0
    vmax = max(vmax, 1.0e-6)
    extent = (
        float(domain.x_coordinates_um[0]),
        float(domain.x_coordinates_um[-1]),
        float(domain.z_coordinates_um[0]),
        float(domain.z_coordinates_um[-1]),
    )

    fig, ax = plt.subplots(figsize=(11.0, 8.0), constrained_layout=True)
    ax.imshow(
        np.ma.array(raster.vessel_id, mask=~raster.lumen_mask).T,
        origin="lower",
        extent=extent,
        cmap="tab20",
        alpha=0.18,
        interpolation="nearest",
        aspect="equal",
    )
    image = ax.imshow(
        divergence.T,
        origin="lower",
        extent=extent,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        alpha=0.78,
        interpolation="nearest",
        aspect="equal",
    )
    ax.contour(
        domain.x_coordinates_um,
        domain.z_coordinates_um,
        raster.lumen_mask.T.astype(float),
        levels=[0.5],
        colors="white",
        linewidths=0.5,
    )
    ax.contour(
        domain.x_coordinates_um,
        domain.z_coordinates_um,
        masks["wall"].T.astype(float),
        levels=[0.5],
        colors="yellow",
        linewidths=0.4,
    )
    if flow.inlet_label is not None:
        ax.contour(
            domain.x_coordinates_um,
            domain.z_coordinates_um,
            (flow.inlet_label > 0).T.astype(float),
            levels=[0.5],
            colors="lime",
            linewidths=0.8,
        )
    if flow.outlet_label is not None:
        ax.contour(
            domain.x_coordinates_um,
            domain.z_coordinates_um,
            (flow.outlet_label > 0).T.astype(float),
            levels=[0.5],
            colors="cyan",
            linewidths=0.8,
        )
    ax.set_title("Divergence heatmap over lumen / inlet / outlet / wall / segment ID")
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Z (um)")
    fig.colorbar(image, ax=ax, label="divergence (1/s)", fraction=0.046, pad=0.03)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _scalar_metadata(metadata: dict[str, object]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in metadata.items():
        if isinstance(value, (str, bool, int, float)):
            result[key] = value
    return result


def _print_diagnostics(summary: dict[str, Any]) -> None:
    print_section("Flow diagnostics")
    print_key_values(
        [(f"File / {name}", path) for name, path in summary["files"].items()]
    )
    print_key_values(
        [
            (f"Divergence / {key}", f"{value:.6g}")
            for key, value in summary["stage_divergence"].items()
        ]
    )
    split_rows: list[tuple[str, object]] = []
    momentum = summary["momentum_balance"]
    split_rows.append(
        (
            "Boundary-fitted FEM solve",
            f"backend={momentum['linear_solver_backend']}, "
            f"reason={momentum['ksp_converged_reason']}, "
            f"iterations={momentum['ksp_iterations']}, "
            f"cells={momentum['cell_count']}, "
            f"converged={momentum['converged']}",
        )
    )
    for key in ("core", "boundary", "lumen"):
        stats = summary["divergence_stats"][key]
        split_rows.append(
            (
                f"Divergence split / {key}",
                f"rms={stats['rms']:.6g}, max_abs={stats['max_abs']:.6g}, "
                f"p99_abs={stats['p99_abs']:.6g}",
            )
        )
    mask = summary["mask_quality"]
    split_rows.extend(
        [
            (
                "Mask quality",
                f"components={mask['n_components']}, holes={mask['hole_cells']}, "
                f"narrow={mask['narrow_cells']} ({100.0 * mask['narrow_lumen_fraction']:.2f}%)",
            ),
            (
                "Narrow-cell locations",
                f"open={mask['narrow_in_inlet_or_outlet']}, "
                f"junction={mask['narrow_in_junction_core']}",
            ),
            (
                "Outlet max relative flux error",
                f"{summary['outlet_flux']['max_relative_error']:.6g}",
            ),
        ]
    )
    print_key_values(split_rows)
