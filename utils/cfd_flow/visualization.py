"""Human-readable production CFD acceptance figures and offline HTML."""

from __future__ import annotations

import html
import math
import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from PIL import Image
from scipy.spatial import cKDTree

from .config import FlowConfig
from .io import FlowError, sha256_file, write_json
from .validated_contract import TARGET_VOLUME_FLOW_M3_S


VISUAL_FAILURE = "CFD_FLOW_PRODUCTION_VISUALIZATION_FAILED"
FIGURES = (
    "01_velocity_overview.png",
    "02_gauge_pressure_overview.png",
    "03_velocity_streamlines.png",
    "04_port_flow_balance.png",
    "05_outlet_flow_split.png",
    "06_pressure_drop_summary.png",
    "07_steady_qc_summary.png",
)


def _figure(config: FlowConfig, *, projection: str | None = None):
    size = (
        config.visualization.width_px / config.visualization.dpi,
        config.visualization.height_px / config.visualization.dpi,
    )
    figure = plt.figure(figsize=size, dpi=config.visualization.dpi, constrained_layout=True)
    axis = figure.add_subplot(111, projection=projection)
    return figure, axis


def _ports(plane_contract: dict[str, Any]) -> dict[str, np.ndarray]:
    return {
        label: np.asarray(record["planes"]["central"]["origin_m"], dtype=float)
        for label, record in plane_contract["ports"].items()
    }


def _label_ports(axis: Any, ports: dict[str, np.ndarray]) -> None:
    labels = {
        "inlet": "INLET", "outlet_01": "OUTLET 01",
        "outlet_02": "OUTLET 02", "outlet_03": "OUTLET 03",
    }
    for key, point in ports.items():
        value = point * 1.0e6
        axis.scatter(*value, c="crimson", s=35, marker="o", depthshade=False)
        axis.text(*value, f"  {labels[key]}", color="black", fontsize=9, weight="bold")


def _field_overview(
    grid: pv.UnstructuredGrid,
    values: np.ndarray,
    path: Path,
    config: FlowConfig,
    ports: dict[str, np.ndarray],
    *,
    title: str,
    color_label: str,
    cmap: str,
    note: str,
) -> tuple[float, float, float, float]:
    centers = np.asarray(grid.cell_centers().points)
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise FlowError(VISUAL_FAILURE, f"non-finite values for {path.name}")
    display_min, display_max = np.percentile(values, [1, 99])
    stride = max(1, len(centers) // 65_000)
    figure, axis = _figure(config, projection="3d")
    artist = axis.scatter(
        centers[::stride, 0] * 1.0e6,
        centers[::stride, 1] * 1.0e6,
        centers[::stride, 2] * 1.0e6,
        c=np.clip(values[::stride], display_min, display_max),
        s=1.0, cmap=cmap, vmin=display_min, vmax=display_max,
        alpha=0.85, linewidths=0,
    )
    _label_ports(axis, ports)
    figure.colorbar(artist, ax=axis, shrink=0.72, pad=0.08, label=color_label)
    axis.set_title(title, fontsize=16, weight="bold")
    axis.set_xlabel("x (µm)")
    axis.set_ylabel("y (µm)")
    axis.set_zlabel("z (µm)")
    axis.view_init(elev=24, azim=-58)
    axis.text2D(0.01, 0.01, note, transform=axis.transAxes, fontsize=9)
    figure.savefig(path)
    plt.close(figure)
    return float(np.min(values)), float(np.max(values)), float(display_min), float(display_max)


def _streamline_seeds(record: dict[str, Any]) -> np.ndarray:
    origin = np.asarray(record["origin_m"], dtype=float)
    basis_u = np.asarray(record["basis_u"], dtype=float)
    basis_v = np.asarray(record["basis_v"], dtype=float)
    ring = np.asarray(record["physical_aperture_contour_uv_m"], dtype=float)
    samples = ring[np.linspace(0, len(ring) - 1, 12, dtype=int)] * 0.40
    return np.vstack(
        [origin, origin + samples[:, :1] * basis_u + samples[:, 1:] * basis_v]
    )


def _nearest_cell_streamlines(
    points_m: np.ndarray,
    velocity_m_s: np.ndarray,
    seeds_m: np.ndarray,
    *,
    dx_m: float,
) -> list[np.ndarray]:
    """Integrate actual vectors with nearest-cell sampling and stop at the wall."""

    tree = cKDTree(points_m)
    lines: list[np.ndarray] = []
    step = 0.55 * dx_m
    for seed in seeds_m:
        current = np.asarray(seed, dtype=float).copy()
        line = [current.copy()]
        for _ in range(1800):
            distance, index = tree.query(current, k=1)
            if distance > 1.05 * dx_m:
                break
            vector = np.asarray(velocity_m_s[index], dtype=float)
            speed = float(np.linalg.norm(vector))
            if not math.isfinite(speed) or speed <= np.finfo(float).tiny:
                break
            following = current + step * vector / speed
            next_distance, _ = tree.query(following, k=1)
            if next_distance > 1.05 * dx_m:
                break
            current = following
            line.append(current.copy())
        if len(line) >= 8:
            lines.append(np.asarray(line))
    return lines


def _streamline_figure(
    grid: pv.UnstructuredGrid,
    points_m: np.ndarray,
    path: Path,
    config: FlowConfig,
    ports: dict[str, np.ndarray],
    inlet_record: dict[str, Any],
) -> tuple[int, float, float]:
    velocity = np.asarray(grid.cell_data["velocity_phy"])
    seeds = _streamline_seeds(inlet_record)
    lines = _nearest_cell_streamlines(
        points_m, velocity, seeds, dx_m=config.mesh.dx_m
    )
    if not lines:
        raise FlowError(VISUAL_FAILURE, "no valid inlet-seeded streamlines")
    all_speed = np.linalg.norm(velocity, axis=1) * 1.0e3
    display_min, display_max = np.percentile(all_speed, [1, 99])
    tree = cKDTree(points_m)
    figure, axis = _figure(config, projection="3d")
    cmap = plt.get_cmap("plasma")
    for line in lines:
        _, indices = tree.query(line, k=1)
        representative = float(np.mean(all_speed[indices]))
        fraction = np.clip(
            (representative - display_min) / max(display_max - display_min, 1.0e-30),
            0.0, 1.0,
        )
        scaled = line * 1.0e6
        axis.plot(scaled[:, 0], scaled[:, 1], scaled[:, 2], color=cmap(fraction), linewidth=1.4)
    _label_ports(axis, ports)
    axis.set_title("Inlet-seeded velocity streamlines", fontsize=16, weight="bold")
    axis.set_xlabel("x (µm)")
    axis.set_ylabel("y (µm)")
    axis.set_zlabel("z (µm)")
    axis.view_init(elev=24, azim=-58)
    axis.text2D(
        0.01, 0.01,
        f"{len(lines)}/{len(seeds)} valid nearest-cell streamlines; actual velocity field, no smoothing",
        transform=axis.transAxes, fontsize=9,
    )
    figure.savefig(path)
    plt.close(figure)
    return len(lines), float(np.min(all_speed)), float(np.max(all_speed))


def _bar_flow(metrics: dict[str, Any], path: Path, config: FlowConfig) -> tuple[float, float]:
    factor = 60.0 * 1.0e12  # m3/s -> nL/min
    names = ("Q target", "Q inlet", "Q1", "Q2", "Q3", "Qout total")
    values = np.asarray(
        [
            TARGET_VOLUME_FLOW_M3_S,
            metrics["Qin_m3_s"], metrics["Q1_m3_s"], metrics["Q2_m3_s"],
            metrics["Q3_m3_s"], metrics["Qout_m3_s"],
        ]
    ) * factor
    figure, axis = _figure(config)
    colors = ("#666666", "#1f77b4", "#4c78a8", "#f58518", "#54a24b", "#2a9d8f")
    bars = axis.bar(names, values, color=colors)
    axis.bar_label(bars, fmt="%.5g", padding=4)
    axis.set_ylim(0.0, float(values.max()) * 1.16)
    axis.set_ylabel("Volumetric flow (nL/min)")
    axis.set_title("Physical port-flow balance", fontsize=16, weight="bold")
    axis.text(
        0.02, 0.90,
        f"Physical interior-plane closure = {100.0 * metrics['physical_volume_closure']:.6f}%",
        transform=axis.transAxes, va="top", fontsize=11,
    )
    axis.grid(axis="y", alpha=0.25)
    figure.savefig(path)
    plt.close(figure)
    return float(values.min()), float(values.max())


def _flow_split(metrics: dict[str, Any], path: Path, config: FlowConfig) -> tuple[float, float]:
    fractions = np.asarray(
        [metrics["flow_fractions"][f"outlet_{index:02d}"] for index in range(1, 4)]
    )
    figure, axis = _figure(config)
    colors = ("#4c78a8", "#f58518", "#54a24b")
    labels = ("Outlet 01", "Outlet 02", "Outlet 03")
    wedges, _ = axis.pie(
        fractions, colors=colors, startangle=90,
        wedgeprops={"width": 0.42, "edgecolor": "white"},
    )
    axis.legend(
        wedges,
        [f"{label}: {100.0 * value:.6f}%" for label, value in zip(labels, fractions)],
        loc="center left", bbox_to_anchor=(0.78, 0.5), fontsize=11,
    )
    axis.text(0, 0, "OUTLET\nFLOW SPLIT", ha="center", va="center", weight="bold")
    axis.set_title("Accepted Base outlet flow distribution", fontsize=16, weight="bold")
    figure.savefig(path)
    plt.close(figure)
    return float(fractions.min()), float(fractions.max())


def _pressure_summary(metrics: dict[str, Any], path: Path, config: FlowConfig) -> tuple[float, float]:
    names = ("Inlet gauge", "Outlet 01", "Outlet 02", "Outlet 03")
    outlet_gauges = np.asarray(
        [
            metrics["inlet_gauge_pressure_pa"],
            metrics["inlet_gauge_pressure_pa"] - metrics["pressure_drops_pa"]["outlet_01"],
            metrics["inlet_gauge_pressure_pa"] - metrics["pressure_drops_pa"]["outlet_02"],
            metrics["inlet_gauge_pressure_pa"] - metrics["pressure_drops_pa"]["outlet_03"],
        ]
    )
    figure, axis = _figure(config)
    bars = axis.bar(names, outlet_gauges, color=("#264653", "#2a9d8f", "#e9c46a", "#e76f51"))
    axis.bar_label(bars, fmt="%.3f Pa", padding=4)
    axis.set_ylim(float(outlet_gauges.min()) * 2.0, float(outlet_gauges.max()) * 1.14)
    axis.set_ylabel("Gauge pressure (Pa)")
    axis.set_title("Gauge-pressure boundary summary", fontsize=16, weight="bold")
    drops = metrics["pressure_drops_pa"]
    axis.text(
        0.02, 0.965,
        "Pressure drops: " + ", ".join(
            f"ΔP0{index}={drops[f'outlet_{index:02d}']:.3f} Pa" for index in range(1, 4)
        ),
        transform=axis.transAxes, va="top", fontsize=11,
    )
    axis.text(
        0.02, 0.03,
        "The ~3.39 MPa reference is an LBM numerical offset; it is not physiological pressure.",
        transform=axis.transAxes, fontsize=10,
    )
    axis.grid(axis="y", alpha=0.25)
    figure.savefig(path)
    plt.close(figure)
    return float(outlet_gauges.min()), float(outlet_gauges.max())


def _qc_summary(
    metrics: dict[str, Any], steady: dict[str, Any], full_v2: dict[str, Any],
    path: Path, config: FlowConfig,
) -> tuple[float, float]:
    rows = (
        ("rho mean", metrics["rho_mean"], "0.9–1.1", 0.9 <= metrics["rho_mean"] <= 1.1),
        (
            "Qin target error",
            abs(metrics["Qin_m3_s"] - TARGET_VOLUME_FLOW_M3_S)
            / TARGET_VOLUME_FLOW_M3_S,
            "≤ 0.01",
            abs(metrics["Qin_m3_s"] - TARGET_VOLUME_FLOW_M3_S)
            / TARGET_VOLUME_FLOW_M3_S
            <= 0.01,
        ),
        ("physical closure", metrics["physical_volume_closure"], "≤ 0.01", metrics["physical_volume_closure"] <= 0.01),
        ("R_mass_short", steady["R_mass_short"], "≤ 0.01", steady["R_mass_short"] <= 0.01),
        ("R_mass_long", steady["R_mass_long"], "≤ 0.01", steady["R_mass_long"] <= 0.01),
        ("R_velocity", steady["R_velocity"], "≤ 0.01", steady["R_velocity"] <= 0.01),
        ("R_pressure", steady["R_pressure"], "≤ 0.005", steady["R_pressure"] <= 0.005),
        ("R_inlet", steady["R_inlet"], "≤ 0.01", steady["R_inlet"] <= 0.01),
        ("minimum PDF", metrics["minimum_pdf"], "> 0", metrics["minimum_pdf"] > 0),
        ("max lattice speed", metrics["maximum_lattice_speed"], "< 0.05", metrics["maximum_lattice_speed"] < 0.05),
        ("Full timestep V2", full_v2["residual"], "≤ 1e-8", full_v2["status"] == "PASS"),
    )
    figure, axis = _figure(config)
    axis.axis("off")
    table_rows = [
        [name, f"{value:.10g}", threshold, "PASS" if passed else "FAIL"]
        for name, value, threshold, passed in rows
    ]
    table = axis.table(
        cellText=table_rows,
        colLabels=("Metric", "Value", "Threshold", "Status"),
        cellLoc="left", colLoc="left", loc="center",
        colWidths=(0.34, 0.24, 0.20, 0.14),
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.65)
    for row_index, row in enumerate(rows, start=1):
        color = "#d8f3dc" if row[3] else "#ffd6d6"
        table[(row_index, 3)].set_facecolor(color)
    axis.set_title("Steady CFD numerical acceptance dashboard", fontsize=16, weight="bold", pad=18)
    figure.savefig(path)
    plt.close(figure)
    values = [float(row[1]) for row in rows]
    return min(values), max(values)


def _write_html(
    path: Path,
    config: FlowConfig,
    metrics: dict[str, Any],
    steady: dict[str, Any],
    full_v2: dict[str, Any],
    coarse_base: dict[str, Any],
    status: str,
) -> None:
    fractions = metrics["flow_fractions"]
    images = "\n".join(
        f'<figure><img src="{name}" alt="{html.escape(name)}"><figcaption>{html.escape(name)}</figcaption></figure>'
        for name in FIGURES
    )
    pressure_reference = config.physics.density_kg_m3 * config.physics.lattice_cs_squared * (
        config.mesh.dx_m / (config.mesh.dx_m**2 / (6 * config.physics.kinematic_viscosity_m2_s))
    ) ** 2
    content = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Production CFD Review</title>
<style>
body{{font-family:Arial,sans-serif;margin:0;background:#f4f7fa;color:#17212b}}main{{max-width:1280px;margin:auto;padding:28px}}
.hero{{background:#12344d;color:white;padding:26px;border-radius:12px}}.pass{{color:#6ee7a2;font-weight:bold}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px;margin:18px 0}}.card{{background:white;padding:15px;border-radius:9px;box-shadow:0 2px 8px #0002}}
figure{{background:white;padding:14px;margin:18px 0;border-radius:9px}}img{{width:100%;height:auto}}table{{border-collapse:collapse;width:100%;background:white}}th,td{{padding:9px;border:1px solid #ccd5df;text-align:left}}
.warning{{background:#fff3cd;border-left:5px solid #d99b00;padding:14px;margin:18px 0}}
</style></head><body><main>
<section class="hero"><h1>PRODUCTION CFD REVIEW</h1><h2>STATUS: <span class="pass">{html.escape(status)}</span></h2>
<p>Pressure colors are gauge pressure; the ~3.39 MPa reference is an LBM numerical offset.</p></section>
<section class="cards">
<div class="card"><b>Q target</b><br>{config.boundary.target_volume_flow_m3_s:.12g} m³/s</div>
<div class="card"><b>Q measured</b><br>{metrics['Qin_m3_s']:.12g} m³/s</div>
<div class="card"><b>Closure</b><br>{100*metrics['physical_volume_closure']:.8f}%</div>
<div class="card"><b>Inlet gauge pressure</b><br>{metrics['inlet_gauge_pressure_pa']:.6f} Pa</div>
<div class="card"><b>Outlet fractions</b><br>{100*fractions['outlet_01']:.3f}% / {100*fractions['outlet_02']:.3f}% / {100*fractions['outlet_03']:.3f}%</div>
<div class="card"><b>Tau / rho mean</b><br>1 / {metrics['rho_mean']:.12g}</div>
</section>
<h2>Contract and provenance</h2><table>
<tr><th>Method</th><td>{html.escape(config.method)}</td></tr>
<tr><th>Source solution</th><td>VALIDATED_RESEARCH_BASE_ACCEPTED_RESTART; fresh full steady production solve = false</td></tr>
<tr><th>Resolution</th><td>Base, dx = {config.mesh.dx_m:.12g} m (0.20 µm)</td></tr>
<tr><th>Scaling</th><td>dt = {config.mesh.dx_m**2/(6*config.physics.kinematic_viscosity_m2_s):.16g} s; tau = 1; omega = 1</td></tr>
<tr><th>Numerical pressure reference</th><td>{pressure_reference:.12f} Pa; not physiological absolute blood pressure</td></tr>
<tr><th>Physical flow</th><td>PHYSICAL_INTERIOR_CROSS_SECTION_VELOCITY_FLUX</td></tr>
<tr><th>Full V2</th><td>{full_v2['status']}, residual {full_v2['residual']:.12g}, gate {full_v2['gate']:.1e}</td></tr>
<tr><th>WSS</th><td>NOT YET FORMALLY VALIDATED — DEFERRED_TO_POST_GRID_PRODUCTION_VALIDATION</td></tr>
</table>
<div class="warning"><b>Resolution scope:</b> TWO-GRID RESOLUTION SENSITIVITY only. Maximum accepted Coarse-to-Base difference = {coarse_base['maximum_absolute_percent_difference']:.6f}%.<br>
Formal three-grid GCI was not completed because Fine steady computation was terminated under resource-budget constraints.</div>
<p>Fine mesh: PASS. Fine controller fix: PASS. Fine 5000-step safety: PASS. Fine steady: NOT_COMPLETED_RESOURCE_BUDGET. No grid-independent claim is made.</p>
<h2>Steady acceptance</h2><p>R_mass_short={steady['R_mass_short']:.12g}; R_mass_long={steady['R_mass_long']:.12g}; R_velocity={steady['R_velocity']:.12g}; R_pressure={steady['R_pressure']:.12g}; R_inlet={steady['R_inlet']:.12g}.</p>
<h2>Visual acceptance package</h2>{images}
</main></body></html>"""
    path.write_text(content, encoding="utf-8")


def _validate_png(path: Path, config: FlowConfig) -> dict[str, Any]:
    with Image.open(path) as image:
        width, height = image.size
        pixels = np.asarray(image.convert("RGB"), dtype=np.float64)
    checks = {
        "size_bytes": path.stat().st_size > 10_000,
        "width": width >= config.visualization.width_px,
        "height": height >= config.visualization.height_px,
        "finite": bool(np.all(np.isfinite(pixels))),
        "nonblank": float(np.std(pixels)) > 1.0,
        "not_all_white": float(np.mean(pixels)) < 254.5,
        "not_all_black": float(np.mean(pixels)) > 0.5,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "width_px": width, "height_px": height,
        "pixel_standard_deviation": float(np.std(pixels)),
        "size_bytes": path.stat().st_size, "checks": checks,
    }


def create_production_visuals(
    *,
    grid: pv.UnstructuredGrid,
    points_m: np.ndarray,
    metrics: dict[str, Any],
    steady_qc: dict[str, Any],
    full_v2: dict[str, Any],
    plane_contract: dict[str, Any],
    coarse_base: dict[str, Any],
    output: Path,
    config: FlowConfig,
    status: str,
) -> dict[str, Any]:
    """Write seven Base-only figures, an offline report, and hard-gate them."""

    output.mkdir(parents=True, exist_ok=True)
    ports = _ports(plane_contract)
    records: list[dict[str, Any]] = []

    speed = np.asarray(grid.cell_data["velocity_magnitude_mm_s"])
    raw_min, raw_max, display_min, display_max = _field_overview(
        grid, speed, output / FIGURES[0], config, ports,
        title="Accepted Base velocity overview", color_label="Velocity magnitude (mm/s)",
        cmap="viridis", note="Display clipped to p1/p99; raw min/max retained in visual_manifest.json",
    )
    records.append({"filename": FIGURES[0], "purpose": "3D physical-velocity overview", "source_data": "accepted Base decoded restart", "units": "mm/s", "raw_min": raw_min, "raw_max": raw_max, "display_min": display_min, "display_max": display_max})

    pressure = np.asarray(grid.cell_data["pressure_gauge_pa"])
    raw_min, raw_max, display_min, display_max = _field_overview(
        grid, pressure, output / FIGURES[1], config, ports,
        title="Accepted Base gauge-pressure overview", color_label="Gauge pressure (Pa)",
        cmap="coolwarm", note="p_gauge = p_solver - P_ref; P_ref is a numerical LBM offset, not physiological pressure",
    )
    records.append({"filename": FIGURES[1], "purpose": "3D gauge-pressure overview", "source_data": "accepted Base decoded restart", "units": "Pa gauge", "raw_min": raw_min, "raw_max": raw_max, "display_min": display_min, "display_max": display_max})

    valid_lines, raw_min, raw_max = _streamline_figure(
        grid, points_m, output / FIGURES[2], config, ports,
        plane_contract["ports"]["inlet"]["planes"]["central"],
    )
    records.append({"filename": FIGURES[2], "purpose": f"inlet-seeded streamlines ({valid_lines} valid)", "source_data": "accepted Base physical velocity; nearest-cell integration", "units": "mm/s", "raw_min": raw_min, "raw_max": raw_max, "display_min": float(np.percentile(speed, 1)), "display_max": float(np.percentile(speed, 99))})

    raw_min, raw_max = _bar_flow(metrics, output / FIGURES[3], config)
    records.append({"filename": FIGURES[3], "purpose": "physical port-flow balance", "source_data": "V3 continuous aperture physical flux", "units": "nL/min", "raw_min": raw_min, "raw_max": raw_max, "display_min": raw_min, "display_max": raw_max})
    raw_min, raw_max = _flow_split(metrics, output / FIGURES[4], config)
    records.append({"filename": FIGURES[4], "purpose": "outlet flow fractions", "source_data": "V3 continuous aperture physical flux", "units": "fraction", "raw_min": raw_min, "raw_max": raw_max, "display_min": raw_min, "display_max": raw_max})
    raw_min, raw_max = _pressure_summary(metrics, output / FIGURES[5], config)
    records.append({"filename": FIGURES[5], "purpose": "gauge pressure and pressure-drop summary", "source_data": "accepted controller pressure and configured outlet gauges", "units": "Pa gauge", "raw_min": raw_min, "raw_max": raw_max, "display_min": raw_min, "display_max": raw_max})
    raw_min, raw_max = _qc_summary(metrics, steady_qc, full_v2, output / FIGURES[6], config)
    records.append({"filename": FIGURES[6], "purpose": "steady numerical QC dashboard", "source_data": "production replay QC and accepted physical-time audit", "units": "mixed; row-labelled", "raw_min": raw_min, "raw_max": raw_max, "display_min": raw_min, "display_max": raw_max})

    html_path = output / "production_review.html"
    _write_html(html_path, config, metrics, steady_qc, full_v2, coarse_base, status)
    html_text = html_path.read_text(encoding="utf-8")
    referenced = re.findall(r'<img src="([^"]+)"', html_text)
    html_checks = {
        "exists": html_path.is_file(),
        "offline": "http://" not in html_text and "https://" not in html_text and "cdn" not in html_text.lower(),
        "all_images_referenced": set(referenced) == set(FIGURES),
        "all_references_exist": all((output / name).is_file() for name in referenced),
        "gci_disclaimer_exact": "Formal three-grid GCI was not completed because Fine steady computation was terminated under resource-budget constraints." in html_text,
        "gauge_pressure_explanation": "Pressure colors are gauge pressure; the ~3.39 MPa reference is an LBM numerical offset." in html_text,
    }
    validations: dict[str, Any] = {}
    for record in records:
        path = output / record["filename"]
        validation = _validate_png(path, config)
        record.update(
            {
                "sha256": sha256_file(path),
                "width_px": validation["width_px"],
                "height_px": validation["height_px"],
                "status": validation["status"],
            }
        )
        validations[record["filename"]] = validation
    passed = (
        all(value["status"] == "PASS" for value in validations.values())
        and all(html_checks.values())
    )
    result = {
        "status": "PASS" if passed else "FAIL",
        "source_solution": "VALIDATED_RESEARCH_BASE_ACCEPTED_RESTART",
        "fine_transient_used": False,
        "items": records,
        "html": {
            "filename": html_path.name, "sha256": sha256_file(html_path),
            "size_bytes": html_path.stat().st_size,
            "status": "PASS" if all(html_checks.values()) else "FAIL",
            "checks": html_checks,
        },
        "validation": validations,
    }
    write_json(output / "visual_manifest.json", result)
    if not passed:
        raise FlowError(VISUAL_FAILURE, "one or more visual/HTML hard gates failed")
    return result


def create_flow_figures(
    grid: pv.UnstructuredGrid,
    partition: Any,
    output: Path,
) -> list[Path]:
    """Retired compatibility entry point; production uses the acceptance package."""

    raise RuntimeError(
        "create_flow_figures is retired; use create_production_visuals with accepted Base provenance"
    )
