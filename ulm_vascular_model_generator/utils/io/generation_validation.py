"""Write compact validation artifacts for a generated three-dimensional tree."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from ..hemodynamics.dimensional_models import haemodynamics_for_config


def generation_validation_summary(vessels, cfg):
    if not vessels:
        raise ValueError("Cannot validate an empty vessel tree.")
    roots = [v for v in vessels if v.parent_id < 0]
    terminals = [v for v in vessels if not v.children]
    coordinates = np.vstack(
        [np.asarray(v.x_p, dtype=float) for v in vessels]
        + [np.asarray(v.x_d, dtype=float) for v in vessels]
    )
    terminal_sum = float(np.sum([v.flow_rate for v in terminals]))
    model = haemodynamics_for_config(cfg)
    inlet_flow = model.inlet_value
    summary = {
        "geometry_mode": str(cfg.geometry_mode),
        "flow_quantity": model.flow_quantity,
        "flow_unit": model.flow_unit,
        "vessel_count": len(vessels),
        "terminal_count": len(terminals),
        "root_count": len(roots),
        "maximum_flow_conservation_residual": float(
            max(v.flow_conservation_residual for v in vessels)
        ),
        "maximum_murray_residual": float(
            max(v.murray_residual for v in vessels)
        ),
        "minimum_radius_um": float(min(v.radius for v in vessels)),
        "maximum_radius_um": float(max(v.radius for v in vessels)),
        "minimum_mean_velocity_mm_s": float(
            min(v.mean_velocity for v in vessels) / 1000.0
        ),
        "maximum_mean_velocity_mm_s": float(
            max(v.mean_velocity for v in vessels) / 1000.0
        ),
        "coordinate_span_x_um": float(np.ptp(coordinates[:, 0])),
        "coordinate_span_y_um": float(np.ptp(coordinates[:, 1])),
        "coordinate_span_z_um": float(np.ptp(coordinates[:, 2])),
    }
    if cfg.geometry_mode == "planar_2d":
        summary.update(
            {
                "configured_inlet_flux_um2_s": inlet_flow,
                "terminal_flux_sum_um2_s": terminal_sum,
                "relative_inlet_terminal_flux_error": abs(
                    terminal_sum - inlet_flow
                )
                / max(inlet_flow, 1.0e-30),
                "minimum_half_width_um": summary.pop("minimum_radius_um"),
                "maximum_half_width_um": summary.pop("maximum_radius_um"),
            }
        )
    else:
        summary.update(
            {
                "configured_inlet_flow_um3_s": inlet_flow,
                "terminal_flow_sum_um3_s": terminal_sum,
                "relative_inlet_terminal_flow_error": abs(
                    terminal_sum - inlet_flow
                )
                / max(inlet_flow, 1.0e-30),
            }
        )
    return summary


def write_generation_validation_outputs(vessels, cfg, output_directory):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    summary = generation_validation_summary(vessels, cfg)
    json_path = output_directory / "generation_validation.json"
    csv_path = output_directory / "generation_validation.csv"
    json_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["metric", "value"])
        writer.writerows(summary.items())
    return json_path, csv_path
