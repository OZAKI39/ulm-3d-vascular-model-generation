"""Finite field-only control selection, before running any mature-cluster branch."""

import numpy as np
from mevo_mb.fast_sim.models import Basis
from .common import BASE, ROOT, PROTOCOL, write
import json


def diagnose(out):
    basis = Basis(BASE / "inputs/basis.h5", 2.5e6)
    base = json.loads((ROOT / "configs/fast_pzt/scene_c.json").read_text())["events"][0]
    center = np.array([0, 80e-6, 40e-6])
    axis = np.array([-30e-6, 0, 30e-6])
    grid = np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T + center

    def sample(event):
        field, _, _ = basis.combine(event["amplitude_Pa"], event["phase_rad"], event["enabled"])
        return field.interpolate_complex_pressure(grid), field.interpolate_complex_gradient(grid)

    p0, g0 = sample(base)
    rms0 = np.sqrt(np.mean(abs(p0) ** 2))
    records = []
    arrays = {}
    chosen = {}
    for family, patterns in [
        (
            "phase",
            [[np.pi / 6, -np.pi / 6, 0, 0], [-np.pi / 6, np.pi / 6, 0, 0], [0, 0, np.pi / 6, -np.pi / 6]],
        ),
        ("amplitude", [[2, 0, 1, 1], [1, 1, 2, 0], [0, 2, 0, 2]]),
    ]:
        valid = []
        for i, pattern in enumerate(patterns):
            event = dict(
                time_s=0,
                amplitude_Pa=list(base["amplitude_Pa"]),
                phase_rad=[0.0, 0.0, 0.0, 0.0],
                enabled=[True] * 4,
            )
            if family == "phase":
                event["phase_rad"] = pattern
            else:
                event["amplitude_Pa"] = (np.array(base["amplitude_Pa"]) * pattern).tolist()
                event["enabled"] = (np.array(pattern) > 0).tolist()
                temp, _ = sample(event)
                scale = rms0 / np.sqrt(np.mean(abs(temp) ** 2))
                event["amplitude_Pa"] = (np.array(event["amplitude_Pa"]) * scale).tolist()
            p, g = sample(event)
            rms = float(np.sqrt(np.mean(abs(p) ** 2)))
            score = float(np.linalg.norm(g - g0) / np.linalg.norm(g0))
            relative = abs(rms / rms0 - 1)
            eligible = (
                relative <= PROTOCOL["field_control_selection"]["local_rms_pressure_relative_tolerance"]
                and min(abs(p)) >= 0.1
                and max(abs(p)) <= 40000
            )
            record = dict(
                name=f"{family}_{i}",
                family=family,
                event=event,
                local_rms_pressure_Pa=rms,
                pressure_relative_difference=float(relative),
                gradient_contrast=score,
                eligible=bool(eligible),
                pressure_min_Pa=float(min(abs(p))),
                pressure_max_Pa=float(max(abs(p))),
            )
            records.append(record)
            arrays[record["name"] + "_pressure_Pa"] = p
            arrays[record["name"] + "_gradient_Pa_m"] = g
            if eligible:
                valid.append(record)
        chosen[family] = max(valid, key=lambda r: r["gradient_contrast"]) if valid else None
    np.savez_compressed(
        out / "FIELD_DIAGNOSTICS.npz",
        points_m=grid,
        baseline_pressure_Pa=p0,
        baseline_gradient_Pa_m=g0,
        **arrays,
    )
    data = dict(
        baseline=base,
        baseline_rms_pressure_Pa=float(rms0),
        candidates=records,
        selected=chosen,
        basis=basis.metadata,
        selection_scope="27 fixed ROI points, before observing post-control trajectories; phase-only amplitudes/enables unchanged",
    )
    write(out / "FIELD_SELECTION.json", data)
    return data
