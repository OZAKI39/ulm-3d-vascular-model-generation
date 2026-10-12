"""Conservative cache scope. Rejection does not reject the full M0 solver."""

import numpy as np
from .metrics import frame_metrics
from .consistent_field import query_numpy


def admit_cache(x, event, objects, mobility="f5_local"):
    p, b, r, g, f = objects
    x = np.asarray(x)
    reasons = []
    if mobility != "f5_local":
        reasons.append("M0_ONLY")
    wall = g.query(x)
    flow = f.sample(x)
    metric = frame_metrics(x, p.radius_m)
    gap = float(np.min(wall["distance_m"] - p.radius_m))
    grad = flow["gradient_s_inv"]
    d = 0.5 * (grad + grad.transpose(0, 2, 1))
    shear = float(np.sqrt(2 * np.sum(d * d, axis=(1, 2))).max())
    if not wall["inside"].all() or gap < 1e-4:
        reasons.append("NEAR_WALL_OR_GEOMETRY_BOUNDARY")
    if flow["error"] or not flow["valid"].all():
        reasons.append("CFD_DOMAIN")
    if shear > 2:
        reasons.append("HIGH_SHEAR")
    if metric["largest_size"] != len(x):
        reasons.append("DISPERSED_COMPONENTS")
    acoustic = bool(np.any(np.asarray(event["amplitude_Pa"]) * np.asarray(event["enabled"])))
    pressure = None
    if acoustic:
        try:
            field, _, _ = b.combine(event["amplitude_Pa"], event["phase_rad"], event["enabled"])
            value, _ = query_numpy(field.pressure_Pa, b.origin, b.spacing, x)
            r.query(value)
            pressure = float(np.abs(value).max())
            if pressure > 1500:
                reasons.append("UNVALIDATED_HIGH_DRIVE")
        except ValueError as exc:
            reasons.append("ACOUSTIC_DOMAIN_OR_RESPONSE_TABLE:" + str(exc))
    return dict(
        status="MULTIRATE_NOT_ADMITTED" if reasons else "ADMITTED_FOR_NUMERICAL_TEST",
        reasons=reasons,
        min_wall_gap_m=gap,
        max_shear_s_inv=shear,
        largest_size=metric["largest_size"],
        max_pressure_Pa=pressure,
        training_ready=False,
    )
