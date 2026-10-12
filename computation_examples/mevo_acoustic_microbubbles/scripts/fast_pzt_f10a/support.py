"""Real CFD plotting samples and independent nearfield diagnostics; no integration."""

import argparse
import time
import numpy as np
from mevo_mb.fast_sim.consistent_field import query_numpy
from mevo_mb.fast_sim.nearfield import coupled_fundamental, feedback_weight
from .common import BASE, objects, js, write, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    t = time.perf_counter()
    p, b, r, g, flow = objects()
    s = js(BASE / "static_v1/RESULTS.json")
    x0 = np.array(s["states"]["near"]["positions_m"]).mean(0)
    axis = np.linspace(-0.0014, 0.0014, 81)
    yy, zz = np.meshgrid(axis, axis)
    xyz = np.c_[np.full(yy.size, x0[0]), yy.ravel(), zz.ravel()]
    cf = flow.sample(xyz)
    geometry = g.query(xyz)
    np.savez_compressed(
        out / "CFD_cross_section.npz",
        positions_m=xyz,
        velocity_m_s=cf["velocity_m_s"],
        valid=cf["valid"],
        inside=geometry["inside"],
        grid_shape=yy.shape,
        wall_gap_m=geometry["distance_m"] - p.radius_m,
    )
    rows = []
    chi = r.query(np.array([r.knots[0] + 0j]))[0, 0] / r.knots[0]
    for row in s["rows"]:
        if row.get("status") != "VALID_STATIC":
            continue
        x = np.array(s["states"][row["state"]]["positions_m"])
        ev = row["event"]
        field, _, _ = b.combine(ev["amplitude_Pa"], ev["phase_rad"], ev["enabled"])
        pp, _ = query_numpy(field.pressure_Pa, b.origin, b.spacing, x)
        h = r.query(pp)
        dets = []
        effective = []
        for i in range(32):
            for j in range(i + 1, 32):
                d = np.linalg.norm(x[j] - x[i])
                if feedback_weight(d, p.radius_m, p.feedback_end_over_sum) > 0:
                    _, _, pi, pj, det = coupled_fundamental(h[i, 0], h[j, 0], pp[i], pp[j], d, p, chi)
                    dets.append(abs(det))
                    effective.extend([abs(pi), abs(pj)])
        rows.append(
            dict(
                state=row["state"],
                mode=row["mode"],
                pressure_Pa=row["pressure_Pa"],
                min_closure_determinant=min(dets),
                min_effective_pair_pressure_Pa=min(effective),
                max_effective_pair_pressure_Pa=max(effective),
            )
        )
    write(
        out / "CHECKS.json",
        dict(
            identity=identity(),
            rows=rows,
            process_wall_s=time.perf_counter() - t,
            static_raw_near_stats_limitation="Inherited predictor initialized min reduction buffers to zero and accumulated maxima across candidates; do not interpret those raw diagnostics as selected-candidate minima. This independent pair-by-pair audit supplies valid selected diagnostics. GPU error-bit legality checks were unaffected.",
        ),
    )


if __name__ == "__main__":
    main()
