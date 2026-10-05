"""Generate dimensional analytic resistance tables, with exclusive output paths."""

from datetime import datetime, timezone
import csv

import numpy as np

from mevo_mb.bubbles.translation_backends import LubricationCoupledTranslation
from scripts.stage4l.dense_probe import BASE, write


def main():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    radius, mu = 1.78e-6, 0.001
    model = LubricationCoupledTranslation(mu)
    for kind in ["pair", "wall"]:
        out = BASE / ("analytic_" + kind) / stamp
        out.mkdir(parents=True, exist_ok=False)
        rows = []
        for eps in [0.2, 0.1, 0.05, 0.01, 0.001, 0.0005, 0.0]:
            e = max(eps, 0.001)
            if kind == "pair":

                def phi(v):
                    return 1 / (2 * v) - 9 / 20 * np.log(v) - 3 / 56 * v * np.log(v)

                expected = 3 * np.pi * mu * radius * (phi(e) - phi(0.1)) if eps < 0.1 else 0.0
            else:

                def phi(v):
                    return 1 / v - 1 / 5 * np.log(v) - 1 / 21 * v * np.log(v)

                expected = 6 * np.pi * mu * radius * (phi(e) - phi(0.1)) if eps < 0.1 else 0.0
            actual = float(model.coefficient(eps * radius, radius, kind))
            if not np.isclose(actual, expected, rtol=1e-13, atol=1e-22):
                raise AssertionError((kind, eps, actual, expected))
            rows.append(
                dict(
                    gap_over_radius=eps,
                    gap_m=eps * radius,
                    coefficient_kg_s=actual,
                    analytic_coefficient_kg_s=expected,
                    absolute_error_kg_s=abs(actual - expected),
                )
            )
        with (out / "coefficients.csv").open("x", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        write(
            out / "SUMMARY.json",
            dict(
                status="PASS",
                radius_m=radius,
                viscosity_Pa_s=mu,
                velocity_convention="relative pair normal speed"
                if kind == "pair"
                else "particle wall-normal speed",
                max_absolute_error_kg_s=max(row["absolute_error_kg_s"] for row in rows),
                rows=rows,
            ),
        )


if __name__ == "__main__":
    main()
