"""Evidence-scoped maps: a sampled uncoated benchmark is not a SonoVue certificate."""

import json
from pathlib import Path

import numpy as np

REGIMES = (
    "FAR_FIELD_SAFE",
    "TRANSITION",
    "NEAR_FIELD_REFERENCE_REQUIRED",
    "CONTACT",
    "OUTSIDE_VALIDATED_RANGE",
)


def equal(a, b):
    return bool(np.isclose(a, b, rtol=1e-10, atol=0))


class PairPhysicsRegimeMap:
    def __init__(self, path, safe_error=0.05, transition_error=0.20, max_kd=0.3):
        self.path = Path(path)
        self.rows = json.loads(self.path.read_text())["rows"]
        if not 0 < safe_error < transition_error:
            raise ValueError("Ordered positive error budgets required")
        self.safe, self.transition, self.max_kd = safe_error, transition_error, max_kd

    def classify(self, R1, R2, separation, frequency, pressure, shell_model, sound_speed=1482.0):
        values = [R1, R2, separation, frequency, pressure, sound_speed]
        if not np.isfinite(values).all() or min(values[:4]) <= 0 or pressure < 0 or sound_speed <= 0:
            raise ValueError("Invalid pair inputs")
        if separation <= R1 + R2:
            return "CONTACT"
        # No monotone/error-smoothness assumption in frequency, pressure or shell.
        rows = [
            r
            for r in self.rows
            if r["shell_model"] == shell_model
            and equal(r["sound_speed_m_s"], sound_speed)
            and equal(r["radius1_m"], min(R1, R2))
            and equal(r["radius2_m"], max(R1, R2))
            and equal(r["frequency_Hz"], frequency)
            and equal(r["pressure_Pa"], pressure)
        ]
        if not rows or 2 * np.pi * frequency * separation / sound_speed > self.max_kd:
            return "OUTSIDE_VALIDATED_RANGE"
        ds = sorted({r["separation_m"] for r in rows})
        if separation < ds[0]:
            return "NEAR_FIELD_REFERENCE_REQUIRED"
        if separation > ds[-1]:
            return "OUTSIDE_VALIDATED_RANGE"
        # Only measured points are certified; untested intervals remain outside.
        selected = [r for r in rows if equal(r["separation_m"], separation)]
        if {r["model"] for r in selected} != {"zhang", "doinikov"}:
            return "OUTSIDE_VALIDATED_RANGE"
        if any(r["max_radial_fraction"] > 0.02 or r["relative_distance_change"] > 0.01 for r in selected):
            return "OUTSIDE_VALIDATED_RANGE"
        if any(r["calculated_sign"] != r["classical_sign"] for r in selected):
            return "NEAR_FIELD_REFERENCE_REQUIRED"
        error = max(max(r["coefficient_relative_error"], r["hydrodynamic_relative_error"]) for r in selected)
        return (
            "FAR_FIELD_SAFE"
            if error <= self.safe
            else ("TRANSITION" if error <= self.transition else "NEAR_FIELD_REFERENCE_REQUIRED")
        )

    def population(self, positions, radii, frequency, pressures, shell_model, sound_speed):
        x, r = np.asarray(positions), np.asarray(radii)
        p = np.broadcast_to(np.asarray(pressures), r.shape)
        counts = {key: 0 for key in REGIMES}
        for i in range(len(r)):
            for j in range(i + 1, len(r)):
                if not equal(p[i], p[j]):
                    key = "OUTSIDE_VALIDATED_RANGE"  # benchmark was equal drive, not pressure/phase gradients
                else:
                    key = self.classify(
                        r[i], r[j], np.linalg.norm(x[i] - x[j]), frequency, p[i], shell_model, sound_speed
                    )
                counts[key] += 1
        n = sum(counts.values())

        def fraction(k):
            return counts[k] / n if n else 0.0

        return dict(
            far_field_fraction=fraction(REGIMES[0]) if n else 1.0,
            transition_fraction=fraction(REGIMES[1]),
            near_field_fraction=fraction(REGIMES[2]),
            contact_fraction=fraction(REGIMES[3]),
            outside_validated_fraction=fraction(REGIMES[4]),
            pair_count=n,
            evidence=str(self.path),
        )


class TranslationValidityMap:
    def __init__(self, path):
        self.data = json.loads(Path(path).read_text())

    def classify(self, radius, pressure, frequency, *, shell_model, cycles=20):
        if (
            not np.isfinite([radius, pressure, frequency]).all()
            or radius <= 0
            or pressure < 0
            or frequency <= 0
        ):
            raise ValueError("Invalid translation inputs")
        for r in self.data["rows"]:
            if (
                r.get("shell_model") == shell_model
                and r.get("cycles") == cycles
                and equal(r["radius_m"], radius)
                and equal(r["pressure_Pa"], pressure)
                and equal(r["frequency_Hz"], frequency)
            ):
                return r["status"]
        return "OUTSIDE_VALIDATED_RANGE"
