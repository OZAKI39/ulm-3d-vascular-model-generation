"""Provenance-bearing SI shell parameters; no implicit extrapolation or calibration."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class CitedParameter:
    value: float
    unit: str
    source_doi: str
    source_table: str
    source_condition: str
    method: str
    uncertainty: tuple[float, float] | None = None
    uncertainty_kind: str = "not reported"

    def __post_init__(self):
        if not np.isfinite(self.value) or not self.source_doi or not self.method:
            raise ValueError("Finite value and explicit provenance required")
        if self.uncertainty is not None:
            lo, hi = self.uncertainty
            if not np.isfinite([lo, hi]).all() or not lo <= self.value <= hi:
                raise ValueError("Uncertainty bounds must contain the value")


@dataclass(frozen=True)
class SonoVueParameterSet:
    """A cited case, not a universal material constant; agent identity is explicit."""

    name: str
    agent: str
    radius_m: CitedParameter
    shell_elasticity_N_m: CitedParameter
    shell_viscosity_kg_s: CitedParameter
    frequency_Hz: CitedParameter
    pressure_Pa: CitedParameter
    fitted_model: str

    def __post_init__(self):
        for name, unit in [
            ("radius_m", "m"),
            ("shell_elasticity_N_m", "N/m"),
            ("shell_viscosity_kg_s", "kg/s"),
            ("frequency_Hz", "Hz"),
            ("pressure_Pa", "Pa"),
        ]:
            p = getattr(self, name)
            if p.unit != unit or p.value <= 0:
                raise ValueError(f"{name} must be positive in {unit}")

    @classmethod
    def from_dict(cls, data):
        data = dict(data)
        for k, v in data.items():
            if isinstance(v, dict):
                data[k] = CitedParameter(**v)
        return cls(**data)

    def to_dict(self):
        return asdict(self)


def load_parameter_cases(path: str | Path) -> dict[str, SonoVueParameterSet]:
    data = json.loads(Path(path).read_text())
    return {c["name"]: SonoVueParameterSet.from_dict(c) for c in data["fixed_cases"]}


class SonoVueParameterDistribution:
    """Fig. 4 interpolation, with digitization uncertainty, not a population fit.

    Radius itself must be supplied by a separate measured size distribution.
    Sampling here varies shell parameters conditional on radius only. No claim
    of biological confidence intervals or independent material variability.
    """

    def __init__(self, path: str | Path):
        self.data = json.loads(Path(path).read_text())["radius_dependence"]
        self._axes = {}
        for key in ("elasticity", "viscosity"):
            a = np.asarray(self.data[key]["samples_SI"], float)
            if a.ndim != 2 or a.shape[1] != 2 or len(a) < 2:
                raise ValueError("At least two radius/value pairs required")
            if not np.isfinite(a).all() or np.any(a <= 0) or np.any(np.diff(a[:, 0]) <= 0):
                raise ValueError("Positive samples and strictly increasing radii required")
            self._axes[key] = a
        self.radius_bounds_m = (
            max(a[0, 0] for a in self._axes.values()),
            min(a[-1, 0] for a in self._axes.values()),
        )

    def evaluate_batch(self, radii_m, *, uncertainty_draw=False, seed=None):
        r = np.asarray(radii_m, dtype=float)
        lo, hi = self.radius_bounds_m
        if not np.isfinite(r).all() or np.any((r < lo) | (r > hi)):
            raise ValueError(f"Outside digitized SonoVue radius range [{lo}, {hi}] m")
        rng = np.random.default_rng(seed)
        out = {}
        for key, a in self._axes.items():
            log = key == "viscosity"
            y = np.log10(a[:, 1]) if log else a[:, 1]
            z = np.interp(r, a[:, 0], y)
            e = self.data[key]["digitization_uncertainty_log10" if log else "digitization_uncertainty_SI"]
            if uncertainty_draw:
                z = z + rng.uniform(-e, e, size=r.shape)
            out[key] = 10**z if log else z
            out[key + "_bounds"] = np.stack(
                (10 ** (z - e), 10 ** (z + e)) if log else (z - e, z + e), axis=-1
            )
        return out

    def parameter_set(self, radius_m):
        p = self.evaluate_batch(radius_m)
        provenance = dict(
            source_doi=self.data["source_doi"],
            source_table="Figure 4",
            source_condition=self.data["source_condition"],
            method="interpolation of manually digitized figure markers",
        )

        def cited(value, unit, bounds=None):
            return CitedParameter(
                float(value),
                unit,
                **provenance,
                uncertainty=None if bounds is None else tuple(map(float, bounds)),
                uncertainty_kind="digitization only; not biological uncertainty",
            )

        return SonoVueParameterSet(
            f"fig4_R{radius_m:.9g}",
            "SonoVue",
            cited(radius_m, "m"),
            cited(p["elasticity"], "N/m", p["elasticity_bounds"]),
            cited(p["viscosity"], "kg/s", p["viscosity_bounds"]),
            cited(2.5e6, "Hz"),
            cited(150e3, "Pa"),
            "Tu2009 Eq6 linearized shell law; not full buckling/rupture calibration",
        )
