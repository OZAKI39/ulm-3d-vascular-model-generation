"""SI/reduced units. Temperature, density and measured viscosity are independent gates."""
from dataclasses import dataclass, asdict
import math
import numpy as np

KB = 1.380649e-23  # exact SI J/K


def positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def kelvin(celsius):
    return positive(float(celsius) + 273.15, "temperature_K")


@dataclass(frozen=True)
class Units:
    L0: float
    M0: float
    t0: float

    def __post_init__(self):
        for key, value in asdict(self).items():
            positive(value, key)

    @property
    def scales(self):
        L, M, t = self.L0, self.M0, self.t0
        return dict(length=L, area=L**2, volume=L**3, time=t, velocity=L/t,
                    acceleration=L/t**2, particle_mass=M, mass_density=M/L**3,
                    number_density=L**-3, force=M*L/t**2, pressure=M/(L*t**2),
                    dynamic_viscosity=M/(L*t), kinematic_viscosity=L**2/t,
                    volume_flow=L**3/t, mass_flow=M/t, energy=M*L**2/t**2)

    def to_si(self, value, quantity):
        a = np.asarray(value, dtype=float)
        if not np.isfinite(a).all():
            raise ValueError("nonfinite quantity")
        out = a * self.scales[quantity]
        if not np.isfinite(out).all():
            raise ValueError("overflow in unit conversion")
        return float(out) if out.ndim == 0 else out

    def to_star(self, value, quantity):
        a = np.asarray(value, dtype=float)
        if not np.isfinite(a).all():
            raise ValueError("nonfinite quantity")
        out = a / self.scales[quantity]
        if not np.isfinite(out).all():
            raise ValueError("overflow in unit conversion")
        return float(out) if out.ndim == 0 else out

    def temperature_K(self, kBT_star):
        return self.to_si(positive(kBT_star, "kBT_star"), "energy") / KB

    def as_dict(self):
        return {**asdict(self), "E0": self.scales["energy"], "si_per_star": self.scales}


def thermal_density_units(L0, rho_si, n_star, m_star, kBT_star, temperature_K):
    for name, value in locals().copy().items():
        positive(value, name)
    M0 = rho_si * L0**3 / (n_star * m_star)
    E0 = KB * temperature_K / kBT_star
    return Units(L0, M0, math.sqrt(M0 * L0**2 / E0))


def consistency(units, *, n_star, m_star, kBT_star, rho_si, temperature_K,
                nu_star=None, nu_si=None, viscosity_tolerance=0.1):
    rho = units.to_si(n_star * m_star, "mass_density")
    T = units.temperature_K(kBT_star)
    result = {"density_si": rho, "temperature_K": T,
              "density_consistent": math.isclose(rho, rho_si, rel_tol=1e-10),
              "thermal_consistent": math.isclose(T, temperature_K, rel_tol=1e-10)}
    if nu_star is not None and nu_si is not None:
        nu = units.to_si(nu_star, "kinematic_viscosity")
        result.update(nu_si=nu, viscosity_relative_error=abs(nu/nu_si-1),
                      viscosity_consistent=abs(nu/nu_si-1) <= viscosity_tolerance)
    result["all_consistent"] = all(v for k, v in result.items() if k.endswith("consistent"))
    return result


def pressure_target_star(p_reference_star, gauge_pa, gauge_reference_pa, units):
    return p_reference_star + units.to_star(gauge_pa - gauge_reference_pa, "pressure")


def viscosity_only_alternative(units, measured_nu_star, target_nu_si, kBT_star):
    """Diagnostic only; NEVER used to remap a sensitivity run or select a candidate."""
    alt = Units(units.L0, units.M0, positive(measured_nu_star, "nu") * units.L0**2 / target_nu_si)
    return {"status": "THERMAL_MAPPING_UNRESOLVED", "selection_allowed": False,
            "t0_s": alt.t0, "implied_temperature_K": alt.temperature_K(kBT_star),
            "reason": "Matching viscosity by changing time units breaks the locked 298.15 K mapping."}
