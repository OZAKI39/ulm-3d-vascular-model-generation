"""Frozen, solver-free scientific contract for validated Tau=1 CFD.

This module owns formulas and input constants only.  Accepted simulation
results and artifact hashes belong to :mod:`validated_evidence`.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping


CONTRACT_REVISION = "VALIDATED_TAU1_CFD_CONTRACT_V1"
FULL_TIMESTEP_REFEREE_REVISION = "MUSUBI_ONE_STEP_DISCRETE_MASS_IDENTITY_V2"
PHYSICAL_PLANE_REVISION = "STANDARDIZED_INTERIOR_PHYSICAL_PORT_PLANES_V3"
PHYSICAL_FLUX_REVISION = "CONTINUOUS_APERTURE_GAUSS_MLS_QUADRATIC_V2"

RHO0_KG_M3 = 1056.0
KINEMATIC_VISCOSITY_M2_S = 3.27e-6
BULK_VISCOSITY_M2_S = 2.18e-6
LATTICE_CS2 = 1.0 / 3.0
BASE_DX_M = 2.0e-7
TARGET_MASS_FLOW_KG_S = 2.890180380479642e-12
TARGET_VOLUME_FLOW_M3_S = 2.7369132390905703e-15
OUTLET_GAUGE_PRESSURES_PA = {
    "outlet_01": 14.544978101274268,
    "outlet_02": 132.20454922317552,
    "outlet_03": -13.700626673311461,
}
PORTS = ("inlet", "outlet_01", "outlet_02", "outlet_03")
OUTLETS = PORTS[1:]

MASS_GATE = 0.01
VELOCITY_GATE = 0.01
PRESSURE_GATE = 0.005
INLET_GATE = 0.01
FLOW_FRACTION_DRIFT_GATE = 0.01
BOUNDARY_WINDOW_CLOSURE_GATE = 0.001
CONTROLLER_GATE = 1.0e-8
FULL_TIMESTEP_IDENTITY_GATE = 1.0e-8
Q_DENSITY_GATE = 0.01
RHO_LATTICE_RANGE = (0.9, 1.1)
MAXIMUM_LATTICE_SPEED = 0.05
SIGNIFICANT_BACKFLOW_FRACTION = 0.05

# Regression-only historical baseline.  It must never be used as the Tau=1
# production reference pressure.
HISTORICAL_DT_S = 2.44140625e-8
HISTORICAL_FIXED_REFERENCE_PRESSURE_PA = 23622.320128


def tau1_time_step_s(
    dx_m: float,
    kinematic_viscosity_m2_s: float = KINEMATIC_VISCOSITY_M2_S,
) -> float:
    """Return the diffusive timestep that gives nu_lattice=1/6 and tau=1."""

    if dx_m <= 0.0 or kinematic_viscosity_m2_s <= 0.0:
        raise ValueError("dx and kinematic viscosity must be positive")
    return float(dx_m) ** 2 / (6.0 * float(kinematic_viscosity_m2_s))


def pressure_reference_pa(
    dx_m: float,
    dt_s: float,
    *,
    rho0_kg_m3: float = RHO0_KG_M3,
    cs2: float = LATTICE_CS2,
) -> float:
    """Return the LBM numerical pressure offset, not physiological pressure."""

    if dx_m <= 0.0 or dt_s <= 0.0 or rho0_kg_m3 <= 0.0 or cs2 <= 0.0:
        raise ValueError("pressure scaling inputs must be positive")
    return float(rho0_kg_m3) * float(cs2) * (float(dx_m) / float(dt_s)) ** 2


def target_lattice_flux(
    target_mass_flow_kg_s: float,
    dx_m: float,
    dt_s: float,
    *,
    rho0_kg_m3: float = RHO0_KG_M3,
) -> float:
    """Convert the physical mass-flow target to Musubi lattice volume flux."""

    if target_mass_flow_kg_s <= 0.0 or dx_m <= 0.0 or dt_s <= 0.0 or rho0_kg_m3 <= 0.0:
        raise ValueError("target scaling inputs must be positive")
    return (
        float(target_mass_flow_kg_s)
        / float(rho0_kg_m3)
        * float(dt_s)
        / float(dx_m) ** 3
    )


def relaxation_from_physical(
    dx_m: float,
    dt_s: float,
    *,
    kinematic_viscosity_m2_s: float = KINEMATIC_VISCOSITY_M2_S,
    cs2: float = LATTICE_CS2,
) -> tuple[float, float, float]:
    """Return ``(nu_lattice, tau, omega)`` for D3Q19 BGK scaling."""

    if min(dx_m, dt_s, kinematic_viscosity_m2_s, cs2) <= 0.0:
        raise ValueError("relaxation scaling inputs must be positive")
    nu_lattice = float(kinematic_viscosity_m2_s) * float(dt_s) / float(dx_m) ** 2
    tau = nu_lattice / float(cs2) + 0.5
    return nu_lattice, tau, 1.0 / tau


def gauge_pressure_pa(physical_pressure_pa: float, reference_pressure_pa: float) -> float:
    return float(physical_pressure_pa) - float(reference_pressure_pa)


def absolute_pressure_pa(gauge_pa: float, reference_pressure_pa: float) -> float:
    return float(reference_pressure_pa) + float(gauge_pa)


@dataclass(frozen=True, slots=True)
class ValidatedTau1Contract:
    """Resolved validated inputs for one spacing, defaulting to Base."""

    dx_m: float = BASE_DX_M
    rho0_kg_m3: float = RHO0_KG_M3
    kinematic_viscosity_m2_s: float = KINEMATIC_VISCOSITY_M2_S
    bulk_viscosity_m2_s: float = BULK_VISCOSITY_M2_S
    cs2: float = LATTICE_CS2
    target_mass_flow_kg_s: float = TARGET_MASS_FLOW_KG_S

    @property
    def dt_s(self) -> float:
        return tau1_time_step_s(self.dx_m, self.kinematic_viscosity_m2_s)

    @property
    def pressure_reference_pa(self) -> float:
        return pressure_reference_pa(
            self.dx_m, self.dt_s, rho0_kg_m3=self.rho0_kg_m3, cs2=self.cs2
        )

    @property
    def nu_lattice(self) -> float:
        return relaxation_from_physical(
            self.dx_m,
            self.dt_s,
            kinematic_viscosity_m2_s=self.kinematic_viscosity_m2_s,
            cs2=self.cs2,
        )[0]

    @property
    def tau(self) -> float:
        return relaxation_from_physical(
            self.dx_m,
            self.dt_s,
            kinematic_viscosity_m2_s=self.kinematic_viscosity_m2_s,
            cs2=self.cs2,
        )[1]

    @property
    def omega(self) -> float:
        return relaxation_from_physical(
            self.dx_m,
            self.dt_s,
            kinematic_viscosity_m2_s=self.kinematic_viscosity_m2_s,
            cs2=self.cs2,
        )[2]

    @property
    def target_lattice_flux(self) -> float:
        return target_lattice_flux(
            self.target_mass_flow_kg_s,
            self.dx_m,
            self.dt_s,
            rho0_kg_m3=self.rho0_kg_m3,
        )

    @property
    def outlet_absolute_pressures_pa(self) -> dict[str, float]:
        return {
            label: absolute_pressure_pa(gauge, self.pressure_reference_pa)
            for label, gauge in OUTLET_GAUGE_PRESSURES_PA.items()
        }

    def lattice_density(self, physical_pressure_pa: float) -> float:
        unit_density_pressure = self.rho0_kg_m3 * self.cs2 * (self.dx_m / self.dt_s) ** 2
        return float(physical_pressure_pa) / unit_density_pressure

    def as_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value.update(
            {
                "revision": CONTRACT_REVISION,
                "dt_s": self.dt_s,
                "pressure_reference_pa": self.pressure_reference_pa,
                "nu_lattice": self.nu_lattice,
                "tau": self.tau,
                "omega": self.omega,
                "target_lattice_flux": self.target_lattice_flux,
                "outlet_gauge_pressures_pa": dict(OUTLET_GAUGE_PRESSURES_PA),
                "outlet_absolute_pressures_pa": self.outlet_absolute_pressures_pa,
            }
        )
        return value


def restart_compatibility(
    saved: Mapping[str, Any], expected: Mapping[str, Any]
) -> dict[str, Any]:
    """Compare fields that make a Tau=1 restart scientifically compatible."""

    fields = (
        "mesh_hashes",
        "dx_m",
        "dt_s",
        "rho0_kg_m3",
        "nu_m2_s",
        "bulk_nu_m2_s",
        "tau",
        "omega",
        "pressure_reference_pa",
        "boundary_contract",
        "outlet_gauge_pressure_pa",
        "outlet_absolute_pressure_pa",
        "target_mass_flow_kg_s",
        "binary_sha256",
        "physical_plane_contract_sha256",
        "layout",
        "collision",
    )
    checks = {field: saved.get(field) == expected.get(field) for field in fields}
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}

