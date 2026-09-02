"""Predeclared dimensionless scenarios for molecular-binding sensitivity studies.

The association rate is an intrinsic input of the molecular model.  It must not
change merely because a particular transport run happened to produce more,
less, or no target contact.  This module therefore constructs every scenario
from axes fixed before the transport run, including an explicit reference time
used only to define the chosen dimensionless Damkoehler number.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
from typing import Sequence

import numpy as np


MOLECULES_PER_UM2_TO_MOLECULES_PER_M2 = 1.0e12


@dataclass(frozen=True)
class DaOnScenario:
    """One predeclared dimensionless association scenario."""

    scenario_index: int
    da_on: float
    da_on_reference_time_s: float
    target_density_molecules_per_um2: float
    target_density_molecules_per_m2: float
    ligand_density_molecules_per_um2: float | None
    ligand_density_molecules_per_m2: float | None
    capture_distance_to_rest_length_ratio: float | None
    rest_length_um: float | None
    capture_distance_um: float | None
    association_rate_m2_per_molecule_s: float


def molecules_per_um2_to_per_m2(value: float) -> float:
    """Convert a surface density from molecule/um^2 to molecule/m^2."""

    density = _finite_positive_float(value, "surface density in molecule/um^2")
    return density * MOLECULES_PER_UM2_TO_MOLECULES_PER_M2


def build_da_on_scenarios(
    *,
    da_on_reference_time_s: float,
    da_on_levels: Sequence[float],
    target_density_molecules_per_um2_levels: Sequence[float],
    ligand_density_molecules_per_um2_levels: Sequence[float] = (),
    capture_distance_to_rest_length_ratios: Sequence[float] = (),
    rest_length_um: float | None = None,
) -> tuple[DaOnScenario, ...]:
    """Build a Cartesian scenario table without reading trajectory outcomes.

    The effective two-surface association rate is defined by

    ``k_on = Da_on / (rho_T * t_ref)``.

    Here ``t_ref`` is the explicit ``da_on_reference_time_s`` supplied before
    simulation.  It is not an observed contact duration from the same run.
    """

    reference_time_s = _finite_positive_float(
        da_on_reference_time_s,
        "da_on_reference_time_s",
    )
    da_levels = _positive_level_tuple(da_on_levels, "da_on_levels")
    target_levels_um2 = _positive_level_tuple(
        target_density_molecules_per_um2_levels,
        "target_density_molecules_per_um2_levels",
    )
    ligand_levels_um2 = _optional_positive_level_tuple(
        ligand_density_molecules_per_um2_levels,
        "ligand_density_molecules_per_um2_levels",
    )
    capture_ratios = _optional_positive_level_tuple(
        capture_distance_to_rest_length_ratios,
        "capture_distance_to_rest_length_ratios",
    )
    rest_length_value = (
        None
        if rest_length_um is None
        else _finite_positive_float(rest_length_um, "rest_length_um")
    )
    if capture_ratios != (None,) and rest_length_value is None:
        raise ValueError(
            "rest_length_um is required when capture-distance/rest-length ratios are provided."
        )

    scenarios: list[DaOnScenario] = []
    axes = product(da_levels, target_levels_um2, ligand_levels_um2, capture_ratios)
    for scenario_index, (da_on, rho_target_um2, rho_ligand_um2, capture_ratio) in enumerate(axes):
        rho_target_m2 = molecules_per_um2_to_per_m2(rho_target_um2)
        rho_ligand_m2 = (
            None
            if rho_ligand_um2 is None
            else molecules_per_um2_to_per_m2(rho_ligand_um2)
        )
        capture_distance_um = (
            None
            if capture_ratio is None or rest_length_value is None
            else capture_ratio * rest_length_value
        )
        scenarios.append(
            DaOnScenario(
                scenario_index=scenario_index,
                da_on=da_on,
                da_on_reference_time_s=reference_time_s,
                target_density_molecules_per_um2=rho_target_um2,
                target_density_molecules_per_m2=rho_target_m2,
                ligand_density_molecules_per_um2=rho_ligand_um2,
                ligand_density_molecules_per_m2=rho_ligand_m2,
                capture_distance_to_rest_length_ratio=capture_ratio,
                rest_length_um=rest_length_value,
                capture_distance_um=capture_distance_um,
                association_rate_m2_per_molecule_s=(
                    da_on / (rho_target_m2 * reference_time_s)
                ),
            )
        )
    return tuple(scenarios)


def _finite_positive_float(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be finite and positive.")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite and positive.") from exc
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")
    return number


def _positive_level_tuple(values: Sequence[float], name: str) -> tuple[float, ...]:
    levels = _optional_positive_level_tuple(values, name)
    if levels == (None,):
        raise ValueError(f"{name} must not be empty.")
    return tuple(float(level) for level in levels)


def _optional_positive_level_tuple(
    values: Sequence[float],
    name: str,
) -> tuple[float | None, ...]:
    if isinstance(values, (str, bytes)):
        raise ValueError(f"{name} must be a sequence of finite positive values.")
    levels = tuple(values)
    if not levels:
        return (None,)
    return tuple(_finite_positive_float(value, name) for value in levels)
