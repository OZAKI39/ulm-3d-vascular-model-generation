from __future__ import annotations

import math

import numpy as np

from utils.cfd_flow.musubi_wall_geometry_benchmark import (
    early_stop_decision,
    effective_radius_ratio,
    observed_order_three,
    poiseuille_pressure_gradient_pa_m,
    poiseuille_velocity_profile,
    source_force_vector,
    wall_stage_decision,
)


def test_periodic_force_conversion_and_poiseuille_profile() -> None:
    gradient = poiseuille_pressure_gradient_pa_m(
        mean_velocity_m_s=0.35e-3, radius_m=1.5e-6, rho_kg_m3=1056.0, nu_m2_s=3.27e-6
    )
    force = source_force_vector(gradient, (0.0, 0.0, 2.0))
    assert np.allclose(force, (0.0, 0.0, gradient))
    velocity = poiseuille_velocity_profile(
        np.asarray((0.0, 1.5e-6)), radius_m=1.5e-6, mean_velocity_m_s=0.35e-3
    )
    assert np.allclose(velocity, (0.70e-3, 0.0))
    assert math.isclose(velocity[0] / 0.35e-3, 2.0)


def test_effective_radius_and_observed_order() -> None:
    assert math.isclose(effective_radius_ratio(16.0, 1.0), 2.0)
    assert math.isclose(observed_order_three((0.4, 0.2, 0.1), (0.16, 0.04, 0.01)), 2.0)


def test_early_stop_and_wall_decision_tree() -> None:
    samples = [
        {
            "iteration": float(iteration),
            "mean_axial_velocity": 1.0 + 1.0e-7 * iteration,
            "profile_l2_error": 0.1 + 1.0e-8 * iteration,
            "minimum_pdf": 0.01,
            "maximum_lattice_speed": 0.001,
            "all_finite": True,
        }
        for iteration in range(0, 221, 20)
    ]
    assert early_stop_decision(samples)["stop"] is True
    assert wall_stage_decision(0.01, oblique_n27_error=0.015) == "WALL_GEOMETRY_ISOLATED_PASS"
    assert wall_stage_decision(0.06) == "CFD_FLOW_WALL_GEOMETRY_ERROR_IDENTIFIED"


def test_module_has_no_vascular_launcher() -> None:
    import utils.cfd_flow.musubi_wall_geometry_benchmark as module

    assert not hasattr(module, "run_vascular_cfd")
