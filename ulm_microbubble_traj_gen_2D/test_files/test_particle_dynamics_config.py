from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import yaml

from ulm_microbubble_traj_gen_2D.utils.core.config import load_config


class ParticleDynamicsConfigTests(unittest.TestCase):
    def test_missing_section_uses_the_only_supported_mobility_defaults(self) -> None:
        config = self._load_config(None)

        self.assertEqual(config.particle_dynamics.time_integrator, "euler")
        self.assertEqual(config.particle_dynamics.integration_substeps, 1)
        self.assertEqual(config.particles.inlet_number_concentration_mb_per_ml, 500_000.0)
        self.assertFalse(config.cardiac_pulsatility.enabled)

    def test_explicit_cardiac_pulsatility_is_loaded_and_validated(self) -> None:
        config = self._load_config(
            None,
            cardiac_pulsatility={
                "enabled": True,
                "waveform": "SYNTHETIC_ECG_ENVELOPE",
                "bpm": 300.0,
                "pulse_propagation_velocity_um_s": 25_000.0,
                "initial_phase_fraction": 0.25,
                "waveform_samples_per_cycle": 1024,
                "preserve_cycle_mean_flow": True,
                "modulation_strength": 0.25,
            },
        )

        cardiac = config.cardiac_pulsatility
        self.assertTrue(cardiac.enabled)
        self.assertEqual(cardiac.waveform, "synthetic_ecg_envelope")
        self.assertEqual(cardiac.bpm, 300.0)
        self.assertEqual(cardiac.pulse_propagation_velocity_um_s, 25_000.0)
        self.assertEqual(cardiac.initial_phase_fraction, 0.25)
        self.assertEqual(cardiac.waveform_samples_per_cycle, 1024)
        self.assertEqual(cardiac.modulation_strength, 0.25)

    def test_invalid_cardiac_inputs_are_rejected(self) -> None:
        invalid_cases = (
            ({"enabled": True, "bpm": 0.0}, "bpm"),
            ({"enabled": True, "pulse_propagation_velocity_um_s": 0.0}, "propagation"),
            ({"enabled": True, "initial_phase_fraction": 1.0}, "phase"),
            ({"enabled": True, "waveform_samples_per_cycle": 16}, "at least 32"),
            ({"enabled": True, "modulation_strength": -0.01}, "modulation_strength"),
            ({"enabled": True, "modulation_strength": 1.01}, "modulation_strength"),
        )
        for cardiac, message in invalid_cases:
            with self.subTest(cardiac=cardiac):
                with self.assertRaisesRegex(ValueError, message):
                    self._load_config(None, cardiac_pulsatility=cardiac)

    def test_explicit_mobility_section_is_loaded_and_normalized(self) -> None:
        config = self._load_config(
            {
                "time_integrator": "HEUN",
                "integration_substeps": 8,
                "near_wall_enabled": False,
                "collisions_enabled": False,
                "xi_min": 0.002,
                "xi_near": 0.2,
                "xi_far": 1.5,
                "collision_layer_um": 0.075,
                "collision_relaxation_time_s": 0.02,
                "neighbor_search": "AUTO",
                "store_full_diagnostics": False,
                "two_wall_warning_gap_ratio": 0.75,
            }
        )

        dynamics = config.particle_dynamics
        self.assertEqual(dynamics.time_integrator, "heun")
        self.assertEqual(dynamics.integration_substeps, 8)
        self.assertFalse(dynamics.near_wall_enabled)
        self.assertFalse(dynamics.collisions_enabled)
        self.assertEqual(dynamics.xi_min, 0.002)
        self.assertEqual(dynamics.xi_near, 0.2)
        self.assertEqual(dynamics.xi_far, 1.5)
        self.assertEqual(dynamics.collision_layer_um, 0.075)
        self.assertEqual(dynamics.collision_relaxation_time_s, 0.02)
        self.assertEqual(dynamics.neighbor_search, "auto")
        self.assertFalse(dynamics.store_full_diagnostics)
        self.assertEqual(dynamics.two_wall_warning_gap_ratio, 0.75)

    def test_invalid_gap_transition_parameters_are_rejected(self) -> None:
        invalid_cases = (
            ({"xi_min": 0.0}, "0 < xi_min"),
            ({"xi_min": 0.2, "xi_near": 0.1}, "xi_near"),
            ({"xi_near": 1.0, "xi_far": 1.0}, "xi_far"),
        )
        for override, message in invalid_cases:
            with self.subTest(override=override):
                dynamics = self._valid_dynamics()
                dynamics.update(override)
                with self.assertRaisesRegex(ValueError, message):
                    self._load_config(dynamics)

    def test_invalid_integrator_is_rejected(self) -> None:
        dynamics = self._valid_dynamics()
        dynamics["time_integrator"] = "rk4"

        with self.assertRaisesRegex(ValueError, "time_integrator"):
            self._load_config(dynamics)

    def test_invalid_integration_substeps_are_rejected(self) -> None:
        invalid_cases = (
            (0, "at least one"),
            (-2, "at least one"),
            (1.5, "integer"),
            ("4", "integer"),
            (True, "integer"),
        )
        for value, message in invalid_cases:
            with self.subTest(value=value):
                dynamics = self._valid_dynamics()
                dynamics["integration_substeps"] = value
                with self.assertRaisesRegex(ValueError, message):
                    self._load_config(dynamics)

    def test_quick_test_can_reduce_the_integration_substep_count(self) -> None:
        dynamics = self._valid_dynamics()
        dynamics["integration_substeps"] = 40

        config = self._load_config(
            dynamics,
            quick_test=True,
            quick_overrides={"integration_substeps": 2},
        )

        self.assertTrue(config.quick_test)
        self.assertEqual(config.particle_dynamics.integration_substeps, 2)

    def test_invalid_collision_parameters_are_rejected(self) -> None:
        invalid_cases = (
            ({"collision_layer_um": -0.01}, "collision_layer_um"),
            ({"collision_layer_um": float("nan")}, "collision_layer_um"),
            ({"collision_relaxation_time_s": 0.0}, "collision_relaxation_time_s"),
            ({"collision_relaxation_time_s": float("inf")}, "collision_relaxation_time_s"),
            ({"two_wall_warning_gap_ratio": float("nan")}, "two_wall_warning_gap_ratio"),
        )
        for override, message in invalid_cases:
            with self.subTest(override=override):
                dynamics = self._valid_dynamics()
                dynamics.update(override)
                with self.assertRaisesRegex(ValueError, message):
                    self._load_config(dynamics)

    def test_cell_list_collision_search_is_accepted(self) -> None:
        dynamics = self._valid_dynamics()
        dynamics["neighbor_search"] = "CELL_LIST"

        config = self._load_config(dynamics)

        self.assertEqual(config.particle_dynamics.neighbor_search, "cell_list")

    def test_contact_controls_are_loaded(self) -> None:
        config = self._load_config(
            None,
            particle_overrides={
                "contact_geometry_tolerance_um": 0.0025,
                "contact_max_time_refinements": 9,
            },
        )

        self.assertEqual(config.particles.contact_geometry_tolerance_um, 0.0025)
        self.assertEqual(config.particles.contact_max_time_refinements, 9)

    def test_invalid_contact_controls_are_rejected(self) -> None:
        invalid_cases = (
            ({"contact_geometry_tolerance_um": 0.0}, "finite and positive"),
            ({"contact_geometry_tolerance_um": float("nan")}, "finite and positive"),
            ({"contact_max_time_refinements": -1}, "non-negative"),
            ({"contact_max_time_refinements": 1.5}, "integer"),
            ({"contact_max_time_refinements": "4"}, "integer"),
            ({"contact_max_time_refinements": True}, "integer"),
        )
        for particle_overrides, message in invalid_cases:
            with self.subTest(particle_overrides=particle_overrides):
                with self.assertRaisesRegex(ValueError, message):
                    self._load_config(
                        None,
                        particle_overrides=particle_overrides,
                    )

    @staticmethod
    def _valid_dynamics() -> dict[str, object]:
        return {
            "time_integrator": "euler",
            "xi_min": 0.001,
            "xi_near": 0.1,
            "xi_far": 1.0,
            "collision_layer_um": 0.05,
            "collision_relaxation_time_s": 0.01,
        }

    @staticmethod
    def _load_config(
        particle_dynamics: dict[str, object] | None,
        *,
        quick_test: bool = False,
        quick_overrides: dict[str, object] | None = None,
        cardiac_pulsatility: dict[str, object] | None = None,
        particle_overrides: dict[str, object] | None = None,
    ):
        document: dict[str, object] = {
            "input": {
                "model_dir": ".",
            },
            "output": {
                "results_dir": "results",
                "timestamp_format": "%Y%m%d_%H%M%S",
                "save_run_config": False,
                "save_npz": False,
            },
            "domain": {
                "grid_spacing_um": 1.0,
                "padding_um": 2.0,
                "min_lumen_radius_cells": 0.5,
                "max_grid_cells": 1000,
            },
            "field": {
                "hybrid_finite_element_distance_um": 2.0,
                "hybrid_transition_width_um": 2.0,
            },
            "particles": {
                "inlet_number_concentration_mb_per_ml": 500_000.0,
                "n_steps": 3,
                "dt_s": 0.001,
                "bubble_diameter_min_um": 1.5,
                "bubble_diameter_max_um": 2.5,
            },
        }
        if particle_dynamics is not None:
            document["particle_dynamics"] = particle_dynamics
        if particle_overrides is not None:
            document["particles"].update(particle_overrides)
        if cardiac_pulsatility is not None:
            document["cardiac_pulsatility"] = cardiac_pulsatility
        if quick_overrides is not None:
            document["quick_test"] = quick_overrides

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.yaml"
            path.write_text(yaml.safe_dump(document), encoding="utf-8")
            return load_config(path, quick_test=quick_test)


if __name__ == "__main__":
    unittest.main()
