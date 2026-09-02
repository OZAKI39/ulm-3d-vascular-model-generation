from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import yaml

from ulm_microbubble_traj_gen_2D.utils.core.config import load_config


class MolecularConfigTests(unittest.TestCase):
    def test_global_and_target_selection_random_seeds_are_independent(self) -> None:
        config = self._load_config(
            random_seed=314159,
            molecular_target_selection={"random_seed": 271828},
        )

        self.assertEqual(config.random_seed, 314159)
        self.assertEqual(config.molecular_target_selection.random_seed, 271828)

    def test_missing_optional_sections_preserve_disabled_defaults(self) -> None:
        config = self._load_config()

        self.assertFalse(config.molecular_target.enabled)
        self.assertEqual(config.molecular_target.region_mode, "none")
        self.assertFalse(config.molecular_binding.enabled)
        self.assertEqual(config.molecular_binding.model, "deterministic_mean_field")
        self.assertFalse(config.binding_scenario_sweep.enabled)
        self.assertEqual(config.binding_scenario_sweep.da_on_reference_time_s, 0.0)
        self.assertEqual(config.molecular_target_selection.default_mode, "manual")
        self.assertIsNone(
            config.molecular_target_selection
            .influence_region_endothelial_wall_area_fraction
        )

    def test_automatic_target_selection_requires_and_loads_v10_spatial_parameters(self) -> None:
        config = self._load_config(
            molecular_target_selection={
                "default_mode": "AUTOMATIC",
                "influence_region_endothelial_wall_area_fraction": 0.10,
                "target_positive_wall_fraction_within_influence": 0.50,
                "target_correlation_length_um": 30.0,
                "random_seed": 7,
                "random_field_modes": 256,
            }
        )

        self.assertEqual(config.molecular_target_selection.default_mode, "automatic")
        self.assertEqual(
            config.molecular_target_selection
            .influence_region_endothelial_wall_area_fraction,
            0.10,
        )
        self.assertEqual(
            config.molecular_target_selection
            .target_positive_wall_fraction_within_influence,
            0.50,
        )
        self.assertEqual(config.molecular_target_selection.random_seed, 7)

        with self.assertRaisesRegex(ValueError, "requires"):
            self._load_config(
                molecular_target_selection={"default_mode": "automatic"}
            )
        with self.assertRaisesRegex(ValueError, "strictly between 0 and 1"):
            self._load_config(
                molecular_target_selection={
                    "default_mode": "automatic",
                    "influence_region_endothelial_wall_area_fraction": 1.0,
                    "target_positive_wall_fraction_within_influence": 0.5,
                    "target_correlation_length_um": 30.0,
                }
            )
    def test_mask_target_and_dimensioned_binding_parameters_are_loaded(self) -> None:
        config = self._load_config(
            molecular_target={
                "enabled": True,
                "region_mode": "MASK_NPZ",
                "mask_npz_path": "selected_target.npz",
                "target_density_molecules_per_m2": 2.7e14,
            },
            molecular_binding={
                "enabled": True,
                "model": "DETERMINISTIC_MEAN_FIELD",
                "ligand_density_molecules_per_m2": 1.0e14,
                "capture_distance_um": 0.04,
                # Capture distance and rest length represent different physics;
                # their relative order is intentionally not hard constrained.
                "rest_length_um": 0.05,
                "association_rate_m2_per_molecule_s": 2.0e-16,
                "zero_force_dissociation_rate_s": 0.06,
                "bond_stiffness_pn_per_um": 10_000.0,
                "reactive_compliance_nm": 0.039,
                "temperature_k": 310.0,
                "mean_field_warning_count": 10.0,
                "bell_exponent_limit": 75.0,
            },
        )

        self.assertEqual(config.molecular_target.region_mode, "mask_npz")
        self.assertEqual(config.molecular_target.target_density_molecules_per_m2, 2.7e14)
        self.assertTrue(config.molecular_binding.enabled)
        self.assertEqual(config.molecular_binding.capture_distance_um, 0.04)
        self.assertEqual(config.molecular_binding.rest_length_um, 0.05)

    def test_mask_path_is_resolved_but_existence_is_deferred_to_field_construction(self) -> None:
        config, config_path = self._load_config(
            molecular_target={
                "enabled": True,
                "region_mode": "mask_npz",
                "mask_npz_path": "target_inputs/not_created_yet.npz",
                "mask_array_key": "roi",
                "x_coordinates_key": "physical_x_um",
                "z_coordinates_key": "physical_z_um",
                "target_density_molecules_per_m2": 0.0,
            },
            return_source_path=True,
        )

        self.assertEqual(
            config.molecular_target.mask_npz_path,
            (config_path.parent / "target_inputs/not_created_yet.npz").resolve(),
        )
        self.assertEqual(config.molecular_target.mask_array_key, "roi")

    def test_binding_requires_a_real_target_region(self) -> None:
        with self.assertRaisesRegex(ValueError, "requires molecular_target.enabled"):
            self._load_config(
                molecular_binding={
                    "enabled": True,
                    "capture_distance_um": 0.04,
                    "rest_length_um": 0.03,
                    "bond_stiffness_pn_per_um": 10_000.0,
                    "reactive_compliance_nm": 0.03,
                }
            )

    def test_invalid_target_inputs_are_rejected(self) -> None:
        invalid_targets = (
            ({"enabled": True, "region_mode": "none"}, "disease ROI"),
            ({"enabled": True, "region_mode": "mask_npz"}, "mask_npz_path"),
            (
                {
                    "enabled": True,
                    "region_mode": "physical_polygon",
                    "polygon_vertices_um": [[0.0, 0.0], [1.0, 0.0]],
                },
                "none.*mask_npz",
            ),
        )
        for target, message in invalid_targets:
            with self.subTest(target=target):
                with self.assertRaisesRegex(ValueError, message):
                    self._load_config(molecular_target=target)

    def test_dimensionally_ambiguous_paper_fields_are_rejected(self) -> None:
        for field in ("paper_forward_rate_s_inv", "paper_encounter_radius_nm"):
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, "dimensionally incompatible"):
                    self._load_config(molecular_binding={field: 1.0})

    def test_enabled_binding_rejects_nonphysical_lengths_and_stiffness(self) -> None:
        target = self._valid_target()
        base_binding = {
            "enabled": True,
            "ligand_density_molecules_per_m2": 1.0e14,
            "capture_distance_um": 0.04,
            "rest_length_um": 0.03,
            "association_rate_m2_per_molecule_s": 1.0e-16,
            "zero_force_dissociation_rate_s": 1.0,
            "bond_stiffness_pn_per_um": 10_000.0,
            "reactive_compliance_nm": 0.03,
        }
        for field in (
            "capture_distance_um",
            "rest_length_um",
            "bond_stiffness_pn_per_um",
            "reactive_compliance_nm",
        ):
            binding = dict(base_binding)
            binding[field] = 0.0
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, field):
                    self._load_config(
                        molecular_target=target,
                        molecular_binding=binding,
                    )

    def test_dimensionless_sweep_requires_complete_positive_axes(self) -> None:
        target = self._valid_target()
        valid_sweep = {
            "enabled": True,
            "da_on_reference_time_s": 0.25,
            "da_on_levels": [0.01, 1.0, 100.0],
            "capture_distance_to_rest_length_ratios": [1.0, 1.5, 2.0],
            "target_density_molecules_per_um2_levels": [27.0, 270.0],
            "ligand_density_molecules_per_um2_levels": [10.0, 100.0],
        }
        config = self._load_config(
            molecular_target=target,
            molecular_binding={"rest_length_um": 0.03},
            binding_scenario_sweep=valid_sweep,
        )
        self.assertEqual(config.binding_scenario_sweep.da_on_levels, (0.01, 1.0, 100.0))
        self.assertEqual(config.binding_scenario_sweep.da_on_reference_time_s, 0.25)

        invalid_sweep = dict(valid_sweep)
        invalid_sweep["da_on_levels"] = [0.1, 0.0]
        with self.assertRaisesRegex(ValueError, "positive values"):
            self._load_config(
                molecular_target=target,
                molecular_binding={"rest_length_um": 0.03},
                binding_scenario_sweep=invalid_sweep,
            )

        with self.assertRaisesRegex(ValueError, "store_full_diagnostics"):
            self._load_config(
                molecular_target=target,
                molecular_binding={"rest_length_um": 0.03},
                binding_scenario_sweep=valid_sweep,
                particle_dynamics={"store_full_diagnostics": False},
            )

        missing_reference_time = dict(valid_sweep)
        missing_reference_time.pop("da_on_reference_time_s")
        with self.assertRaisesRegex(ValueError, "fixed before the transport run"):
            self._load_config(
                molecular_target=target,
                molecular_binding={"rest_length_um": 0.03},
                binding_scenario_sweep=missing_reference_time,
            )

        invalid_reference_time = dict(valid_sweep)
        invalid_reference_time["da_on_reference_time_s"] = float("nan")
        with self.assertRaisesRegex(ValueError, "finite and non-negative"):
            self._load_config(
                molecular_target=target,
                molecular_binding={"rest_length_um": 0.03},
                binding_scenario_sweep=invalid_reference_time,
            )

    @staticmethod
    def _valid_target() -> dict[str, object]:
        return {
            "enabled": True,
            "region_mode": "mask_npz",
            "mask_npz_path": "selected_target.npz",
            "target_density_molecules_per_m2": 2.7e14,
        }

    @staticmethod
    def _load_config(
        *,
        random_seed: int | None = None,
        molecular_target: dict[str, object] | None = None,
        molecular_target_selection: dict[str, object] | None = None,
        molecular_binding: dict[str, object] | None = None,
        binding_scenario_sweep: dict[str, object] | None = None,
        particle_dynamics: dict[str, object] | None = None,
        return_source_path: bool = False,
    ):
        document: dict[str, object] = {
            "input": {"model_dir": "."},
            "output": {
                "results_dir": "results",
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
        if random_seed is not None:
            document["random_seed"] = random_seed
        if molecular_target is not None:
            document["molecular_target"] = molecular_target
        if molecular_target_selection is not None:
            document["molecular_target_selection"] = molecular_target_selection
        if molecular_binding is not None:
            document["molecular_binding"] = molecular_binding
        if binding_scenario_sweep is not None:
            document["binding_scenario_sweep"] = binding_scenario_sweep
        if particle_dynamics is not None:
            document["particle_dynamics"] = particle_dynamics

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.yaml"
            path.write_text(yaml.safe_dump(document), encoding="utf-8")
            config = load_config(path)
            if return_source_path:
                return config, path
            return config


if __name__ == "__main__":
    unittest.main()
