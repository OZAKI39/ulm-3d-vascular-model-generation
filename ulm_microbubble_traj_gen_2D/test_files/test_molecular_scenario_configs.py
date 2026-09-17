from __future__ import annotations

import unittest
from pathlib import Path

import yaml


SCENARIO_DIR = (
    Path(__file__).resolve().parents[1] / "configs" / "molecular_binding_scenarios"
)


class MolecularScenarioConfigTests(unittest.TestCase):
    def test_no_experiment_sweep_uses_a_predeclared_reference_time(self) -> None:
        scenario = self._read("no_experiment_dimensionless_sweep.yaml")
        sweep = scenario["binding_scenario_sweep"]

        self.assertTrue(sweep["enabled"])
        self.assertEqual(sweep["da_on_reference_time_s"], 1.0)
        self.assertEqual(sweep["da_on_levels"], [0.01, 0.1, 1.0, 10.0, 100.0])
        self.assertEqual(
            sweep["capture_distance_to_rest_length_ratios"], [1.0, 1.5, 2.0]
        )
        self.assertEqual(
            scenario["molecular_binding"]["association_rate_m2_per_molecule_s"], 0.0
        )

    def test_literature_presets_cannot_silently_become_production_kon(self) -> None:
        expected = {
            "literature_slex_eselectin.yaml": (0.029, 0.06, 10_000.0, 0.039),
            "literature_psa_pselectin.yaml": (0.05, 0.000034, 50_000.0, 0.021),
        }
        for name, values in expected.items():
            with self.subTest(name=name):
                scenario = self._read(name)
                binding = scenario["molecular_binding"]
                cautions = " ".join(
                    scenario["scenario_metadata"]["paper_parameter_cautions"]
                )

                self.assertFalse(binding["enabled"])
                self.assertEqual(binding["association_rate_m2_per_molecule_s"], 0.0)
                self.assertEqual(binding["capture_distance_um"], 0.0)
                self.assertNotIn("forward_rate_s_inv", binding)
                self.assertNotIn("encounter_radius_nm", binding)
                self.assertIn("not the v7 effective 2D k_on", cautions)
                self.assertIn("not the v7 capture distance", cautions)
                self.assertEqual(
                    (
                        binding["rest_length_um"],
                        binding["zero_force_dissociation_rate_s"],
                        binding["bond_stiffness_pn_per_um"],
                        binding["reactive_compliance_nm"],
                    ),
                    values,
                )

    @staticmethod
    def _read(name: str) -> dict[str, object]:
        document = yaml.safe_load((SCENARIO_DIR / name).read_text(encoding="utf-8"))
        if not isinstance(document, dict):
            raise AssertionError(f"Scenario {name} must contain a YAML mapping.")
        return document


if __name__ == "__main__":
    unittest.main()
