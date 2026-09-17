from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import yaml

from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_binding_scenarios import (
    MOLECULES_PER_UM2_TO_MOLECULES_PER_M2,
    build_da_on_scenarios,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_contact_pilot import (
    analyze_molecular_contact_pilot,
    contact_pilot_report_mapping,
    render_contact_pilot_yaml,
    save_contact_pilot_yaml,
)


class MolecularContactPilotTests(unittest.TestCase):
    def test_ragged_records_report_observed_exposure_without_kon_inference(self) -> None:
        pilot = self._contact_pilot()

        self.assertEqual(pilot.n_frames, 4)
        self.assertEqual(pilot.n_records, 10)
        self.assertEqual(pilot.n_unique_bubbles, 3)
        self.assertEqual(pilot.n_contacting_bubbles, 3)
        self.assertEqual(pilot.n_contact_records, 5)
        self.assertEqual(
            pilot.sampling_convention,
            "saved_frame_rectangle_each_positive_record_times_output_dt",
        )
        self.assertAlmostEqual(pilot.total_contact_time_s, 0.5)
        self.assertAlmostEqual(pilot.median_positive_bubble_contact_time_s, 0.2)
        self.assertAlmostEqual(
            pilot.contact_weighted_mean_reaction_area_um2,
            2.3,
        )
        self.assertAlmostEqual(
            pilot.contact_weighted_mean_abs_tangential_slip_um_s,
            6.0,
        )
        self.assertIn("not used to infer", pilot.exposure_summary.lower())
        self.assertEqual(pilot.right_censored_event_count, 1)
        self.assertAlmostEqual(pilot.right_censored_contact_time_s, 0.1)

        summaries = {summary.bubble_id: summary for summary in pilot.bubble_summaries}
        self.assertAlmostEqual(summaries[0].cumulative_contact_time_s, 0.2)
        self.assertAlmostEqual(summaries[1].cumulative_contact_time_s, 0.1)
        self.assertAlmostEqual(summaries[2].cumulative_contact_time_s, 0.2)
        self.assertEqual([event.bubble_id for event in pilot.events], [0, 2, 1])
        self.assertEqual(
            [(event.first_frame, event.last_frame) for event in pilot.events],
            [(1, 2), (1, 2), (3, 3)],
        )

    def test_nonconsecutive_positive_records_form_separate_events(self) -> None:
        pilot = analyze_molecular_contact_pilot(
            frame_offsets=np.array([0, 1, 2, 3], dtype=np.int64),
            bubble_id=np.array([7, 7, 7], dtype=np.int64),
            reaction_area_um2=np.array([1.0, 0.0, 2.0]),
            tangential_slip_um_s=np.array([3.0, 4.0, 5.0]),
            output_dt_s=0.01,
        )

        self.assertEqual(len(pilot.events), 2)
        self.assertEqual(pilot.bubble_summaries[0].event_count, 2)
        self.assertAlmostEqual(pilot.bubble_summaries[0].cumulative_contact_time_s, 0.02)

    def test_zero_exposure_still_allows_the_predeclared_scenario_table(self) -> None:
        pilot = analyze_molecular_contact_pilot(
            frame_offsets=np.array([0, 1, 2], dtype=np.int64),
            bubble_id=np.array([3, 3], dtype=np.int64),
            reaction_area_um2=np.zeros(2),
            tangential_slip_um_s=np.array([1.0, 2.0]),
            output_dt_s=0.001,
        )
        scenarios = build_da_on_scenarios(
            da_on_reference_time_s=0.25,
            da_on_levels=[1.0],
            target_density_molecules_per_um2_levels=[27.0],
        )

        self.assertEqual(pilot.total_contact_time_s, 0.0)
        self.assertIsNone(pilot.median_positive_bubble_contact_time_s)
        self.assertIn("scenario table remains unchanged", pilot.exposure_summary)
        self.assertEqual(len(scenarios), 1)
        report = contact_pilot_report_mapping(pilot, scenarios)
        self.assertEqual(report["pilot"]["total_contact_time_s"], 0.0)
        self.assertEqual(len(report["da_on_scenarios"]), 1)
        self.assertNotIn("association_rate_identifiable", report["pilot"])
        self.assertNotIn("contact_time_usable_for_kon", report["pilot"])

    def test_numerical_lock_marks_exposure_but_does_not_redefine_scenarios(self) -> None:
        pilot = analyze_molecular_contact_pilot(
            frame_offsets=np.array([0, 1, 2], dtype=np.int64),
            bubble_id=np.array([3, 3], dtype=np.int64),
            reaction_area_um2=np.array([1.0, 1.0]),
            tangential_slip_um_s=np.array([2.0, 0.0]),
            output_dt_s=0.01,
            numerical_wall_lock=np.array([False, True]),
        )
        scenarios = build_da_on_scenarios(
            da_on_reference_time_s=0.25,
            da_on_levels=[1.0],
            target_density_molecules_per_um2_levels=[27.0],
        )

        self.assertEqual(pilot.numerical_wall_lock_contact_records, 1)
        self.assertAlmostEqual(pilot.numerical_wall_lock_contact_time_s, 0.01)
        self.assertIn("numerically suspect", pilot.exposure_summary)
        self.assertEqual(len(scenarios), 1)

    def test_da_on_table_uses_fixed_reference_time_and_density_conversion(self) -> None:
        scenarios = build_da_on_scenarios(
            da_on_reference_time_s=0.25,
            da_on_levels=[0.1, 1.0],
            target_density_molecules_per_um2_levels=[27.0],
            ligand_density_molecules_per_um2_levels=[10.0, 100.0],
            capture_distance_to_rest_length_ratios=[1.0, 2.0],
            rest_length_um=0.05,
        )

        self.assertEqual(len(scenarios), 8)
        first = scenarios[0]
        expected_target_density_m2 = 27.0 * MOLECULES_PER_UM2_TO_MOLECULES_PER_M2
        self.assertEqual(first.target_density_molecules_per_m2, expected_target_density_m2)
        self.assertEqual(first.da_on_reference_time_s, 0.25)
        self.assertEqual(
            first.ligand_density_molecules_per_m2,
            10.0 * MOLECULES_PER_UM2_TO_MOLECULES_PER_M2,
        )
        self.assertAlmostEqual(first.capture_distance_um, 0.05)
        self.assertAlmostEqual(
            first.association_rate_m2_per_molecule_s,
            0.1 / (expected_target_density_m2 * 0.25),
        )

    def test_yaml_report_can_be_rendered_and_saved(self) -> None:
        pilot = self._contact_pilot()
        scenarios = build_da_on_scenarios(
            da_on_reference_time_s=1.0,
            da_on_levels=[1.0],
            target_density_molecules_per_um2_levels=[27.0],
        )
        study_context = {
            "capture_distance_to_rest_length_ratio": 1.5,
            "scenario_execution": "definitions_only",
        }

        rendered = render_contact_pilot_yaml(
            pilot,
            scenarios,
            study_context=study_context,
        )
        document = yaml.safe_load(rendered)
        self.assertEqual(document["da_on_scenarios"][0]["da_on_reference_time_s"], 1.0)
        self.assertEqual(document["study_context"], study_context)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nested" / "contact_pilot.yaml"
            saved_path = save_contact_pilot_yaml(
                path,
                pilot,
                scenarios,
                study_context=study_context,
            )
            self.assertEqual(saved_path, path.resolve())
            self.assertEqual(path.read_text(encoding="utf-8"), rendered)

    def test_invalid_scenario_reference_time_is_rejected(self) -> None:
        for value in (0.0, -1.0, float("nan"), True):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "finite and positive"):
                    build_da_on_scenarios(
                        da_on_reference_time_s=value,
                        da_on_levels=[1.0],
                        target_density_molecules_per_um2_levels=[27.0],
                    )

    def test_invalid_ragged_input_is_rejected(self) -> None:
        invalid_cases = (
            (
                np.array([0, 2], dtype=np.int64),
                np.array([1], dtype=np.int64),
                np.array([0.0]),
                "frame_offsets",
            ),
            (
                np.array([0, 2], dtype=np.int64),
                np.array([1, 1], dtype=np.int64),
                np.array([0.0, 0.0]),
                "at most once per frame",
            ),
            (
                np.array([0, 1], dtype=np.int64),
                np.array([1], dtype=np.int64),
                np.array([-1.0]),
                "non-negative",
            ),
        )
        for offsets, ids, area, message in invalid_cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    analyze_molecular_contact_pilot(
                        frame_offsets=offsets,
                        bubble_id=ids,
                        reaction_area_um2=area,
                        tangential_slip_um_s=np.zeros(ids.size),
                        output_dt_s=0.01,
                    )

    @staticmethod
    def _contact_pilot():
        return analyze_molecular_contact_pilot(
            frame_offsets=np.array([0, 2, 5, 7, 10], dtype=np.int64),
            bubble_id=np.array([0, 1, 0, 1, 2, 0, 2, 0, 1, 2], dtype=np.int64),
            reaction_area_um2=np.array([0.0, 0.0, 1.0, 0.0, 2.0, 1.5, 3.0, 0.0, 4.0, 0.0]),
            tangential_slip_um_s=np.array([0.0, 0.0, 4.0, 0.0, -6.0, -2.0, 8.0, 0.0, -10.0, 0.0]),
            output_dt_s=0.1,
        )


if __name__ == "__main__":
    unittest.main()
