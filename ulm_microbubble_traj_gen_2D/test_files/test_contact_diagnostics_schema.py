from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.io.field_io import save_trajectories_npz
from ulm_microbubble_traj_gen_2D.utils.core.types import ParticleTrajectories
from ulm_microbubble_traj_gen_2D.utils.visualization.results.result_loader import _trajectory_arrays


class ContactDiagnosticsSchemaTests(unittest.TestCase):
    def test_new_writer_saves_contact_constraint_diagnostics_only(self) -> None:
        trajectories = _trajectory(
            metadata={
                "dt_s": 0.1,
                "output_dt_s": 0.1,
                "contact_geometry_tolerance_um": 0.001,
                "wall_contact_integrator": "revised_v16_continuous_geometry_predictive_mobility_unilateral_single_wall",
            }
        )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.npz"
            save_trajectories_npz(path, trajectories)
            with np.load(path, allow_pickle=False) as saved:
                self.assertIn("record_contact_constraint_active", saved.files)
                self.assertIn("record_contact_reaction_force_pn", saved.files)
                self.assertIn(
                    "record_contact_free_normal_velocity_um_s",
                    saved.files,
                )
                self.assertIn(
                    "record_contact_constrained_normal_velocity_um_s",
                    saved.files,
                )
                self.assertNotIn("record_wall_slide_failure_duration_s", saved.files)
                display = _trajectory_arrays(saved)

        self.assertTrue(display["strict_nonnegative_wall_gap"])
        self.assertEqual(display["wall_gap_invalid_below_um"], 0.0)
        self.assertAlmostEqual(display["contact_geometry_tolerance_um"], 0.001)
        np.testing.assert_array_equal(
            display["contact_constraint_active"],
            np.asarray([[True], [False]]),
        )
        np.testing.assert_allclose(
            display["contact_reaction_force_pn"],
            np.asarray([[2.0], [0.0]]),
        )
        np.testing.assert_allclose(
            display["contact_free_normal_velocity_um_s"],
            np.asarray([[-3.0], [1.0]]),
        )
        np.testing.assert_allclose(
            display["contact_constrained_normal_velocity_um_s"],
            np.asarray([[0.0], [1.0]]),
        )
        self.assertEqual(
            display["numerical_wall_lock_source"],
            "v16_accepted_internal_velocity_vs_same_id_realized_velocity",
        )

    def test_constraint_records_do_not_branch_on_metadata_version(self) -> None:
        trajectories = _trajectory(
            metadata={
                "dt_s": 0.1,
                "output_dt_s": 0.1,
                "contact_geometry_tolerance_um": 0.001,
            }
        )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "legacy_constraint_trajectory.npz"
            save_trajectories_npz(path, trajectories)
            with np.load(path, allow_pickle=False) as saved:
                display = _trajectory_arrays(saved)

        self.assertEqual(
            display["numerical_wall_lock_source"],
            "v16_accepted_internal_velocity_vs_same_id_realized_velocity",
        )

    def test_legacy_reader_retains_old_penetration_warning_threshold(self) -> None:
        trajectories = _trajectory(
            metadata={
                "dt_s": 0.1,
                "output_dt_s": 0.1,
                "penetration_tolerance_um": 0.05,
            },
            include_contact_diagnostics=False,
        )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "legacy_trajectory.npz"
            save_trajectories_npz(path, trajectories)
            with np.load(path, allow_pickle=False) as saved:
                payload = {key: np.asarray(saved[key]) for key in saved.files}
            payload["record_wall_slide_failure_duration_s"] = np.asarray([0.1, 0.0])
            np.savez_compressed(path, **payload)
            with np.load(path, allow_pickle=False) as saved:
                display = _trajectory_arrays(saved)

        self.assertFalse(display["strict_nonnegative_wall_gap"])
        self.assertAlmostEqual(display["wall_gap_invalid_below_um"], -0.05)
        self.assertAlmostEqual(display["contact_geometry_tolerance_um"], 0.05)
        np.testing.assert_array_equal(
            display["contact_constraint_active"],
            np.zeros((2, 1), dtype=bool),
        )
        self.assertTrue(
            np.all(np.isnan(display["contact_free_normal_velocity_um_s"]))
        )
        np.testing.assert_array_equal(
            display["numerical_wall_lock"],
            np.asarray([[True], [False]]),
        )
        self.assertEqual(
            display["numerical_wall_lock_source"],
            "legacy_wall_slide_failure_duration",
        )


def _trajectory(
    *,
    metadata: dict[str, object],
    include_contact_diagnostics: bool = True,
) -> ParticleTrajectories:
    return ParticleTrajectories(
        frame_offsets=np.asarray([0, 1, 2], dtype=np.int64),
        bubble_id=np.asarray([0, 0], dtype=np.int64),
        positions_um=np.asarray([[1.0, 0.0, 1.0], [1.1, 0.0, 1.0]]),
        velocities_um_s=np.asarray([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        wall_shear_stress_pa=np.zeros(2),
        vessel_id=np.ones(2, dtype=np.int32),
        active=np.ones(2, dtype=bool),
        diameter_um=np.full(2, 2.0),
        wall_gap_um=np.asarray([0.0, 0.2]),
        wall_contact=np.asarray([True, False]),
        wall_normal_xz=np.asarray([[0.0, 1.0], [0.0, 1.0]]),
        registry_bubble_id=np.asarray([0], dtype=np.int64),
        registry_diameter_um=np.asarray([2.0]),
        birth_frame=np.asarray([0], dtype=np.int32),
        death_frame=np.asarray([-1], dtype=np.int32),
        termination_reason=np.asarray([0], dtype=np.uint8),
        active_count_per_frame=np.asarray([1, 1], dtype=np.int32),
        injected_count_per_frame=np.asarray([1, 0], dtype=np.int32),
        terminated_count_per_frame=np.asarray([0, 0], dtype=np.int32),
        metadata=metadata,
        contact_constraint_active=(
            np.asarray([True, False]) if include_contact_diagnostics else None
        ),
        contact_reaction_force_pn=(
            np.asarray([2.0, 0.0]) if include_contact_diagnostics else None
        ),
        contact_free_normal_velocity_um_s=(
            np.asarray([-3.0, 1.0]) if include_contact_diagnostics else None
        ),
        contact_constrained_normal_velocity_um_s=(
            np.asarray([0.0, 1.0]) if include_contact_diagnostics else None
        ),
    )


if __name__ == "__main__":
    unittest.main()
