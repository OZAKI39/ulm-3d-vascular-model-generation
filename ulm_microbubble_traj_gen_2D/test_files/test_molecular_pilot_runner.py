from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import yaml

from ulm_microbubble_traj_gen_2D.utils.core.config import (
    BindingScenarioSweepConfig,
    MolecularBindingConfig,
)
from ulm_microbubble_traj_gen_2D.utils.workflows.molecular_pilot_runner import (
    _numerical_wall_lock_mask,
    run_configured_contact_pilot,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_binding import (
    reaction_disk_radius_um,
    surface_slip_velocity_um_s,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_contact_pilot import (
    analyze_molecular_contact_pilot,
    contact_pilot_report_mapping,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_trajectory_kinematics import (
    realized_center_velocities_um_s,
)
from ulm_microbubble_traj_gen_2D.utils.core.types import ParticleTrajectories


class _FullyPositiveTarget:
    def __init__(self) -> None:
        self.last_points_xz_um = None
        self.point_batches: list[np.ndarray] = []

    def reaction_area_um2(
        self,
        points_xz_um,
        tangents_xz,
        reaction_radius_um,
    ):
        self.last_points_xz_um = np.asarray(points_xz_um, dtype=np.float64).copy()
        self.point_batches.append(self.last_points_xz_um)
        del tangents_xz
        radius = np.asarray(reaction_radius_um, dtype=np.float64)
        return np.pi * radius * radius


class MolecularPilotRunnerTests(unittest.TestCase):
    def test_matching_capture_distance_reports_internal_substep_authority(self) -> None:
        trajectories = replace(
            _two_frame_trajectory(),
            registry_target_exposure_time_s=np.asarray([0.003]),
            registry_target_exposure_event_count=np.asarray([2], dtype=np.int32),
            registry_target_reaction_area_time_um2_s=np.asarray([0.004]),
            registry_target_exposure_right_censored=np.asarray([True]),
            registry_target_exposure_quantitative_applicability_fraction=np.asarray(
                [0.75]
            ),
        )
        binding_cfg = MolecularBindingConfig(
            enabled=True,
            capture_distance_um=0.2,
            rest_length_um=0.1,
        )
        sweep_cfg = BindingScenarioSweepConfig(
            enabled=True,
            da_on_reference_time_s=0.25,
            da_on_levels=(1.0,),
            capture_distance_to_rest_length_ratios=(2.0,),
            target_density_molecules_per_um2_levels=(100.0,),
            ligand_density_molecules_per_um2_levels=(10.0,),
        )

        with tempfile.TemporaryDirectory() as directory:
            path = run_configured_contact_pilot(
                Path(directory),
                trajectories,
                _FullyPositiveTarget(),
                binding_cfg,
                sweep_cfg,
            )[0]
            report = yaml.safe_load(path.read_text(encoding="utf-8"))

        context = report["study_context"]
        self.assertEqual(context["exposure_authority"], "internal_physical_substep_registry")
        internal = context["internal_physical_substep_exposure"]
        self.assertEqual(internal["encountered_bubble_count_N_enc"], 1)
        self.assertEqual(internal["independent_exposure_event_count"], 2)
        self.assertEqual(internal["total_exposure_time_s"], 0.003)
        self.assertEqual(internal["total_reaction_area_time_E_T_um2_s"], 0.004)
        self.assertAlmostEqual(
            internal["exposure_weighted_rbc_quantitative_applicability_fraction"],
            0.75,
        )

    def test_lock_detection_compares_accepted_velocity_with_realized_motion(self) -> None:
        trajectories = replace(
            _two_frame_trajectory(),
            velocities_um_s=np.asarray(
                [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                dtype=np.float32,
            ),
            contact_constraint_active=np.asarray([True, True]),
        )
        realized = np.zeros((2, 3), dtype=np.float64)

        locked, source = _numerical_wall_lock_mask(
            trajectories,
            realized,
        )

        np.testing.assert_array_equal(locked, np.asarray([True, False]))
        self.assertEqual(
            source,
            "v16_accepted_internal_velocity_vs_same_id_realized_velocity",
        )

    def test_lock_detection_does_not_branch_on_metadata_version(self) -> None:
        trajectories = replace(
            _two_frame_trajectory(),
            metadata={"output_dt_s": 0.01},
            contact_constraint_active=np.asarray([True, True]),
        )

        _, source = _numerical_wall_lock_mask(
            trajectories,
            np.zeros((2, 3), dtype=np.float64),
        )

        self.assertEqual(
            source,
            "v16_accepted_internal_velocity_vs_same_id_realized_velocity",
        )

    def test_each_ratio_reports_exposure_but_uses_the_same_predeclared_rate(self) -> None:
        trajectories = _two_frame_trajectory()
        binding_cfg = MolecularBindingConfig(
            enabled=False,
            rest_length_um=0.1,
        )
        sweep_cfg = BindingScenarioSweepConfig(
            enabled=True,
            da_on_reference_time_s=0.25,
            da_on_levels=(1.0,),
            capture_distance_to_rest_length_ratios=(0.5, 2.0),
            target_density_molecules_per_um2_levels=(100.0,),
            ligand_density_molecules_per_um2_levels=(10.0,),
        )

        target = _FullyPositiveTarget()
        original_lexsort = np.lexsort
        with patch(
            "ulm_microbubble_traj_gen.utils.molecular.molecular_contact_pilot.np.lexsort",
            wraps=original_lexsort,
        ) as lexsort, tempfile.TemporaryDirectory() as directory:
            paths = run_configured_contact_pilot(
                Path(directory), trajectories, target, binding_cfg, sweep_cfg
            )
            self.assertEqual(len(paths), 2)
            first = yaml.safe_load(paths[0].read_text(encoding="utf-8"))
            second = yaml.safe_load(paths[1].read_text(encoding="utf-8"))

        self.assertEqual(first["pilot"]["total_contact_time_s"], 0.0)
        self.assertEqual(len(first["da_on_scenarios"]), 1)
        self.assertEqual(
            first["study_context"]["capture_distance_to_rest_length_ratio"],
            0.5,
        )
        self.assertEqual(first["study_context"]["records_without_unique_wall_normal"], 1)
        self.assertIn(
            "not automatically rerun",
            first["study_context"]["scenario_execution"],
        )
        self.assertGreater(second["pilot"]["total_contact_time_s"], 0.0)
        self.assertEqual(len(second["da_on_scenarios"]), 1)
        self.assertEqual(
            first["da_on_scenarios"][0]["association_rate_m2_per_molecule_s"],
            second["da_on_scenarios"][0]["association_rate_m2_per_molecule_s"],
        )
        self.assertEqual(
            first["da_on_scenarios"][0]["da_on_reference_time_s"],
            0.25,
        )
        self.assertIn(
            "predeclared_dimensionless_scenario",
            first["study_context"]["association_rate_parameter_source"],
        )
        np.testing.assert_allclose(
            target.last_points_xz_um,
            trajectories.positions_um[[1]][:, (0, 2)],
        )
        self.assertEqual([batch.shape[0] for batch in target.point_batches], [1])
        self.assertEqual(
            lexsort.call_count,
            1,
            "The ratio sweep must sort ratio-independent trajectory records once.",
        )

        # Reconstruct the former dense fully-positive target calculation and
        # require every reported exposure field to remain exactly identical.
        normals = np.asarray(trajectories.wall_normal_xz, dtype=np.float64)
        normal_norm = np.linalg.norm(normals, axis=1)
        valid_normals = np.isfinite(normal_norm) & (normal_norm > 1.0e-12)
        safe_normals = np.zeros_like(normals)
        safe_normals[valid_normals] = (
            normals[valid_normals] / normal_norm[valid_normals, None]
        )
        safe_normals[~valid_normals, 1] = 1.0
        tangents = np.column_stack((-safe_normals[:, 1], safe_normals[:, 0]))
        realized = realized_center_velocities_um_s(
            trajectories.frame_offsets,
            trajectories.bubble_id,
            trajectories.positions_um,
            0.01,
        )
        slip = surface_slip_velocity_um_s(
            np.sum(realized[:, (0, 2)] * tangents, axis=1),
            0.5 * np.asarray(trajectories.diameter_um, dtype=np.float64),
            np.asarray(trajectories.angular_velocity_rad_s, dtype=np.float64),
        )
        wall_lock, _ = _numerical_wall_lock_mask(trajectories, realized)
        for ratio, actual in zip((0.5, 2.0), (first, second), strict=True):
            reaction_radius = reaction_disk_radius_um(
                0.5 * np.asarray(trajectories.diameter_um, dtype=np.float64),
                np.asarray(trajectories.wall_gap_um, dtype=np.float64),
                ratio * binding_cfg.rest_length_um,
            )
            reaction_radius[~valid_normals] = 0.0
            dense_reference = analyze_molecular_contact_pilot(
                trajectories.frame_offsets,
                trajectories.bubble_id,
                np.pi * reaction_radius * reaction_radius,
                slip,
                0.01,
                numerical_wall_lock=wall_lock,
            )
            reference = contact_pilot_report_mapping(dense_reference)
            self.assertEqual(actual["pilot"], reference["pilot"])
            self.assertEqual(
                actual["per_bubble_contact"], reference["per_bubble_contact"]
            )
            self.assertEqual(actual["contact_events"], reference["contact_events"])


def _two_frame_trajectory() -> ParticleTrajectories:
    return ParticleTrajectories(
        frame_offsets=np.asarray([0, 1, 2], dtype=np.int64),
        bubble_id=np.asarray([0, 0], dtype=np.int64),
        positions_um=np.asarray([[1.0, 0.0, 1.1], [1.1, 0.0, 1.1]], dtype=np.float32),
        velocities_um_s=np.asarray([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32),
        wall_shear_stress_pa=np.zeros(2, dtype=np.float32),
        vessel_id=np.ones(2, dtype=np.int32),
        active=np.ones(2, dtype=bool),
        diameter_um=np.full(2, 2.0, dtype=np.float32),
        wall_gap_um=np.full(2, 0.1, dtype=np.float32),
        wall_contact=np.zeros(2, dtype=bool),
        wall_normal_xz=np.asarray([[0.0, 0.0], [0.0, 1.0]], dtype=np.float32),
        registry_bubble_id=np.asarray([0], dtype=np.int64),
        registry_diameter_um=np.asarray([2.0], dtype=np.float32),
        birth_frame=np.asarray([0], dtype=np.int32),
        death_frame=np.asarray([-1], dtype=np.int32),
        termination_reason=np.asarray([0], dtype=np.uint8),
        active_count_per_frame=np.asarray([1, 1], dtype=np.int32),
        injected_count_per_frame=np.asarray([1, 0], dtype=np.int32),
        terminated_count_per_frame=np.asarray([0, 0], dtype=np.int32),
        metadata={
            "output_dt_s": 0.01,
            "wall_contact_integrator": "revised_v16_continuous_geometry_predictive_mobility_unilateral_single_wall",
        },
        angular_velocity_rad_s=np.zeros(2, dtype=np.float32),
    )


if __name__ == "__main__":
    unittest.main()
