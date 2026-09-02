from __future__ import annotations

import unittest
from dataclasses import replace

import numpy as np
from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.core.config import ParticleDynamicsConfig
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_binding import (
    MolecularBindingParameters,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_hydrodynamic_fields import (
    build_particle_hydrodynamic_fields,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_mobility_transport import (
    _EvaluationContext,
    _evaluate_rhs,
    ParticleWallGapInvariantError,
)
from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import straight_channel_case


class ParticleMobilityTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.domain, cls.raster, cls.flow = straight_channel_case()

    def test_disabled_wall_and_collision_terms_reproduce_background_velocity(self) -> None:
        context = _evaluation_context(
            self.domain,
            self.raster,
            self.flow,
            _dynamics_config(near_wall_enabled=False, collisions_enabled=False),
            radii_um=np.full(2, 0.5, dtype=np.float64),
        )
        evaluation = _evaluate_rhs(
            np.asarray([[8.0, 2.5], [9.0, 3.5]], dtype=np.float64),
            np.asarray([0, 1], dtype=np.int64),
            np.ones(2, dtype=bool),
            context,
        )

        np.testing.assert_allclose(
            evaluation.particle_velocity_xz_um_s,
            evaluation.fluid_velocity_xz_um_s,
            rtol=0.0,
            atol=1.0e-12,
        )

    def test_synchronous_collision_rhs_is_invariant_to_lane_permutation(self) -> None:
        context = _evaluation_context(
            self.domain,
            self.raster,
            self.flow,
            _dynamics_config(
                near_wall_enabled=False,
                collisions_enabled=True,
                collision_layer_um=0.1,
                collision_relaxation_time_s=0.01,
            ),
            radii_um=np.full(3, 0.5, dtype=np.float64),
        )
        positions = np.asarray(
            [[8.0, 3.0], [8.8, 3.0], [9.6, 3.0]],
            dtype=np.float64,
        )
        bubble_ids = np.asarray([0, 1, 2], dtype=np.int64)
        active = np.ones(3, dtype=bool)
        baseline = _evaluate_rhs(positions, bubble_ids, active, context)

        permutation = np.asarray([2, 0, 1], dtype=np.int64)
        permuted_ids = bubble_ids[permutation]
        permuted = _evaluate_rhs(
            positions[permutation],
            permuted_ids,
            active[permutation],
            context,
        )
        order_by_id = np.argsort(permuted_ids)

        np.testing.assert_allclose(
            permuted.particle_velocity_xz_um_s[order_by_id],
            baseline.particle_velocity_xz_um_s,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            permuted.collision_force_xz_pn[order_by_id],
            baseline.collision_force_xz_pn,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_array_equal(
            permuted.collision_neighbor_count[order_by_id],
            baseline.collision_neighbor_count,
        )
        np.testing.assert_allclose(
            np.sum(baseline.collision_force_xz_pn, axis=0),
            np.zeros(2),
            rtol=0.0,
            atol=1.0e-12,
        )

    def test_numba_batch_rhs_matches_python_reference(self) -> None:
        python_context = _evaluation_context(
            self.domain,
            self.raster,
            self.flow,
            _dynamics_config(
                near_wall_enabled=True,
                collisions_enabled=True,
                collision_layer_um=0.1,
                collision_relaxation_time_s=0.01,
                neighbor_search="all_pairs",
            ),
            radii_um=np.asarray([0.5, 0.6, 0.45], dtype=np.float64),
        )
        numba_context = replace(python_context, use_numba=True)
        for field_name in (
            "velocity_xz_um_s",
            "wall_shear_stress_pa",
            "local_vessel_radius_um",
            "velocity_gradient_s_inv",
            "dynamic_viscosity_pa_s",
        ):
            self.assertEqual(getattr(numba_context, field_name).dtype, np.float32)
        positions = np.asarray(
            [[8.0, 1.65], [8.8, 2.7], [9.6, 4.25]],
            dtype=np.float64,
        )
        bubble_ids = np.asarray([0, 1, 2], dtype=np.int64)
        active = np.ones(3, dtype=bool)

        reference = _evaluate_rhs(positions, bubble_ids, active, python_context)
        accelerated = _evaluate_rhs(positions, bubble_ids, active, numba_context)

        for field_name in (
            "particle_velocity_xz_um_s",
            "fluid_velocity_xz_um_s",
            "angular_velocity_rad_s",
            "collision_force_xz_pn",
            "wall_gap_um",
            "wall_normal_xz",
            "gap_ratio",
            "near_wall_weight",
            "generalized_mobility",
        ):
            np.testing.assert_allclose(
                getattr(accelerated, field_name),
                getattr(reference, field_name),
                rtol=2.0e-13,
                atol=2.0e-13,
                equal_nan=True,
            )
        np.testing.assert_array_equal(
            accelerated.collision_neighbor_count,
            reference.collision_neighbor_count,
        )
        np.testing.assert_array_equal(
            accelerated.two_wall_warning,
            reference.two_wall_warning,
        )
        np.testing.assert_array_equal(
            accelerated.wall_normal_valid,
            reference.wall_normal_valid,
        )
        self.assertAlmostEqual(
            accelerated.maximum_collision_speed_um_s,
            reference.maximum_collision_speed_um_s,
        )

    def test_rhs_wall_normal_comes_from_exact_authoritative_face(self) -> None:
        context = _evaluation_context(
            self.domain,
            self.raster,
            self.flow,
            _dynamics_config(near_wall_enabled=True, collisions_enabled=False),
            radii_um=np.asarray([0.5], dtype=np.float64),
        )
        position = np.asarray([[2.25, 1.4]], dtype=np.float64)

        evaluation = _evaluate_rhs(
            position,
            np.asarray([0], dtype=np.int64),
            np.asarray([True]),
            context,
        )
        expected_normal = context.boundary_geometry.exact_solid_wall_state_xz_um(
            context.boundary_geometry.grid_to_world_xz(position)
        ).inward_normal_xz[0]

        np.testing.assert_allclose(
            evaluation.wall_normal_xz[0],
            expected_normal,
            rtol=0.0,
            atol=1.0e-14,
        )

    def test_generalized_mobility_is_available_before_first_contact(self) -> None:
        positions = np.asarray([[8.0, 1.65], [8.0, 3.0]], dtype=np.float64)
        context = _evaluation_context(
            self.domain,
            self.raster,
            self.flow,
            _dynamics_config(near_wall_enabled=True, collisions_enabled=False),
            radii_um=np.asarray([0.5, 0.5], dtype=np.float64),
        )
        probe = _evaluate_rhs(
            positions,
            np.asarray([0, 1], dtype=np.int64),
            np.ones(2, dtype=bool),
            context,
        )
        tolerance = float(np.mean(probe.wall_gap_um))
        evaluation = _evaluate_rhs(
            positions,
            np.asarray([0, 1], dtype=np.int64),
            np.ones(2, dtype=bool),
            replace(context, contact_geometry_tolerance_um=tolerance),
        )
        contact = evaluation.wall_gap_um <= tolerance

        self.assertEqual(int(np.count_nonzero(contact)), 1)
        self.assertGreater(np.linalg.norm(evaluation.generalized_mobility[contact]), 0.0)
        self.assertGreater(np.linalg.norm(evaluation.generalized_mobility[~contact]), 0.0)

    def test_true_zero_gap_is_not_replaced_by_hydrodynamic_regularization(self) -> None:
        radius_um = 0.5
        context = _evaluation_context(
            self.domain,
            self.raster,
            self.flow,
            _dynamics_config(near_wall_enabled=True, collisions_enabled=False),
            radii_um=np.asarray([radius_um], dtype=np.float64),
        )
        evaluation = _evaluate_rhs(
            np.asarray([[8.0, 0.5]], dtype=np.float64),
            np.asarray([0], dtype=np.int64),
            np.asarray([True]),
            context,
        )

        self.assertEqual(float(evaluation.wall_gap_um[0]), 0.0)
        self.assertEqual(
            float(evaluation.hydrodynamic_gap_um[0]),
            radius_um * float(context.dynamics.xi_min),
        )

    def test_molecular_target_query_receives_particle_centre_and_exact_gap_radius(self) -> None:
        radius_um = 0.5
        capture_distance_um = 0.1
        context = _evaluation_context(
            self.domain,
            self.raster,
            self.flow,
            _dynamics_config(near_wall_enabled=True, collisions_enabled=False),
            radii_um=np.asarray([radius_um], dtype=np.float64),
        )
        target_spy = _TargetReactionAreaSpy()
        parameters = MolecularBindingParameters(
            ligand_density_molecules_m2=0.0,
            target_density_molecules_m2=0.0,
            association_rate_m2_per_molecule_s=0.0,
            zero_force_dissociation_rate_s=0.0,
            rest_length_um=0.0,
            spring_stiffness_pn_per_um=0.0,
            force_sensitivity_length_nm=0.0,
            temperature_k=310.0,
        )
        context = replace(
            context,
            origin_xz_um=(
                float(self.domain.origin_um[0]),
                float(self.domain.origin_um[2]),
            ),
            molecular_target_field=target_spy,
            molecular_binding_parameters=parameters,
            molecular_capture_distance_um=capture_distance_um,
        )
        position_grid = np.asarray([[8.0, 0.5]], dtype=np.float64)

        evaluation = _evaluate_rhs(
            position_grid,
            np.asarray([0], dtype=np.int64),
            np.asarray([True]),
            context,
        )

        expected_centre = np.asarray(
            [
                float(self.domain.origin_um[0])
                + float(self.domain.spacing_um) * position_grid[0, 0],
                float(self.domain.origin_um[2])
                + float(self.domain.spacing_um) * position_grid[0, 1],
            ]
        )
        expected_reaction_radius = np.sqrt(
            2.0 * radius_um * capture_distance_um - capture_distance_um**2
        )
        np.testing.assert_allclose(target_spy.points_xz_um, expected_centre[None, :])
        np.testing.assert_allclose(
            target_spy.reaction_radius_um,
            np.asarray([expected_reaction_radius]),
        )
        self.assertEqual(float(evaluation.wall_gap_um[0]), 0.0)

    def test_roundoff_negative_gap_is_normalized_before_molecular_capture(self) -> None:
        radius_um = 0.5
        capture_distance_um = 0.1
        context = _evaluation_context(
            self.domain,
            self.raster,
            self.flow,
            _dynamics_config(near_wall_enabled=True, collisions_enabled=False),
            radii_um=np.asarray([radius_um], dtype=np.float64),
        )
        target_spy = _TargetReactionAreaSpy()
        parameters = MolecularBindingParameters(
            ligand_density_molecules_m2=0.0,
            target_density_molecules_m2=0.0,
            association_rate_m2_per_molecule_s=0.0,
            zero_force_dissociation_rate_s=0.0,
            rest_length_um=0.0,
            spring_stiffness_pn_per_um=0.0,
            force_sensitivity_length_nm=0.0,
            temperature_k=310.0,
        )
        context = replace(
            context,
            origin_xz_um=(
                float(self.domain.origin_um[0]),
                float(self.domain.origin_um[2]),
            ),
            molecular_target_field=target_spy,
            molecular_binding_parameters=parameters,
            molecular_capture_distance_um=capture_distance_um,
        )

        evaluation = _evaluate_rhs(
            np.asarray([[8.0, radius_um - 1.0e-13]], dtype=np.float64),
            np.asarray([0], dtype=np.int64),
            np.asarray([True]),
            context,
            time_s=0.59,
        )

        self.assertEqual(float(evaluation.wall_gap_um[0]), 0.0)
        expected = np.sqrt(
            2.0 * radius_um * capture_distance_um
            - capture_distance_um**2
        )
        np.testing.assert_allclose(target_spy.reaction_radius_um, [expected])

    def test_material_wall_gap_reports_particle_identity_and_time(self) -> None:
        context = _evaluation_context(
            self.domain,
            self.raster,
            self.flow,
            _dynamics_config(near_wall_enabled=True, collisions_enabled=False),
            radii_um=np.asarray([0.5], dtype=np.float64),
        )

        with self.assertRaises(ParticleWallGapInvariantError) as captured:
            _evaluate_rhs(
                np.asarray([[8.0, 0.5 - 1.0e-6]], dtype=np.float64),
                np.asarray([0], dtype=np.int64),
                np.asarray([True]),
                context,
                time_s=0.59,
            )

        message = str(captured.exception)
        self.assertIn("permanent_microbubble_id=0", message)
        self.assertIn("physical_time_s=0.58999999999999997", message)
        self.assertIn("position_xz_um=", message)
        self.assertIn("wall_gap_um=", message)
        self.assertIn("geometry_roundoff_tolerance_um=", message)


def _dynamics_config(**overrides: object) -> ParticleDynamicsConfig:
    return replace(
        ParticleDynamicsConfig(
            time_integrator="euler",
            near_wall_enabled=True,
            collisions_enabled=True,
            collision_layer_um=0.05,
            collision_relaxation_time_s=0.01,
            neighbor_search="all_pairs",
            store_full_diagnostics=True,
        ),
        **overrides,
    )


def _evaluation_context(
    domain,
    raster,
    flow,
    dynamics_cfg: ParticleDynamicsConfig,
    *,
    radii_um: np.ndarray,
) -> _EvaluationContext:
    vessel = Vessel(vid=0, parent_id=-1, children=[])
    vessel.x_p = np.asarray([0.0, 0.0, 3.0])
    vessel.x_d = np.asarray([19.0, 0.0, 3.0])
    vessel.radius = 3.0
    geometry = build_continuous_vessel_geometry([vessel], domain)
    fields = build_particle_hydrodynamic_fields(
        domain, raster, flow, continuous_geometry=geometry
    )
    return _EvaluationContext(
        velocity_xz_um_s=np.ascontiguousarray(flow.velocity_xz_um_s),
        wall_shear_stress_pa=np.ascontiguousarray(flow.wall_shear_stress_pa),
        vessel_id=np.ascontiguousarray(raster.vessel_id, dtype=np.int32),
        local_vessel_radius_um=np.ascontiguousarray(raster.radius_um),
        velocity_gradient_s_inv=np.ascontiguousarray(fields.velocity_gradient_s_inv),
        dynamic_viscosity_pa_s=np.ascontiguousarray(fields.dynamic_viscosity_pa_s),
        all_bubble_radii_um=np.ascontiguousarray(radii_um, dtype=np.float64),
        spacing_um=float(domain.spacing_um),
        dynamics=dynamics_cfg,
        contact_geometry_tolerance_um=1.0e-6,
        use_numba=False,
        boundary_geometry=fields.boundary_geometry,
        hybrid_velocity=fields.hybrid_velocity,
    )


class _TargetReactionAreaSpy:
    def __init__(self) -> None:
        self.points_xz_um: np.ndarray | None = None
        self.reaction_radius_um: np.ndarray | None = None

    def reaction_area_um2(
        self,
        points_xz_um: np.ndarray,
        tangents_xz: np.ndarray,
        reaction_radius_um: np.ndarray,
    ) -> np.ndarray:
        self.points_xz_um = np.asarray(points_xz_um, dtype=np.float64).copy()
        self.reaction_radius_um = np.asarray(
            reaction_radius_um, dtype=np.float64
        ).copy()
        return np.zeros(self.points_xz_um.shape[0], dtype=np.float64)


if __name__ == "__main__":
    unittest.main()
