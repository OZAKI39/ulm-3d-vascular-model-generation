"""Dependency-light tests for the optional boundary-fitted flow backend."""

from __future__ import annotations

import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from ulm_vascular_model_generator.utils.core.models import Vessel

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from ulm_microbubble_traj_gen_2D.utils.core.config import (
    DomainConfig,
    FieldConfig,
    _validate_field_config,
    load_config,
)
from ulm_microbubble_traj_gen_2D.utils.flow.dolfinx_gmsh_solver import (
    INLET_PHYSICAL_TAG_OFFSET,
    OUTLET_PHYSICAL_TAG_OFFSET,
    WALL_PHYSICAL_TAG,
    _export_finite_element_velocity,
    _local_shear_stress_on_grid,
    _maximum_labeled_flux_relative_error,
    build_boundary_mesh_blueprint,
    solve_dolfinx_stokes_gmsh_2d,
)
from ulm_microbubble_traj_gen_2D.utils.flow.hybrid_velocity import (
    sample_finite_element_velocity,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.grid_domain import (
    build_domain_from_vessels,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.vessel_rasterizer import (
    rasterize_vessels,
)
from ulm_microbubble_traj_gen_2D.utils.io.field_io import (
    FieldReuseValidationError,
    _validate_saved_solver_metadata,
)


class DolfinxGmshBackendTests(unittest.TestCase):
    def test_removed_solver_options_are_rejected_at_config_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            config_path = Path(directory) / "legacy.yaml"
            config_path.write_text(
                "field:\n  solver_mode: phiflow_viscous_fv_projection_2d\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "Unsupported field options"):
                load_config(config_path)

    def test_boundary_blueprint_separates_wall_inlet_and_outlet(self) -> None:
        vessel, _, _, geometry, _ = _straight_case()
        blueprint = build_boundary_mesh_blueprint(geometry)
        tags = blueprint.segment_physical_tag

        self.assertEqual(blueprint.segment_count, geometry.full_boundary_segment_count)
        self.assertEqual(blueprint.points_xz_um.shape[1], 2)
        self.assertIn(WALL_PHYSICAL_TAG, tags)
        self.assertIn(INLET_PHYSICAL_TAG_OFFSET + 1, tags)
        self.assertIn(OUTLET_PHYSICAL_TAG_OFFSET + 1, tags)
        self.assertEqual(
            set(blueprint.physical_tag_to_section_index),
            {INLET_PHYSICAL_TAG_OFFSET + 1, OUTLET_PHYSICAL_TAG_OFFSET + 1},
        )

        segment_length = np.linalg.norm(
            geometry.full_boundary_end_xz_um - geometry.full_boundary_start_xz_um,
            axis=1,
        )
        expected_width = 2.0 * float(vessel.radius)
        for tag in (INLET_PHYSICAL_TAG_OFFSET + 1, OUTLET_PHYSICAL_TAG_OFFSET + 1):
            represented_width = float(np.sum(segment_length[tags == tag]))
            self.assertAlmostEqual(
                represented_width,
                expected_width,
                delta=max(0.05, 0.02 * expected_width),
            )

    def test_solver_requires_continuous_geometry(self) -> None:
        vessel, _, domain, _, raster = _straight_case()
        with self.assertRaisesRegex(ValueError, "ContinuousVesselGeometry"):
            solve_dolfinx_stokes_gmsh_2d(
                domain, raster, _field_config(), [vessel], None
            )

    def test_dolfinx_specific_configuration_is_validated(self) -> None:
        cfg = _field_config()
        for changes in (
            {"gmsh_bulk_mesh_size_um": -1.0},
            {"gmsh_wall_mesh_size_um": float("nan")},
            {"gmsh_element_order": 2},
            {"dolfinx_velocity_degree": 1},
            {"dolfinx_pressure_degree": 0},
            {"dolfinx_ksp_rtol": 0.0},
            {"blood_density_kg_m3": 0.0},
        ):
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError):
                    _validate_field_config(replace(cfg, **changes))

    def test_reuse_accepts_only_converged_nonempty_dolfinx_mesh(self) -> None:
        metadata = {
            "physical_converged": True,
            "solver_mode": "dolfinx_stokes_gmsh_2d",
            "physical_acceptance_schema": "dolfinx_fixed_flow_petsc_stokes_v3",
            "dolfinx_ksp_converged_reason": 2,
            "mesh_cell_count": 100,
        }
        _validate_saved_solver_metadata(metadata)
        _validate_saved_solver_metadata(
            {
                **metadata,
                "physical_acceptance_schema": "dolfinx_fixed_flow_scipy_stokes_v3",
                "linear_solver_converged_reason": 1,
                "linear_solver_relative_residual": 1.0e-12,
            }
        )
        for changes in (
            {"dolfinx_ksp_converged_reason": -3},
            {"mesh_cell_count": 0},
        ):
            with self.subTest(changes=changes):
                with self.assertRaises(FieldReuseValidationError):
                    _validate_saved_solver_metadata({**metadata, **changes})

    def test_labeled_flux_error_checks_each_terminal(self) -> None:
        target = np.asarray([0.0, 10.0, 10.0])
        actual = np.asarray([0.0, 11.0, 9.0])
        self.assertAlmostEqual(
            _maximum_labeled_flux_relative_error(actual, target),
            0.1,
        )

    def test_labeled_flux_error_rejects_missing_targets(self) -> None:
        self.assertEqual(
            _maximum_labeled_flux_relative_error(
                np.asarray([0.0]),
                np.asarray([0.0]),
            ),
            float("inf"),
        )

    def test_dolfinx_velocity_is_exported_without_using_dof_ordering(self) -> None:
        mesh = SimpleNamespace(
            geometry=SimpleNamespace(
                dofmap=np.asarray([[0, 1, 2]], dtype=np.int32),
                x=np.asarray(
                    [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]]
                ),
            )
        )

        class Velocity:
            @staticmethod
            def eval(points, cells):
                del cells
                x = points[:, 0]
                z = points[:, 1]
                return np.column_stack(
                    (1.0 + x + z + x * z, 2.0 * x - z * z)
                )

        field = _export_finite_element_velocity(
            Velocity(),
            mesh,
            degree=2,
            preferred_bin_size_um=1.0,
        )
        points = np.asarray([[0.5, 0.75]])
        velocity, gradient, cells = sample_finite_element_velocity(
            field,
            points,
            np.asarray([True]),
            use_numba=False,
        )
        np.testing.assert_allclose(
            velocity, [[1.0 + 0.5 + 0.75 + 0.5 * 0.75, 1.0 - 0.75**2]]
        )
        np.testing.assert_allclose(
            gradient, [[[1.75, 1.5], [2.0, -1.5]]]
        )
        np.testing.assert_array_equal(cells, [0])

    def test_local_shear_matches_simple_shear_and_rejects_rigid_rotation(self) -> None:
        _, _, domain, _, raster = _straight_case()
        gradient = np.zeros((*domain.shape, 2, 2), dtype=np.float64)

        # u_x = gamma*z gives tau_local = mu*|gamma|.
        gradient[..., 0, 1] = 10.0
        simple_shear = _local_shear_stress_on_grid(
            gradient=gradient,
            raster=raster,
            domain=domain,
        )
        np.testing.assert_allclose(
            simple_shear[raster.lumen_mask],
            3.0e-3 * 10.0,
            rtol=1.0e-12,
            atol=1.0e-12,
        )

        # Rigid-body rotation has an antisymmetric gradient and no viscous strain.
        gradient.fill(0.0)
        gradient[..., 0, 1] = -7.0
        gradient[..., 1, 0] = 7.0
        rigid_rotation = _local_shear_stress_on_grid(
            gradient=gradient,
            raster=raster,
            domain=domain,
        )
        np.testing.assert_allclose(rigid_rotation[raster.lumen_mask], 0.0)


def _straight_case():
    vessel = Vessel(
        vid=0,
        parent_id=-1,
        children=[],
        x_p=np.asarray([0.0, 0.0, 0.0]),
        x_d=np.asarray([20.0, 0.0, 0.0]),
        radius=2.0,
        flow_rate=40.0,
    )
    domain_cfg = DomainConfig(
        grid_spacing_um=0.5,
        padding_um=2.0,
        min_lumen_radius_cells=0.5,
        min_resolved_diameter_cells=4.0,
        max_grid_cells=100_000,
    )
    domain = build_domain_from_vessels([vessel], domain_cfg)
    geometry = build_continuous_vessel_geometry(
        [vessel],
        domain,
        maximum_boundary_element_length_um=0.25,
    )
    raster = rasterize_vessels(
        [vessel],
        domain,
        domain_cfg,
        effective_thickness_um=10.0,
        continuous_geometry=geometry,
    )
    return vessel, domain_cfg, domain, geometry, raster


def _field_config() -> FieldConfig:
    return FieldConfig(
        effective_thickness_um=10.0,
        boundary_depth_cells=1.5,
        flux_tolerance=1.0e-3,
        kinematic_viscosity_um2_s=3.0e6,
        hybrid_finite_element_distance_um=2.0,
        hybrid_transition_width_um=2.0,
        gmsh_bulk_mesh_size_um=1.0,
        gmsh_wall_mesh_size_um=0.5,
        gmsh_wall_refinement_distance_um=3.0,
        gmsh_element_order=1,
        dolfinx_velocity_degree=2,
        dolfinx_pressure_degree=1,
        dolfinx_ksp_rtol=1.0e-8,
    )


if __name__ == "__main__":
    unittest.main()
