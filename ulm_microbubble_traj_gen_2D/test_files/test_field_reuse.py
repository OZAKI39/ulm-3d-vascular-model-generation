from __future__ import annotations

import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import (
    straight_channel_case,
)
from ulm_microbubble_traj_gen_2D.utils.core.config import (
    DomainConfig,
    FieldConfig,
)
from ulm_microbubble_traj_gen_2D.utils.io.field_io import (
    FieldReuseValidationError,
    build_field_reuse_config_contract,
    load_reusable_flow_field,
    save_field_npz,
    save_run_config,
)

_FIELD_NAME = "velocity_and_wall_shear_field.npz"


def _continuous_geometry() -> SimpleNamespace:
    section_points = np.asarray([[0.0, 3.0], [19.0, 3.0]], dtype=np.float64)
    section_normals = np.asarray([[-1.0, 0.0], [1.0, 0.0]], dtype=np.float64)
    return SimpleNamespace(
        geometry_schema="v16_continuous_swept_vessel_boundary",
        geometry_hash_sha256="0123456789abcdef" * 4,
        solid_face_start_xz_um=np.asarray([[0.0, 0.5], [19.0, 5.5]], dtype=np.float64),
        solid_face_end_xz_um=np.asarray([[19.0, 0.5], [0.0, 5.5]], dtype=np.float64),
        solid_face_inward_normal_xz=np.asarray([[0.0, 1.0], [0.0, -1.0]], dtype=np.float64),
        solid_face_length_um=np.asarray([19.0, 19.0], dtype=np.float64),
        solid_face_ring_index=np.asarray([0, 0], dtype=np.int32),
        solid_face_arclength_start_um=np.asarray([0.0, 19.0], dtype=np.float64),
        solid_face_arclength_end_um=np.asarray([19.0, 38.0], dtype=np.float64),
        boundary_ring_length_um=38.0,
        open_section_point_xz_um=section_points,
        open_section_outward_normal_xz=section_normals,
        open_section_tangent_xz=np.asarray([[0.0, 1.0], [0.0, 1.0]], dtype=np.float64),
        open_section_half_width_um=np.asarray([2.5, 2.5], dtype=np.float64),
        open_section_label=np.asarray([1, 1], dtype=np.int32),
        open_section_kind=np.asarray([-1, 1], dtype=np.int8),
        open_section_vessel_id=np.asarray([0, 0], dtype=np.int32),
    )


def _complete_flow(domain, flow):
    nx, nz = domain.shape
    zeros = np.zeros(domain.shape, dtype=np.float32)
    zeros_xz = np.zeros((*domain.shape, 2), dtype=np.float32)
    open_face_cell_ij = np.asarray([[0, 3], [nx - 1, 3]], dtype=np.int32)
    open_face_index_ij = np.asarray([[0, 3], [nx, 3]], dtype=np.int32)
    open_normals = np.asarray([[-1.0, 0.0], [1.0, 0.0]], dtype=np.float32)
    section_points = np.asarray([[0.0, 3.0], [19.0, 3.0]], dtype=np.float64)
    return replace(
        flow,
        initial_velocity_xz_um_s=flow.velocity_xz_um_s.copy(),
        initial_speed_um_s=flow.speed_um_s.copy(),
        divergence_s_inv=zeros.copy(),
        wall_penetration_um_s=zeros.copy(),
        pressure=zeros.copy(),
        boundary_velocity_xz_um_s=zeros_xz.copy(),
        boundary_normal_xz=zeros_xz.copy(),
        boundary_weight=zeros.copy(),
        boundary_edge_length_um=zeros.copy(),
        open_boundary_flux_um2_s=zeros.copy(),
        face_flux_x_um2_s=np.zeros((nx + 1, nz), dtype=np.float32),
        face_flux_z_um2_s=np.zeros((nx, nz + 1), dtype=np.float32),
        inlet_target_by_label_um2_s=np.asarray([40.0], dtype=np.float64),
        outlet_target_by_label_um2_s=np.asarray([40.0], dtype=np.float64),
        inlet_actual_by_label_um2_s=np.asarray([40.0], dtype=np.float64),
        outlet_actual_by_label_um2_s=np.asarray([40.0], dtype=np.float64),
        open_face_cell_ij=open_face_cell_ij,
        open_face_index_ij=open_face_index_ij,
        open_face_axis=np.asarray([0, 0], dtype=np.int8),
        open_face_normal_xz=open_normals,
        open_face_center_xz_um=np.asarray([[-0.5, 3.0], [19.5, 3.0]], dtype=np.float64),
        open_face_length_um=np.asarray([1.0, 1.0], dtype=np.float64),
        open_face_label=np.asarray([1, 1], dtype=np.int32),
        open_face_kind=np.asarray([-1, 1], dtype=np.int8),
        open_section_point_xz_um=section_points,
        open_section_outward_normal_xz=open_normals.astype(np.float64),
        open_section_tangent_xz=np.asarray([[0.0, 1.0], [0.0, 1.0]], dtype=np.float64),
        open_section_half_width_um=np.asarray([2.5, 2.5], dtype=np.float64),
        open_section_label=np.asarray([1, 1], dtype=np.int32),
        open_section_kind=np.asarray([-1, 1], dtype=np.int8),
        solver_metadata={
            "solver_mode": "dolfinx_stokes_gmsh_2d",
            "physical_converged": True,
            "physical_acceptance_schema": "dolfinx_fixed_flow_scipy_stokes_v3",
            "linear_solver_converged_reason": 1,
            "mesh_cell_count": 100,
        },
    )


def _config(root: Path):
    return SimpleNamespace(
        model_dir=root.resolve(),
        domain=DomainConfig(
            grid_spacing_um=1.0,
            padding_um=2.0,
            min_lumen_radius_cells=0.5,
            min_resolved_diameter_cells=5.0,
            max_grid_cells=10_000,
        ),
        field=FieldConfig(
            effective_thickness_um=10.0,
            boundary_depth_cells=1.5,
            flux_tolerance=1.0e-3,
            kinematic_viscosity_um2_s=3.0e6,
            hybrid_finite_element_distance_um=1.0,
            hybrid_transition_width_um=1.0,
        ),
        quick_test=False,
    )


def _write_reusable_result(root: Path):
    domain, raster, basic_flow = straight_channel_case()
    flow = _complete_flow(domain, basic_flow)
    geometry = _continuous_geometry()
    cfg = _config(root)

    field_path = root / _FIELD_NAME
    save_field_npz(
        field_path,
        domain,
        raster,
        flow,
        continuous_geometry=geometry,
    )
    save_run_config(
        root / "run_config.yaml",
        {
            "input": {"model_dir": str(root)},
            "domain": asdict(cfg.domain),
            "field": asdict(cfg.field),
            "quick_test": {"enabled": False},
            "_resolved_field_reuse_contract": build_field_reuse_config_contract(cfg),
        },
    )
    return cfg, domain, raster, flow, geometry, field_path


class FieldReuseTests(unittest.TestCase):
    def test_round_trip_accepts_directory_and_direct_npz_path(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, expected, geometry, field_path = _write_reusable_result(root)

            from_directory = load_reusable_flow_field(
                root,
                cfg=cfg,
                domain=domain,
                raster=raster,
                continuous_geometry=geometry,
            )
            from_file = load_reusable_flow_field(
                field_path,
                cfg=cfg,
                domain=domain,
                raster=raster,
                continuous_geometry=geometry,
            )

            self.assertEqual(from_directory.field_npz_path, field_path.resolve())
            self.assertEqual(from_file.field_npz_path, field_path.resolve())
            self.assertEqual(
                from_directory.run_config_path,
                (root / "run_config.yaml").resolve(),
            )
            np.testing.assert_array_equal(
                from_directory.flow.velocity_xz_um_s,
                expected.velocity_xz_um_s,
            )
            np.testing.assert_array_equal(
                from_directory.flow.face_flux_x_um2_s,
                expected.face_flux_x_um2_s,
            )
            np.testing.assert_array_equal(
                from_directory.flow.open_face_kind,
                expected.open_face_kind,
            )
            np.testing.assert_array_equal(
                from_directory.flow.hybrid_velocity.finite_element
                .velocity_coefficients_um_s,
                expected.hybrid_velocity.finite_element
                .velocity_coefficients_um_s,
            )
            self.assertEqual(
                from_directory.flow.hybrid_velocity.regular_grid_distance_um,
                expected.hybrid_velocity.regular_grid_distance_um,
            )
            self.assertTrue(from_directory.flow.solver_metadata["physical_converged"])

    def test_v18_cache_reconstructs_missing_local_shear_from_fem(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, _, geometry, field_path = _write_reusable_result(root)
            with np.load(field_path, allow_pickle=False) as archive:
                legacy = {
                    key: np.asarray(archive[key])
                    for key in archive.files
                    if key != "local_shear_stress_pa"
                }
            legacy["field_schema_version"] = np.asarray(
                "v18_hybrid_fem_grid_continuous_wall"
            )
            np.savez_compressed(field_path, **legacy)

            reused = load_reusable_flow_field(
                field_path,
                cfg=cfg,
                domain=domain,
                raster=raster,
                continuous_geometry=geometry,
            )

            self.assertIsNotNone(reused.flow.local_shear_stress_pa)
            self.assertEqual(reused.flow.local_shear_stress_pa.shape, domain.shape)
            self.assertTrue(
                np.isfinite(reused.flow.local_shear_stress_pa).all()
            )
            self.assertEqual(
                reused.flow.solver_metadata["local_shear_gradient_sampling"],
                "reconstructed_from_exported_fem_velocity_at_cartesian_lumen_centres",
            )

    def test_rejects_current_raster_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, _, geometry, _ = _write_reusable_result(root)
            changed_radius = raster.radius_um.copy()
            changed_radius[0, 0] += 0.25
            mismatched_raster = replace(raster, radius_um=changed_radius)

            with self.assertRaises(FieldReuseValidationError):
                load_reusable_flow_field(
                    root,
                    cfg=cfg,
                    domain=domain,
                    raster=mismatched_raster,
                    continuous_geometry=geometry,
                )

    def test_rejects_continuous_geometry_hash_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, _, geometry, _ = _write_reusable_result(root)
            mismatched_geometry = SimpleNamespace(
                **{
                    **vars(geometry),
                    "geometry_hash_sha256": "f" * 64,
                }
            )

            with self.assertRaises(FieldReuseValidationError):
                load_reusable_flow_field(
                    root,
                    cfg=cfg,
                    domain=domain,
                    raster=raster,
                    continuous_geometry=mismatched_geometry,
                )

    def test_rejects_effective_field_config_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, _, geometry, _ = _write_reusable_result(root)
            mismatched_cfg = SimpleNamespace(
                model_dir=cfg.model_dir,
                domain=cfg.domain,
                field=replace(
                    cfg.field,
                    effective_thickness_um=cfg.field.effective_thickness_um + 1.0,
                ),
                quick_test=False,
            )

            with self.assertRaises(FieldReuseValidationError):
                load_reusable_flow_field(
                    root,
                    cfg=mismatched_cfg,
                    domain=domain,
                    raster=raster,
                    continuous_geometry=geometry,
                )

    def test_rejects_resolved_input_path_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, _, geometry, _ = _write_reusable_result(root)
            mismatched_cfg = SimpleNamespace(
                model_dir=(root / "different_model").resolve(),
                domain=cfg.domain,
                field=cfg.field,
                quick_test=False,
            )

            with self.assertRaises(FieldReuseValidationError):
                load_reusable_flow_field(
                    root,
                    cfg=mismatched_cfg,
                    domain=domain,
                    raster=raster,
                    continuous_geometry=geometry,
                )

    def test_rejects_field_not_marked_physically_converged(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, _, geometry, field_path = _write_reusable_result(root)
            with np.load(field_path, allow_pickle=False) as saved:
                payload = {name: saved[name].copy() for name in saved.files}
            payload["solver_physical_converged"] = np.asarray([False])
            np.savez_compressed(field_path, **payload)

            with self.assertRaises(FieldReuseValidationError):
                load_reusable_flow_field(
                    root,
                    cfg=cfg,
                    domain=domain,
                    raster=raster,
                    continuous_geometry=geometry,
                )

    def test_rejects_field_without_current_momentum_acceptance_schema(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, _, geometry, field_path = _write_reusable_result(root)
            with np.load(field_path, allow_pickle=False) as saved:
                payload = {name: saved[name].copy() for name in saved.files}
            payload.pop("solver_physical_acceptance_schema")
            np.savez_compressed(field_path, **payload)

            with self.assertRaises(FieldReuseValidationError):
                load_reusable_flow_field(
                    root,
                    cfg=cfg,
                    domain=domain,
                    raster=raster,
                    continuous_geometry=geometry,
                )

    def test_rejects_field_with_failed_linear_solver(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, _, geometry, field_path = _write_reusable_result(root)
            with np.load(field_path, allow_pickle=False) as saved:
                payload = {name: saved[name].copy() for name in saved.files}
            payload["solver_linear_solver_converged_reason"] = np.asarray([-3])
            np.savez_compressed(field_path, **payload)

            with self.assertRaises(FieldReuseValidationError):
                load_reusable_flow_field(
                    root,
                    cfg=cfg,
                    domain=domain,
                    raster=raster,
                    continuous_geometry=geometry,
                )

    def test_rejects_missing_source_run_config(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cfg, domain, raster, _, geometry, _ = _write_reusable_result(root)
            (root / "run_config.yaml").unlink()

            with self.assertRaises(FieldReuseValidationError):
                load_reusable_flow_field(
                    root,
                    cfg=cfg,
                    domain=domain,
                    raster=raster,
                    continuous_geometry=geometry,
                )


if __name__ == "__main__":
    unittest.main()
