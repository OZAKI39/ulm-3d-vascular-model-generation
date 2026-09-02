from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.io.field_io import save_field_npz
from ulm_microbubble_traj_gen_2D.utils.particles.particle_hydrodynamic_fields import (
    build_particle_hydrodynamic_fields,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_inlet_flux import (
    build_inlet_flux_model,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_geometry_tolerance import (
    project_roundoff_negative_wall_gap_xz_um,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_perfusion_schedule import (
    build_perfusion_schedule,
)
from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import (
    particle_config,
    straight_channel_case,
)


class ParticleInletFluxV15Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.domain, cls.raster, cls.flow = straight_channel_case()
        cls.root = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 3.0]),
            x_d=np.asarray([19.0, 0.0, 3.0]),
            radius=3.0,
            flow_rate=50.0,
        )
        cls.geometry = build_continuous_vessel_geometry(
            [cls.root], cls.domain
        )
        cls.fields = build_particle_hydrodynamic_fields(
            cls.domain,
            cls.raster,
            cls.flow,
            continuous_geometry=cls.geometry,
        )

    def _model(self, diameter_um: float = 1.0):
        geometry = self.fields.boundary_geometry
        self.assertIsNotNone(geometry)
        return build_inlet_flux_model(
            self.domain,
            self.flow,
            [self.root],
            geometry,
            replace(
                particle_config(n_steps=1),
                bubble_diameter_min_um=diameter_um,
                bubble_diameter_max_um=diameter_um,
            ),
            effective_thickness_um=1.0,
            boundary_depth_cells=1.5,
        )

    def test_every_sampled_inlet_centre_is_inside_clear_and_inlet_connected(self) -> None:
        model = self._model(1.0)
        geometry = self.fields.boundary_geometry
        radius_um = 0.5
        section = model.sections[0]
        self.assertTrue(np.all(section.inside_lumen))
        self.assertEqual(float(section.wall_distance_um[0]), 0.0)
        self.assertEqual(float(section.wall_distance_um[-1]), 0.0)

        for quantile in np.linspace(0.0, 1.0 - 1.0e-12, 101):
            position_grid = model.sample_position_grid(radius_um, float(quantile))
            position_xz_um = geometry.grid_to_world_xz(position_grid)
            self.assertTrue(bool(geometry.contains_xz_um(position_xz_um)))
            self.assertGreaterEqual(
                float(geometry.true_gap_at_xz_um(position_xz_um, radius_um)),
                0.0,
            )
            self.assertTrue(bool(geometry.is_accessible_grid(position_grid, radius_um)))

        self.assertAlmostEqual(model.raw_section_flow_um2_s, 47.5)
        self.assertAlmostEqual(model.size_accessible_section_flow_um2_s, 47.5)
        self.assertAlmostEqual(model.injection_rate_per_s, 2.375)

    def test_radius_larger_than_inlet_connected_bottleneck_is_rejected(self) -> None:
        with self.assertRaisesRegex(
            ValueError,
            "finite-size bubble distribution cannot pass|zero MB/s flux",
        ):
            self._model(6.1)

    def test_roundoff_negative_inlet_gap_is_projected_strictly_feasible(self) -> None:
        corrected, gap_um, tolerance_um = (
            project_roundoff_negative_wall_gap_xz_um(
                np.asarray([5.0, 5.0 + 3.7e-13]),
                1.0,
                self.geometry,
            )
        )

        self.assertGreater(tolerance_um, 3.7e-13)
        self.assertGreaterEqual(gap_um, 0.0)
        self.assertTrue(bool(self.geometry.contains_xz_um(corrected)))
        self.assertTrue(bool(self.geometry.is_accessible_xz_um(corrected, 1.0)))

    def test_material_inlet_penetration_is_never_hidden_by_roundoff_repair(self) -> None:
        original = np.asarray([5.0, 5.0 + 1.0e-4])
        corrected, gap_um, tolerance_um = (
            project_roundoff_negative_wall_gap_xz_um(
                original,
                1.0,
                self.geometry,
            )
        )

        np.testing.assert_array_equal(corrected, original)
        self.assertLess(gap_um, -tolerance_um)

    def test_corrected_flux_does_not_create_an_event_by_point_two_seconds(self) -> None:
        model = self._model(1.0)
        short_schedule = build_perfusion_schedule(model, 0.2)
        full_schedule = build_perfusion_schedule(model, 1.0)

        self.assertEqual(short_schedule.count, 0)
        half_interval = 0.5 * model.mean_injection_interval_s
        np.testing.assert_allclose(
            full_schedule.planned_time_s,
            [half_interval, 3.0 * half_interval],
        )

    def test_saved_field_contains_authoritative_open_boundary_schema(self) -> None:
        geometry = self.fields.boundary_geometry
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "field.npz"
            save_field_npz(
                path,
                self.domain,
                self.raster,
                self.flow,
                continuous_geometry=geometry,
            )
            with np.load(path, allow_pickle=False) as saved:
                self.assertEqual(
                    str(saved["field_schema_version"].item()),
                    "v19_hybrid_fem_grid_local_shear",
                )
                self.assertIn("local_shear_stress_pa", saved.files)
                np.testing.assert_array_equal(
                    saved["lumen_fraction"],
                    self.raster.lumen_fraction,
                )
                for name in (
                    "continuous_wall_start_xz_um",
                    "continuous_wall_end_xz_um",
                    "continuous_wall_inward_normal_xz",
                    "continuous_open_section_point_xz_um",
                    "continuous_open_section_outward_normal_xz",
                    "continuous_open_section_tangent_xz",
                ):
                    self.assertIn(name, saved.files)


if __name__ == "__main__":
    unittest.main()
