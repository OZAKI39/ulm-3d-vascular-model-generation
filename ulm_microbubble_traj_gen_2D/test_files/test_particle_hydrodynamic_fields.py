from __future__ import annotations

import unittest

import numpy as np
from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.core.types import (
    FlowField,
    GridDomain,
    RasterizedVessels,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_hydrodynamic_fields import (
    build_particle_hydrodynamic_fields,
)
from ulm_microbubble_traj_gen_2D.test_files.hybrid_velocity_fixture import (
    rectangular_hybrid_velocity,
)


class ParticleHydrodynamicFieldTests(unittest.TestCase):
    def test_linear_velocity_field_has_exact_gradient_and_vorticity(self) -> None:
        domain, raster, flow, geometry = _channel_case()
        x = domain.x_coordinates_um[:, None]
        z = domain.z_coordinates_um[None, :]
        velocity = np.empty((*domain.shape, 2), dtype=np.float32)
        velocity[..., 0] = 1.0 + 2.0 * x + 3.0 * z
        velocity[..., 1] = -4.0 - 5.0 * x + 7.0 * z
        flow = _flow_with_velocity(flow, velocity)

        fields = build_particle_hydrodynamic_fields(
            domain, raster, flow, continuous_geometry=geometry
        )

        expected_gradient = np.broadcast_to(
            np.asarray([[2.0, 3.0], [-5.0, 7.0]], dtype=np.float32),
            (*domain.shape, 2, 2),
        )
        np.testing.assert_allclose(fields.velocity_gradient_s_inv, expected_gradient, rtol=0.0, atol=1.0e-6)

    def test_zero_solid_viscosity_is_filled_from_nearest_fluid(self) -> None:
        domain, raster, flow, geometry = _channel_case(viscosity_mpas=4.25)

        fields = build_particle_hydrodynamic_fields(
            domain, raster, flow, continuous_geometry=geometry
        )

        np.testing.assert_allclose(fields.dynamic_viscosity_pa_s, 4.25e-3, rtol=0.0, atol=1.0e-9)
        self.assertTrue(np.all(fields.dynamic_viscosity_pa_s[~raster.lumen_mask] > 0.0))

    def test_wall_normal_points_from_each_solid_wall_into_lumen(self) -> None:
        domain, raster, flow, geometry = _channel_case()

        fields = build_particle_hydrodynamic_fields(
            domain, raster, flow, continuous_geometry=geometry
        )
        points = np.asarray([[2.0, 0.5], [2.0, 2.5]])
        state = fields.boundary_geometry.exact_solid_wall_state_xz_um(points)
        np.testing.assert_allclose(state.distance_um, [0.25, 0.25], atol=1.0e-6)
        np.testing.assert_allclose(
            state.inward_normal_xz,
            [[0.0, 1.0], [0.0, -1.0]],
            atol=1.0e-6,
        )

    def test_shapes_dtypes_and_contiguous_storage_are_stable(self) -> None:
        domain, raster, flow, geometry = _channel_case()

        fields = build_particle_hydrodynamic_fields(
            domain, raster, flow, continuous_geometry=geometry
        )

        expected_float_fields = {
            "velocity_gradient_s_inv": ((*domain.shape, 2, 2), np.float32),
            "dynamic_viscosity_pa_s": (domain.shape, np.float32),
        }
        for name, (shape, dtype) in expected_float_fields.items():
            values = getattr(fields, name)
            self.assertEqual(values.shape, shape)
            self.assertEqual(values.dtype, dtype)
            self.assertTrue(values.flags.c_contiguous)
        for name in ("solid_site_mask", "open_boundary_mask"):
            values = getattr(fields, name)
            self.assertEqual(values.shape, domain.shape)
            self.assertEqual(values.dtype, np.bool_)
            self.assertTrue(values.flags.c_contiguous)


def _channel_case(viscosity_mpas: float = 3.5):
    nx, nz = 9, 7
    shape = (nx, nz)
    spacing_um = 0.5
    lumen = np.zeros(shape, dtype=bool)
    lumen[:, 1:-1] = True
    velocity = np.zeros((*shape, 2), dtype=np.float32)
    viscosity = np.zeros(shape, dtype=np.float32)
    viscosity[lumen] = np.float32(viscosity_mpas)
    inlet = np.zeros(shape, dtype=np.int32)
    outlet = np.zeros(shape, dtype=np.int32)
    inlet[0, 2:-2] = 1
    outlet[-1, 2:-2] = 1

    domain = GridDomain(
        origin_um=np.asarray([0.0, 0.0, 0.0]),
        spacing_um=spacing_um,
        shape=shape,
        fixed_y_um=0.0,
        x_coordinates_um=np.arange(nx, dtype=float) * spacing_um,
        z_coordinates_um=np.arange(nz, dtype=float) * spacing_um,
    )
    raster = RasterizedVessels(
        lumen_mask=lumen,
        wall_mask=lumen.copy(),
        vessel_id=np.zeros(shape, dtype=np.int32),
        radius_um=np.full(shape, 1.5, dtype=np.float32),
        flow_rate_um3_s=np.ones(shape, dtype=np.float32),
        q2d_flow_um2_s=np.ones(shape, dtype=np.float32),
        viscosity_mpas=viscosity,
        direction_xz=np.zeros((*shape, 2), dtype=np.float32),
        distance_to_centerline_um=np.zeros(shape, dtype=np.float32),
        distance_to_wall_um=np.zeros(shape, dtype=np.float32),
        wall_normal_xz=np.zeros((*shape, 2), dtype=np.float32),
        lumen_fraction=lumen.astype(np.float32),
    )
    flow = FlowField(
        velocity_xz_um_s=velocity,
        speed_um_s=np.zeros(shape, dtype=np.float32),
        wall_shear_stress_pa=np.zeros(shape, dtype=np.float32),
        hybrid_velocity=rectangular_hybrid_velocity(domain),
        inlet_label=inlet,
        outlet_label=outlet,
        solver_metadata={"physical_converged": True},
    )
    vessel = Vessel(vid=0, parent_id=-1, children=[])
    vessel.x_p = np.asarray([0.0, 0.0, 1.5])
    vessel.x_d = np.asarray([4.0, 0.0, 1.5])
    vessel.radius = 1.25
    geometry = build_continuous_vessel_geometry([vessel], domain)
    return domain, raster, flow, geometry


def _flow_with_velocity(flow: FlowField, velocity: np.ndarray) -> FlowField:
    return FlowField(
        velocity_xz_um_s=velocity,
        speed_um_s=np.linalg.norm(velocity, axis=-1).astype(np.float32),
        wall_shear_stress_pa=flow.wall_shear_stress_pa,
        hybrid_velocity=flow.hybrid_velocity,
        inlet_label=flow.inlet_label,
        outlet_label=flow.outlet_label,
        solver_metadata=flow.solver_metadata,
    )


if __name__ == "__main__":
    unittest.main()
