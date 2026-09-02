from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace
import unittest

import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.cardiac.cardiac_waveform import build_periodic_ecg_flow_waveform
from ulm_microbubble_traj_gen_2D.utils.cardiac.cardiac_pulsatility import (
    CardiacPulsatility,
    build_cardiac_pulsatility,
)
from ulm_microbubble_traj_gen_2D.utils.core.config import (
    CardiacPulsatilityConfig,
    ParticleDynamicsConfig,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_mobility_transport import _evaluate_rhs
from ulm_microbubble_traj_gen_2D.utils.particles.particle_perfusion_schedule import (
    build_perfusion_schedule,
)
from ulm_microbubble_traj_gen_2D.utils.core.types import GridDomain, RasterizedVessels
from ulm_microbubble_traj_gen_2D.test_files.test_particle_mobility_transport import (
    _evaluation_context,
)
from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import (
    straight_channel_case as _straight_channel_case,
)


class CardiacWaveformTests(unittest.TestCase):
    def test_waveform_is_positive_periodic_and_cycle_mean_preserving(self) -> None:
        waveform = build_periodic_ecg_flow_waveform(300.0, 2048)

        self.assertAlmostEqual(waveform.period_s, 0.2, places=14)
        self.assertGreater(float(np.min(waveform.multiplier)), 0.0)
        self.assertAlmostEqual(waveform.cycle_mean, 1.0, places=13)
        np.testing.assert_allclose(
            waveform.evaluate(np.asarray([-0.037, 0.013, 0.267])),
            waveform.evaluate(np.asarray([-0.037, 0.013, 0.267]) + waveform.period_s),
            rtol=0.0,
            atol=2.0e-13,
        )
        self.assertAlmostEqual(
            float(waveform.integrate_s(0.0, waveform.period_s)),
            waveform.period_s,
            places=13,
        )

    def test_modulation_strength_scales_only_deviation_from_mean(self) -> None:
        full = build_periodic_ecg_flow_waveform(300.0, 2048, modulation_strength=1.0)
        reduced = build_periodic_ecg_flow_waveform(
            300.0, 2048, modulation_strength=0.25
        )
        steady = build_periodic_ecg_flow_waveform(300.0, 2048, modulation_strength=0.0)

        np.testing.assert_allclose(
            reduced.multiplier,
            1.0 + 0.25 * (full.multiplier - 1.0),
            rtol=0.0,
            atol=2.0e-15,
        )
        np.testing.assert_allclose(steady.multiplier, 1.0, rtol=0.0, atol=0.0)
        self.assertAlmostEqual(reduced.cycle_mean, 1.0, places=13)
        self.assertAlmostEqual(steady.cycle_mean, 1.0, places=13)

    def test_invalid_modulation_strength_is_rejected(self) -> None:
        for strength in (-0.01, 1.01, float("nan")):
            with self.subTest(strength=strength):
                with self.assertRaisesRegex(ValueError, "modulation_strength"):
                    build_periodic_ecg_flow_waveform(
                        300.0, 256, modulation_strength=strength
                    )

    def test_cardiac_schedule_inverts_integrated_number_flux(self) -> None:
        waveform = build_periodic_ecg_flow_waveform(300.0, 1024)
        cardiac = CardiacPulsatility(
            waveform=waveform,
            phase_offset_s=0.0,
            path_distance_um=np.zeros((1, 1), dtype=np.float32),
            delay_s=np.zeros((1, 1), dtype=np.float64),
            delay_gradient_s_per_um=np.zeros((1, 1, 2), dtype=np.float64),
            waveform_name="synthetic_ecg_envelope",
            preserve_cycle_mean_flow=True,
            modulation_strength=1.0,
            pulse_propagation_velocity_um_s=25_000.0,
        )
        model = SimpleNamespace(
            injection_rate_per_s=100.0,
            sample_radius_um=lambda quantile: 1.0,
            sample_position_grid=lambda radius, quantile: np.asarray([quantile, radius]),
        )

        schedule = build_perfusion_schedule(model, waveform.period_s, cardiac)

        self.assertEqual(schedule.count, 20)
        cumulative = 100.0 * np.asarray(
            cardiac.integrate_inlet_multiplier_s(0.0, schedule.planned_time_s)
        )
        np.testing.assert_allclose(
            cumulative,
            np.arange(schedule.count, dtype=float) + 0.5,
            rtol=0.0,
            atol=2.0e-12,
        )
        self.assertGreater(
            float(np.max(np.abs(np.diff(schedule.planned_time_s) - 0.01))),
            1.0e-5,
        )

class CardiacPropagationTests(unittest.TestCase):
    def test_parent_child_path_distance_is_continuous_at_the_junction(self) -> None:
        nx, nz = 21, 3
        shape = (nx, nz)
        domain = GridDomain(
            origin_um=np.zeros(3, dtype=float),
            spacing_um=1.0,
            shape=shape,
            fixed_y_um=0.0,
            x_coordinates_um=np.arange(nx, dtype=float),
            z_coordinates_um=np.arange(nz, dtype=float),
        )
        vessel_id = np.zeros(shape, dtype=np.int32)
        vessel_id[10:, :] = 1
        scalar = np.ones(shape, dtype=np.float32)
        vector = np.zeros((*shape, 2), dtype=np.float32)
        raster = RasterizedVessels(
            lumen_mask=np.ones(shape, dtype=bool),
            wall_mask=np.zeros(shape, dtype=bool),
            vessel_id=vessel_id,
            radius_um=scalar,
            flow_rate_um3_s=scalar,
            q2d_flow_um2_s=scalar,
            viscosity_mpas=scalar,
            direction_xz=vector,
            distance_to_centerline_um=scalar,
            distance_to_wall_um=scalar,
            wall_normal_xz=vector,
            lumen_fraction=np.ones(shape, dtype=np.float32),
        )
        root = Vessel(
            vid=0,
            parent_id=-1,
            children=[1],
            x_p=np.asarray([0.0, 0.0, 1.0]),
            x_d=np.asarray([10.0, 0.0, 1.0]),
            radius=1.0,
        )
        child = Vessel(
            vid=1,
            parent_id=0,
            children=[],
            x_p=np.asarray([10.0, 0.0, 1.0]),
            x_d=np.asarray([20.0, 0.0, 1.0]),
            radius=1.0,
        )

        cardiac = build_cardiac_pulsatility(
            domain,
            raster,
            [root, child],
            CardiacPulsatilityConfig(enabled=True, waveform_samples_per_cycle=256),
        )

        assert cardiac is not None
        np.testing.assert_allclose(
            cardiac.path_distance_um[:, 1], np.arange(nx, dtype=float), rtol=0.0, atol=1.0e-7
        )

    def test_root_path_delay_uses_physical_seconds_not_period_sample_count(self) -> None:
        nx, nz = 1001, 3
        shape = (nx, nz)
        domain = GridDomain(
            origin_um=np.zeros(3, dtype=float),
            spacing_um=1.0,
            shape=shape,
            fixed_y_um=0.0,
            x_coordinates_um=np.arange(nx, dtype=float),
            z_coordinates_um=np.arange(nz, dtype=float),
        )
        lumen = np.ones(shape, dtype=bool)
        scalar = np.ones(shape, dtype=np.float32)
        vector = np.zeros((*shape, 2), dtype=np.float32)
        raster = RasterizedVessels(
            lumen_mask=lumen,
            wall_mask=np.zeros(shape, dtype=bool),
            vessel_id=np.zeros(shape, dtype=np.int32),
            radius_um=scalar,
            flow_rate_um3_s=scalar,
            q2d_flow_um2_s=scalar,
            viscosity_mpas=scalar,
            direction_xz=vector,
            distance_to_centerline_um=scalar,
            distance_to_wall_um=scalar,
            wall_normal_xz=vector,
            lumen_fraction=np.ones(shape, dtype=np.float32),
        )
        root = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 1.0]),
            x_d=np.asarray([1000.0, 0.0, 1.0]),
            radius=1.0,
        )
        cfg = CardiacPulsatilityConfig(
            enabled=True,
            bpm=300.0,
            pulse_propagation_velocity_um_s=25_000.0,
            waveform_samples_per_cycle=1024,
        )

        cardiac = build_cardiac_pulsatility(domain, raster, [root], cfg)

        self.assertIsNotNone(cardiac)
        assert cardiac is not None
        self.assertAlmostEqual(float(cardiac.path_distance_um[-1, 1]), 1000.0, places=6)
        self.assertAlmostEqual(float(cardiac.delay_s[-1, 1]), 0.04, places=13)
        sample = cardiac.sample(np.asarray([[1000.0, 1.0]]), time_s=0.07)
        expected = cardiac.waveform.evaluate(0.07 - 0.04)
        self.assertAlmostEqual(float(sample.multiplier[0]), float(expected), places=13)

        positions = np.asarray([[125.25, 0.75], [500.5, 1.0], [999.0, 1.25]])
        reference = cardiac.sample(positions, time_s=0.071, use_numba=False)
        accelerated = cardiac.sample(positions, time_s=0.071, use_numba=True)
        np.testing.assert_allclose(
            accelerated.multiplier, reference.multiplier, rtol=0.0, atol=2.0e-15
        )
        np.testing.assert_allclose(
            accelerated.gradient_per_um,
            reference.gradient_per_um,
            rtol=0.0,
            atol=2.0e-15,
        )

    def test_rhs_scales_carrier_velocity_and_sampled_shear_at_internal_time(self) -> None:
        domain, raster, flow = _straight_channel_case()
        root = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 3.0]),
            x_d=np.asarray([19.0, 0.0, 3.0]),
            radius=3.0,
        )
        cardiac = build_cardiac_pulsatility(
            domain,
            raster,
            [root],
            CardiacPulsatilityConfig(enabled=True, waveform_samples_per_cycle=1024),
        )
        assert cardiac is not None
        context = _evaluation_context(
            domain,
            raster,
            flow,
            ParticleDynamicsConfig(near_wall_enabled=False, collisions_enabled=False),
            radii_um=np.asarray([0.5]),
        )
        context = replace(
            context,
            wall_shear_stress_pa=np.full(domain.shape, 2.0, dtype=np.float64),
            cardiac=cardiac,
        )
        position = np.asarray([[10.0, 3.0]])
        time_s = 0.061

        evaluation = _evaluate_rhs(
            position,
            np.asarray([0]),
            np.asarray([True]),
            context,
            time_s,
        )

        multiplier = float(cardiac.sample(position, time_s).multiplier[0])
        self.assertAlmostEqual(float(evaluation.cardiac_multiplier[0]), multiplier, places=13)
        self.assertAlmostEqual(float(evaluation.fluid_velocity_xz_um_s[0, 0]), 10.0 * multiplier)
        self.assertAlmostEqual(float(evaluation.particle_velocity_xz_um_s[0, 0]), 10.0 * multiplier)
        self.assertAlmostEqual(
            float(evaluation.sampled_wall_shear_stress_pa[0]), 2.0 * multiplier
        )

    def test_fused_numba_pulsatile_rhs_matches_python_reference(self) -> None:
        domain, raster, flow = _straight_channel_case()
        root = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 3.0]),
            x_d=np.asarray([19.0, 0.0, 3.0]),
            radius=3.0,
        )
        cardiac = build_cardiac_pulsatility(
            domain,
            raster,
            [root],
            CardiacPulsatilityConfig(enabled=True, waveform_samples_per_cycle=1024),
        )
        assert cardiac is not None
        reference_context = _evaluation_context(
            domain,
            raster,
            flow,
            ParticleDynamicsConfig(near_wall_enabled=False, collisions_enabled=False),
            radii_um=np.asarray([0.5, 0.6, 0.45]),
        )
        reference_context = replace(reference_context, cardiac=cardiac)
        accelerated_context = replace(reference_context, use_numba=True)
        positions = np.asarray(
            [[4.25, 2.1], [10.5, 3.0], [15.75, 3.8]],
            dtype=np.float64,
        )
        bubble_ids = np.arange(positions.shape[0], dtype=np.int64)
        active = np.ones(positions.shape[0], dtype=bool)
        time_s = 0.061

        reference = _evaluate_rhs(
            positions, bubble_ids, active, reference_context, time_s
        )
        accelerated = _evaluate_rhs(
            positions, bubble_ids, active, accelerated_context, time_s
        )

        for field_name in (
            "particle_velocity_xz_um_s",
            "fluid_velocity_xz_um_s",
            "angular_velocity_rad_s",
            "sampled_wall_shear_stress_pa",
            "wall_gap_um",
            "wall_normal_xz",
            "cardiac_multiplier",
            "generalized_mobility",
        ):
            np.testing.assert_allclose(
                getattr(accelerated, field_name),
                getattr(reference, field_name),
                rtol=2.0e-13,
                atol=2.0e-13,
                equal_nan=True,
            )


if __name__ == "__main__":
    unittest.main()
