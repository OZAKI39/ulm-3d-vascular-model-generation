from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.visualization.results.bubble_style import (
    BubbleGlowStyle,
    bubble_highlight_geometry,
    build_bubble_glow_style,
    glow_layer_diameters,
    glow_layer_diameters_points,
)
from ulm_microbubble_traj_gen_2D.utils.visualization.results.result_loader import (
    _realized_center_velocities_um_s,
)
from ulm_microbubble_traj_gen_2D.utils.visualization.results.plotting import (
    _active_frame_selection,
    _bubble_color_values,
    _bubble_contact_ring_colors,
    _bubble_edge_colors,
    _bubble_physical_diameters,
    _colorbar_label,
    _select_bubbles,
    _tail_segments,
)


class ParticleVisualizationIdentityTests(unittest.TestCase):
    def test_realized_speed_uses_post_constraint_motion_and_never_crosses_replacement_ids(self) -> None:
        positions = np.zeros((5, 1, 3), dtype=float)
        positions[:, 0, 0] = np.asarray([0.0, 1.0, 1.0, 10.0, 10.0])
        bubble_id = np.asarray([[3], [3], [3], [8], [8]], dtype=np.int64)
        active = np.ones((5, 1), dtype=bool)

        realized = _realized_center_velocities_um_s(positions, bubble_id, active, dt_s=0.5)

        np.testing.assert_allclose(realized[:, 0, 0], np.asarray([2.0, 0.0, 0.0, 0.0, 0.0]))
        self.assertEqual(float(realized[2, 0, 0]), 0.0)
        self.assertEqual(float(realized[3, 0, 0]), 0.0)

    def test_speed_color_prefers_realized_center_velocity_over_pre_constraint_rhs(self) -> None:
        data = SimpleNamespace(
            positions_um=np.zeros((2, 1, 3), dtype=float),
            velocities_um_s=np.asarray([[[500.0, 0.0, 0.0]], [[500.0, 0.0, 0.0]]]),
            realized_velocities_um_s=np.zeros((2, 1, 3), dtype=float),
        )

        values = _bubble_color_values(data, np.asarray([0]), frame=1, color_mode="speed")

        np.testing.assert_array_equal(values, np.asarray([0.0]))
        self.assertEqual(_colorbar_label("speed"), "realized center speed (um/s)")

    def test_tail_starts_again_when_a_display_lane_receives_a_new_unique_id(self) -> None:
        positions = np.zeros((4, 1, 3), dtype=float)
        positions[:, 0, 0] = np.asarray([0.0, 1.0, 10.0, 11.0])
        data = SimpleNamespace(
            positions_um=positions,
            bubble_id=np.asarray([[3], [3], [8], [8]], dtype=np.int64),
            active=np.ones((4, 1), dtype=bool),
        )
        segments = _tail_segments(data, np.asarray([0]), frame=3, tail_length=10)
        self.assertEqual(len(segments), 1)
        np.testing.assert_array_equal(segments[0][:, 0], np.asarray([10.0, 11.0]))

    def test_zero_display_limit_selects_the_full_population(self) -> None:
        data = SimpleNamespace(positions_um=np.zeros((2, 7, 3), dtype=float))
        np.testing.assert_array_equal(_select_bubbles(data, 0), np.arange(7))

    def test_frame_selection_draws_only_real_active_bubbles(self) -> None:
        data = SimpleNamespace(
            positions_um=np.asarray([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [np.nan, 0.0, 0.0]]]),
            active=np.asarray([[True, False, True]]),
            bubble_id=np.asarray([[10, 11, -1]], dtype=np.int64),
        )
        selected = _active_frame_selection(data, np.arange(3), frame=0)
        np.testing.assert_array_equal(selected, np.asarray([0]))

    def test_marker_area_reflects_each_bubble_diameter(self) -> None:
        data = SimpleNamespace(bubble_diameter_um=np.asarray([[1.5, 2.0, 2.5]], dtype=float))
        diameters = _bubble_physical_diameters(data, np.arange(3), frame=0)
        np.testing.assert_array_equal(diameters, np.asarray([1.5, 2.0, 2.5]))

    def test_glow_enlarges_only_visual_layers_and_keeps_physical_diameters_unchanged(self) -> None:
        physical = np.asarray([1.5, 2.0, 2.5], dtype=float)
        original = physical.copy()
        style = BubbleGlowStyle(outer_diameter_scale=3.0)
        layers = glow_layer_diameters(physical, style)

        np.testing.assert_array_equal(physical, original)
        np.testing.assert_allclose(style.layer_diameter_scales, (3.0, 2.2, 1.5))
        np.testing.assert_allclose(layers[0], physical * 3.0)
        np.testing.assert_allclose(layers[-1], physical * 1.5)

    def test_spherical_highlight_remains_inside_the_real_bubble_disk(self) -> None:
        xy = np.asarray([[10.0, 20.0], [30.0, 40.0]])
        diameters = np.asarray([1.5, 2.5])
        style = BubbleGlowStyle()
        offsets, highlight_diameters = bubble_highlight_geometry(xy, diameters, style)

        center_shift = np.linalg.norm(offsets - xy, axis=1)
        self.assertTrue(np.all(center_shift + 0.5 * highlight_diameters < 0.5 * diameters))

    def test_highlight_follows_saved_positive_y_rotation_angle(self) -> None:
        xy = np.asarray([[10.0, 20.0], [10.0, 20.0]])
        diameters = np.asarray([2.0, 2.0])
        style = BubbleGlowStyle()

        offsets, _ = bubble_highlight_geometry(
            xy,
            diameters,
            style,
            rotation_angle_rad=np.asarray([0.0, 0.5 * np.pi]),
        )

        np.testing.assert_allclose(offsets[0] - xy[0], np.asarray([-0.3, 0.3]), atol=1.0e-12)
        np.testing.assert_allclose(offsets[1] - xy[1], np.asarray([0.3, 0.3]), atol=1.0e-12)

    def test_contact_uses_outer_orange_ring_without_recoloring_physical_core(self) -> None:
        data = SimpleNamespace(
            bubble_wall_gap_um=np.asarray([[0.01, -1.0e-12]], dtype=float),
            bubble_wall_contact=np.asarray([[True, True]], dtype=bool),
            wall_gap_invalid_below_um=0.0,
        )
        selected = np.asarray([0, 1])

        core = _bubble_edge_colors(data, selected, frame=0)
        ring = _bubble_contact_ring_colors(data, selected, frame=0)

        np.testing.assert_array_equal(core[0], np.asarray([0.0, 0.0, 0.0, 1.0]))
        np.testing.assert_array_equal(ring[0], np.asarray([1.0, 0.65, 0.0, 1.0]))
        np.testing.assert_array_equal(core[1], np.asarray([1.0, 0.0, 0.0, 1.0]))
        np.testing.assert_array_equal(ring[1], np.asarray([1.0, 0.0, 0.0, 1.0]))

    def test_near_wall_and_target_adhesion_use_different_ring_colors(self) -> None:
        data = SimpleNamespace(
            bubble_wall_gap_um=np.asarray([[0.2, 0.1, 2.0]], dtype=float),
            bubble_wall_contact=np.zeros((1, 3), dtype=bool),
            bubble_near_wall_weight=np.asarray([[0.5, 0.8, 0.0]], dtype=float),
            bubble_bond_count_expected=np.asarray([[0.0, 2.5, 0.0]], dtype=float),
            wall_gap_invalid_below_um=0.0,
        )

        ring = _bubble_contact_ring_colors(data, np.arange(3), frame=0)

        np.testing.assert_array_equal(ring[0], np.asarray([1.0, 0.65, 0.0, 1.0]))
        np.testing.assert_array_equal(ring[1], np.asarray([1.0, 0.1, 0.65, 1.0]))
        np.testing.assert_array_equal(ring[2], np.zeros(4))

    def test_legacy_contact_colors_keep_the_saved_penetration_allowance(self) -> None:
        data = SimpleNamespace(
            bubble_wall_gap_um=np.asarray([[-0.01, -0.10]], dtype=float),
            bubble_wall_contact=np.asarray([[True, True]], dtype=bool),
            wall_gap_invalid_below_um=-0.05,
        )
        selected = np.asarray([0, 1])

        core = _bubble_edge_colors(data, selected, frame=0)

        np.testing.assert_array_equal(core[0], np.asarray([0.0, 0.0, 0.0, 1.0]))
        np.testing.assert_array_equal(core[1], np.asarray([1.0, 0.0, 0.0, 1.0]))

    def test_overview_halo_has_a_screen_size_floor_without_changing_the_core(self) -> None:
        physical = np.asarray([1.5, 2.5], dtype=float)
        style = BubbleGlowStyle(outer_diameter_scale=3.0, minimum_outer_diameter_points=6.0)
        layers = glow_layer_diameters_points(physical, points_per_um=0.1, style=style)

        np.testing.assert_allclose(layers[0], np.asarray([6.0, 6.0]))
        np.testing.assert_allclose(layers[1], np.asarray([6.0 * 2.2 / 3.0] * 2))
        np.testing.assert_array_equal(physical, np.asarray([1.5, 2.5]))

    def test_glow_scale_controls_the_overview_size_as_well_as_the_physical_ratio(self) -> None:
        style = build_bubble_glow_style(enabled=True, outer_diameter_scale=4.0)
        self.assertEqual(style.outer_diameter_scale, 4.0)
        self.assertEqual(style.minimum_outer_diameter_points, 8.0)


if __name__ == "__main__":
    unittest.main()
