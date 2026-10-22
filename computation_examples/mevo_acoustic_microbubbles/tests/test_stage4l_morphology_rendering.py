"""Rendering checks for physical size, occlusion and unchanged input geometry."""

import numpy as np
from PIL import Image
import pytest

from scripts.stage4l import render_morphology as r


@pytest.fixture
def scene():
    fig = r.style.base.plt.figure(figsize=(19.2, 10.8), dpi=200)
    ax = fig.add_subplot(projection="3d")
    r.style.original.axes_style(ax, np.array([[-1.0, 1.0]] * 3), ["x", "y", "z"])
    fig.canvas.draw()
    yield ax, Image.new("RGB", (r.WIDTH, r.HEIGHT))
    r.style.base.plt.close(fig)


def test_depth_occlusion_is_independent_of_submission_order(scene):
    ax, background = scene
    _, _, _, _, eye = r.projection(ax, np.zeros((1, 3)))
    positions = np.array([-eye * 20e-6, eye * 20e-6])
    radii = np.full(2, 4e-6)
    ids = np.array([11, 29])
    before = positions.copy()
    a, evidence = r.analytic_spheres(background, ax, positions, radii, ids)
    b, reverse = r.analytic_spheres(background, ax, positions[::-1], radii[::-1], ids[::-1])
    assert evidence["submitted_ids"] == [11, 29]
    assert evidence["covered_ids"] == reverse["covered_ids"] == [29]
    assert np.array_equal(np.asarray(a), np.asarray(b))
    assert np.array_equal(positions, before)


def test_subpixel_sphere_is_not_enlarged_to_a_marker(scene):
    ax, background = scene
    x = np.zeros((1, 3))
    _, scale, _, _, _ = r.projection(ax, x)
    physical_radii = np.array([0.20 / scale])
    _, evidence = r.analytic_spheres(background, ax, x, physical_radii, np.array([7]), supersampling=64)
    assert np.isclose(evidence["radii_px"][0], 0.20)
    assert abs(evidence["covered_pixel_area"] - np.pi * 0.20**2) < 0.003
    assert not evidence["minimum_glyph_size"]


def test_separated_bubbles_remain_distinct_and_use_common_scale(scene):
    ax, background = scene
    _, scale, u, _, _ = r.projection(ax, np.zeros((1, 3)))
    x = np.array([-u * 15e-6, u * 15e-6])
    radii = np.array([2e-6, 4e-6])
    _, evidence = r.analytic_spheres(background, ax, x, radii, np.array([4, 8]))
    assert evidence["covered_ids"] == [4, 8]
    np.testing.assert_allclose(evidence["radii_px"], radii * scale, rtol=1e-12)
    assert np.isclose(np.linalg.norm(np.diff(evidence["centers_px"], axis=0)), 30e-6 * scale)
