import numpy as np
import pytest
from particle_3d.geometry import TetraGeometry
from particle_3d.validation_cases import AFFINE_ATOL, AFFINE_RTOL


def test_known_weights_all_point_kinds(affine):
    field = affine["field"]
    actual = field.geometry.weights(affine["positions"], np.zeros(len(affine["positions"]), dtype=int))
    np.testing.assert_allclose(actual, affine["weights"], atol=AFFINE_ATOL, rtol=AFFINE_RTOL)
    np.testing.assert_allclose(actual.sum(axis=1), 1., atol=AFFINE_ATOL, rtol=0)
    assert affine["batch"].inside_lumen.all()
    assert set(affine["kinds"]) == {"centroid", "random_interior", "face", "edge", "vertex"}
    print("max weight error", np.abs(actual - affine["weights"]).max())


@pytest.mark.parametrize("scale,offset", [(1e-6, 1e-4), (1., 0.), (1e3, 1e4)])
def test_scale_and_translation_roundoff_budget(scale, offset):
    vertices = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]) * scale + offset
    geometry = TetraGeometry(vertices, np.array([[0, 1, 2, 3]]))
    expected = np.array([.1, .2, .3, .4])
    actual = geometry.weights(expected @ vertices, 0)
    np.testing.assert_allclose(actual, expected, atol=geometry.weight_tolerance[0], rtol=0)


def test_degenerate_tetra_rejected():
    with pytest.raises(ValueError, match="degenerate"):
        TetraGeometry(np.zeros((4, 3)), np.array([[0, 1, 2, 3]]))
