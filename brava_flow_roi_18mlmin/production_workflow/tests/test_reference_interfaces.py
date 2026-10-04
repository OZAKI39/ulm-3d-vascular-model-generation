"""Small in-memory interface regressions; no mesher, CFD or trajectory is run."""
import sys
from pathlib import Path
import numpy as np
import pyvista as pv
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from pipe_probe_geometry import strict_cells
from poiseuille import reference


@pytest.fixture
def tetra_fixture():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    tetra = np.array([[0, 1, 2, 3]])
    grid = pv.UnstructuredGrid(np.array([4, 0, 1, 2, 3]), np.array([10], dtype=np.uint8), points)
    return grid, points, tetra


def test_strict_contains_interior_face_and_rejects_outside(tetra_fixture):
    grid, points, tetra = tetra_fixture
    probes = np.array([[.1, .2, .3], [0., .2, .3], [-1e-8, .2, .3], [.5, .5, .5]])
    cells, _ = strict_cells(grid, points, tetra, probes)
    assert cells.tolist() == [0, 0, -1, -1]


def test_strict_rejects_locator_false_positive(tetra_fixture, monkeypatch):
    grid, points, tetra = tetra_fixture
    monkeypatch.setattr(type(grid), 'find_containing_cell', lambda self, q: np.zeros(len(q), dtype=int))
    cells, original = strict_cells(grid, points, tetra, np.array([[-1e-8, .2, .3]]))
    assert original.tolist() == [0] and cells.tolist() == [-1]


def test_strict_recovers_locator_false_negative(tetra_fixture, monkeypatch):
    grid, points, tetra = tetra_fixture
    monkeypatch.setattr(type(grid), 'find_containing_cell', lambda self, q: np.full(len(q), -1, dtype=int))
    cells, original = strict_cells(grid, points, tetra, np.array([[.1, .2, .3]]))
    assert original.tolist() == [-1] and cells.tolist() == [0]


def test_analytic_axis_wall_and_incompressibility():
    radius, length, mean, mu = 1e-3, 12e-3, .1, .00345312
    result = reference(np.array([[0., 0., 0.], [radius, 0., length]]), radius, length, mean, mu)
    np.testing.assert_allclose(result['u_m_s'][:, 2], [2 * mean, 0.], atol=1e-16)
    np.testing.assert_allclose(np.trace(result['G_s_inv'], axis1=1, axis2=2), 0.)
    assert result['p_Pa'][0] > 0 and result['p_Pa'][1] == 0
    assert result['G_s_inv'][1, 2, 0] < 0


def test_analytic_rejects_outside_cylinder():
    with pytest.raises(ValueError, match='outside'):
        reference(np.array([[1.001e-3, 0., 0.]]), 1e-3, .012, .1, .00345312)
