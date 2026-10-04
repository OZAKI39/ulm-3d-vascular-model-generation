"""Analytical counterexamples for the finite-facet geometry acceptance bound."""
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'gmsh_reference'))
from audit_reference import finite_facet_deviation


def test_cylindrical_chord_uses_edge_minimum_not_triangle_centroid():
    radius, angle = .0005967792299824433, np.pi/272
    x, y = radius*np.cos(angle), radius*np.sin(angle)
    # The most inward point is on an edge; the centroid underestimates sagitta.
    triangle = np.array([[[x,-y,0.], [x,y,0.], [x,y,.001]]])
    actual = finite_facet_deviation(triangle, radius)
    np.testing.assert_allclose(actual, radius*(1-np.cos(angle)), rtol=1e-10)
    assert actual[0] > radius-np.linalg.norm(triangle.mean(1)[0,:2])


def test_projected_triangle_enclosing_axis_cannot_pass_from_good_vertices():
    theta = np.arange(3)*2*np.pi/3
    triangle = np.array([np.c_[np.cos(theta),np.sin(theta),np.zeros(3)]])
    np.testing.assert_allclose(finite_facet_deviation(triangle, 1.), [1.])
