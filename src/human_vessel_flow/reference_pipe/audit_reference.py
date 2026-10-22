"""Exact finite-facet radial error, retained from Gmsh V2 QA."""
import numpy as np

def finite_facet_deviation(xyz, radius):
    """Exact max radial deviation on each planar triangle, not centroid sagitta."""
    xy = xyz[..., :2]
    minima = np.full(len(xy), np.inf)
    for i, j in [(0, 1), (1, 2), (2, 0)]:
        a, delta = xy[:, i], xy[:, j]-xy[:, i]
        denom = np.sum(delta**2, axis=1)
        t = np.divide(-np.sum(a*delta, axis=1), denom, out=np.zeros(len(a)), where=denom > 0)
        minima = np.minimum(minima, np.linalg.norm(a + np.clip(t, 0, 1)[:, None]*delta, axis=1))
    cross = lambda a, b: a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]
    signs = np.stack([cross(xy[:, j]-xy[:, i], -xy[:, i]) for i,j in [(0,1),(1,2),(2,0)]], axis=1)
    inside = (np.all(signs > 0, axis=1) | np.all(signs < 0, axis=1))
    minima[inside] = 0
    return np.maximum(np.abs(radius-minima), np.max(np.abs(np.linalg.norm(xy, axis=2)-radius), axis=1))
