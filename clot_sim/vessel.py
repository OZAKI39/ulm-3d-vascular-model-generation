"""Frozen native tetrahedral lumen; exact strict membership queries used by V4."""

from pathlib import Path
import numpy as np
import pyvista as pv
import vtk
import yaml
from .config import INPUTS


class VesselSources:
    def __init__(self):
        root = INPUTS / "vessel"
        self.config = yaml.safe_load((root / "bifurcation.yaml").read_text())
        self.mesh = pv.read(root / "mesh-complete.mesh.vtu")
        self.flow = pv.read(root / "result_121.vtu")
        self.faces = {
            name: pv.read(root / (name + ".vtp")).triangulate()
            for name in ["WALL", "INLET", "OUTLET_A", "OUTLET_B"]
        }
        self.tet = self.flow.cells_dict[pv.CellType.TETRA]
        vertices = np.asarray(self.flow.points[self.tet], float)
        self.origins = vertices[:, 0]
        self.inverse = np.linalg.inv(
            np.transpose(vertices[:, 1:] - vertices[:, :1], (0, 2, 1))
        )
        self.locator = vtk.vtkStaticCellLocator()
        self.locator.SetDataSet(self.flow)
        self.locator.BuildLocator()


def locator(surface):
    loc = vtk.vtkStaticCellLocator()
    loc.SetDataSet(surface)
    loc.BuildLocator()
    return loc


def hit(loc, a, b):
    t = vtk.mutable(0.0)
    x = [0.0, 0.0, 0.0]
    p = [0.0, 0.0, 0.0]
    sub = vtk.mutable(0)
    cell = vtk.mutable(0)
    ok = loc.IntersectWithLine(a, b, 1e-13, t, x, p, sub, cell)
    return bool(ok), np.array(x), int(cell)


class LumenSampler:

    def __init__(self, sources=None):
        self.sources = sources or VesselSources()
        self.queries = 0
        self.rejected = 0
        self.fallbacks = 0

    def locate(self, positions):
        s = self.sources
        X = np.atleast_2d(np.asarray(positions, float))
        self.queries += len(X)
        ids = np.array([s.locator.FindCell(p) for p in X], int)

        def bary(p, ii):
            w = np.einsum("nij,nj->ni", s.inverse[ii], p - s.origins[ii])
            return np.c_[1 - w.sum(1), w]

        weights = bary(X, np.maximum(ids, 0))
        bad = (ids < 0) | (weights.min(1) < -1e-10) | (weights.max(1) > 1 + 1e-10)
        for i in np.flatnonzero(bad):
            found = False
            candidate = vtk.vtkIdList()
            eps = 1e-13
            s.locator.FindCellsWithinBounds(
                np.c_[X[i] - eps, X[i] + eps].ravel(), candidate
            )
            cells = np.array(
                [candidate.GetId(j) for j in range(candidate.GetNumberOfIds())], int
            )
            if len(cells):
                w = bary(X[i : i + 1], cells)
                good = (w.min(1) >= -1e-10) & (w.max(1) <= 1 + 1e-10)
                if good.any():
                    j = np.flatnonzero(good)[0]
                    ids[i] = cells[j]
                    weights[i] = w[j]
                    found = True
                    self.fallbacks += 1
            if not found:
                ids[i] = -1
                weights[i] = np.nan
        self.rejected += int((ids < 0).sum())
        return (ids, weights)
