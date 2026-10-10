"""Strict P1 interpolation: verify the containing cell returned by a VTK locator.

VTK FindCell/probe can accept a neighbouring tetrahedron with small negative
parametric coordinates even when the probe tolerance is zero. This adapter
checks barycentric membership and searches intersecting cell bounds when needed.
No nearest-neighbour velocity or outside-domain extrapolation is used.
"""
import numpy as np
import vtk


class TetraP1Sampler:
    def __init__(self, grid, values):
        import pyvista as pv
        if set(grid.cells_dict) != {pv.CellType.TETRA}:
            raise ValueError("Strict tetra_p1 interpolation requires linear tetrahedra")
        self.cells = np.asarray(grid.cells_dict[pv.CellType.TETRA],int)
        points = np.asarray(grid.points,float)[self.cells]
        self.origins = points[:,0]
        self.inverse = np.linalg.inv(np.swapaxes(points[:,1:]-points[:,:1],1,2))
        self.values = np.asarray(values,float)[self.cells]
        self.locator = vtk.vtkStaticCellLocator()
        self.locator.SetDataSet(grid)
        self.locator.BuildLocator()
        self.roundoff = 1e-10
        self.bbox_epsilon = max(float(np.ptp(grid.points,axis=0).max())*1e-12,1e-15)
        self.fallback_queries = 0
        self.queries = 0

    def weights(self, positions, cells):
        local = np.einsum("nij,nj->ni",self.inverse[cells],positions-self.origins[cells])
        return np.column_stack([1-local.sum(axis=1),local])

    def sample(self, positions):
        x = np.asarray(positions,float)
        if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
            raise ValueError("Finite N by 3 SI coordinates required")
        self.queries += len(x)
        if not len(x):
            return np.empty((0,)+self.values.shape[2:])
        cells = np.array([self.locator.FindCell(p) for p in x],int)
        w = self.weights(x,np.maximum(cells,0))
        bad = (cells < 0) | (w.min(axis=1) < -self.roundoff) | (w.max(axis=1) > 1+self.roundoff)
        for i in np.flatnonzero(bad):
            self.fallback_queries += 1
            ids = vtk.vtkIdList()
            bounds = np.column_stack([x[i]-self.bbox_epsilon,x[i]+self.bbox_epsilon]).ravel()
            self.locator.FindCellsWithinBounds(bounds,ids)
            candidates = np.array([ids.GetId(j) for j in range(ids.GetNumberOfIds())],int)
            if len(candidates):
                weights = self.weights(x[i:i+1],candidates)
                valid = (weights.min(axis=1) >= -self.roundoff) & (weights.max(axis=1) <= 1+self.roundoff)
                if valid.any():
                    j = int(np.flatnonzero(valid)[0])
                    cells[i],w[i] = candidates[j],weights[j]
                    continue
            raise ValueError("Flow query outside VTK mesh; strict tetrahedral membership failed; no extrapolation")
        return np.einsum("ni,ni...->n...",w,self.values[cells])
