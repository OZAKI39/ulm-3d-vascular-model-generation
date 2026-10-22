import numpy as np
import pyvista as pv
import pytest
from mevo_mb.flow.fields import VTKFlowField


def test_strict_p1_rejects_wrong_neighbour_and_outside(tmp_path):
    # Continuous P1 field with different gradients on either side of z=0.
    points = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]],np.float32)*.001
    cells = np.array([[4,0,1,2,3],[4,0,2,1,4]])
    grid = pv.UnstructuredGrid(cells,np.full(2,pv.CellType.TETRA),points)
    grid["Velocity"] = np.column_stack([points[:,0]+2*points[:,1]+3*np.abs(points[:,2]),
                                        np.zeros(5),np.zeros(5)])
    path = tmp_path / "piecewise.vtu"
    grid.save(path)
    flow = VTKFlowField(path,spatial_method="tetra_p1")
    sampler = flow.fields[0].sampler
    original = sampler.locator

    class WrongNeighbour:
        def FindCell(self,point):
            return 0

        def FindCellsWithinBounds(self,bounds,ids):
            return original.FindCellsWithinBounds(bounds,ids)

    sampler.locator = WrongNeighbour()
    xyz = np.array([[.0002,.0002,-1e-7],[.0002,.0002,.0001],[0,0,-.001]])
    expected = xyz[:,0]+2*xyz[:,1]+3*np.abs(xyz[:,2])
    np.testing.assert_allclose(flow.velocity(xyz,0.)[:,0],expected,rtol=1e-7,atol=1e-15)
    assert sampler.fallback_queries == 2
    with pytest.raises(ValueError,match="outside VTK mesh"):
        flow.velocity([[.002,.002,.002]],0.)
