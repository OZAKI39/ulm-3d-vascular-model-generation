"""Reject geometry changes while accepting an exact float32 writer permutation."""
from pathlib import Path
import sys
import numpy as np
import pyvista as pv
import pytest

ROOT=Path(__file__).resolve().parents[1]
from human_vessel_flow.reference_pipe.pipe_measurements import serialized_node_mapping


@pytest.fixture
def mesh():
    p=np.array([[.000123456789,0,0],[.000523456789,0,0],
                [.000123456789,.0003,0],[.000123456789,0,.0003],
                [.000523456789,.0003,.0003]])
    t=np.array([[0,1,2,3],[1,2,3,4]])
    permutation=np.array([3,1,4,0,2]);inverse=np.argsort(permutation)
    grid=pv.UnstructuredGrid(np.c_[np.full(2,4),inverse[t]].ravel(),
                             np.full(2,10,np.uint8),p[permutation].astype(np.float32))
    return p,t,grid,inverse


def test_float32_permutation_preserves_node_identity(mesh):
    p,t,grid,inverse=mesh
    assert np.array_equal(serialized_node_mapping(p,t,grid),inverse)


def test_one_ulp_coordinate_change_is_rejected(mesh):
    p,t,grid,_=mesh
    grid.points[0,0]=np.nextafter(grid.points[0,0],np.float32(np.inf))
    with pytest.raises(ValueError,match='exact mesh serialization'):
        serialized_node_mapping(p,t,grid)


def test_changed_tetra_connectivity_is_rejected(mesh):
    p,_,grid,_=mesh
    with pytest.raises(ValueError,match='connectivity'):
        serialized_node_mapping(p,np.array([[0,1,2,3],[0,2,3,4]]),grid)


def test_rounding_collision_is_rejected(mesh):
    p,t,grid,_=mesh
    p[4]=p[0]+np.array([1e-16,0,0])
    with pytest.raises(ValueError,match='collision'):
        serialized_node_mapping(p,t,grid)
