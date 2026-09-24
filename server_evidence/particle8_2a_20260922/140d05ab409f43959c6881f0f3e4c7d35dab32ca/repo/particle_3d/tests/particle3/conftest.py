from pathlib import Path
import sys
import pytest
import numpy as np
PACKAGE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PACKAGE/'src'))


@pytest.fixture(scope='session')
def p3_repo():return PACKAGE.parent


@pytest.fixture(scope='session')
def plane_triangle():return np.array([[-20,-20,0],[0,20,0],[20,-20,0]],dtype=float)*1e-6


@pytest.fixture(scope='session')
def plane_wall(plane_triangle):
    from particle_3d.wall_geometry import WallGeometry
    return WallGeometry(plane_triangle[None,:,:])


@pytest.fixture(scope='session')
def real_wall(p3_repo):
    from particle_3d.wall_geometry import WallGeometry
    return WallGeometry.from_frozen(p3_repo/'formal_3D_flow_solver/FEM_SimVascular')


@pytest.fixture(scope='session')
def p3_population():
    from particle_3d.rbc_distribution import sample_rbc_geometries
    return sample_rbc_geometries(100000,2026092002)


@pytest.fixture(scope='session')
def p3_geometries(p3_population):
    from particle_3d.rbc import RBCGeometry
    from particle_3d.rbc_distribution import quantile_indices
    return [RBCGeometry.from_population(p3_population,i) for i in quantile_indices(p3_population.samples)]
