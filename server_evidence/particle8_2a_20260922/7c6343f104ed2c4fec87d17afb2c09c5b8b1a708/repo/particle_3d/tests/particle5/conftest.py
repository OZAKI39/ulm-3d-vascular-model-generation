from pathlib import Path
import sys
import pytest
PACKAGE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PACKAGE/'src'))

@pytest.fixture(scope='session')
def repo():
    return PACKAGE.parent

@pytest.fixture(scope='session')
def mu(repo):
    from particle_3d.hydrodynamic_resistance import viscosity_from_frozen
    return viscosity_from_frozen(repo/'formal_3D_flow_solver/FEM_SimVascular')[0]
