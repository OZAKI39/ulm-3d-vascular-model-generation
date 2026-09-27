from pathlib import Path
import sys,json
import pytest
PACKAGE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PACKAGE/'src'))
@pytest.fixture(scope='session')
def repo():return PACKAGE.parent
@pytest.fixture(scope='session')
def mu(repo):
    from particle_3d.hydrodynamic_resistance import viscosity_from_frozen
    return viscosity_from_frozen(repo/'formal_3D_flow_solver/FEM_SimVascular')[0]
@pytest.fixture(scope='session')
def provenance(repo):
    from particle_3d.particle6_checkpoint import frozen_provenance
    return frozen_provenance(repo/'formal_3D_flow_solver/FEM_SimVascular')
@pytest.fixture(scope='session')
def parity_cases(mu):
    from particle_3d.particle6_validation import run_parity
    return [run_parity(mu,mixed) for mixed in [False,True]]
@pytest.fixture(scope='session')
def restart_cases(mu,repo,provenance,tmp_path_factory):
    from particle_3d.particle6_validation import restart_case
    base=tmp_path_factory.mktemp('actual_lammps_restart')
    return [restart_case(mu,base/str(mixed),repo,provenance,mixed) for mixed in [False,True]]
@pytest.fixture(scope='session')
def real_result(mu,repo):
    from particle_3d.particle6_validation import real_case
    return real_case(repo,mu)
