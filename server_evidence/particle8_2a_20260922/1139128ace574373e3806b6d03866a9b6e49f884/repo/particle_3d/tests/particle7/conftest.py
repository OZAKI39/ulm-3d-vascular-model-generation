from pathlib import Path
import sys,json
import pytest
PACKAGE=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(PACKAGE/'src'))
@pytest.fixture(scope='session')
def repo(): return PACKAGE.parent
@pytest.fixture(scope='session')
def report(): return PACKAGE/'reports/particle7'
@pytest.fixture(scope='session')
def contract(): return json.loads((PACKAGE/'contracts/PARTICLE7_INLET_POPULATION_V0.json').read_text())
@pytest.fixture
def source():
 from particle_3d.injection_population import PopulationSource
 from particle_3d.particle7_cases import SONOVUE
 return PopulationSource(SONOVUE)
@pytest.fixture
def scheduler(source):
 from particle_3d.injection_population import InjectionScheduler,LinearProfile
 return InjectionScheduler(source,LinearProfile([0],[1e-12]))
@pytest.fixture
def engine():
 from particle_3d.particle7_cases import synthetic_engine
 e=synthetic_engine()
 yield e
 e.bridge.close()
@pytest.fixture(scope='session')
def real():
 from particle_3d.particle7_cases import real_inlet
 return real_inlet()
