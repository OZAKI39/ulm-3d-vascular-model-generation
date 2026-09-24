from pathlib import Path
import sys,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from particle_3d.particle81_simulation import OUTPUT
from particle_3d.particle81_replay import Scene
@pytest.fixture(scope='session')
def root():return OUTPUT
@pytest.fixture(scope='session')
def scene():return Scene()
