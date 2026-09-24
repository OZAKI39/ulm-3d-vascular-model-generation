from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from particle_3d.particle8_replay import read
@pytest.fixture(scope='session')
def output(): return ROOT/'reports/particle8'
@pytest.fixture(scope='session')
def scenes(output):
    return {k:read(output/'data'/(k+'_scene.json')) for k in ['synthetic','restarted','real_mixed','real_single_mb']}
