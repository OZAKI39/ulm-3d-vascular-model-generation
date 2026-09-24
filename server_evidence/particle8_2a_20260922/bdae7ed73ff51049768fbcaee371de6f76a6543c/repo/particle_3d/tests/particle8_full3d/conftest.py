from pathlib import Path
import sys,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from particle_3d.particle8_full3d_data import OUTPUT,scene_from_name
@pytest.fixture(scope='session')
def root():return OUTPUT
@pytest.fixture(scope='session')
def scenes(root):return {n:scene_from_name(root,n) for n in ['real_single_mb_entry','real_mixed_inlet_smoke','synthetic_open_section','synthetic_restarted']}
