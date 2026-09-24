from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pytest
from common import *

@pytest.fixture(scope='session')
def elevated():
    with np.load(CASE/'SV_MESH/tet10_arrays_si.npz') as a:return {k:a[k] for k in a.files}

@pytest.fixture
def synthetic():
    from p2 import EDGES
    corners=np.array([[1.,0,0],[0,1.,0],[0,0,1.],[0.,0,0]])
    return corners,np.vstack([corners,corners[EDGES].mean(axis=1)])
