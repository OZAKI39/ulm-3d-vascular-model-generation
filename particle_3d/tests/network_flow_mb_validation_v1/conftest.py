from pathlib import Path
import json,sys
import numpy as np
import pytest
R=Path(__file__).resolve().parents[2]/"reports/network_derived_flow_mb_validation_v1"
sys.path.insert(0,str(R/"scripts"))
sys.path.insert(0,str(R.parents[1]/"src"))
def read(p):return json.loads(Path(p).read_text())
@pytest.fixture(scope="session")
def report():return R
@pytest.fixture(scope="session")
def events():return read(R/"data/paired_events.json")
@pytest.fixture(scope="session")
def newenv():
 from runner import make_environment
 return make_environment(R/"server_bundle","NEW")
@pytest.fixture(scope="session")
def oldenv():
 from runner import make_environment
 return make_environment(R/"server_bundle","OLD")
@pytest.fixture(scope="session")
def newpaths(events):
 return [(e,np.load(R/"outputs/NEW/trajectories"/f"mb_{e['particle_id']:06d}.npz")["samples"]) for e in events]
