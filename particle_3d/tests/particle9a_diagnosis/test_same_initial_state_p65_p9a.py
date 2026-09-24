from pathlib import Path
from copy import deepcopy
import json
import numpy as np
from particle_3d.particle81_simulation import environment, SavedTrajectoryStepper
from particle_3d.particle9a_motion import Particle9AStepper
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.lammps_neighbors import ValidationNeighborPolicy
from particle_3d.particle8_replay import canonical_hash

def test_all_twelve_actual_admissions_identical():
    root=Path(__file__).resolve().parents[2]
    events=json.loads((root/'reports/particle9a_2mmps_diagnosis/data/original_smoke_audit.json').read_text())['events']
    events=list(events.values()) if isinstance(events,dict) else events
    env=environment(); checker=FiniteSizeAdmission(None,None,wall=env.wall,field=env.field)
    def provider(i,s,t):
        f=env.field.sample(s.center_m); return f.velocity_m_s,.5*f.vorticity_s_inv
    for e in events:
        h=canonical_hash(e); p1=checker.check(deepcopy(e),np.array(e['birth_center_m']),{})[0]
        p2=checker.check(deepcopy(e),np.array(e['birth_center_m']),{})[0]
        query=ValidationNeighborPolicy(20e-6,.5e-6,'DIAG_SAME_INPUT')
        a=SavedTrajectoryStepper([p1],query,provider,env.mu,wall=env.wall,boundary_classifier=env.classifier)
        b=Particle9AStepper([p2],query,provider,env.mu,wall=env.wall,boundary_classifier=env.classifier,
                          gradient_provider=lambda x:env.field.sample(x).velocity_gradient_s_inv)
        assert np.array_equal(a.samples,b.samples)
        assert canonical_hash(e)==h
    assert len(events)==12
