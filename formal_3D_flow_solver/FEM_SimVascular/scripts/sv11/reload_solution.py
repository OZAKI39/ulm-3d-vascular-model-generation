#!/usr/bin/env python3
"""Fresh process independently reopens the selected new-stage VTU."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv11 import load,REPORT
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.provenance import sha256,write_json
name=sys.argv[1] if len(sys.argv)>1 else 'petsc_short'
state=load(name+'_qc')['states'][-1]
policy=json.loads((ROOT/'configs/time_policy.json').read_text())
measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',state['Q_target_m3_s'],policy['Umean_m_s'])
path=ROOT/state['path'];u,p=measure.read(path);again=measure.measure(u,p)
keys=['velocity_L2','pressure_range_pa','Q_in_m3_s','outlet_flows_m3_s','Q_out_total_m3_s','epsilon_mass','epsilon_Q',
      'wall_velocity_max_m_s','wall_velocity_P95_m_s','wall_noslip_pass','velocity_finite','pressure_finite']
equal={key:again[key]==state[key] for key in keys}
write_json(REPORT/'solution_reload.json',{'status':'PASS' if all(equal.values()) and sha256(path)==state['sha256'] else 'FAIL',
           'fresh_process':True,'artifact_scope':'accepted solution' if (REPORT/'accepted_solution.json').exists() else 'unaccepted transient diagnostic',
           'path':state['path'],'sha256':sha256(path),'equal':equal,'recomputed':{key:again[key] for key in keys}})
print('Reload',all(equal.values()))
