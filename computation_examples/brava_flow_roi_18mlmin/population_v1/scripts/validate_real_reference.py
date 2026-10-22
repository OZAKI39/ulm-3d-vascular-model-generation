"""Phase 2: old/new integration parity at actual BraVa birth and near-wall states."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--source-root',type=Path,required=True)
args=parser.parse_args()
root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(args.source_root/'particle_3d/src'))
sys.path.insert(0,str(args.source_root/'particle_3d/tests/population_v1'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.provenance import digest,write_json
from test_phase2 import particle,compare

active=json.loads((root/'reports/ACTIVE_FLOW.json').read_text())
boundaries={name:str(root/'mesh/SV_MESH/mesh-surfaces'/f'{name}.vtp')
            for name in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03','WALL']}
spec=dict(mesh=str(root/'mesh/SV_MESH/mesh-complete.mesh.vtu'),
          flow=str(root/active['case']/'frozen_flow/steady_flow.vtu'),
          boundaries=boundaries,mu_Pa_s=active['mu_Pa_s'],
          sonovue_root=str(args.source_root/'sonovue_size_distribution_v0'))
spec['sha256']=dict(mesh=active['mesh_sha256'],flow=active['flow_sha256'],
                    **{name:digest(path) for name,path in boundaries.items()})
env=load_environment(spec)
write_json(root/'population_v1/data/environment.json',spec)
cohort=json.loads((root/'microbubble/data/cohort.json').read_text())
event=cohort['events'][0];rows=[]
rows.append(dict(location='ACTUAL_FIRST_BIRTH',**compare(env,[particle(env,1,event['birth_center_m'],event['radius_m'])],40)))
saved=np.load(root/'microbubble/pilot/reference/mb_000001/trajectory.npz')['samples']
k=int(np.argmin(np.abs(saved[:,0]-3.)))
rows.append(dict(location='ACTUAL_NEAR_WALL_AT_PILOT_3S',original_age_s=float(saved[k,0]),
                **compare(env,[particle(env,1,saved[k,1:4],event['radius_m'])],20)))
record=dict(PASS=all(r['PASS'] for r in rows),phase=2,flow_sha256=env.flow_sha256,
            reference='Unmodified Particle9AStepper driven on the same physical-time partition',results=rows)
write_json(root/'population_v1/validation/phase2_real_flow.json',record)
print(json.dumps(record,indent=2),flush=True)
