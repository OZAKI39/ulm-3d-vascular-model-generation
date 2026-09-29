"""Actual-flow complete integrator parity, including a saved near-wall state."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import torch
ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);args=ap.parse_args()
sys.path.insert(0,str(args.source_root/'particle_3d/src'));sys.path.insert(0,str(args.source_root/'particle_3d/tests/population_v1'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
from particle_3d.population_v1.gpu_field import GPUField
from particle_3d.population_v1.gpu_wall import GPUWall
from particle_3d.population_v1.provenance import write_json
from test_phase2 import particle,event,reference
root=Path(__file__).resolve().parents[1];base=root.parent
env=load_environment(json.loads((root/'data/environment.json').read_text()))
gf=GPUField(env.field,cache_path=root/'data/tetra_acceleration.npz');gw=GPUWall(env.wall)
cohort=json.loads((base/'microbubble/data/cohort.json').read_text())['events'][0]
saved=np.load(base/'microbubble/pilot/reference/mb_000001/trajectory.npz')['samples'];k=int(np.argmin(abs(saved[:,0]-3.)))
rows=[]
for location,position,steps in [('ACTUAL_FIRST_BIRTH',cohort['birth_center_m'],400),('ACTUAL_NEAR_WALL_AT_PILOT_3S',saved[k,1:4],200)]:
    p=particle(env,1,position,cohort['radius_m']);old=reference(env,p)
    gpu=GPUPopulationEngine(env,pairs=True,gpu_field=gf,gpu_wall=gw);gpu.insert_accepted(p,event(p))
    errors=dict(position_m=0.,velocity_m_s=0.,omega_s_inv=0.,wall_gap_m=0.);started=time.perf_counter()
    for step in range(1,steps+1):
        t=step*.0005;old.step_to(t);gpu.advance(t);ref=old.read()[0]
        for name,left,right in [('position_m',gpu.pool.x[0].cpu().numpy(),ref.position),('velocity_m_s',gpu.pool.v[0].cpu().numpy(),ref.velocity),
            ('omega_s_inv',gpu.pool.omega[0].cpu().numpy(),ref.omega),('wall_gap_m',float(gpu.pool.wall_gap[0]),old.samples[-1][14])]:
            errors[name]=max(errors[name],float(np.max(abs(np.asarray(left)-right))))
        if old.boundary_event.startswith('OUTLET_'):break
    passed=errors['position_m']<8*env.wall.roundoff_m and errors['velocity_m_s']<1e-9 and errors['omega_s_inv']<1e-3
    row=dict(PASS=bool(passed),location=location,steps=step,errors=errors,wall_roundoff_m=env.wall.roundoff_m,
        runtime_s=time.perf_counter()-started,backend_counts=dict(gpu.backend_counts),reference_terminal=old.boundary_event,
        gpu_terminals=gpu.terminals)
    rows.append(row);print(json.dumps(row,indent=2),flush=True)
record=dict(PASS=all(r['PASS'] for r in rows),phase=7,flow_sha256=env.flow_sha256,results=rows,
    reference='UNMODIFIED_PARTICLE9A_STEPPER_AT_IDENTICAL_NOMINAL_PHYSICAL_TIME_GRID',
    gpu_role='ACTUAL_STATE_UPDATES_AND_ISOLATED_RESISTANCE_CONTACT_NOT_RENDERING_ONLY')
write_json(root/'validation/phase7_real_gpu_dynamics.json',record)
if not record['PASS']:raise SystemExit(1)
