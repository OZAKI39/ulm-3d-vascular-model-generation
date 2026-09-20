#!/usr/bin/env python3
"""Replay actual SV1.2 fields with the same monitor used by live production."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13 import *
from sv_validation.provenance import write_json,sha256
from sv_validation.postprocess import SolutionMeasurements
check_reference();policy=frozen_policy();ref=load('reference_freeze')
measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',ref['accepted_solution']['Q_target_m3_s'],policy['Umean_m_s'])
monitor=SteadyStopMonitor(measure,policy);rows=[]
history=load('solver_history','sv1_2')
assert history['linear_failures']==history['nonlinear_failures']==0
assert all(r['converged'] for r in history['nonlinear'])
assert all(r['linear_converged'] and r['residual_finite'] for r in history['linear'])
for original in load('saved_state_qc','sv1_2')['states']:
    path=ROOT/original['path'];assert sha256(path)==original['sha256']
    u,p=measure.read(path);state=measure.measure(u,p)
    state.update(step=original['step'],time_s=original['time_s'],path=str(path.relative_to(ROOT)),sha256=sha256(path))
    stop=monitor.observe(state,u,history['linear_failures'],history['nonlinear_failures'])
    rows.append({'step':state['step'],'stop':stop,'epsilon_mass':state['epsilon_mass'],
                 'interval':monitor.intervals[-1] if monitor.intervals else None})
first=next((r['step'] for r in rows if r['stop']),None)
expected=ref['baseline_steady']['first_five_joint_intervals_end_step']
result={'status':'PASS' if first==expected else 'FAIL','reason':None if first==expected else 'STEADY_STOP_LOGIC_MISMATCH',
        'first_stop_step':first,'expected_first_stop_step':expected,'actual_field_count':len(rows),'rows':rows,
        'source':'independent reread of every actual SV1.2 saved VTU; same live monitor implementation'}
write_json(REPORT/'steady_stop_replay.json',result)
assert first==expected, 'STEADY_STOP_LOGIC_MISMATCH'
print(f'Replay PASS: first qualifying saved state {first}, actual fields {len(rows)}',flush=True)
