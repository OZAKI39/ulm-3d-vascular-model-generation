#!/usr/bin/env python3
"""Validate an exited CPU run; never starts or reruns the native solver."""
import json,sys
from pathlib import Path
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13 import *
from sv_validation.sv11 import parse_solver_log,linear_gate,nonlinear_gate,parse_boundary
from sv_validation.sv12 import checkpoint_audit
from sv_validation.provenance import write_json
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.validation import integral_agreement
check_reference();policy=frozen_policy();ref=load('reference_freeze');execution=load('cpu_execution')
require(load('running')['status']=='EXITED','Native solver must already be exited')
case=OUTPUT/'cpu_early_stop';failure=execution['failure'];code=execution['exit_code']
stop_request=execution['stop_request'];wall=execution['GNU_wall_s'];elapsed=execution['elapsed_monotonic_s']
initial_step=execution['initial_step'];rows=execution['linear'];nl={r['step']:r for r in rows}
history=parse_solver_log((ROOT/execution['log']).read_text(),policy['dt_s'])
require(len(history['linear_solves'])==len(rows),'Incomplete linear history')
qc=load('cpu_saved_states');monitor=SimpleNamespace(states=qc['states'],intervals=qc['intervals'])
measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',ref['accepted_solution']['Q_target_m3_s'],policy['Umean_m_s'])
try:
    require(failure is None,'EARLY_STOP_INVALID: live failure')
    linear_gate({'exit_code':code,'history':history});nonlinear_gate(history)
    require(stop_request and stop_request['kind']=='STEADY','EARLY_STOP_INVALID: no automatic steady stop')
    final=monitor.states[-1];require(final['step']==max(nl)==stop_request['requested_final_step'],'Final field/checkpoint does not match safe stop')
    require(stop_gate(monitor.intervals,final,0,0,policy),'EARLY_STOP_INVALID: final gates')
    cp=checkpoint_audit(case/'4-procs/stFile_last.bin',final['step'],policy['dt_s'])
    write_json(REPORT/'cpu_final_checkpoint.json',cp)
    eq=equivalence(measure,ROOT/final['path'],ROOT/ref['accepted_solution']['path'],policy)
    write_json(REPORT/'cpu_reference_equivalence.json',eq);require(eq['status']=='PASS','EARLY_STOP_INVALID: field equivalence')
    native=parse_boundary(case/'4-procs/B_NS_Velocity_flux.txt',measure.Q)
    last=next(r for r in native if r['step']==final['step'])
    final['native_flux_difference_over_Q']=integral_agreement({r:last['signed_outward_flows_m3_s'][r] for r in final['signed_outward_boundary_flows_m3_s']},final['signed_outward_boundary_flows_m3_s'],measure.Q)
    final['accepted']=True;write_json(REPORT/'cpu_accepted_solution.json',final)
    baseline=ref['baseline_resources'];baseline_wall=baseline['total_solver_wall_s']
    inherited=baseline['inherited_first_ten_python_wall_s']
    performance={'status':'PASS','reference_measured_wall_s':baseline_wall,'candidate_measured_wall_s':wall,
       'measured_continuation_speedup':baseline_wall/wall,'saved_measured_seconds':baseline_wall-wall,
       'same_initial_checkpoint':True,'reference_new_steps':ref['accepted_solution']['step']-initial_step,'candidate_new_steps':len(nl),
       'reference_endpoint':ref['accepted_solution']['step'],'candidate_endpoint':final['step'],
       'inherited_first10_wall_s':inherited,'reference_accounted_total_s':baseline['total_python_monotonic_s']+inherited,
       'candidate_accounted_total_s':elapsed+inherited,'accounted_total_speedup':(baseline['total_python_monotonic_s']+inherited)/(elapsed+inherited),
       'timing_scope':'fresh same-machine native continuations from identical step10; total includes historical seed cost, not a newly measured t=0 run'}
    write_json(REPORT/'cpu_performance.json',performance)
    write_json(REPORT/'cpu_validation.json',{'status':'PASS','reason':None,'stop_step':final['step'],
       'first_qualifying_step':stop_request['detail']['qualifying_saved_step'],'science_equivalent':True,'checkpoint_complete':True,
       'linear_failures':0,'nonlinear_failures':0,'mass_pass':True,'steady_pass':True,'native_exit':code})
    print(json.dumps(performance,indent=2),flush=True)
except Exception as exc:
    write_json(REPORT/'cpu_validation.json',{'status':'FAIL','reason':'EARLY_STOP_INVALID','detail':str(exc)})
    raise
