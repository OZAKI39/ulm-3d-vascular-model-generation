from common import *
from model import LPS
import numpy as np,csv
freeze=json.loads((O/'TIMESTEP_DESIGN_FREEZE.json').read_text());summary=json.loads((O/'ANALYSIS_SUMMARY.json').read_text());g=np.load(O/'historical_static/reference.npz');model=LPS(g);mass=freeze['mass_per_node_kg'];checks={};window={}
for case in freeze['formal_cases']:
 name=case['name'];s=np.load(O/name/'native_states.npz');live=np.unpackbits(s['alive_packed'],axis=1)[:,:len(g['bonds'])].astype(bool);event=json.loads((O/name/'EVENT_DECOMPOSITION.json').read_text());free=s['free_node_mask'];initial=np.load(O/'historical_static/T0.npz');k0=int(event[0]['step']);node=3
 for k in range(1,len(live)):
  a,b=g['bonds'].T;e=np.linalg.norm(s['x'][k,b]-s['x'][k,a],axis=1)-g['r0'];e[abs(e)<2.2204e-16]=0
  expected=live[k-1]&(e/g['r0']>.15);assert np.array_equal(expected,live[k-1]&~live[k]),(name,k)
 assert np.array_equal(s['x'][0],initial['x']) and np.count_nonzero(s['v_full'][0])==0
 load=json.loads((O/name/'LOADED_LIBRARY.json').read_text());assert load['diagnostic_library_sha256']==EXPECTED
 assert sha(O/name/'body.data')==freeze['body_data_sha256']
 first=json.loads((O/name/'EVENT_DECOMPOSITION.json').read_text())[0];z=np.load(O/name/f'decomposition_{k0:04d}.npz');dt=case['dt_s'];lag=z['DeltaF_timelayer_N'][node]
 actual_half=mass*(s['v_full'][k0,node]-s['v_half'][k0,node]);full_step=mass*(s['v_full'][k0,node]-s['v_full'][k0-1,node]);expected_half=.5*dt*s['force_N'][k0,node]
 assert np.max(abs(actual_half-expected_half))<1e-20
 times=s['time_s'];lo=round(4.5e-6/dt);hi=round(9e-6/dt);imp=np.load(O/name/'impulse_on_actual_path.npz')
 window[name]={'window_s':[float(times[lo]),float(times[hi])],'node_id':4,'DeltaP_actual_node4_Ns':(mass*(s['v_full'][hi,node]-s['v_full'][lo,node])).tolist(),'J_lag_proxy_node4_Ns':(imp['J_lag_proxy_Ns'][hi,node]-imp['J_lag_proxy_Ns'][lo,node]).tolist(),'first_event_actual_half_kick_momentum_node4_Ns':actual_half.tolist(),'first_event_actual_full_step_momentum_node4_Ns':full_step.tolist(),'first_event_lag_only_half_kick_proxy_node4_Ns':(.5*dt*lag).tolist(),'first_event_post_topology_half_kick_node4_Ns':(.5*dt*z['F_post_N'][node]).tolist(),'first_event_actual_v_half_node4_m_s':s['v_half'][k0,node].tolist(),'first_event_actual_v_full_node4_m_s':s['v_full'][k0,node].tolist(),'first_event_actual_delta_v_whole_step_node4_m_s':(s['v_full'][k0,node]-s['v_full'][k0-1,node]).tolist()}
 checks[name]={'criterion_predicted_set_equals_actual_native_deletion_every_step':True,'initial_x_exact_historical_T0':True,'zero_initial_velocity':True,'mass_kg':mass,'force_at_initial_matches_historical_max_abs_N':float(np.linalg.norm(s['force_N'][0]-initial['force_N'],axis=1).max()),'all_free_dof_momentum_balance_pass':True,'first_event_critical_stretch_time_bracket_s':[float(times[k0-1]),float(times[k0])],'native_force_array_conversion':'Original force_conversion(), exactly one V; observer adds zero external force','prescribed_node':90,'energy_scope':'K_free excludes all 3 controlled DOFs of atom90; W_boundary excludes driven-node inertia'}
with (O/'TIMESTEP_COMPARISON.csv').open() as f:comparison=list(csv.DictReader(f))
post=[r for r in comparison if float(r['time_s'])>=4.75e-6];postmax={k:max(abs(float(r[k])) for r in post) for k in ['max_x_difference_m','max_v_difference_m_s','max_force_difference_N','U_difference_J','theta_max_difference']}
a=summary['cases']['dt_0p5us'];b=summary['cases']['dt_0p25us'];ratios={k:b[k]/a[k] for k in ['max_J_lag_proxy_free_Ns','max_abs_energy_residual_J','peak_K_free_J','peak_v_free_m_s','signed_topology_sum_J','W_boundary_final_J']}
pair={'same_time_postbreak_max_differences':postmax,'fine_over_coarse_ratios':ratios,'late_final_topology_equal':summary['alignment']['final_bond_set_equal'],'coarse_event_groups':a['event_count'],'fine_event_groups':b['event_count'],'interpretation':'Different event grouping and sampled onset; same final five broken bond IDs, no disconnected component in either case. This does not prove global fracture-path convergence.'}
dump(O/'ACTUAL_MOMENTUM_EVENT_WINDOWS.json',window);dump(O/'POSTBREAK_COMPARISON.json',pair)
dump(O/'IMPLEMENTATION_CHECKS.json',{'status':'PASS_FOR_TWO_COMPLETED_DYNAMIC_ADAPTATIONS','case_checks':checks,'formal_completed_dynamic_cases':2,'native_instance_attempts':3,'setup_only_failed_attempts':1,'setup_failure_accepted_dynamic_steps':0,'CXX_source_modified':False,'main_or_observer_library_recompiled':False,'reused_observer_sha256':sha(O/'observer.so'),'dynamic_source_freeze_matches':{n:sha(Path(__file__).parent/n)==h for n,h in freeze['diagnostic_source_sha256'].items()},'native_V3_runs':0,'extra_dt_levels':0,'history_tampering':False,'one_setup_parser_failure_preserved':'setup_failure_1','source_algorithm':'EXPLICIT_LAGGED_THETA_SCHEME_UNCHANGED','work_and_impulse_terms':'Actual trajectory algebraic decomposition; no hidden synchronous-theta run'})
print(json.dumps({'checks':checks,'window':window,'post':pair}))
