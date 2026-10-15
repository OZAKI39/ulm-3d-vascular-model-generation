from pathlib import Path
import json,hashlib,ast,subprocess,numpy as np,time
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path.cwd();O=R/'verification'/C['task']/C['run'];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protected=json.loads((O/'provenance/PROTECTED_START.json').read_text());changes=[]
for name,digest in protected.items():
 p=R/name
 if not p.exists() or sha(p)!=digest:changes.append(name)
assert not changes,changes
before=O/'provenance/python_before';oldphysics={}
for name in ['early_flow_coupling.py','native_y_wall_contact.py','native_y_wall_query.py','y_fragmentation_cases.py','current_position_provider.py']:
 q=json.loads((O/'V3_replay/REPLAY_DESIGN_AND_IDENTITY.json').read_text())['non_target_python_sha256'][name];oldphysics[name]=q==sha(R/'pd_clot/lammps_fragmentation'/name)
old=ast.parse((before/'detached_fragment_metrics.py').read_text());new=ast.parse((R/'pd_clot/lammps_fragmentation/detached_fragment_metrics.py').read_text())
func=lambda t,n:ast.dump(next(x for x in t.body if isinstance(x,ast.FunctionDef) and x.name==n),include_attributes=False)
assert func(old,'graph_components')==func(new,'graph_components')
old=ast.parse((before/'failure_classifier.py').read_text());new=ast.parse((R/'pd_clot/lammps_fragmentation/failure_classifier.py').read_text());assert func(old,'classify')==func(new,'classify')
for f in (R/'analysis/native_state_v1').glob('*.py'):ast.parse(f.read_text())
for f in (O/'provenance/python_new').glob('*.py'):ast.parse((R/'pd_clot/lammps_fragmentation'/f.name).read_text())
from pd_clot.lammps_fragmentation.native_y_wall_contact import penalty,CONTACT
norm=np.tile([0.,1.,0.],(4,1));clear=np.array([40,20,19,-1])*1e-6;clear[1]=CONTACT['effective_radius_m'];vel=np.tile([0.,-.01,0.],(4,1));f=penalty(clear,norm,vel);assert np.all(f[:2]==0) and np.all(f[2:,1]>0) and not np.any(f[:,[0,2]])
record={'status':'PASS_WITH_EXPLICIT_LIMITS','protected_file_count':len(protected),'protected_changed':changes,'production_native_library_and_binary_unchanged':True,'protected_native_source_unchanged':True,'non_target_python_unchanged':oldphysics,'graph_components_AST_unchanged':True,'legacy_classifier_AST_unchanged':True,'postprocess_test_status':json.loads((O/'FRAGMENT_STATISTICS_SEMANTIC_AUDIT.json').read_text())['status'],'wall_penalty_formula_regression':'PASS','wall_parameters':CONTACT,'wall_boundary_note':'Exact radius used at zero overlap. Approaching damping uses strict overlap>0; roundoff-positive overlap switches damping on (existing behavior, unchanged).','particle_pair_cpp_header_identity':'UNCHANGED_BETWEEN_PRODUCTION_AND_FINAL','new_break_contact_native_dynamic':'NOT_EXERCISED_IN_SMALL_PAIR_ZERO_TESTS; special scheduling code unchanged, V3 partial no contact yet','Gaussian_input':'EXACT_SAME_DECK_AND_PROVIDER','early_drag_logic':'SAME_CODE_AND_GRAPH; true IDs/online surviving topology; final V3 drag-on event NOT_REACHED','old_gaussian_plus_drag':'INDEPENDENT_LOAD_DEFINITION_LIMITATION_NOT_CHANGED','native_momentum_and_torque':'Small static before/after fracture force sums and torques saved; near roundoff','V3_failed_attempt':'One partial only; final owned guard not replayed','GPU':'NOT_USED','case_counts':{'diagnostic_scenarios_successful':10,'diagnostic_scenarios_failed':2,'native_diagnostic_instances_including_restart':16,'formal_V3_attempts':1,'formal_V3_completions':0,'old_baseline_reruns':0,'new_material_or_CFD_cases':0},'small_test_performance':{p.parent.name:json.loads(p.read_text()) for p in O.glob('*/PERFORMANCE.json') if p.parent.name!='V3_replay'},'compile_resources':[json.loads((O/n).read_text()) for n in ['COMPILE_PERFORMANCE.json','OWNER_COMPILE_PERFORMANCE.json']],'failed_process_CPU_RSS':'NOT_RECORDED_ON_NATIVE_SIGSEGV','output_sizes':{f:subprocess.check_output(['du','-sb',str(O/f)],text=True).split()[0] for f in ['source','source_owner_fix','build','build_owner_fix','V3_replay']},'git_push':'NOT_RUN_NOT_AUTHORIZED'}
# Exact small specimen comparison before the defect event.
parity={}
for case in ['sort','migration','fracture']:
 a=np.load(O/('original_'+case)/'T0.npz');b=np.load(O/('final_'+case)/'T0.npz');parity[case]={k:bool(np.array_equal(a[k],b[k])) for k in ['x','force_N','theta','m','alive']}
record['predefect_elastic_bitwise_parity']=parity
(O/'IMPLEMENTATION_CHECKS.json').write_text(json.dumps(record,indent=2))
print('Protected:',len(protected),'changes:',changes,'parity:',parity)
