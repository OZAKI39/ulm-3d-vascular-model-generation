"""Read-only scoped identity review and one-case preregistration. No simulation."""
from common import *
import shutil,subprocess,datetime,platform
O.mkdir(parents=True,exist_ok=False);F.mkdir(parents=True,exist_ok=False)
(O/'provenance').mkdir();(O/'keyframes').mkdir();(O/'source_snapshot').mkdir()
oldmanifest=json.loads((OLD/'PAYLOAD_MANIFEST.json').read_text())
protected={str(R/p):d['sha256'] for p,d in oldmanifest['files'].items()}
extra=[LIB,R/'analysis/softening_feasibility_v1/model.py',Path('/workspace/y_clot/env.sh'),R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so',PREV/'build/liblammps.so',R/'verification/lps_native_state_identity_and_fracture_bugfix_v1/20261011T203015Z/build_owner_fix/liblammps.so',OLD/'PAYLOAD_MANIFEST.json']
identity=json.loads((OLD/'FINAL_LIBRARY_IDENTITY.json').read_text())
for n,h in identity['source_files'].items():
 for p in [R/'analysis/dynamic_softening_v1'/n,OLD/'source/src/BPM'/n]:
  assert sha(p)==h,(p,'source mismatch');extra.append(p)
for n,h in json.loads((OLD/'DYNAMIC_DAMAGE_COMMIT_DESIGN_FREEZE.json').read_text())['source_sha256'].items():
 if 'variational' not in n:
  p=OLD/'source/src'/n;assert sha(p)==h;extra.append(p)
for p in extra:protected[str(p)]=sha(p)
assert all(sha(p)==h for p,h in protected.items())
assert sha(LIB)=='866065394c48ebc8a6d32f147bde992a325c2a503a49d50392243eee10fea9bd'
assert sha(R/'analysis/softening_feasibility_v1/model.py')=='3eb1396b376cf30e40586c26ba649909ab9382d9eac9607666359258273acf74'
dump(O/'PROTECTED_FILES_BEFORE.json',dict(scope='Prior sealed 195-file payload plus directly used source/library/env assets; NO recursive 31001-file rescan',files=protected))
gitdiff=subprocess.check_output(['git','diff','--binary'],cwd=R);(O/'tracked_diff.before.patch').write_bytes(gitdiff)
(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md').write_bytes((R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md').read_bytes())
inventory=[]
for p in sorted((R/'analysis/dynamic_softening_v1').iterdir()):
 if p.is_file():inventory.append(dict(path=str(p),sha256=sha(p),disposition='RETAIN_SEALED_REPRODUCTION_EVIDENCE'))
dump(O/'ACTIVE_ENTRYPOINT_AUDIT.json',dict(status='PASS_MINIMAL_SCOPED_REVIEW',candidates=inventory,safe_deletions=0,deleted_files=[],reason='All prior dynamic files are sealed source or replay/audit evidence. No proven disposable temporary source.',historical_recursive_scan=False,prior_manifest_sha256=sha(OLD/'PAYLOAD_MANIFEST.json')))
for n in ['reference.npz','reference.native','body.data','original_bond_ids.npy','341_boundary_ids.npz']:
 shutil.copy2(OLD/'provenance'/n,O/'provenance'/n)
for n in identity['source_files']:shutil.copy2(OLD/'source/src/BPM'/n,O/'source_snapshot'/n)
shutil.copy2(R/'analysis/softening_feasibility_v1/model.py',O/'source_snapshot/model.py')
old=json.loads((OLD/'DYNAMIC_DAMAGE_COMMIT_DESIGN_FREEZE.json').read_text())
ref=np.load(O/'provenance/reference.npz');ids=np.load(O/'provenance/341_boundary_ids.npz')
assert len(ref['X'])==341 and len(ref['bonds'])==12047
for name in ['left_ids','right_ids','free_ids']:assert np.array_equal(ids[name],old['boundary'][name])
assert np.array_equal(np.load(O/'provenance/original_bond_ids.npy'),np.arange(1,12048))
direct_evidence=['FINAL_STATUS.json','DYNAMIC_341_NODE_AUDIT.json','DYNAMIC_ENERGY_LEDGER.json','DT_HALF_WINDOW_COMPARISON.json','FIRST_SOFTENING_AND_FAILURE_EVENTS.json','STOP_STATE_TOPOLOGY_GRADIENT_AUDIT.json','RUN0_AND_TRIAL_STATE_SAFETY.json','INSTRUMENTATION_RECOVERY.json','DYNAMIC_INTEGRATOR_IMPLEMENTATION_AUDIT.md','REPRODUCE.txt']
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',dict(status='PASS',library=str(LIB),library_sha256=sha(LIB),source_files=identity['source_files'],source_root=str(OLD/'source'),model_sha256=sha(R/'analysis/softening_feasibility_v1/model.py'),CPP_changed=False,compiled=False,prior_certificates={n:sha(OLD/n) for n in direct_evidence},native_wrapper_sha256=sha(A/'native.py'),wrapper_identical_to_prior=sha(A/'native.py')==sha(R/'analysis/dynamic_softening_v1/native.py')))
freeze=dict(task=TASK,run_id=RUN,frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),library=str(LIB),library_sha256=sha(LIB),source_identity_sha256=sha(O/'SOURCE_AND_LIBRARY_IDENTITY.json'),reference_sha256={n:sha(O/'provenance'/n) for n in ['reference.npz','reference.native','body.data','original_bond_ids.npy','341_boundary_ids.npz']},original_bonds=ref['bonds']+1,original_bond_ids=np.arange(1,12048),material=old['material'],boundary=old['boundary'],actual_analytic_expression='ql=ternary(time<T,time/T,1); qu=ternary(time<T+H,0,ternary(time<2*T+H,(time-T-H)/T,1)); u=A/2*(cos(pi*qu)-cos(pi*ql)); v=A*pi/(2*T)*(sin(pi*ql)-sin(pi*qu)); left=-u,-v; right=+u,+v; grip y,z fixed',dt_s=.25e-6,steps=2800,initial_state='Original body.data, zero velocity, kappa=d=D=0, no restart/import',dynamic_hook=old['dynamic_hook'],time_layer=old['time_layer'],gates=old['gates'],contact=old['contact'],contact_path_guard='Additional geometric safety guard: stop if active pair reaches r/dx<0.10 (near-coincidence), or relative pair step exceeds inherited 0.10dx. This is a conservative numerical guard, not a calibrated admissible physical overlap. Record maximum penetration without interpreting shallow overlap as a hard-wall violation.',ledger=old['ledger'],normalization='running max over [1e-15 J, abs(W), abs(U-U0)+abs(K-K0)+(D-D0)+abs(Uc-Uc0)]; exact prior definition, never reset',external_work='sum over actual steps: -0.5*(F_grip_old+F_grip_new) dot dx_grip + exact DeltaK_grip',D='At identical new x: U(x_new,d_old)-U(x_new,d_new); signed_topology=-D for display only',independent_audit='Every actual step: frozen model U,force,kappa,d,D; native VV equations, momentum, IDs, contact; explicit VV/spatial/switch/contact defect decomposition',save_strategy='Every step scalar CSV, FP64 x/v/physical total and material forces once in compressed trajectory. Sparse changes of d save original ID,kappa,d; all kappa can be reconstructed exactly from x and frozen model and is checked online. Full keyframes at initial/onset/first failure/new record maxDeltaD>=0.08/phase transitions/contact activation/stop or end. Stop also saves previous full state and restart. No duplicated full time-by-all-bond matrices.',budget=dict(main_cases=1,max_steps=2800,main_wall_s=300,MPI_ranks=1,threads=1,RSS_stop_GiB=4,payload_limit_MiB=600,compile_calls=0,additional_dynamic_cases=0),stop_handling='First violation stops before another run 1; persist violating state, previous state, event IDs, energy and exact source SHA. No dt/4 continuation. No P4 specimen.',limitations=['No global time convergence certificate','No natural fracture area/Gc/macroscopic strength certificate','No active dynamic contact PASS unless observed'],driver_sha256={p.name:sha(p) for p in A.glob('*.py')})
dump(O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json',freeze)
(O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.sha256').write_text(sha(O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json')+'\n')
env={p:Path(p).read_text().strip() for p in ['/sys/fs/cgroup/cpu/cpu.cfs_quota_us','/sys/fs/cgroup/cpu/cpu.cfs_period_us','/sys/fs/cgroup/memory/memory.limit_in_bytes','/sys/fs/cgroup/memory/memory.usage_in_bytes']}
dump(O/'ENVIRONMENT_BEFORE.json',dict(python=platform.python_version(),numpy=np.__version__,cgroup=env,disk=shutil.disk_usage(R)._asdict(),git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),processes=subprocess.check_output(['ps','-eo','pid,pcpu,rss,comm','--sort=-pcpu'],text=True),GPU_used=False))
print('FROZEN',O,sha(O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json'),flush=True)
