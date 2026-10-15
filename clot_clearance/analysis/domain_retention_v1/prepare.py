from common import *
import shutil,difflib
assert not O.exists(),'Never overwrite previous results'
for folder in [O,F,O/'provenance',O/'preflight',O/'fault_test',O/'full',O/'source_identity'] :folder.mkdir(parents=True,exist_ok=True)
previous=json.loads((PREV/'PAYLOAD_MANIFEST.json').read_text());bad=[]
for name,m in previous['files'].items():
 if sha(R/name)!=m['sha256']:bad.append(name)
assert not bad,bad
shutil.copy2(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
for name in ['body.data','reference.npz','reference.native','original_bond_ids.npy','341_boundary_ids.npz']:shutil.copy2(PREV/'provenance'/name,O/'provenance'/name)
identity=json.loads((PREV/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());assert sha(LIB)==identity['new_sha256']=='08ccf24110e2cf2d07dbc8c0f31c0f368368b19b5427e6e5de9cc108a5ae115d'
for name in ['bond_bpm_peri_variational.cpp','bond_bpm_peri_variational.h','fix_lps_variational_dynamic.cpp','fix_lps_variational_dynamic.h','fix_move.cpp','fix_nve.cpp','verlet.cpp','update.cpp','event_state_audit.cpp','model.py']:
 shutil.copy2(PREV/'source_identity'/name,O/'source_identity'/name)
for name in ['change_box.cpp','read_restart.cpp']:shutil.copy2(PREV/'source/src'/name,O/'source_identity'/name)
z=np.load(PREV/'full/trajectory.npz');cols=list(z['columns']);arr=z['series'];speed=float(arr[:,cols.index('max_speed_m_s')].max());end=float(arr[-1,cols.index('time_s')]);keys=list((PREV/'full/keyframes').glob('*.npz'));x=np.concatenate([np.load(p)['x'] for p in keys]);saved_abs=float(abs(x).max());remaining=.0007-end;padding=2*DX;required=saved_abs+speed*remaining+padding;rounded=np.ceil(required/DX)*DX
assert abs(BOX-rounded)<1e-15
dump(O/'DOMAIN_CHANGE_AUDIT.json',dict(status='FROZEN_PENDING_INITIAL_PARITY',old_bounds_m=[[-250e-6,250e-6]]*3,new_bounds_m=[[-BOX,BOX]]*3,boundary='f f f',remap=False,numerical_container_not_physical_wall=True,saved_keyframe_count=len(keys),saved_xyz_min_m=x.min(0),saved_xyz_max_m=x.max(0),scope='All 18 saved full states; prior full per-step coordinates not available. Scalar max speed is over every accepted step.',observed_max_speed_m_s=speed,remaining_s=remaining,estimated_additional_travel_m=speed*remaining,blank_padding_m=padding,required_halfwidth_m=required,rounded_halfwidth_m=BOX,not_a_rigorous_future_velocity_bound=True,online_guard='Check candidate native-VV drift before every trial; checkpoint on clearance<dx; stop if predicted outside frozen box',change_box_command=f'change_box all x final {-BOX:.17g} {BOX:.17g} y final {-BOX:.17g} {BOX:.17g} z final {-BOX:.17g} {BOX:.17g} units box',restart_timing='Fresh: after read_data, before fixes. Restart: restore same-ID fixes and run0 before change_box only if stored box differs. New rolling/rejection restarts already contain new box.',source_checked='source_identity/change_box.cpp',old_atom73_cause=json.loads((PREV/'NATIVE_ATOM_LOSS_DIAGNOSIS.json').read_text())))
protected={str(p):sha(p) for p in [LIB,R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so',Path('/workspace/y_clot/env.sh'),R/'analysis/softening_feasibility_v1/model.py',PREV/'PAYLOAD_MANIFEST.json']}
dump(O/'PROTECTED_BASELINES.json',protected)
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',dict(library=str(LIB),library_sha256=sha(LIB),same_as_previous=True,recompiled=False,upstream=identity['upstream'],native_source_changes=[],native_source_sha256={p.name:sha(p) for p in (O/'source_identity').iterdir()},frozen_model_sha256=sha(R/'analysis/softening_feasibility_v1/model.py'),previous_manifest_sha256=sha(PREV/'PAYLOAD_MANIFEST.json'),previous_286_files_verified=True,production_library_sha256=protected[str(R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so')],production_modified=False))
for name in ['common.py','native.py','event_integrator.py','run_protocol.py']:
 old=R/'analysis/event_resolved_v1'/name;new=A/name
 if old.exists() and new.exists():
  with (O/'INPUT_AND_DRIVER_PATCH.diff').open('a') as f:f.write(''.join(difflib.unified_diff(old.read_text().splitlines(True),new.read_text().splitlines(True),fromfile='previous/'+name,tofile='domain_retention/'+name)))
print('PREPARED',O, 'BOX',BOX)
