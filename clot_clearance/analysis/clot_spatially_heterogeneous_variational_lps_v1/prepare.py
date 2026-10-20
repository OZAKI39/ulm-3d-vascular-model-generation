from common import *
import shutil,subprocess,datetime
assert not (O/'P0_SCOPE.json').exists()
for p in [O/'provenance',O/'build',O/'checks',F/'figures',F/'videos']:p.mkdir(parents=True,exist_ok=True)
for name in ['reference.npz','reference.native','body.data']:shutil.copy2(MV/'provenance'/name,O/'provenance'/name)
for name in ['MATERIAL_AND_SPECIMEN_FREEZE.json','CURRENT_EXPOSED_SURFACE_SPEC.json','FRAGMENT_FLUID_TRACTION.json','MULTIFRAGMENT_OBJECTIVE.json']:
 shutil.copy2(MV/name,O/name)
shutil.copy2(MV/'provenance/auxiliary_surface_geometry.npz',O/'provenance/auxiliary_surface_geometry.npz')
import pandas as pd
h=pd.read_csv(MV/'cases/checkerboard/history.csv');t=np.r_[0.,h.time_s];amps=np.r_[3200.,h.nominal_amplitude_Pa]
assert np.all(np.diff(t)>0) and np.isfinite(amps).all() and np.min(amps)>=0 and np.max(amps)<=6000
np.savez_compressed(O/'SHEAR_DRIVE.npz',time_s=t,amplitude_Pa=amps)
# Verify the extraction against independently saved native callback amplitudes.
errors=[]
for p in (MV/'cases/checkerboard/keys').glob('*.npz'):
 s=load_npz(p);errors.append(abs(float(s['field_amplitude_Pa'])-np.interp(float(s['time']),t,amps)))
assert max(errors)<1e-8
shutil.copy2(MV/'cases/checkerboard/history.csv',O/'provenance/prior_checkerboard_history.csv')
dump(O/'SHEAR_DRIVE_FREEZE.json',dict(status='FROZEN_OPEN_LOOP_BEFORE_NEW_DYNAMICS',source=str(MV/'cases/checkerboard/history.csv'),source_sha256=sha(MV/'cases/checkerboard/history.csv'),drive_sha256=sha(O/'SHEAR_DRIVE.npz'),knots=len(t),interpolation='Piecewise linear interpolation of prior accepted amplitude endpoints including initial 3200 Pa; prior half-cosine envelope applied separately unchanged',max_native_key_amplitude_error_Pa=max(errors),both_cases_identical=True,feedback=False,time_s=dict(rise=RISE,hold=HOLD,off=OFF,end=END),spatial='Unchanged 250 um checkerboard surface tangent streamfunction and current-space C2 window',surface='Previous current supported surface and exterior-connected crack proxy, unchanged',label=LABEL,acoustic_realizability='UNESTABLISHED'))
ref=load_npz(O/'provenance/reference.npz')
dump(O/'P0_SCOPE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),git_sha=subprocess.check_output(['git','-c',f'safe.directory={R}','rev-parse','HEAD'],cwd=R,text=True).strip(),deleted_files=[],cleanup_reason='No verified disposable unreferenced code; isolate new work, preserve historical scientific evidence',baseline_library_sha256=sha(OLDLIB),original_points=len(ref['X']),material_nodes=int(ref['clot_count']),original_bonds=len(ref['bonds']),precut=int(ref['initial_cut'].sum()),baseline_files={str(p):sha(p) for p in [OLDLIB,MV/'FINAL_STATUS.json',MV/'DELIVERY_MANIFEST.json',HA/'bond_bpm_peri_variational.cpp',HA/'bond_bpm_peri_variational.h']}))
(O/'provenance/git_existing.diff').write_bytes(subprocess.check_output(['git','-c',f'safe.directory={R}','diff'],cwd=R))
print('PREPARED',O,flush=True)
