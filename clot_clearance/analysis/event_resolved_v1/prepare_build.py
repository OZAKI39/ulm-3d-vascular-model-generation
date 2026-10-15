from common import *
import shutil,subprocess,time,difflib
assert not LIB.exists() and not (O/'replay/STARTED.json').exists(), 'Never overwrite a computed research run'
O.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True)
for d in ['provenance','source_identity','replay','local','full']:(O/d).mkdir(exist_ok=True)
shutil.copytree(PREV/'provenance',O/'provenance',dirs_exist_ok=True)
shutil.copy2(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
paths=[p for p in PREV.rglob('*') if p.is_file()]+[OLD/'build/liblammps.so',R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so',Path('/workspace/y_clot/env.sh'),R/'analysis/softening_feasibility_v1/model.py']
dump(O/'PROTECTED_BASELINES.json',{str(p):sha(p) for p in paths})
assert sha(OLD/'build/liblammps.so')=='866065394c48ebc8a6d32f147bde992a325c2a503a49d50392243eee10fea9bd'
if not (O/'source').exists():shutil.copytree(OLD/'source',O/'source',symlinks=True)
p=O/'source/src/fix_move.cpp';before=p.read_text()
old='''  error->all(FLERR, Error::NOLASTLINE, "Resetting timestep size is not allowed with fix move");'''
new='''  // Isolated event integrator: only fully prescribed equal-variable motion.
  // These six expressions use native cumulative `time`; no elapsed-step*dt
  // displacement or partially integrated component is admitted here.
  if (mstyle != VARIABLE || !xvarstr || !yvarstr || !zvarstr ||
      !vxvarstr || !vyvarstr || !vzvarstr ||
      xvarstyle != EQUAL || yvarstyle != EQUAL || zvarstyle != EQUAL ||
      vxvarstyle != EQUAL || vyvarstyle != EQUAL || vzvarstyle != EQUAL)
    error->all(FLERR, Error::NOLASTLINE,
               "Variable dt requires six prescribed equal variables in isolated fix move");
  dt = dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;'''
assert before.count(old)==1
after=before.replace(old,new);p.write_text(after)
(O/'NATIVE_SOURCE_PATCH.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='certified/src/fix_move.cpp',tofile='isolated/src/fix_move.cpp')))
(O/'source_identity/fix_move.before.cpp').write_text(before)
shutil.copy2(p,O/'source_identity/fix_move.cpp')
shutil.copy2(A/'event_state_audit.cpp',O/'source/src/event_state_audit.cpp')
shutil.copy2(A/'event_state_audit.cpp',O/'source_identity/event_state_audit.cpp')
with (O/'NATIVE_SOURCE_PATCH.diff').open('a') as f:
 lines=(A/'event_state_audit.cpp').read_text().splitlines();f.write('\n--- /dev/null\n+++ isolated/src/event_state_audit.cpp\n@@ -0,0 +1,%d @@\n'%len(lines));f.write('\n'.join('+'+line for line in lines)+'\n')
core=['src/BPM/bond_bpm_peri_variational.cpp','src/BPM/bond_bpm_peri_variational.h','src/fix_lps_variational_dynamic.cpp','src/fix_lps_variational_dynamic.h','src/fix_nve.cpp','src/verlet.cpp','src/update.cpp','src/input.cpp','src/read_restart.cpp','src/BPM/pair_bpm_peri.cpp']
identity={}
for name in core:
 p=O/'source'/name
 if not p.exists():
  found=list((O/'source/src').rglob(Path(name).name));assert len(found)==1; p=found[0];name=str(p.relative_to(O/'source'))
 assert sha(p)==sha(OLD/'source'/name)
 shutil.copy2(p,O/'source_identity'/p.name);identity[name]=dict(old_sha256=sha(OLD/'source'/name),new_sha256=sha(p),unchanged=True)
start=time.perf_counter()
with (O/'build.log').open('w') as f:
 subprocess.run(['cmake','-S',str(O/'source/cmake'),'-B',str(O/'build'),'-DCMAKE_BUILD_TYPE=Release','-DBUILD_MPI=ON','-DBUILD_OMP=ON','-DBUILD_SHARED_LIBS=ON','-DPKG_BPM=ON','-DPKG_MOLECULE=ON','-DMPI_CXX_COMPILER=/workspace/y_clot/deps/openmpi/bin/mpicxx'],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=120)
 subprocess.run(['cmake','--build',str(O/'build'),'-j12'],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=900)
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',dict(old_library=str(OLD/'build/liblammps.so'),old_sha256=sha(OLD/'build/liblammps.so'),new_library=str(LIB),new_sha256=sha(LIB),only_modified_cpp='src/fix_move.cpp',additional_readonly_cpp='event_state_audit.cpp',observer_sha256=sha(A/'event_state_audit.cpp'),fix_move_before_sha256=sha(O/'source_identity/fix_move.before.cpp'),fix_move_after_sha256=sha(p:=O/'source/src/fix_move.cpp'),unchanged_sources=identity,upstream='8de817dd79bfe4525d5d39246a212d833e6dee07',compile_wall_s=time.perf_counter()-start,compile_jobs=12,production_replaced=False))
print('BUILD_DONE',sha(LIB),flush=True)
