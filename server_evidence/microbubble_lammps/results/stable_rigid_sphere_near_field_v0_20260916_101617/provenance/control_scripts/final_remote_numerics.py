from pathlib import Path
import subprocess,os,json,hashlib,shutil
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617');W=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617');PY='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
with open(R/'logs/REMOTE_NUMERICAL_FINALIZER.log','w') as f:p=subprocess.run([PY,'-B',str(R/'scripts/finalize_rigid_v0.py')],cwd=R,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',HWLOC_COMPONENTS='-gl'),stdout=f,stderr=subprocess.STDOUT)
print('FINALIZER_RC',p.returncode);print((R/'logs/REMOTE_NUMERICAL_FINALIZER.log').read_text()[-4800:]);p.check_returncode()
# Stop only this task's completed supervisor, preserving every previous failure/run record.
q=subprocess.run(['supervisorctl','-c',str(W/'supervisord.conf'),'shutdown'],capture_output=True,text=True);print(q.stdout,q.stderr)
for n in ['supervisord.conf','supervisor.log']:
 shutil.copy2(W/n,R/'provenance'/n)
build={'active_release':'22Jul2025 Update 6','active_commit':'9c5ab448c78a14fd534619622162ba418d6a1fb1','upstream_modified':False,'source_root':str(R/'src'),'build_root':str(W/'build'),'binaries':{}}
for n in ['rigid_lmp_cpu','rigid_lmp_gpu','librigid_math.so']:
 p=W/'build'/n;ld=subprocess.run(['ldd',str(p)],capture_output=True,text=True);assert 'lammps_migration' not in ld.stdout
 build['binaries'][n]={'path':str(p),'size_bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'ldd':ld.stdout}
for n in ['CMakeCache.txt','build.ninja']:shutil.copy2(W/'build'/n,R/'provenance'/n)
(R/'provenance/PROJECT_BUILD_PROVENANCE.json').write_text(json.dumps(build,indent=2)+'\n')
