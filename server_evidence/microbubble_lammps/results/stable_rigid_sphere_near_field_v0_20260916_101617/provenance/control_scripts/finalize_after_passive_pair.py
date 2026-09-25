from pathlib import Path
import subprocess,os,json,shutil
R=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617');W=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617');PY='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
assert json.loads((R/'PHASE_PASSIVE_PAIR_STATE.json').read_text())['status']=='PASS'
q=subprocess.run(['supervisorctl','-c',str(W/'j2_supervisord.conf'),'shutdown'],capture_output=True,text=True);print(q.stdout,q.stderr,flush=True);q.check_returncode()
for n in ['j2_supervisord.conf','j2_supervisor.log']:shutil.copy2(W/n,R/'provenance'/n)
with open(R/'logs/REMOTE_NUMERICAL_FINALIZER.log','w') as f:p=subprocess.run([PY,'-B',str(R/'scripts/finalize_rigid_v0.py')],cwd=R,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',HWLOC_COMPONENTS='-gl'),stdout=f,stderr=subprocess.STDOUT)
print('FINALIZER_RC',p.returncode);print((R/'logs/REMOTE_NUMERICAL_FINALIZER.log').read_text()[-4800:]);p.check_returncode()
