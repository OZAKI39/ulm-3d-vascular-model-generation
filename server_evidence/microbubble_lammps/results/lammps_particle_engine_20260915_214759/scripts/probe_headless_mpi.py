from pathlib import Path
import subprocess,json,os,time,datetime
w=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759');p=w/'build_provenance'
env=os.environ.copy();display=env.pop('DISPLAY',None)
cmd=['timeout','-k','2','20','mpirun','--allow-run-as-root','--bind-to','none','-np','1',str(w/'build_gpu/lmp'),'-h'];t=time.time()
with (p/'gpu_headless_mpi1_help.txt').open('w') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
text=(p/'gpu_headless_mpi1_help.txt').read_text()
info={'status':'PASS' if r.returncode==0 and 'nve/sphere/kk' in text and 'gran/hooke/history/kk' in text else 'FAIL','returncode':r.returncode,'elapsed_seconds':time.time()-t,'inherited_DISPLAY':display,'change':'Remove DISPLAY only from child MPI execution environment','evidence':'Strace of mpirun hostname blocked in X11 handshake poll on localhost TCP 6006, before launching an application. No numerical solver case had started.','command':cmd,'solver_timesteps':0,'physical_contract_changed':False,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(p/'HEADLESS_MPI_STARTUP_VERIFICATION.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info));assert info['status']=='PASS'
