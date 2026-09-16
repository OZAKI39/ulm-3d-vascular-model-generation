from pathlib import Path
import subprocess,os,signal,json,time
w=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759');p=w/'build_provenance'
for pid in [52014]:
 path=Path(f'/proc/{pid}/cmdline')
 if path.exists():
  cmd=path.read_bytes().replace(b'\x00',b' ').decode();assert 'mpirun' in cmd and 'build_gpu/lmp -h' in cmd
  os.kill(pid,signal.SIGKILL);print('STOPPED_TIMED_OUT_HELP_PROBE',pid)
cmd=['timeout','-k','2','15','strace','-f','-tt','-s','128','-o',str(p/'mpi_startup_strace.txt'),'mpirun','--allow-run-as-root','--bind-to','none','-np','1','/bin/hostname']
t=time.time()
with (p/'mpi_hostname_probe.txt').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
info={'command':cmd,'returncode':r.returncode,'elapsed':time.time()-t,'solver_timesteps':0}
(p/'MPI_STARTUP_DIAGNOSTIC.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info));print((p/'mpi_hostname_probe.txt').read_text())
print('\n'.join((p/'mpi_startup_strace.txt').read_text().splitlines()[-45:]))
