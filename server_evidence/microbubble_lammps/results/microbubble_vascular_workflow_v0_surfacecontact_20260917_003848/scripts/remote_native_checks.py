from pathlib import Path
import subprocess,json,os,hashlib,datetime,platform
S=Path(__file__).resolve().parents[1];W=Path(str(S).replace('/results/','/work/'));B=W/'build';assert str(S).startswith('/workspace/') and 'surfacecontact_' in str(S)
env={**os.environ,'HWLOC_COMPONENTS':'-gl','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
checks=[]
for name,cmd in [('injection',[str(B/'test_injection')])]+[(f'ports_lifecycle_mpi{n}',['mpirun','--allow-run-as-root','--bind-to','none','--oversubscribe','-n',str(n),str(B/'test_ports_lifecycle'),str(S)]) for n in [1,4]]:
 p=subprocess.run(cmd,cwd=W,env=env,text=True,capture_output=True,timeout=120);(S/'validation'/f'NATIVE_{name}.log').write_text(p.stdout+p.stderr);assert p.returncode==0,(name,p.stderr);checks.append(dict(name=name,command=cmd,status='PASS',stdout=p.stdout))
(S/'validation/NATIVE_INHERITED_CHECKS.json').write_text(json.dumps(dict(status='PASS',checks=checks),indent=2)+'\n')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
j=dict(status='PASS',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),hostname=platform.node(),source_hashes={str(p.relative_to(S)):sha(p) for p in sorted((S/'src').glob('*')) if p.is_file()},binaries={n:dict(path=str(B/n),sha256=sha(B/n)) for n in ['workflow_lmp','librigid_math.so','surface_contact_probe','test_surface_trajectories','test_injection','test_ports_lifecycle']},CMake_sha256=sha(S/'CMakeLists.txt'),compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0],MPI=subprocess.check_output(['mpirun','--version'],text=True).splitlines()[0],LAMMPS_release='22Jul2025 Update6',LAMMPS_commit='9c5ab448c78a14fd534619622162ba418d6a1fb1',LAMMPS_static_library=dict(path='/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759/build_cpu/liblammps.a',sha256=sha(Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759/build_cpu/liblammps.a'))),CPU_only=True,upstream_modified=False,raw_resistance_model_changed=False,wall_hydrodynamics_main_used=False,geometry_modified=False)
(S/'provenance/PROJECT_BUILD_PROVENANCE.json').write_text(json.dumps(j,indent=2)+'\n');print('NATIVE_INHERITED_AND_BUILD_PROVENANCE_PASS')
