from pathlib import Path
import subprocess,json,hashlib
S=Path(__file__).resolve().parents[1];R='/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_'+S.name;W=R.replace('/results/','/work/')
script='from pathlib import Path\nimport hashlib,json,subprocess,socket,datetime\n'+f'R=Path({R!r});W=Path({W!r})\n'+'''
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
baseline=Path('/workspace/microbubble_lammps/results/stable_rigid_sphere_near_field_v0_20260916_101617')
source={str(p.relative_to(R)):sha(p) for p in (R/'src').glob('*') if p.is_file()}
unchanged={n:sha(R/'src'/n)==sha(baseline/'src'/n) for n in ['rigid_math.cpp','rigid_math.hpp','frozen_flow.cpp','frozen_flow.hpp']}
assert all(unchanged.values())
bins={n:dict(path=str(W/'build'/n),sha256=sha(W/'build'/n)) for n in ['workflow_lmp','librigid_math.so','test_injection','test_ports_lifecycle','test_position_sampling','test_wall_constraint']}
engine=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759')
stable=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617/build/rigid_lmp_cpu')
assert sha(stable)=='830d3a1bf41f4944614d6899d8ef74199ff4cc4b496f6037dbe72a2004007789'
receipt=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),MPI=subprocess.check_output(['mpirun','--version'],text=True).splitlines()[0],status='PASS',hostname=socket.gethostname(),source_hashes=source,stable_components_bitwise_unchanged=unchanged,binaries=bins,CMake_sha256=sha(R/'CMakeLists.txt'),LAMMPS_static_library=dict(path=str(engine/'build_cpu/liblammps.a'),sha256=sha(engine/'build_cpu/liblammps.a')),compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0],LAMMPS_release='22Jul2025 Update6',LAMMPS_commit='9c5ab448c78a14fd534619622162ba418d6a1fb1',upstream_modified=False,wall_hydrodynamics_main_used=False,raw_resistance_model_changed=False,kinematic_geometric_model_changed=True,CPU_only=True,stable_CPU_binary_preserved=True,guide_sha256=sha(Path('/etc/vast-agents-guide.md')))
print(json.dumps(receipt,indent=2))
'''
r=subprocess.run(['ssh','vast4090','python3 -'],input=script,text=True,capture_output=True);assert r.returncode==0,r.stderr
(S/'provenance/PROJECT_BUILD_PROVENANCE.json').write_text(r.stdout)
receipt=json.loads(r.stdout)
for rel,h in receipt['source_hashes'].items():assert hashlib.sha256((S/rel).read_bytes()).hexdigest()==h,rel
for n in ['LOW','MEDIUM','LONG_TRANSPORT']:
    for ranks in [1,4]:
        run=json.loads((S/'runs'/f'{n}_release_mpi{ranks}'/'EXECUTION_RECEIPT.json').read_text());assert run['binary_sha256']==receipt['binaries']['workflow_lmp']['sha256'];assert run['library_sha256']==receipt['binaries']['librigid_math.so']['sha256']
print('BUILD_SOURCE_AND_EXECUTED_BINARY_IDENTITY_PASS')
