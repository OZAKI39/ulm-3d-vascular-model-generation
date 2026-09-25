from pathlib import Path
import subprocess,shlex,json,hashlib
p={'source': '/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_adaptive_flux_20260916_222104/src', 'work': '/workspace/microbubble_lammps/work/microbubble_vascular_workflow_v0_adaptive_flux_20260916_222104', 'field': '/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_20260916_214122/fields/FROZEN_FLOW_FIELD_V0.h5'}
flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','hdf5'],text=True))
exe=Path(p['work'])/'flow_sampler_probe'
cmd=['c++','-std=c++17','-O2',p['source']+'/frozen_flow.cpp',p['source']+'/flow_sampler_probe.cpp','-o',str(exe)]+flags
r=subprocess.run(cmd,text=True,capture_output=True,check=True)
print(json.dumps({'command':cmd,'stdout':r.stdout,'stderr':r.stderr,'binary':str(exe),'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'binary_bytes':exe.stat().st_size,'field':p['field'],'runtime_scope':'Read-only velocity queries; no LAMMPS dynamics'},indent=2))
