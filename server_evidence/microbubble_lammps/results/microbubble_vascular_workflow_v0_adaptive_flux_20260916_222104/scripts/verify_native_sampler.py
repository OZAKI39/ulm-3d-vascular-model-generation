"""Compile and execute the unchanged native frozen sampler on Vast, for audit only."""
from pathlib import Path
import subprocess,shlex,json,hashlib,io
import numpy as np
S=Path(__file__).resolve().parents[1];p=json.loads((S/'provenance/STAGE_PATHS.json').read_text())
files=[S/'src'/n for n in ['frozen_flow.cpp','frozen_flow.hpp','flow_sampler_probe.cpp']]
r=subprocess.run(['rsync','-a',*[str(x) for x in files],p['ssh_alias']+':'+p['remote_stage']+'/src/'],capture_output=True,text=True,check=True)
log=r.stdout+r.stderr
payload={'source':p['remote_stage']+'/src','work':p['remote_work'],'field':p['donor_remote']+'/fields/FROZEN_FLOW_FIELD_V0.h5'}
code='''from pathlib import Path
import subprocess,shlex,json,hashlib
p=PAYLOAD
flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','hdf5'],text=True))
exe=Path(p['work'])/'flow_sampler_probe'
cmd=['c++','-std=c++17','-O2',p['source']+'/frozen_flow.cpp',p['source']+'/flow_sampler_probe.cpp','-o',str(exe)]+flags
r=subprocess.run(cmd,text=True,capture_output=True,check=True)
print(json.dumps({'command':cmd,'stdout':r.stdout,'stderr':r.stderr,'binary':str(exe),'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'binary_bytes':exe.stat().st_size,'field':p['field'],'runtime_scope':'Read-only velocity queries; no LAMMPS dynamics'},indent=2))
'''.replace('PAYLOAD',repr(payload))
(S/'provenance/NATIVE_PROBE_BUILD_SCRIPT.py').write_text(code)
r=subprocess.run(['ssh',p['ssh_alias'],'python3','-'],input=code,capture_output=True,text=True,check=True);log+=r.stderr
build=json.loads(r.stdout);(S/'provenance/NATIVE_PROBE_BUILD.json').write_text(json.dumps(build,indent=2)+'\n')
cmd=shlex.join([build['binary'],build['field']])
r=subprocess.run(['ssh',p['ssh_alias'],cmd],input=(S/'inputs/SAMPLER_COMPARISON_POINTS.txt').read_text(),capture_output=True,text=True,check=True);log+=r.stderr
(S/'raw/NATIVE_SAMPLER_OUTPUT.txt').write_text(r.stdout);(S/'provenance/NATIVE_PROBE_SSH.log').write_text(log)
native=np.loadtxt(io.StringIO(r.stdout));ref=np.load(S/'raw/SAMPLER_PYTHON_REFERENCE.npz')
assert native.shape==(len(ref['status']),4)
status_equal=bool(np.array_equal(native[:,0].astype(int),ref['status']))
error=float(np.abs(native[:,1:]-ref['velocity']).max());scale=float(np.abs(ref['velocity']).max())
assert status_equal and error<1e-12*max(scale,1e-10),(status_equal,error,scale)
res={'status':'PASS','scope':'native/Python sampler status and velocity equality on actual quadrature positions','points':len(native),'status_equal':status_equal,'max_velocity_absolute_difference_m_s':error,'velocity_scale_m_s':scale,'status_counts':{str(i):int((ref['status']==i).sum()) for i in range(5)},'full_cross_section_flux_is_not_validated_by_this_test':True}
(S/'validation/NATIVE_SAMPLER_EQUIVALENCE.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
