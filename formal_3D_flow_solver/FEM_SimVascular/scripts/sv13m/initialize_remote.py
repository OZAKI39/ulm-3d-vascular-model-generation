"""Freeze remote Stage L and prepare exactly one same-input baseline run."""
import json,shutil
from pathlib import Path
from runner_remote import *
ref=json.loads((BASE/'configs/reference_manifest.json').read_text())
stack=ref['stack'];old=Path(stack['svmp_gpu_build']['build']).parents[2]
assert old.name=='sv1_3l'
for n in ('svmp_gpu_build','compatibility_winner'):
 write('baseline_L_'+n,stack[n])
for n in ('flow_input_manifest.json','petsc_options.json','mpi_resolution.json'):
 shutil.copyfile(old/'configs'/n,BASE/'configs'/n)
external=BASE/'external';external.mkdir(exist_ok=True)
shutil.copyfile(old/'external/flow_inputs.tar.gz',external/'flow_inputs.tar.gz')
(BASE/'reports/remote_agent_guide.txt').write_text(Path('/etc/vast-agents-guide.md').read_text())
print('Stage L stack and input descriptors copied; original files remain read-only.',flush=True)
