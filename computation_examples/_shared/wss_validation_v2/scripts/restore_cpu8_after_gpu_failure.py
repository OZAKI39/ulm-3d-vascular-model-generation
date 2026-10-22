"""Preserve failed GPU1 medium run and restore the previously registered CPU8 inputs."""
from pathlib import Path
import json,shutil,time,hashlib,argparse
root=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
case=root/'stage3/vessel_medium';oldlog=case/'run/solver.log';assert 'FRESH_RETRY_FAILED' in oldlog.read_text()
failed=root/'stage3/failed_attempts/vessel_medium_gpu1';assert not failed.exists();failed.parent.mkdir(exist_ok=True);case.rename(failed)
dump(failed/'reports/execution.json',dict(status='FAIL',reason='First fresh GPU1 solve DIVERGED_BREAKDOWN; another fresh retry failed. Process stopped; no qualified CFD state.',exit_code=None,solver_log_sha256=sha(failed/'run/solver.log'),relocated_from=str(case),relocated_to=str(failed)))
case.mkdir();shutil.copytree(failed/'SV_MESH',case/'SV_MESH');shutil.copytree(failed/'reports',case/'reports')
# Execution records belong solely to preserved failed case, never to fresh case.
for name in ['launch.json','progress.json','execution.json','linear_attempt_acceptance.json']:
 p=case/'reports'/name
 if p.exists():p.unlink()
(case/'run').mkdir();shutil.copy2(failed/'run/solver.xml',case/'run/solver.xml')
for name in ['mesh_request.json','mesh_launch.json']:
 if (failed/name).exists():shutil.copy2(failed/name,case/name)
for name in ['vessel_medium','vessel_fine']:
 c=root/'stage3'/name;archive=c/'reports/unexecuted_cpu8_inputs';assert archive.exists() and not (c/'run/solver.log').exists()
 preserved=c/'reports/gpu1_input_variant';preserved.mkdir(exist_ok=True)
 source=failed if name=='vessel_medium' else c
 for rel in ['policy.json','input_hashes.json','run/PETSC_OPTIONS.txt']:
  target=preserved/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source/rel,target)
  target=c/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(archive/rel,target)
 hashes=json.loads((c/'input_hashes.json').read_text());assert all(sha(c/k)==v for k,v in hashes.items())
 assert json.loads((c/'policy.json').read_text())['MPI_ranks']==8
 dump(c/'reports/CPU8_after_GPU1_failure.json',dict(unix=time.time(),reason='Actual medium GPU1 first fresh solve and fresh retry failed; previous CPU8 option restored with byte-identical registered input hashes. No physical, mesh, timestep or tolerance change.',failed_case=str(failed),original_CPU8_manifest_restored=True,initial_state='zero'))
 print(name,'original CPU8 inputs restored; mesh/XML/criteria unchanged')
