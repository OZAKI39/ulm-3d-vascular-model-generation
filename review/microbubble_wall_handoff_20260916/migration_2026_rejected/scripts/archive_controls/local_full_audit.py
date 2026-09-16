from pathlib import Path
import json,hashlib,subprocess,sys
C=Path(__file__).parent;S=json.loads((C/'TASK_PATHS.json').read_text());R=Path(S['local_result'])
v=subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=R,capture_output=True,text=True);assert v.returncode==0,v.stderr
(C/'WSL_SHA256_INITIAL.log').write_text(v.stdout+v.stderr)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
oldbase=Path('/home/lzy/projects/compre_output');old={'engine':oldbase/'lammps_particle_engine/20260915_214759','coupling':oldbase/'palabos_lammps_coupling_v0/20260915_224019','passive':oldbase/'passive_transport_v0/20260915_233516'}
inputs={}
for name,root in old.items():
 mismatch=[];allowed=[];count=0
 for line in (root/'FROZEN_INPUT_SHA256SUMS').read_text().splitlines():
  h,n=line.split('  ',1);count+=1
  if sha(R/'migration'/name/n)!=h:
   (allowed if n in ['src/coupled_lammps_main.cpp','src/passive_main.cpp'] else mismatch).append(n)
 assert not mismatch,(name,mismatch)
 inputs[name]={'files_compared':count,'numerical_input_and_math_mismatches':mismatch,'explicit_API_compatibility_changes':allowed}
(C/'LOCAL_BASELINE_INPUT_COMPARISON.json').write_text(json.dumps({'status':'PASS','stages':inputs},indent=2)+'\n')
out=C/'independent_finalizer';out.mkdir(exist_ok=True)
p=subprocess.run(['/usr/bin/python3','-B',str(R/'scripts/finalize_migration_candidate.py'),'--root',str(R),'--output-dir',str(out)],capture_output=True,text=True);(C/'LOCAL_INDEPENDENT_FINALIZER.log').write_text(p.stdout+p.stderr);print(p.stdout);assert p.returncode==0,p.stderr
result=json.loads((out/'INDEPENDENT_MIGRATION_EVALUATION.json').read_text());assert result['status']=='FAIL' and result['source_correction_gate']=='FAIL' and all(x['status']=='PASS'for x in result['regressions'].values()),result
remote=json.loads((R/'validation/remote_independent/SOURCE_CORRECTION_RECOMPUTED.json').read_text());local=json.loads((out/'SOURCE_CORRECTION_RECOMPUTED.json').read_text());assert remote==local
# Preserve output outside the downloaded immutable snapshot until the final metadata seal.
receipt={'status':'PASS','scientific_migration_status':'FAIL','all_three_numerical_regressions':'PASS','source_correction_gate':'FAIL_INDEPENDENTLY_REPRODUCED','all_frozen_math_and_inputs_unchanged':'YES','source_correction_recomputed_identical_to_remote':True,'manifest_sha256':sha(R/'SHA256SUMS'),'manifest_entries':len((R/'SHA256SUMS').read_text().splitlines())}
(C/'LOCAL_FINAL_AUDIT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
