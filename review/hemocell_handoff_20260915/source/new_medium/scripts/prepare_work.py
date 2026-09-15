from pathlib import Path
import os,json,hashlib,shutil,subprocess,time
N=Path(__file__).resolve().parents[1];B=Path('/workspace/hemocell_restore/work/stage4_restore_20260915_150938');A=Path('/workspace/hemocell_restore/archive/minimal_restore_20260915_143631');P=A/'payload/hemocell_rbc_stage1/20260915_010138'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
assert json.loads((B/'FINALIZATION_TERMINAL.json').read_text())['restore_status']=='PASS'
assert not (N/'provenance/PREPARATION_COMPLETE.json').exists()
freeze=json.loads((B/'provenance/WORKING_COPY_SHA256.json').read_text());bad=[p for p,h in freeze.items() if sha(B/p)!=h];assert not bad,bad
oldbin=json.loads((B/'provenance/STAGE4_BINARY.json').read_text());assert sha(Path(oldbin['path']))==oldbin['sha256']
protected={str(B/p):h for p,h in freeze.items()}
for p in [Path(oldbin['path']),B/'build/official/libpalabos.a',B/'BUILD_PROVENANCE.json',B/'OFFICIAL_GPU_SMOKE.json']:
 protected[str(p)]=sha(p)
for p in B.glob('stage4_*/RUN_TERMINAL.json'):protected[str(p)]=sha(p)
save(N/'provenance/PROTECTED_BASELINE_SHA256.json',protected)
for p in (B/'source').iterdir():
 if p.is_file():shutil.copy2(p,N/'source'/p.name)
shutil.copy2(P/'scripts/prepare_numerics.py',N/'provenance/ORIGINAL_PREPARE_NUMERICS.py')
shutil.copy2(P/'inputs/NUMERICS_INPUT.json',N/'provenance/ORIGINAL_NUMERICS_INPUT.json')
shutil.copy2(B/'frozen_bundle/frozen_inputs/step3_baseline_source/prepare_numerics.py',N/'provenance/STEP3_PREPARE_NUMERICS.py')
(N/'frozen_bundle').symlink_to(B/'frozen_bundle',target_is_directory=True)
for name in ['multiplane_quadrature.tsv','control_volume_indices.txt']:
 shutil.copy2(B/'stage4_200/contracts'/name,N/'contracts'/name)
shutil.copy2(B/'stage4_200/contracts/solver_parameters.txt',N/'inputs/OLD_SOLVER_PARAMETERS.txt')
x=json.loads((P/'inputs/NUMERICS_INPUT.json').read_text());x['scope']='PURE_FLUID_NEW_MEDIUM_SMOKE provisional development medium; 500 then5000, no RBC';x['baseline_work']=str(B);x['old_solver_parameters_sha256']=sha(N/'inputs/OLD_SOLVER_PARAMETERS.txt');x['new_medium_inlet_multiplier_validation']='NOT_PERFORMED'
x['geometry_contract_source']=str(B/'stage4_200/contracts/lattice_unit_contract.json')
if Path(x['geometry_contract_source']).is_file():x['geometry_contract_sha256']=sha(Path(x['geometry_contract_source']))
assert x['dx_m']==float((N/'inputs/OLD_SOLVER_PARAMETERS.txt').read_text().splitlines()[2].split()[6])
save(N/'inputs/NUMERICS_INPUT.json',x)
for n in (500,5000):
 d=N/f'run_{n}'
 for sub in ['diagnostics/field_samples','logs','provenance']:(d/sub).mkdir(parents=True,exist_ok=True)
 (d/'contracts').symlink_to(N/'contracts',target_is_directory=True)
save(N/'provenance/PREPARATION_COMPLETE.json',dict(status='PASS',protected_baseline_files=len(protected),old_binary_only_hashed_not_executed=True,geometry_shared_read_only=str(N/'frozen_bundle'),GPU_core_source_sha256={p.name:sha(p) for p in (N/'source').glob('*.hpp')},new_dt_not_yet_generated=True,unix=time.time()))
print('WORK_PREPARATION_PASS',flush=True)
