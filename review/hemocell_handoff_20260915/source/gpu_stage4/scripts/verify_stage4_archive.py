#!/usr/bin/env python3
"""Read-only audit (unless --seal) of the completed Stage4 archive."""
from pathlib import Path
import hashlib,json,gzip,sys,csv,collections,subprocess,os
R=Path(__file__).resolve().parents[1];B=R.parent/'20260914_stage3_sustained'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def j(p):return json.loads(p.read_text())
def verify_manifest(root,m):
 count=0
 for line in m.read_text().splitlines():
  digest,rel=line.split('  ',1);p=root/rel
  assert p.is_file() and sha(p)==digest,rel;count+=1
 return count
def check():
 terminal=j(R/'provenance/SCHEDULE_TERMINAL.json');selection=j(R/'FINAL_GPU_SELECTION.json')
 for rel,d in j(R/'provenance/execution_hashes.json').items():assert sha(R/rel)==d,rel
 for rel,d in j(R/'provenance/BASELINE_SOURCE_SHA256.json').items():assert sha(R/'baseline'/rel)==d,rel
 assert sha(B/'SHA256SUMS')==j(R/'provenance/STAGE3_PROTECTED_SEAL.json')['sha256sums']
 old_count=verify_manifest(B,B/'SHA256SUMS')
 assert sha(R/'verification/GPU_NUMERICAL_COMPARISON_CONTRACT.json')=='5e78c61974be7ca2a52b3add288858330804c8d74d47ae909fced717d1f02456'
 code=(R/'source/gpu_adapter.hpp').read_text();orig=(R/'baseline/source/gpu_adapter.hpp').read_text()
 assert code[code.index('    Stats inspect'):]==orig[orig.index('    Stats inspect'):]
 alo=orig.index('// Scratch MPI1');ahi=orig.index('struct GPUState {')
 blo=code.index('// Scratch MPI1');bhi=code.index('#include "stage4_lbm_batch.hpp"')
 assert orig[alo:ahi]==code[blo:bhi]
 guo='        {Stage2Scope p("GUO_FUSED_WALL_INLET_OUTLETS");'
 assert orig[orig.index(guo):orig.index('    Stats inspect')]==code[code.index(guo):code.index('    Stats inspect')]
 native=(R/'baseline/native/palabos/src/atomicBlock/atomicAcceleratedLattice3D.hh').read_text()
 start=native.index('void AtomicAcceleratedLattice3D<T, Descriptor>::collideAndStream(CollFun const &collFun)')
 tail=native[start:native.index('template <typename T, template <typename U> class Descriptor>',start+5)]
 body=tail.split('            size_t i = &f0 - populations;\n',1)[1].split('\n        });',1)[0]
 adapted=(R/'upstream/palabos/src/atomicBlock/atomicAcceleratedLattice3D.hh').read_text()
 assert body in adapted[adapted.index('::Stage4BatchView::collideCell('):]
 assert hashlib.sha256(body.encode()).hexdigest()==j(R/'provenance/STATIC_MATH_REUSE_PROOF.json')['body_extraction_sha256']
 assert j(R/'provenance/binary_hashes.json')['gpu']==sha(R/'build_gpu/vascular/vascularPoC')
 counts=[];decoded=0;maximum=0
 for d in sorted((R/'runs').iterdir()):
  start=j(d/'RUN_STARTED.json');assert 0<start['total_steps']<=5000
  rec=j(d/'RUN_TERMINAL.json') if (d/'RUN_TERMINAL.json').exists() else None
  assert rec is not None,'Incomplete process receipt '+d.name
  maximum=max(maximum,rec['completed_steps'])
  if rec['status']=='PASS':
   assert rec['completed_steps']==start['total_steps'] and rec['returncode']==0 and rec['binding_status']=='PASS'
   stored=j(d/'provenance/LOSSLESS_STORAGE.json')
   for entry in stored['entries']:
    p=d/(entry['path']+'.gz');assert sha(p)==entry['gzip_sha256']
    h=hashlib.sha256();n=0
    with gzip.open(p,'rb') as f:
     for b in iter(lambda:f.read(2**20),b''):h.update(b);n+=len(b)
    assert h.hexdigest()==entry['raw_sha256'] and n==entry['raw_bytes'];decoded+=1
  count=verify_manifest(d,d/'SHA256SUMS');counts.append(dict(run=d.name,files=count,status=rec['status']))
 assert maximum<=5000
 # Validate realized schedule limits; no new runs were silently substituted.
 names={x['run'] for x in counts}
 assert all(n in {'base_correctness','batch_correctness','cpu_sanity','batch_profile'} or
            any(n==f'{v}_{h}_r{i}' for v in ['baseline','batch'] for h,rs in [(200,3),(1000,2),(5000,2)] for i in range(1,rs+1)) for n in names)
 if any('_1000_' in n for n in names):assert j(R/'verification/SHORT_GATE.json')['status']=='PASS'
 if any('_5000_' in n for n in names):assert j(R/'verification/HORIZON_1000_GATE.json')['status']=='PASS'
 if selection['STAGE4_OPTIMIZATION_ACCEPTED']=='YES':
  assert terminal['status']=='ACCEPTED' and selection['PRODUCTION_GPU_CANDIDATE']=='STAGE4'
  for key in ['INITIAL_CORRECTNESS_GATE','SHORT_GATE','HORIZON_1000_GATE','HORIZON_5000_GATE','MIGRATION_GATE','STABILITY_GATE']:assert j(R/'verification'/f'{key}.json')['status']=='PASS'
 else:assert selection['PRODUCTION_GPU_CANDIDATE']=='STAGE3'
 required=['GPU_STAGE4_REPORT.md','LBM_DISPATCH_INVENTORY.csv','STAGE4_CORRECTNESS.csv','STAGE4_SHORT_BENCHMARK.csv','STAGE4_1000_BENCHMARK.csv','STAGE4_5000_BENCHMARK.csv','PRE_POST_KERNEL_PROFILE.csv','PRE_POST_MEMORY_MIGRATION.json','PRE_POST_TIMESTEP_BREAKDOWN.csv','GPU_RESOURCE_TRACE.csv','STAGE4_PATCH.diff','FINAL_GPU_SELECTION.json']
 for n in required:assert (R/n).is_file(),n
 return dict(archive_audit='PASS',stage4_result=terminal['status'],protected_stage3_files=old_count,source_frozen_files=len(j(R/'provenance/execution_hashes.json')),independent_baseline_files=len(j(R/'provenance/BASELINE_SOURCE_SHA256.json')),runs=counts,lossless_snapshots_verified=decoded,max_actual_steps=maximum,physics_changed=False,numerical_contract_changed=False,formal_source_modified=False)
if __name__=='__main__':
 result=check()
 if '--seal' in sys.argv:
  (R/'verification/FINAL_INTEGRITY_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
  # Supervisor must already be stopped, so its logs cannot change after sealing.
  assert not (R/'stage4_supervisor.pid').exists(),'Stop task-local supervisor before sealing'
  paths=[p for p in sorted(R.rglob('*')) if p.is_file() and p.relative_to(R).parts[0]!='tmp' and p!=R/'SHA256SUMS']
  (R/'SHA256SUMS').write_text(''.join(sha(p)+'  '+str(p.relative_to(R))+'\n' for p in paths))
  result['sealed_files']=verify_manifest(R,R/'SHA256SUMS')
 elif (R/'SHA256SUMS').exists():result['sealed_files']=verify_manifest(R,R/'SHA256SUMS')
 print(json.dumps(result,indent=2))
