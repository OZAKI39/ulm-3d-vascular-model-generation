"""Read-only identity verification and final receipt; never launches computation."""
from common import *
import datetime,subprocess

def main():
 assert not (HERE/'TASK.json').exists(),'supervisor still active'
 receipts=[]
 for p in sorted((HERE/'monitoring').glob('*/EXECUTION.json')):
  r=read(p);assert r['status']=='EXITED';receipts.append(dict(job=p.parent.name,exit_code=r['exit_code'],stop_reason=r['stop_reason'],elapsed_s=r['elapsed_s'],peak_tree_RSS_GiB=r['peak_tree_RSS_GiB'],initial_memory=r['initial_memory'],sample_count=len(r['samples']),OOM_events_delta=(r['samples'][-1]['oom_events']-r['initial_memory']['oom_events']) if r['samples'] else 0,limits={k:v for k,v in r['task'].items() if k.startswith('max_')},execution_sha256=sha(p),log_sha256=sha(p.parent/'job.log')))
 total=sum(r['elapsed_s'] for r in receipts);assert total<5400
 previous=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json');source=read(IDENTITY/'EXACT_DEPLOYED_SOURCE_MANIFEST.json')
 assert all(sha(RETAINED/p)==h for p,h in source.items());assert sha(previous['actual_binary'])==previous['actual_binary_sha256'];assert sha(previous['PETSc']['path'])==previous['PETSc']['sha256']
 cuda=read(HERE/'CUDA_LIBRARY_IDENTITY.json');libs={}
 def walk(v):
  if isinstance(v,dict):
   if 'path' in v and 'sha256' in v and Path(v['path']).is_file():libs[v['path']]=v['sha256']
   libs.update(v.get('current_libraries_sha256',{}))
   for x in v.values():walk(x)
  elif isinstance(v,list):
   for x in v:walk(x)
 walk(cuda)
 assert all(sha(p)==h for p,h in libs.items())
 inputs=read(HERE/'INPUT_SHA256.json');assert all(sha(BASE/p)==h for p,h in inputs.items())
 runtime=[]
 for name in ['PIPE_CURVED_ELEVATION.npz','dispatch_original','dispatch_trace','native_dispatch.o','set_bc.cpp.o','nn.cpp.o','eq_assem.cpp.o']:
  path=HERE/'runtime'/name;runtime.append(dict(relative_path='runtime/'+name,server_path=str(path),size_bytes=path.stat().st_size,sha256=sha(path),tracked=False))
 put('RUNTIME_ARTIFACTS.json',dict(artifacts=runtime,original_case_meshes='Read-only symlinks to authoritative historical meshes; not newly generated/modified',full_CFD_binary='NOT_BUILT',high_order_CFD_checkpoint='NOT_RUN'))
 put('EXECUTION_RECEIPT.json',dict(status='ALL_SUPERVISED_TASKS_EXITED',jobs=receipts,non_CFD_budget_s=5400,supervised_compute_elapsed_s=total,peak_tree_RSS_GiB=max(r['peak_tree_RSS_GiB'] for r in receipts),CFD_runs=0,CFD_wall_time_s=0,linear_solves=0,nonlinear_steps=0,TASK_removed=True,no_automatic_followup=True,read_report_transfer_time_excluded=True,final_identity_checks=dict(source_files=len(source),input_files=len(inputs),production_binary_unchanged=True,PETSc_unchanged=True,CUDA_libraries_unchanged=len(libs)),initial_build_failure='Missing all_fun.h declaration in new dispatch harness; header added, same native libraries/flags. No CFD failure or restart.',resource_assessment_revision='Initial no-overlap sensitivity marked NEEDS_REVIEW and retained. Final includes measured original ASM duplication; no scientific gates changed.',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 # All delivered source/metadata/figures/data, excluding immutable runtime and logs.
 outputs={str(p.relative_to(HERE)):sha(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name not in ['OUTPUT_SHA256.json','TASK.json']}
 put('OUTPUT_SHA256.json',outputs)
 print(json.dumps(dict(delivered_artifacts=len(outputs),jobs=len(receipts),elapsed_s=total,peak_tree_RSS_GiB=max(r['peak_tree_RSS_GiB'] for r in receipts),all_tasks_exited=True)))
if __name__=='__main__':main()
