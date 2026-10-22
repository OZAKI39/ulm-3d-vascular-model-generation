"""Read-only historical protection and complete artifact hashes; no Git push."""
import time,subprocess,shutil,json,numpy as np
from pathlib import Path
from .common import *
def main():
    start=time.perf_counter();verify();expected=read(BASE/'HISTORICAL_RAW_SHA256.json');history=[dict(path=p,sha256=h,unchanged=sha(Path(p))==h) for p,h in expected.items()]
    # F10AB uses several campaign subdirectories, not only runs/.
    ab=read(AB/'RAW_DATA_MANIFEST.json')['records'];abchecks=[dict(path=r['raw_path'],sha256=r['raw_sha256'],unchanged=sha(AB/r['raw_path'])==r['raw_sha256']) for r in ab]
    assert all(r['unchanged'] for r in history+abchecks)
    write(BASE/'evidence/HISTORICAL_PRESERVATION.json',dict(AC_history=history,AB_history=abchecks,all_unchanged=True,protected_source_and_inputs=True,identity=identity()))
    server=Path('/workspace/mevo_acoustic_microbubbles/f10ae_target_20261019');dst=BASE/'evidence/server_logs';dst.mkdir(parents=True,exist_ok=True)
    for p in (server/'logs').iterdir():
        if p.is_file():shutil.copy2(p,dst/p.name)
    for name in ['setup.sh','install_optimizer.sh']:
        if (server/name).exists():shutil.copy2(server/name,dst/name)
    g=read(BASE/'FINAL_GATE.json');g['gates']['HISTORICAL_TRAJECTORIES_PRESERVED'].update(status='PASS',evidence=['evidence/HISTORICAL_PRESERVATION.json'],failure_reason=None)
    visual=read(BASE/'VISUAL_INSPECTION.json') if (BASE/'VISUAL_INSPECTION.json').exists() else None
    g['gates']['FIGURES_VALIDATED'].update(status='PASS' if visual and visual.get('figures_pass') else 'INCONCLUSIVE',failure_reason=None if visual and visual.get('figures_pass') else 'Awaiting visual inspection')
    vm=read(BASE/'videos/MANIFEST.json');vpass=all(v['fully_decoded'] and v['uniform_PTS'] for v in vm) and visual and visual.get('videos_pass')
    g['gates']['VIDEOS_VALIDATED'].update(status='PASS' if vpass else 'INCONCLUSIVE',failure_reason=None if vpass else 'Decode checks complete; visual inspection pending')
    for v in g['gates'].values():assert v['status'] in ['PASS','FAIL','NOT_RUN','INCONCLUSIVE']
    write(BASE/'FINAL_GATE.json',g)
    stages={}
    for file in ['LOW_PRESSURE_ANALYSIS.json','BASELINE_VALIDATION.json','OPTIMIZER_COMPARISON.json','SELECTED_FROZEN.json','VALIDATION_COMPLETE.json','WORKFLOW_COMPLETE.json','LOW_PRESSURE_ADMISSION.json','AUDIT.json','figures/MANIFEST.json']:
        r=read(BASE/file);stages[file]=r.get('wall_s')
    stages['videos']=sum(v['wall_s'] for v in vm)
    raw=read(BASE/'LOW_PRESSURE_RAW.json');stages['APECSS_sum_process_elapsed_s']=sum(r['wall_s'] for r in raw)
    rows=read(BASE/'ALL_RUNS.json');cost=dict(stages=stages,physical_attempts=len(rows),complete=sum(r['status']=='PASS' for r in rows),solver_sum_wall_s=sum(r.get('wall_s',0) for r in rows),orchestration_task_sum_wall_s=sum(r.get('task_wall_s',0) for r in rows),optimizer_s=read(BASE/'OPTIMIZER_COMPARISON.json')['methods'][1]['optimizer_s'],warning='Parallel summed run times are not campaign elapsed. Workflow includes child stages; never add child stage walls to it. APECSS sum includes the failed 120 s refinement.',identity=identity())
    # Supervisor timestamps are actual process-wall evidence; exclude job waits from summed compute.
    code='''import xmlrpc.client,json
from supervisor.xmlrpc import SupervisorTransport
s=xmlrpc.client.ServerProxy("http://localhost",transport=SupervisorTransport(None,None,"unix:///var/run/supervisor.sock"))
rows=[]
for name in ["mevo_f10ae_env_20261019","mevo_f10ae_optimizer_env","mevo_f10ae_prepare","mevo_f10ae_lowpressure","mevo_f10ae_lowpressure_analyze","mevo_f10ae_science","mevo_f10ae_post","mevo_f10ae_review","mevo_f10ae_supplement"]:
 try:
  r=s.supervisor.getProcessInfo(name);rows.append({k:r[k] for k in ["name","start","stop","statename","exitstatus"]})
 except Exception:pass
print(json.dumps(rows))
'''
    jobs=json.loads(subprocess.check_output(['/usr/bin/python3','-c',code],text=True));cost['supervisor_processes']=jobs
    low_jobs=sum(r['stop']-r['start'] for r in jobs if r['name'] in ['mevo_f10ae_lowpressure','mevo_f10ae_lowpressure_analyze'] and r['stop']>=r['start'])
    cost['nonoverlapping_core_scientific_wall_s']=low_jobs+read(BASE/'WORKFLOW_COMPLETE.json')['wall_s']+read(BASE/'LOW_PRESSURE_ADMISSION.json')['wall_s']+read(BASE/'AUDIT.json')['wall_s']
    if (BASE/'evidence/AUDIT_INITIAL.json').exists():cost['nonoverlapping_core_scientific_wall_s']+=read(BASE/'evidence/AUDIT_INITIAL.json')['wall_s']
    if (BASE/'LOW_DOMAIN_FOLLOWUP.json').exists():cost['nonoverlapping_core_scientific_wall_s']+=read(BASE/'LOW_DOMAIN_FOLLOWUP.json')['wall_s']
    cost['latest_figures_and_videos_wall_s']=read(BASE/'figures/MANIFEST.json')['wall_s']+sum(v['wall_s'] for v in vm)
    cost['time_scope']='Core = original radial jobs incl failed refinement + full scientific workflow + extra sub-0.1 Pa admission + full independent audit + incremental audit + single frozen low-domain followup. Excludes installs, developer time, network copies, idle job waits, visual/manual review, and repeated render versions. Artifact rendering reported separately.'
    complete=[r for r in rows if r['status']=='PASS']
    cost['solver_complete_wall_s']={k:float(v) for k,v in zip(['median','p95','max'],np.percentile([r['wall_s'] for r in complete],[50,95,100]))}
    cost['process_rss_peak_MiB']=max(r.get('rss_peak_MiB',0) for r in rows)
    cost['incomplete_fraction']=1-len(complete)/len(rows)
    cost['independently_confirmed_target_count']=sum(v['confirmed'] for v in read(BASE/'INTERPRETATION.json')['targets'].values())
    cost['confirmed_targets_per_core_hour']=cost['independently_confirmed_target_count']*3600/cost['nonoverlapping_core_scientific_wall_s']
    samples=[json.loads(s) for s in (BASE/'RESOURCES.jsonl').read_text().splitlines()]
    gpu=np.array([[float(x) for x in r['gpu'].split(',')] for r in samples])
    cost['resource_sampling']=dict(samples=len(samples),interval_s=3,gpu_utilization_mean_percent=float(gpu[:,0].mean()),gpu_utilization_max_percent=float(gpu[:,0].max()),gpu_memory_max_MiB=float(gpu[:,1].max()),cpu_count=samples[0]['cpu_count'],worker_processes=4,BLAS_threads_each=1,scope='Whole GPU snapshot, not exclusive per-process attribution; CPU RSS recorded per trajectory')
    write(BASE/'COMPUTE_COST.json',cost)
    # Complete server manifest; local copy is deliberately selective.
    files={}
    for p in BASE.rglob('*'):
        if p.is_file() and p.name!='INPUT_OUTPUT_SHA256.json' and not p.name.endswith('.lock'):files[str(p.relative_to(BASE))]=dict(sha256=sha(p),bytes=p.stat().st_size)
    write(BASE/'INPUT_OUTPUT_SHA256.json',dict(files=files,identity=identity(),excluded=['self manifest','inactive locks']))
    write(server/'DELIVERY.json',dict(identity=identity(),source_execution=str(ROOT),results=str(BASE),manifest_sha256=sha(BASE/'INPUT_OUTPUT_SHA256.json'),scientific_execution=read(BASE/'WORKFLOW_COMPLETE.json')['identity'],training_ready=False))
    print('F10AE_ARTIFACTS',len(files),'finalize wall',time.perf_counter()-start,flush=True)
if __name__=='__main__':main()
