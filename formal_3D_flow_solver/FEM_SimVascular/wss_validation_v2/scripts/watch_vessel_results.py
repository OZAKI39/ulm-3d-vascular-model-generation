"""Retrieve and independently analyze each completed, named stage-three run.

Does not start CFD or stage four, and does not accept mere solver exit success.
"""
import subprocess,sys,time,json,os,hashlib
from pathlib import Path
V=Path(__file__).resolve().parents[1]
R='/workspace/wss_validation_v2_20260927T1230Z'
names=['vessel_baseline_cpu_mpi8','vessel_baseline_gpu_mpi1_halfdt','vessel_medium']
done=[]
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
state=V/'logs/stage3_analysis_watch.json'
if state.exists():
    old=json.loads(state.read_text())
    done=[name for name in old.get('analyzed',[]) if name in names and (V/'stage3'/name/'reports/control_volume_flux.csv').exists() and (V/'evidence/local_wall_snapshots'/(name+'.npz')).exists()]
def save(**kw):
    state.write_text(json.dumps(dict(observed_unix=time.time(),analyzed=done,**kw),indent=2)+'\n')
def run(args,log):
    with (V/'logs'/log).open('w') as f:
        subprocess.run(args,env=env,cwd=V.parent,stdout=f,stderr=subprocess.STDOUT,check=True)
while len(done)<len(names):
    code='''from pathlib import Path
import json
root=Path(%r)
out={}
for name in %r:
 p=root/'stage3'/name/'reports/execution.json'
 out[name]=json.loads(p.read_text())['status'] if p.exists() else 'PENDING'
print(json.dumps(out))
'''%(R,names)
    r=subprocess.run(['ssh','vast4090','python3 -'],input=code,text=True,capture_output=True)
    if r.returncode:
        save(status='WAITING_SSH',error=r.stderr[-1500:]);time.sleep(30);continue
    observed=json.loads(r.stdout)
    for name,status in observed.items():
        if name in done or status=='PENDING':continue
        if status!='PASS':
            save(status='SOLVER_FAILURE_REQUIRES_REVIEW',case=name,remote_status=status)
            raise SystemExit(2)
        save(status='RETRIEVING_AND_ANALYZING',case=name)
        case=V/'stage3'/name
        try:
            run(['rsync','-az','vast4090:'+R+'/stage3/'+name+'/',str(case)+'/'],'retrieve_'+name+'.log')
            for script,extra in [('analyze_vessel.py',['--case',str(case)]),('mesh_metrics.py',[str(case)]),('diagnose_wall_stencil.py',['--case',str(case)]),('control_volume_checks.py',['--case',str(case)]),('export_local_wall_evidence.py',['--case',str(case)])]:
                run([sys.executable,'-B',str(V/'scripts'/script),*extra],script.removesuffix('.py')+'_'+name+'.log')
            if name=='vessel_medium':
                quality=json.loads((case/'reports/flow_quality.json').read_text());assert quality['accepted_final_and_log_checks']
                gate=dict(accepted=True,source_flow_sha256=hashlib.sha256((case/'frozen_flow/flow_arrays_si.npz').read_bytes()).hexdigest(),checked_unix=time.time(),checks=['production_WSS','independent_log_and_field','fixed_regions_and_sections','control_volume_identity','local_subset_export'])
                (case/'reports/independent_acceptance.json').write_text(json.dumps(gate,indent=2)+'\n')
                run(['rsync','-az',str(case/'reports/flow_quality.json'),str(case/'reports/independent_acceptance.json'),'vast4090:'+R+'/stage3/vessel_medium/reports/'],'medium_independent_gate_sync.log')
            done.append(name)
            run([sys.executable,'-B',str(V/'scripts/collect_results.py')],'collect_stage3_latest.log')
            run([sys.executable,'-B',str(V/'scripts/collect_geometry.py')],'geometry_stage3_latest.log')
            run([sys.executable,'-B',str(V/'scripts/compare_vessel_control_metrics.py')],'control_metric_comparisons.log')
        except Exception as e:
            save(status='ANALYSIS_FAILURE_REQUIRES_REVIEW',case=name,error=repr(e));raise
    save(status='COMPLETE' if len(done)==len(names) else 'WAITING_SOLVERS',remote_status=observed)
    if len(done)<len(names):time.sleep(30)
print('All authorized stage-three results retrieved and independently checked; fine CFD SKIPPED_BY_USER.')
