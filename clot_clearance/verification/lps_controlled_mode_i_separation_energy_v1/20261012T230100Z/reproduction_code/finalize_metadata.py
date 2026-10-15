from common import *
import time,resource,shutil,subprocess
start=time.perf_counter();protected=json.loads((O/'HISTORICAL_PROTECTION_BEFORE.json').read_text());missing=[];changed=[]
for p,v in protected.items():
 if not Path(p).is_file():missing.append(p)
 elif sha(p)!=v['sha256']:changed.append(p)
assert not missing and not changed,(missing,changed)
after={'status':'PASS','protected_files':len(protected),'protected_bytes':sum(v['bytes'] for v in protected.values()),'missing':missing,'changed':changed,'verification_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'production_library_sha256':sha(PROD),'final_doublefix_sha256':sha(B/'liblammps.so'),'env_sha256':sha('/workspace/y_clot/env.sh')};dump(O/'HISTORICAL_PROTECTION_AFTER.json',after)
status=json.loads((O/'FINAL_STATUS.json').read_text());doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text())
for key in ['PRODUCTION_LIBRARY_UNCHANGED','HISTORICAL_RESULTS_PROTECTED']:status['statuses'][key]['status']='PASS';status['statuses'][key]['scope']=f'{len(protected)}个受保护文件SHA一致，缺失0、变化0'
status['statuses']['WSL_DOCUMENT_UPDATE']={'status':'PASS_ACTUALLY_UPDATED_AND_READ_BACK','scope':'真实WSL文件并发检查、原子写入、完整旧文保留及读回通过','evidence':'WSL_DOCUMENT_UPDATE.json','after_sha256':doc['after_sha256']};status['completed_UTC']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());dump(O/'FINAL_STATUS.json',status)
repro=O/'reproduction_code';repro.mkdir(exist_ok=True)
for p in Path(__file__).parent.glob('*'):
 if p.is_file():shutil.copy2(p,repro/p.name)
phases={}
for name in ['PREPARATION_RESOURCES','ELASTIC_PRECHECK','NATIVE_STATIC_CHECK','AUDIT_RESOURCES','FIGURE_RESOURCES.initial','FIGURE_RESOURCES']:
 z=json.loads((O/(name+'.json')).read_text());phases[name]={k:z[k] for k in ['wall_s','CPU_user_s','CPU_system_s','peak_RSS_KiB'] if k in z}
ff=[p for base in [O,F,R/'analysis/mode_i_v1'] for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
dump(O/'RUN_RESOURCES.json',{'status':'COMPLETED_WITH_DESIGN_STOP','phases':phases,'measured_stage_wall_sum_s':sum(q.get('wall_s',0.) for q in phases.values()),'measured_CPU_sum_s':sum(q.get('CPU_user_s',0.)+q.get('CPU_system_s',0.) for q in phases.values()),'maximum_measured_process_RSS_MiB':max(q.get('peak_RSS_KiB',0) for q in phases.values())/1024,'whole_elapsed_wall_s':None,'whole_elapsed_reason':'Human/agent design, SSH transfer, file inspection and edits are not inside one continuously timed process. Per-stage timing retained instead.','precheck_count':1,'precheck_evaluations':1132,'formal_fracture_cases':0,'native_static_instances':1,'dynamic_steps':0,'MPI_ranks_used':1,'threads_per_process':1,'GPU_used':False,'cgroup_CPU_cores':15.36,'cgroup_RAM_bytes':45176848384,'precheck_budget':json.loads((O/'PRECHECK_FREEZE.json').read_text())['budget'],'formal_budget_reserved_not_spent':json.loads((O/'PRECHECK_FREEZE.json').read_text())['formal_candidate_if_admitted'],'formal_performance_projection':None,'projection_reason':'Design failed admission before formal performance scheduling; budget was not the reason for stopping. Static force-evaluation timings are not claimed as dynamic step timings.','output_snapshot':{'files':len(ff),'bytes':sum(p.stat().st_size for p in ff),'excludes_historical_inputs_outside_run':True,'final_manifest_is_authoritative':True},'final_protection_check_wall_s':time.perf_counter()-start})
print('METADATA_FINAL',len(protected),'historical files unchanged',flush=True)
