"""Gate formal production and reuse matching, already integrated benchmark tracks."""
from pathlib import Path
import json,subprocess,sys,shlex,xml.etree.ElementTree as ET
from resource_plan import make_plan
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import *
R=ROOT/REL;c=json.loads((R/'data/run_context.json').read_text());remote=c['remote']
def passed(p):
 root=ET.parse(p).getroot();suites=[root] if root.tag=='testsuite' else list(root)
 return all(int(s.get(k,0))==0 for s in suites for k in ['failures','errors','skipped'])
assert passed(R/'logs/baseline_tests.xml') and passed(R/'logs/new_tests.xml')
# Retrieve the completed benchmark first; no provisional throughput selection.
subprocess.run(['scp','-q','vast4090:'+remote+'/'+REL+'/data/worker_scaling_results.json',str(R/'data/worker_scaling_results.json')],check=True)
b=json.loads((R/'data/worker_scaling_results.json').read_text());assert b['benchmark_science_identical'] and b['exact_same_dt_reference_24_trajectory_parity'] and b['numerical_safety_pass'] and b['dt_s']==DT
subprocess.run(['scp','-q','vast4090:'+remote+'/'+REL+'/outputs/benchmark/workers_01/batch_result.json',str(R/'data/benchmark_first_worker_result.json')],check=True)
free_bytes=int(subprocess.check_output(['ssh','-o','BatchMode=yes','vast4090','python3','-c',shlex.quote('import shutil;print(shutil.disk_usage("'+remote+'").free)')],text=True))
plan=make_plan(b,json.loads((R/'data/benchmark_first_worker_result.json').read_text()),free_bytes)
write_new(R/'data/RESOURCE_PLAN.json',plan)
subprocess.run(['scp','-q',str(R/'data/RESOURCE_PLAN.json'),'vast4090:'+remote+'/'+REL+'/data/RESOURCE_PLAN.json'],check=True)
# A controller-only preproduction correction preserves the original deployment
# and benchmark evidence. No trajectory kernel/state/identity changes.
subprocess.run(['scp','-q',str(R/'scripts/run_server.py'),'vast4090:'+remote+'/run_server_resource_revision.pending.py'],check=True)
revision="""from pathlib import Path
import hashlib,json,datetime
root=Path(REMOTE);manifest=root/'deployment_manifest.json';m=json.loads(manifest.read_text());rel=REL+'/scripts/run_server.py';old=root/rel;new=root/'run_server_resource_revision.pending.py'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(old)==m[rel]
(root/'logs/benchmark_controller_original.py').write_bytes(old.read_bytes())
(root/'logs/deployment_manifest_before_resource_plan.json').write_bytes(manifest.read_bytes())
record=dict(reason='Preproduction storage estimate must include adaptive substeps and rejected-trial logs',old_controller_sha256=sha(old),new_controller_sha256=sha(new),scientific_kernel_unchanged=True,trajectory_identity_unchanged=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
new.replace(old);m[rel]=sha(old);manifest.write_text(json.dumps(m,indent=2)+'\\n')
(root/'logs/preproduction_controller_revision.json').write_text(json.dumps(record,indent=2)+'\\n')
""".replace('REMOTE',repr(remote)).replace('REL',repr(REL))
subprocess.run(['ssh','-o','BatchMode=yes','vast4090','python3','-'],input=revision,text=True,check=True)

gate=dict(identity=b['identity'],new_tests_pass=True,baseline_tests_pass=True,
 baseline_tests_sha256=digest(R/'logs/baseline_tests.xml'),new_tests_sha256=digest(R/'logs/new_tests.xml'),
 core_cohort_sha256=digest(R/'data/CORE500_COHORT.json'),resource_policy_sha256=digest(R/'data/RESOURCE_PLAN.json'))
write_new(R/'data/production_gate.json',gate)
subprocess.run(['scp','-q',str(R/'data/production_gate.json'),'vast4090:'+remote+'/'+REL+'/data/production_gate.json'],check=True)
# The point script only runs after the final MB cohort has been fixed.
subprocess.run(['scp','-q',str(R/'scripts/run_points.py'),'vast4090:'+remote+'/'+REL+'/scripts/run_points.py'],check=True)
launch='''from pathlib import Path
import json,os,sys,shutil,subprocess
root=Path(REMOTE);sys.path.insert(0,str(root/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import completion_matches,write_new,digest
r=root/REL;b=json.loads((r/'data/worker_scaling_results.json').read_text());b.update(json.loads((r/'data/RESOURCE_PLAN.json').read_text()));core=json.loads((r/'data/CORE500_COHORT.json').read_text())
reused=[]
for e in core['events'][:24]:
 source=r/'outputs/benchmark'/f"workers_{b['workers_production']:02d}"/'tracks'/f"mb_{e['particle_id']:06d}"
 target=r/'outputs/formal/core500/tracks'/source.name
 if not completion_matches(source,b['identity'],e):raise ValueError('Cannot reuse unmatched benchmark trajectory')
 if not target.exists():target.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(source,target)
 if not completion_matches(target,b['identity'],e):raise ValueError('Copied benchmark trajectory differs')
 reused.append(dict(particle_id=e['particle_id'],source=str(source),target=str(target),completion_sha256=digest(target/'COMPLETE.json')))
write_new(r/'data/benchmark_reuse.json',dict(count=len(reused),reason='Frozen CORE500 preexisted benchmark; exact same source/state/dt/horizon identity; do not recompute valid completed trajectories',rows=reused))
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
script=root/'launch_formal_and_points.sh'
script.write_text('#!/bin/bash\\nset -euo pipefail\\ncd '+str(root)+'\\n/root/particle8_2_runs/env/bin/python '+str(r/'scripts/run_server.py')+' --stage production\\n/root/particle8_2_runs/env/bin/python '+str(r/'scripts/run_points.py')+'\\n')
f=(root/'logs/production.txt').open('xb')
p=subprocess.Popen(['bash',str(script)],cwd=root,env=env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
(root/'logs/production_pid.json').write_text(json.dumps(dict(pid=p.pid,production_workers=b['workers_production'],resource_N_max=b['resource_limited_N_max']))+'\\n');print(json.dumps(dict(pid=p.pid,workers=b['workers_production'],maximum_N=b['resource_limited_N_max'],reused=len(reused))))
'''.replace('REMOTE',repr(remote)).replace('REL',repr(REL))
subprocess.run(['ssh','-o','BatchMode=yes','vast4090','python3','-'],input=launch,text=True,check=True)
