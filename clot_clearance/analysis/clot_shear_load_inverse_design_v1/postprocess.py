"""Run read-only reductions after both native cases finish; no extra dynamics."""
from common import *
import subprocess,sys,time,resource
while not all((O/'cases'/name/'SUMMARY.json').exists() for name in CASES):time.sleep(10)
for name in CASES:assert json.loads((O/'cases'/name/'SUMMARY.json').read_text())['status']=='COMPLETED'
records=[];figures_only='--figures-only' in sys.argv
scripts=['surface_coverage_audit.py','figures.py','report.py'] if figures_only else ['analyze.py','validate_outputs.py','surface_coverage_audit.py','figures.py','report.py']
if '--render-only' in sys.argv:scripts=['figures.py','report.py']
resource_name='POSTPROCESS_RENDER_REVIEW_RESOURCES' if '--render-only' in sys.argv else 'POSTPROCESS_FINALIZATION_RESOURCES' if figures_only else 'POSTPROCESS_RESOURCES'
resource_path=O/(resource_name+'.json');index=1
while resource_path.exists():
    index+=1;resource_path=O/(resource_name+f'_pass{index}.json')
for script in scripts:
    start=time.perf_counter();print('START',script,flush=True)
    result=subprocess.run([sys.executable,str(A/script)],cwd=R)
    records.append(dict(script=script,returncode=result.returncode,wall_s=time.perf_counter()-start))
    usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    dump(resource_path,dict(stages=records,cumulative_child_cpu_s=usage.ru_utime+usage.ru_stime,max_child_RSS_bytes=usage.ru_maxrss*1024,scope='Read-only reductions and figure rendering; no native physics rerun'))
    if result.returncode:raise SystemExit(result.returncode)
    print('DONE',script,records[-1]['wall_s'],flush=True)
print('POSTPROCESS_COMPLETE',flush=True)
