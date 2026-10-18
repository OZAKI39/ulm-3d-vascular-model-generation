from common import *
import time,subprocess,sys,datetime

phase=sys.argv[1]
while not (O/'cases/biphasic/SUMMARY.json').exists():time.sleep(15)
if phase=='prepare':commands=[['analyze.py'],['compare_motion.py'],['coverage_audit.py'],['figures.py'],['render_videos.py','--smoke']]
else:
    while not (O/'SMOKE_VISUAL_REVIEW.json').exists():time.sleep(15)
    review=json.loads((O/'SMOKE_VISUAL_REVIEW.json').read_text());assert review['status']=='PASS'
    commands=[['render_videos.py'],['report.py']]
for args in commands:
    source=A/args[0];start=datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (O/'POST_EXECUTIONS.jsonl').open('a') as f:f.write(json.dumps(dict(start_utc=start,phase=phase,args=args,source_sha256=sha(source)))+'\n')
    print('POST_START',start,args,flush=True);subprocess.run([sys.executable,str(source),*args[1:]],cwd=R,check=True);print('POST_DONE',args,flush=True)
dump(O/f'POST_{phase.upper()}_COMPLETE.json',dict(status='PASS',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),physical_runs_started=0))
