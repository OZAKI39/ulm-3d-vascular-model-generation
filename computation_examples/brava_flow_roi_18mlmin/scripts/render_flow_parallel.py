"""Render the two rigid poses concurrently; each owns a distinct output folder."""
from pathlib import Path
import argparse,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--stills-only',action='store_true');args=parser.parse_args()
active=json.loads((ROOT/'reports/ACTIVE_FLOW.json').read_text());assert active['status']=='PASS'
jobs=[];start=time.time()
for pose in ['0','15']:
    log=(ROOT/'logs'/f"flow_{'preview' if args.stills_only else 'render'}_candidate_{pose}.log").open('a')
    command=[sys.executable,'-B',str(ROOT/'visualization/render_flow.py'),'--pose',pose]
    if args.stills_only:command.append('--stills-only')
    log.write(json.dumps(dict(command=command,source_flow_sha256=active['flow_sha256'],start_unix=start))+'\n');log.flush()
    p=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy())
    jobs.append((pose,p,log))
results=[]
for pose,p,log in jobs:
    code=p.wait();log.close();results.append(dict(pose=pose,pid=p.pid,exit_code=code));print(pose,'exit',code,flush=True)
record=dict(status='PASS' if all(x['exit_code']==0 for x in results) else 'FAIL',jobs=results,elapsed_s=time.time()-start,source_flow_sha256=active['flow_sha256'])
(ROOT/'reports'/('FLOW_PREVIEW_EXECUTION.json' if args.stills_only else 'FLOW_RENDER_EXECUTION.json')).write_text(json.dumps(record,indent=2)+'\n')
assert record['status']=='PASS',record
