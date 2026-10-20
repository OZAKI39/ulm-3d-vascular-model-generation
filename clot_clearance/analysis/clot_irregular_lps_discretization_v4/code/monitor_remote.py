"""Read-only remote progress; retain a small local recovery pointer."""
from pathlib import Path
import datetime,json,subprocess,sys

D=Path(__file__).resolve().parents[1]
SSH=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15','-p','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','root@95.3.33.46']
REMOTE=r'''
import csv,json,pathlib,subprocess
o=pathlib.Path('/workspace/y_clot/project/verification/clot_irregular_lps_discretization_v4/20261019T221000Z')
c=o/'cases/H_IRREGULAR_LPS_V4'
out={}
q=subprocess.run(['supervisorctl','status','clot_v4_formal_20261019T221000Z','clot_v4_post_20261019T221000Z'],capture_output=True,text=True)
out['services']=q.stdout.strip().splitlines()
if (c/'history.csv').exists():
 with (c/'history.csv').open() as f:
  last=None
  for row in csv.DictReader(f):
   if row.get('relative_R') is not None:last=row
 if last:
  keys=['step','time_s','bulk_failed','interface_failed','softened','detached_nodes','largest_nodes','N20_200','relative_R','rejects']
  out['latest_accepted']={k:last.get(k) for k in keys}
for name,where in [('checkpoint',c/'LAST_LEGAL_CHECKPOINT.json'),('native_result',o/'NATIVE_V4_RESULT.json'),('postprocess',o/'POSTPROCESS_STATUS.json')]:
 if where.exists():
  q=json.loads(where.read_text())
  out[name]={k:v for k,v in q.items() if k in ['path','time_s','step','status','stage','accepted_time_s','exception','traceback']}
print(json.dumps(out))
'''

def main():
 q=subprocess.run(SSH+['python3','-'],input=REMOTE,capture_output=True,text=True,timeout=45)
 if q.returncode:raise RuntimeError(q.stderr)
 data=json.loads(q.stdout);stamp=datetime.datetime.now(datetime.timezone.utc).isoformat();data['observed_utc']=stamp
 path=D/'RESUME_STATUS.json';record=json.loads(path.read_text());record['recorded_utc']=stamp
 if 'latest_accepted' in data:record['last_observed_time_s']=float(data['latest_accepted']['time_s'])
 record['latest_monitor']=data
 if 'native_result' in data:record['status']=data['native_result']['status']
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n');tmp.replace(path)
 print(json.dumps(data,ensure_ascii=False))
if __name__=='__main__':main()
