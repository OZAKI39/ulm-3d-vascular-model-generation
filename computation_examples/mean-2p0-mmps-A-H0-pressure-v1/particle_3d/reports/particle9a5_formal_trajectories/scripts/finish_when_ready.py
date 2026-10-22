"""Collect, audit, render and regress the complete formal run; never claim launch as completion."""
from pathlib import Path
import json,subprocess,sys,time,datetime
ROOT=Path(__file__).resolve().parents[4];R=ROOT/'particle_3d/reports/particle9a5_formal_trajectories';c=json.loads((R/'data/run_context.json').read_text())
probe='''from pathlib import Path
import json
root=Path(REMOTE);r=root/'particle_3d/reports/particle9a5_formal_trajectories';done=r/'data/point_complete.json';pidfile=root/'logs/production_pid.json'
state=dict(complete=done.exists(),started=pidfile.exists(),alive=False)
if pidfile.exists():
 pid=json.loads(pidfile.read_text())['pid'];p=Path('/proc')/str(pid)
 state['alive']=p.exists() and str(root).encode() in (p/'cmdline').read_bytes()
print(json.dumps(state))
'''.replace('REMOTE',repr(c['remote']))
errors=0
while True:
 try:
  result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15','vast4090','python3','-'],input=probe,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True,timeout=40)
  state=json.loads(result.stdout);errors=0
 except Exception as e:
  errors+=1;print('CONNECTION_RETRY',errors,str(e),flush=True)
  if errors>=5:raise
  time.sleep(45);continue
 if state['complete']:break
 if state['started'] and not state['alive']:raise RuntimeError('Production stopped before posthoc point completion; inspect preserved records')
 time.sleep(45)
print('FORMAL_AND_POINT_COMPLETE',datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
commands=['collect.py','compare_dt.py','analyze_formal.py','render_formal.py','render_animations.py','final_regression.py','verify_protection.py','analyze_formal.py','build_review.py','verify_delivery.py']
for i,script in enumerate(commands):
 print('RUNNING',script,flush=True)
 # Delivery manifest is written last and must not hash its own changing stdout log.
 if script=='verify_delivery.py':subprocess.run([sys.executable,str(R/'scripts'/script)],cwd=ROOT,check=True)
 else:
  with (R/'logs'/f'finish_{i:02d}_{script[:-3]}.txt').open('x') as log:subprocess.run([sys.executable,str(R/'scripts'/script)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
print('AUTOMATED_DELIVERY_CHECKS_COMPLETE_MANUAL_VISUAL_REVIEW_PENDING',flush=True)
