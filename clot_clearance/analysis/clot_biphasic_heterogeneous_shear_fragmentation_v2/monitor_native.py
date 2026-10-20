"""Read-only WSL progress monitor; never starts/restarts or alters a native case."""
import subprocess,time,json
REMOTE="""import csv,json,pathlib
o=pathlib.Path('/workspace/y_clot/project/verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z')
result={};done=True
for m in ['uniform','biphasic']:
 p=o/'cases'/m;s=p/'SUMMARY.json';h=p/'history.csv'
 if s.exists():
  q=json.loads(s.read_text());r=q['last_row'] or {};result[m]={'status':q['status'],'t_ms':q['time_s']*1000,'steps':q['accepted_steps'],'soft':r.get('softened'),'failed_bulk':r.get('bulk_failed'),'failed_interface':r.get('interface_failed'),'free':r.get('components'),'target':r.get('N20_200'),'largest':r.get('largest_nodes'),'R_percent':q['max_relative_energy_residual']*100,'stop':q['stop_reasons']}
 else:
  done=False
  if not h.exists():result[m]={'status':'NOT_STARTED'};continue
  rows=[r for r in csv.DictReader(h.open()) if r.get('relative_R') and r.get('Re_h_slip_max')];r=rows[-1] if rows else {}
  result[m]={'status':'RUNNING','t_ms':float(r.get('time_s',0))*1000,'steps':r.get('step'),'tau_Pa':r.get('shear_peak_Pa'),'soft':r.get('softened'),'failed_bulk':r.get('bulk_failed'),'failed_interface':r.get('interface_failed'),'free':r.get('components'),'target':r.get('N20_200'),'largest':r.get('largest_nodes'),'R_percent':100*float(r.get('relative_R',0))}
print(json.dumps({'cases':result,'both_done':done}))
"""
while True:
 try:
  p=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15','-p','15283','-i','/home/lzy/.ssh/vast_step3b_ed25519','root@95.3.33.46','python3 -'],input=REMOTE,text=True,capture_output=True,timeout=30)
  if p.returncode:print('MONITOR_CONNECTION_ERROR',p.returncode,p.stderr[-500:],flush=True)
  else:
   row=json.loads(p.stdout.strip().splitlines()[-1]);print(json.dumps(row),flush=True)
   if row['both_done']:break
 except Exception as e:print('MONITOR_READ_ERROR',repr(e),flush=True)
 time.sleep(50)
