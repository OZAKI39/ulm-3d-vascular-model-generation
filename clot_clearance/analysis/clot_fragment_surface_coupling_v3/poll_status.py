"""Cheap one-shot monitoring of complete accepted CSV rows and process health."""
from pathlib import Path
import csv,json,subprocess,datetime,os

ROOT=Path(os.environ.get('Y_CLOT_PROJECT_ROOT','/workspace/y_clot/project'))
OUT=ROOT/'verification/clot_fragment_surface_coupling_v3/20261018T170400Z'
case=OUT/'cases/biphasic';request=json.loads((OUT/'FORMAL_REQUEST.json').read_text())
status=subprocess.run(['supervisorctl','status',request['name']],capture_output=True,text=True).stdout.strip()
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),supervisor=status)
p=case/'history.csv'
if p.exists():
    with p.open('rb') as f:
        header=next(csv.reader([f.readline().decode()]));f.seek(0,2);size=f.tell();start=max(0,size-64000);f.seek(start);tail=f.read()
    lines=tail.split(b'\n')[:-1]
    if start:lines=lines[1:]
    for line in reversed(lines):
        values=next(csv.reader([line.decode()]))
        if len(values)!=len(header):continue
        row=dict(zip(header,values))
        try:
            result.update(accepted_step=int(row['step']),physical_time_ms=float(row['time_s'])*1000,completion_fraction=float(row['time_s'])/.0036,
                target_20_200=int(row['N20_200']),largest_free_nodes=int(row['largest_nodes']),detached_nodes=int(row['detached_nodes']),
                softened_bonds=int(row['softened']),relative_energy_residual=float(row['relative_R']),
                L2_active_patches=int(row['L2_active_patches']),L2_effective_area_m2=float(row['L2_effective_area_m2']))
        except (ValueError,KeyError):continue
        break
pid=subprocess.run(['supervisorctl','pid',request['name']],capture_output=True,text=True).stdout.strip()
if pid.isdigit() and int(pid)>0:
    result['process']=subprocess.run(['ps','-p',pid,'-o','pid=,stat=,etime=,time=,%cpu=,rss='],capture_output=True,text=True).stdout.strip()
summary=case/'SUMMARY.json'
if summary.exists():
    q=json.loads(summary.read_text());result.update(formal_status=q['status'],stop_reasons=q['stop_reasons'],formal_done=True)
else:result['formal_done']=False
result['post_prepare_done']=(OUT/'POST_PREPARE_COMPLETE.json').exists()
result['post_prepare_status']=subprocess.run(['supervisorctl','status','clot_surface_v3_post_prepare_20261018T170400Z'],capture_output=True,text=True).stdout.strip()
print('MONITOR '+json.dumps(result,ensure_ascii=False),flush=True)
