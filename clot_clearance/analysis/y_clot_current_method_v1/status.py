from common import *
import csv,datetime
case=O/'cases/biphasic'
z=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
if (case/'history.csv').exists():
    with (case/'history.csv').open() as f:rows=[q for q in csv.DictReader(f) if None not in q and all(v is not None for v in q.values())]
    if rows:
        last=rows[-1]
        z.update({k:last[k] for k in ['step','time_s','dt_s','bulk_failed','interface_failed','detached_nodes','largest_nodes','N20_200','relative_R','wall_contact_nodes','min_wall_clearance_m','L2_active_patches']})
        z['max_relative_R']=max(float(q['relative_R']) for q in rows)
if (case/'SUMMARY.json').exists():z['formal_status']=json.loads((case/'SUMMARY.json').read_text())['status']
if (O/'FINAL_STATUS.json').exists():z['report_ready']=True
print(json.dumps(z,ensure_ascii=False),flush=True)
