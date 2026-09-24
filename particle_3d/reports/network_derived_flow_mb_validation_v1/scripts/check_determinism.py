from pathlib import Path
import json,gzip,hashlib
import numpy as np
R=Path(__file__).resolve().parents[1]
rows=[]
for pid in [3,15,22]:
    paths=[R/'outputs'/label/'trajectories'/f'mb_{pid:06d}' for label in ['determinism_w1','determinism_w3']]
    a,b=[json.loads(p.with_suffix('.json').read_text()) for p in paths]
    x,y=[np.load(p.with_suffix('.npz'))['samples'] for p in paths]
    logs=[gzip.decompress(p.with_suffix('.audit.jsonl.gz').read_bytes()) for p in paths]
    row=dict(bubble_id=pid,status_equal=a['end_reason']==b['end_reason'] and a['failure_detail']==b['failure_detail'],
        outlet_equal=a['exit_outlet']==b['exit_outlet'],final_position_bitwise_equal=x[-1,1:4].tobytes()==y[-1,1:4].tobytes(),
        transit_time_equal=a['residence_time_s']==b['residence_time_s'],all_samples_bitwise_equal=x.shape==y.shape and x.tobytes()==y.tobytes(),
        event_log_bitwise_equal=logs[0]==logs[1],event_log_sha256=[hashlib.sha256(v).hexdigest() for v in logs],
        samples_shape=list(x.shape))
    row['pass']=all(row[k] for k in ['status_equal','outlet_equal','final_position_bitwise_equal','transit_time_equal','all_samples_bitwise_equal','event_log_bitwise_equal'])
    rows.append(row)
result=dict(status='PASS' if all(r['pass'] for r in rows) else 'FAIL',workers=[1,3],N=3,ids=[3,15,22],
    selection_rule='First historical smoke ID; historical stationary ID 15; ID 22 with historical inlet/handoff review',rows=rows)
(R/'data/worker_count_determinism_subset.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2));assert result['status']=='PASS'
