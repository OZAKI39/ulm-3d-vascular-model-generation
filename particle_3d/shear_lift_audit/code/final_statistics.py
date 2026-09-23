"""Server-only complementary summaries; no force or trajectory changes."""
from pathlib import Path
import sys,json,csv,socket
import numpy as np
from run_audit import METRICS,dump
from audit_math import stats,force_ratio

def main(config):
    cfg=json.loads(Path(config).read_text());out=Path(cfg['output'])
    with np.load(out/'all_scalar_states.npz') as f:d={k:f[k] for k in f.files}
    masks={'whole':np.ones(len(d['state']),bool),'accepted_intervals':~d['initial_diagnostic'],
           'initial_diagnostics':d['initial_diagnostic'],
           **{k[7:]:v for k,v in d.items() if k.startswith('region_')}}
    bins={};guards={}
    for group,m in masks.items():
        guards[group]={key:{str(c):int(np.sum(d[key][m]==c)) for c in range(4)}
                       for key in ['lift_drag_status','lift_lubrication_status']}
        bins[group]={}
        for key in ['lift_drag_ratio','lift_lubrication_ratio']:
            x=d[key][m];valid=np.isfinite(x);x=x[valid]
            bins[group][key]=dict(resolved=int(valid.sum()),undefined=int((~valid).sum()),
              lt_0p01=int(np.sum(x<.01)),from_0p01_to_0p1=int(np.sum((x>=.01)&(x<.1))),
              from_0p1_to_1=int(np.sum((x>=.1)&(x<=1))),gt_1=int(np.sum(x>1)))
    accepted={k:{metric:stats(d[metric][m]) for metric in METRICS} for k,m in masks.items()
              if k in ['accepted_intervals','initial_diagnostics']}
    sensitivity={str(f):[] for f in [.1,1.,10.]}
    receipts=json.loads((out/'trajectory_audit_receipts.json').read_text())
    for r in receipts:
        if not r['count']:continue
        with np.load(r['file']) as z:
            m=z['region_near_wall'];a=z['CANDIDATE_SAFFMAN_MAGNITUDE'][m];b=z['lubrication_N'][m]
            ta=z['lift_numerical_zero_N'][m];tb=z['lubrication_numerical_zero_N'][m]
            for factor in sensitivity:
                ratio,status=force_ratio(a,b,tb*float(factor),ta*float(factor))
                sensitivity[factor].append(ratio)
    sensitivity={factor:stats(np.concatenate(parts)) for factor,parts in sensitivity.items()}
    result=dict(compute_host=socket.gethostname(),label='DESCRIPTIVE_ONLY_NOT_VALIDATED_INCLUSION_CRITERIA',
                status_labels={'0':'RESOLVED','1':'BOTH_FORCE_MAGNITUDES_NEAR_ZERO',
                               '2':'DENOMINATOR_NEAR_ZERO_NOT_EVALUABLE','3':'NO_ACTIVE_LUBRICATION'},
                ratio_bins=bins,ratio_guard_counts=guards,initial_and_accepted_statistics=accepted,
                numerical_guard_multiplier_sensitivity=sensitivity,
                numerical_guard_note='Self-scaled single-wall matrix condition only; not a rigorous error bound for multi-contact constrained solves.')
    dump(out/'supplementary_statistics.json',result)
    with (out/'descriptive_ratio_bins.csv').open('w',newline='') as f:
        writer=csv.writer(f);keys=list(bins['whole']['lift_drag_ratio'])
        writer.writerow(['group','metric','role']+keys)
        for group,rows in bins.items():
            for metric,v in rows.items():writer.writerow([group,metric,result['label']]+[v[k] for k in keys])
    with (out/'accepted_interval_statistics.csv').open('w',newline='') as f:
        writer=csv.writer(f);keys=['count','median','P90','P95','P99','max'];writer.writerow(['group','metric']+keys)
        for group,rows in accepted.items():
            for metric,v in rows.items():writer.writerow([group,metric]+[v[k] for k in keys])
    print(json.dumps(sensitivity,indent=2))

if __name__=='__main__':main(sys.argv[1])
