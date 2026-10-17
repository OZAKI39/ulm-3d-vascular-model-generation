import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *

def score(rows):
    result={'n':len(rows),'all_six_valid_fraction':np.mean([all(r['estimate']['valid'].values()) for r in rows]) if rows else None,'per_parameter':{},'algorithm_ms':stats([r['algorithm_ms'] for r in rows])}
    bad_global=[r for r in rows if r.get('truth') is None and r.get('case')=='local_opposite']
    result['global_nonaffine_correct_rejection_fraction']=np.mean([not any(r['estimate']['valid'].values()) for r in bad_global]) if bad_global else None
    result['global_nonaffine_n']=len(bad_global)
    for j,k in enumerate(NAMES):
        usable=[r for r in rows if r.get('truth') is not None and r['estimate']['valid'][k] and r['estimate']['values'][k] is not None]
        err=np.array([r['estimate']['values'][k]-r['truth'][j] for r in usable]);limit=P['parameter_error_limits'][j]
        nz=[r for r in usable if abs(r['truth'][j])>1e-7]
        ci=[]
        for r in usable:
            u=r['estimate'].get('uncertainty')
            if u and 'intervals' in u:
                low,high=u['intervals'][k];ci.append(low<=r['truth'][j]<=high)
        result['per_parameter'][k]={'valid_n':len(usable),'mae':np.mean(np.abs(err)) if len(err) else None,'bias':np.mean(err) if len(err) else None,'error_std':np.std(err) if len(err) else None,'abs_error_p95':np.percentile(np.abs(err),95) if len(err) else None,
          'sign_correct_fraction_nonzero_truth':np.mean([np.sign(r['estimate']['values'][k])==np.sign(r['truth'][j]) for r in nz]) if nz else None,'nonzero_truth_valid_n':len(nz),
          'false_accept_count':int(np.sum(np.abs(err)>limit)),'false_accept_fraction_of_accepted':np.mean(np.abs(err)>limit) if len(err) else None,
          'nominal95_coverage':np.mean(ci) if ci else None,'CI_n':len(ci),'error_limit':limit,'unit':P['parameter_error_units'][j]}
    accepted=[r for r in rows if r.get('truth') is not None and all(r['estimate']['valid'].values())]
    failures=[r for r in accepted if np.any(np.abs(np.array([r['estimate']['values'][k] for k in NAMES])-r['truth'])>P['parameter_error_limits'])]
    result['false_accept_any_parameter_n']=len(failures);result['false_accept_any_parameter_fraction_of_accepted']=len(failures)/len(accepted) if accepted else None
    result['false_accept_cases']=[{k:r[k] for k in ['case','seed','window','end'] if k in r} for r in failures]
    return result

test=json.loads((ROOT/'results/benchmark_test.json').read_text())['rows'];methods=sorted(set(r['method'] for r in test));summary={}
for m in methods:
    rr=[r for r in test if r['method']==m]
    summary[m]={'overall':score(rr),'windows':{str(n):score([r for r in rr if r['window']==n]) for n in P['windows']},'cases':{c:score([r for r in rr if r['case']==c]) for c in P['cases']}}
ai=ROOT/'results/ai_comparison.json'
if ai.exists():
    ar=json.loads(ai.read_text())
    if ar['status']=='COMPLETE':
        arows=[r for r in ar['rows'] if r['truth'] is not None];summary['SEA_RAFT_S']={'overall':score(arows),'windows':{str(n):score([r for r in arows if r['window']==n]) for n in P['windows']}}
        keys={(r['case'],r['seed'],r['window'],r['end']) for r in arows}
        summary['matched_AI_comparison']={m:score([r for r in test if r['method']==m and (r['case'],r['seed'],r['window'],r['end']) in keys]) for m in methods}
        summary['matched_AI_comparison']['SEA_RAFT_S']=score(arows)
real=json.loads((ROOT/'results/real_observation.json').read_text());rr=real['rows'];real_summary={}
for name in sorted(set(k for r in rr for k in r['methods'])):
    values=[r['methods'][name] for r in rr if name in r['methods']]
    real_summary[name]={'n':len(values),'valid_fraction':np.mean([all(v['estimate']['valid'].values()) for v in values]),'algorithm_ms':stats([v['algorithm_ms'] for v in values]),'per_parameter_unobservable_fraction':{k:1-np.mean([v['estimate']['valid'][k] for v in values]) for k in NAMES}}
response=json.loads((ROOT/'results/window_response.json').read_text());lags={}
for n in P['windows']:
    rrn=[r for r in response['rows'] if r['window']==n and r['end_t']>response['onset_t']]
    hit=[r for r in rrn if r['estimate']['values']['translation_vx'] is not None and r['estimate']['values']['translation_vx']>=27]
    lags[str(n)]={'first_90pct_step_response_s':hit[0]['end_t']-response['onset_t'] if hit else None,'endpoint_span_s':(n-1)*P['synthetic_dt_s'],'midpoint_age_s':(n-1)*P['synthetic_dt_s']/2}
dump(ROOT/'results/metrics.json',{'truth_tests':summary,'real_recording':real_summary,'step_response':lags,'definitions':{'error_limits':P['parameter_error_limits'],'false_accept':'Any accepted parameter exceeds preregistered absolute error limit; numerator/accepted denominator explicit. Not a universal scientific accuracy certificate.','correct_reject':'Only local_opposite cases have no single global affine truth. Other difficult-case rates reported separately.','CI':'Empirical spatial-block bootstrap nominal 95%; coverage measured, not assumed.'}})
print(json.dumps({'PIV':summary['PIV_multipass']['overall'],'CPD':summary['CPD']['overall'],'real':real_summary,'lags':lags},indent=2))
