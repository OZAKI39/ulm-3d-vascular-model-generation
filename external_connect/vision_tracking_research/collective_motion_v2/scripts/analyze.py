import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from vcommon import OUT,DATA,P,dump,identity,LIMITS,NAMES,sha,clean
from metrics import summarize

def main():
    frozen=json.loads((OUT/'results/model_freeze.json').read_text());evaluated=json.loads((OUT/'results/final_evaluation_complete.json').read_text());assert evaluated['freeze_sha256']==sha(OUT/'results/model_freeze.json')
    truth=np.load(DATA/'final_test/truth.npy');label=np.load(DATA/'final_test/label.npy');meta=[json.loads(s) for s in (DATA/'final_test/index.jsonl').read_text().splitlines()]
    outputs={};failfile=OUT/'results/all_failures.jsonl'
    with failfile.open('w') as fw:
        for row in evaluated['evaluations']:
            name=row['model'];key=name+('_'+row['corruption'] if row['corruption'] else '');p=dict(np.load(OUT/'results'/row['file']));r=frozen['models'][name]
            predclass=p['probability'].argmax(1);eligible=p['algebra_valid']&(predclass==0);accepted=eligible&(p['score']<=r['threshold']) if r['threshold'] is not None else np.zeros(len(label),bool)
            outputs[key]=summarize(p['mean'],truth,label,accepted,predclass,p['score'],eligible,meta,p['sigma'],r['quantiles'])
            bad=(label!=0)|~np.isfinite(p['mean']).all(1)|np.any(np.abs(p['mean']-truth)>LIMITS,axis=1)
            for i in np.where(bad)[0]:
                fw.write(json.dumps(clean(dict(method=key,index=int(i),case=meta[i]['case'],source_frame=meta[i]['source_frame'],seed=meta[i]['seed'],true_class=int(label[i]),accepted=bool(accepted[i]),truth=truth[i].tolist(),prediction=p['mean'][i].tolist(),error=(p['mean'][i]-truth[i]).tolist(),score=float(p['score'][i]))),allow_nan=False)+'\n')
        for method in ['PIV_SINGLE','PIV_MULTIPASS','CPD']:
            file=OUT/'results'/f'{method}_final_test.json'
            if not file.exists():outputs[method]={'status':'NOT_RUN'};continue
            raw=json.loads(file.read_text());rows=raw['rows'];idx=np.array([r['index'] for r in rows]);mean=np.array([[np.nan if x is None else x for x in r['mean']] for r in rows]);accepted=np.array([r['accepted'] for r in rows]);predclass=np.array([r['predicted_class'] for r in rows]);score=np.array([r['score'] for r in rows]);eligible=accepted
            outputs[method]=summarize(mean,truth[idx],label[idx],accepted,predclass,score,eligible,[meta[i] for i in idx]);outputs[method]['subset_indices']=idx.tolist();outputs[method]['execution']=dict(wall_s=raw['wall_s'],workers=raw['parallel_workers'],scope=raw['scope'])
    matched={}
    for coverage in [.1,.25,.5,.75]:
        matched[str(coverage)]={}
        for name,r in outputs.items():
            curve=r.get('risk_coverage',{});n=int(round(coverage*len(label)))
            if len(curve.get('risk',[]))>=n and r['n_total']==len(label):matched[str(coverage)][name]=dict(n=n,risk=curve['risk'][n-1])
            else:matched[str(coverage)][name]=dict(status='COVERAGE_NOT_REACHABLE_OR_SUBSET')
    report=dict(identity=identity(),freeze_sha256=evaluated['freeze_sha256'],limits=LIMITS.tolist(),units=P['units'],generalization_scope=P['generalization_scope'],champions=frozen['champions_by_validation_only'],methods=outputs,matched_coverage=matched,
        statistical_scope='Seeds are synthetic realizations from one recording, not independent physical experiments. Four source frame blocks permit descriptive variability only.',real_accuracy='NOT_VERIFIED',training_ready=False)
    dump(OUT/'results/analysis.json',report)
    print({k:{x:v.get(x) for x in ['n_total','n_accepted','accepted_failure_fraction','mae_all_finite_global']} for k,v in outputs.items()})
if __name__=='__main__':main()
