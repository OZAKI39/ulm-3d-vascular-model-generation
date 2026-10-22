"""Immutable final-test metrics. Diagnostic curves do not tune publication gates."""
import sys,json,collections
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from benchmark import LIMITS
from figures_full import plt,save,LABELS
METHODS=['PIV','CPD','POT_partial','POT_unbalanced','LRAFT_ST','CoTracker3','UnLiteFlowNet']
SHORT=['PIV','CPD','Partial OT','Unbalanced OT','LRAFT-ST','CoTracker3','UnLiteFlowNet']
COLORS=['#5064a7','#a0a0a0','#9776b5','#598d79','#b7814c','#447b9b','#b9646b']
def metrics(rows,selector):
    accepted=[r for r in rows if r[selector]];bad=sum(r['bad_if_accepted'] for r in accepted)
    return dict(n_total=len(rows),accepted=len(accepted),bad_accepted=bad,risk=bad/len(accepted) if accepted else None,coverage=len(accepted)/len(rows))
def main():
    data={m:json.loads((OUT/'results'/f'final_{m}.json').read_text()) for m in METHODS};analysis={};curves={};conditions={}
    primary=['PIV','POT_unbalanced','LRAFT_ST','CoTracker3'];common=sorted(set.intersection(*[{r['index'] for r in data[m] if r['class_id']==0 and r['absolute_error'] is not None} for m in primary]));common_mae={}
    for m in METHODS:
        rows=data[m];err=np.asarray([r['absolute_error'] for r in rows if r['class_id']==0 and r['absolute_error'] is not None]);rank=sorted([r for r in rows if r['base_valid']],key=lambda r:-r['quality_score']);risk=np.cumsum([r['bad_if_accepted'] for r in rank])/np.arange(1,len(rank)+1);coverage=np.arange(1,len(rank)+1)/len(rows);curves[m]=dict(coverage=coverage,risk=risk,scope='posthoc descriptive quality ranking; not new thresholds')
        matched={}
        for c in [.1,.25,.5]:
            k=int(np.ceil(c*len(rows)));matched[str(c)]=dict(n=k,risk=float(risk[k-1])) if len(rank)>=k else dict(n=0,risk=None,status='INSUFFICIENT_BASE_COVERAGE')
        translation=[]
        for r in rows:
            if r['case'].startswith('translation_') and r['raw'] is not None:
                j=int('_y_' in r['case']);truth=r['truth'][j];pred=r['raw'][j];translation.append(dict(index=r['index'],case=r['case'],truth=truth,pred=pred,sign_correct=pred*truth>0,amplitude_ratio=abs(pred/truth),near_zero=abs(pred/truth)<.2))
        conditions[m]={}
        for case in sorted(set(r['case'] for r in rows)):
            rr=[r for r in rows if r['case']==case];conditions[m][case]=dict(base=metrics(rr,'base_valid'),gate=metrics(rr,'risk_gate_valid'))
        analysis[m]=dict(base=metrics(rows,'base_valid'),gate=metrics(rows,'risk_gate_valid'),finite_global_mae_n=len(err),finite_global_mae=err.mean(0) if len(err) else None,matched_coverage=matched,translation_diagnostics=translation,adapter_ms=stats([r['timing'].get('adapter_total_ms',np.nan) for r in rows]),model_ms=stats([r['timing'].get('model_ms',np.nan) for r in rows]),global_missing_raw=sum(r['class_id']==0 and r['raw'] is None for r in rows))
        if m in primary:common_mae[m]=np.mean([r['absolute_error'] for r in rows if r['index'] in common],0)
    fig,ax=plt.subplots(2,2,figsize=(11,8))
    for k,m in enumerate(METHODS):
        trs=analysis[m]['translation_diagnostics'];ax[0,0].scatter([r['truth'] for r in trs],[r['pred'] for r in trs],s=15,alpha=.6,color=COLORS[k],label=SHORT[k]);ax[0,1].bar(k,np.mean([r['sign_correct'] for r in trs]) if trs else 0,color=COLORS[k]);ax[1,0].bar(k,np.median([r['amplitude_ratio'] for r in trs]) if trs else 0,color=COLORS[k]);ax[1,1].bar(k,sum(r['near_zero'] for r in trs),color=COLORS[k])
    ax[0,0].plot([-40,40],[-40,40],ls='--',color='black',lw=1);ax[0,0].legend(fontsize=8,ncol=2);ax[0,0].set(xlabel='Known translation (px/s)',ylabel='Diagnostic raw estimate (px/s)')
    for a,label in [(ax[0,1],'Translation sign recovery'),(ax[1,0],'Median |prediction / truth|'),(ax[1,1],'Near-zero (<0.2 × truth) cases')]:a.set(xticks=range(7),xticklabels=SHORT,ylabel=label);a.tick_params(axis='x',rotation=35)
    save(fig,'Fig_F10W_8_zero_shot_recovery','Official zero-shot motion recovery; finite estimates only, failures retained')
    fig,ax=plt.subplots(2,2,figsize=(12,8));mae=np.array([analysis[m]['finite_global_mae'] for m in METHODS])/LIMITS
    im=ax[0,0].imshow(mae,aspect='auto',cmap='magma_r',vmin=0,vmax=3);ax[0,0].set(yticks=range(7),yticklabels=[f'{s} (n={analysis[m]["finite_global_mae_n"]})' for s,m in zip(SHORT,METHODS)],xticks=range(6),xticklabels=['vx','vy','div','vort','s1','s2'],title='Conditional MAE / frozen tolerance');fig.colorbar(im,ax=ax[0,0],label='Error / tolerance')
    for j,m in enumerate(primary):ax[0,1].plot(range(6),common_mae[m]/LIMITS,'o-',label=m)
    ax[0,1].set(xticks=range(6),xticklabels=['vx','vy','div','vort','s1','s2'],ylabel='Common-support MAE / tolerance',title=f'Same {len(common)} global examples');ax[0,1].legend(fontsize=8)
    x=np.arange(7);ax[1,0].bar(x-.16,[analysis[m]['base']['accepted'] for m in METHODS],.32,label='Geometry gate');ax[1,0].bar(x+.16,[analysis[m]['gate']['accepted'] for m in METHODS],.32,label='+ validation risk gate');ax[1,0].legend(fontsize=8);ax[1,0].set(ylabel='Accepted / 138',xticks=x,xticklabels=SHORT);ax[1,0].tick_params(axis='x',rotation=35)
    ax[1,1].bar(x,[analysis[m]['gate']['bad_accepted'] for m in METHODS],color=COLORS);ax[1,1].set(ylabel='Wrong accepted cases',xticks=x,xticklabels=SHORT,title='Zero acceptance is not success');ax[1,1].tick_params(axis='x',rotation=35)
    save(fig,'Fig_F10W_9_six_parameter_comparison','Frozen final-test comparison; synthetic known motion, one recording only')
    fig,ax=plt.subplots(1,2,figsize=(13,6))
    for k,m in enumerate(METHODS):ax[0].plot(curves[m]['coverage'],curves[m]['risk'],label=SHORT[k],color=COLORS[k])
    ax[0].axhline(.05,color='black',ls='--');ax[0].legend(fontsize=8);ax[0].set(xlabel='Accepted fraction of all 138 cases',ylabel='Joint error among accepted',ylim=(-.02,1.02),title='Descriptive final risk–coverage curve')
    cases=['fixed_background','dense','dropout','boundary','symmetric','non_affine','defocus','occluded'];mat=np.array([[conditions[m][c]['base']['bad_accepted'] for c in cases] for m in METHODS]);im=ax[1].imshow(mat,aspect='auto',cmap='OrRd',vmin=0,vmax=6);ax[1].set(yticks=range(7),yticklabels=SHORT,xticks=range(len(cases)),xticklabels=[c.replace('_',' ') for c in cases],title='Wrong geometry-accepted cases / 6 per condition');ax[1].tick_params(axis='x',rotation=55);fig.colorbar(im,ax=ax[1],label='Wrong accepted count')
    save(fig,'Fig_F10W_10_risk_and_failures','Risk, coverage and difficult-case false acceptance')
    replay={}
    for m in ['PIV','LRAFT_ST','CoTracker3','POT_unbalanced']:
        p=OUT/'results'/f'causal_replay_{m}.json'
        if p.exists():r=json.loads(p.read_text());replay[m]={k:v for k,v in r.items() if k not in ['rows','producer','consumer']}
    if replay:
        fig,ax=plt.subplots(2,2,figsize=(11,8));names=list(replay);x=np.arange(len(names));a=[replay[m]['full_pipeline_ms']['p50'] or 0 for m in names];b=[replay[m]['full_pipeline_ms']['p95'] or 0 for m in names];ax[0,0].bar(x-.17,a,.34,label='p50');ax[0,0].bar(x+.17,b,.34,label='p95');ax[0,0].axhline(33.3,ls='--',color='#b86b60');ax[0,0].set(ylabel='Complete ready pipeline (ms)');ax[0,0].legend();ax[0,1].bar(x,[replay[m]['fresh_six_hz'] for m in names]);ax[0,1].set(ylabel='Fresh accepted six-parameter Hz');ax[1,0].bar(x,[replay[m]['n_dropped']/937 for m in names]);ax[1,0].set(ylabel='Dropped input-frame fraction');ax[1,1].bar(x,[replay[m]['first_ready_pts_s'] or 0 for m in names]);ax[1,1].set(ylabel='First full-window source PTS (s)')
        for a in ax.flat:a.set(xticks=x,xticklabels=names);a.tick_params(axis='x',rotation=20)
        save(fig,'Fig_F10W_11_paced_runtime','Full-video paced replay; latest slot = 1, recording is not a live camera')
    dump(OUT/'results/model_analysis.json',dict(identity=identity(),methods=analysis,risk_curves=curves,conditions=conditions,common_support=dict(methods=primary,n=len(common),indices=common,mae=common_mae),replay=replay,source_scope='SINGLE_RECORDING_GENERALIZATION_NOT_VERIFIED',training_ready=False))
if __name__=='__main__':main()
