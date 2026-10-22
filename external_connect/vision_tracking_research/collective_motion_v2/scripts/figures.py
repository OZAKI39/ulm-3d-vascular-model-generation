"""Publication figures from frozen predictions only; no synthetic score fillers."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,cv2
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from vcommon import ROOT,OUT,DATA,P,NAMES,LIMITS,dump,sha,identity
from heatmaps import centers_to_heatmap

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#222222','xtick.color':'#222222','ytick.color':'#222222','pdf.fonttype':42,'savefig.facecolor':'white'})
COLORS=['#344c68','#7a54a0','#2986ac','#d48a30','#3a9174','#c45e66'];SHORT=['vx','vy','divergence','vorticity','shear 1','shear 2'];FILES=[]
def save(fig,name,inputs):
    folder=OUT/'figures';folder.mkdir(exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(folder/f'{name}.{ext}',dpi=260,bbox_inches='tight')
    plt.close(fig);FILES.append(dict(name=name,inputs=inputs,outputs={e:sha(folder/f'{name}.{e}') for e in ['png','pdf']}))
def image_axes(ax,img,title):ax.imshow(img,cmap='gray',vmin=0,vmax=255);ax.set_title(title,fontsize=10);ax.set_xticks([]);ax.set_yticks([])
def pretty(name):
    if name=='PIV_SINGLE':return 'PIV single pass'
    for k,v in [('direct_fused_robust','Direct6 + heat (robust)'),('direct_gray','Direct6 gray'),('direct_fused','Direct6 + heat'),('direct_heat','Direct6 heat only'),('geo_fused','GeoMotion')]:
        if name.startswith(k):return v
    return name

def fig1():
    old=ROOT.parent/'smart_online_tracking_v2';folder=old/'results/heldout_49_50_20261016_v1';freeze=old/'results/expert_freezes/vision2_heldout_49_50_20261016_v1'
    rows=json.loads((folder/'all_rows.json').read_text());r=next(r for r in rows if r['model']=='opencv' and r['domain']=='ROI');fr=json.loads((freeze/'frozen.json').read_text());source=next(x for x in fr['rows'] if x['record']['id']==r['image_id'])
    full=cv2.imread(str(freeze/source['frozen_image_file']));y0,x0,y1,x1=r['box_yxyx'];gray=cv2.cvtColor(full[y0:y1,x0:x1],cv2.COLOR_BGR2GRAY);pred=np.load(folder/'predictions'/(r['artifact']+'.npz'));xy=np.array([o['xy_local'] for o in r['objects']]);heat=centers_to_heatmap(xy,gray.shape)
    fig,ax=plt.subplots(2,3,figsize=(12,8),layout='constrained')
    image_axes(ax[0,0],gray,'(a) Real Chronos reviewed ROI')
    image_axes(ax[0,1],gray,'(b) Frozen OpenCV candidate centers');ax[0,1].scatter(*xy.T,s=15,facecolors='none',edgecolors='#2b9bb9',linewidth=.9)
    ax[0,2].imshow(heat,cmap='magma',vmin=0,vmax=1);ax[0,2].set(title='(c) Center heatmap (not ground truth)',xticks=[],yticks=[])
    image_axes(ax[1,0],gray,'(d) Human review: successes and errors')
    good={m['pred'] for m in r['metrics']['matches']};fp=set(r['metrics']['false_positive_ids'])
    for o in r['objects']:
        if o['label'] in good:ax[1,0].plot(*o['xy_local'],'o',ms=4,mfc='none',color='#1a9aa3')
        if o['label'] in fp:ax[1,0].plot(*o['xy_local'],'x',ms=5,color='#c24538')
    for label in r['metrics']['false_negative_ids']:
        yy,xx=np.where(pred['truth']==label);ax[1,0].plot(xx.mean(),yy.mean(),'+',ms=7,color='#d69a21')
    ax[1,0].legend(handles=[Line2D([],[],marker='o',ls='',mfc='none',color='#1a9aa3',label='Matched'),Line2D([],[],marker='x',ls='',color='#c24538',label='False positive'),Line2D([],[],marker='+',ls='',color='#d69a21',label='Missed')],loc='upper right',fontsize=8)
    image_axes(ax[1,1],gray,'(e) Gray + actual detector heat');ax[1,1].imshow(heat,cmap='magma',vmin=0,vmax=1,alpha=np.minimum(.7,heat*.8))
    image_axes(ax[1,2],cv2.cvtColor(full,cv2.COLOR_BGR2GRAY),'(f) Full image: channel and fixed background');ax[1,2].add_patch(plt.Rectangle((x0,y0),x1-x0,y1-y0,fill=False,color='#2b9bb9',lw=1.5))
    fig.suptitle('Real image evidence | historical frozen detector; partial human review only\nNo new detection accuracy claim; no unreviewed pixels treated as negatives',fontsize=13)
    save(fig,'Fig_F10V_1_detection_heatmaps',[str(folder/'all_rows.json'),str(freeze/source['frozen_image_file'])])

def main():
    assert (OUT/'results/final_evaluation_complete.json').exists()
    analysis=json.loads((OUT/'results/analysis.json').read_text());methods=analysis['methods'];ch=analysis['champions'];keys=['PIV_SINGLE']+[ch[n] for n in ['direct_gray','direct_fused','direct_heat','direct_fused_robust','geo_fused']]
    raw=np.load(DATA/'final_test/raw.npy',mmap_mode='r');truth=np.load(DATA/'final_test/truth.npy');meta=[json.loads(s) for s in (DATA/'final_test/index.jsonl').read_text().splitlines()]
    fig1()
    cases=['translation_x','rotation_cw','contraction','expansion','shear','fixed_background'];fig,ax=plt.subplots(2,6,figsize=(15,6.8),layout='constrained')
    for col,case in enumerate(cases):
        i=next(i for i,m in enumerate(meta) if m['case']==case)
        for row,k in enumerate([0,2]):image_axes(ax[row,col],raw[i,k],f'{case.replace("_"," ")}\nt={meta[i]["times_s"][k]*1000:.1f} ms')
        ax[1,col].set_xlabel(f'vx={truth[i,0]:.1f}, vy={truth[i,1]:.1f} px/s\ndiv={truth[i,2]:.2f}, vort={truth[i,3]:.2f}\ns1={truth[i,4]:.2f}, s2={truth[i,5]:.2f} (1/s)',fontsize=9)
    fig.suptitle('Real-patch synthetic motion | one process per triple; interval-equivalent truth\nFeathered patches and stationary source dots remain visible domain limitations',fontsize=13)
    save(fig,'Fig_F10V_2_synthetic_truth',['data/final_test/raw.npy','data/final_test/index.jsonl','data/final_test/truth.npy'])
    fig,ax=plt.subplots(2,2,figsize=(13,9),layout='constrained');arr=np.array([methods[k]['mae_all_finite_global'] for k in keys])/LIMITS
    im=ax[0,0].imshow(arr,cmap='magma',vmin=0,vmax=max(1.5,np.nanmax(arr)));ax[0,0].set(yticks=range(len(keys)),yticklabels=[pretty(k) for k in keys],xticks=range(6),xticklabels=SHORT,title='(a) MAE / unchanged error tolerance')
    ax[0,0].tick_params(axis='x',rotation=30);fig.colorbar(im,ax=ax[0,0],label='Normalized MAE',shrink=.8)
    for i in range(len(keys)):
        for j in range(6):ax[0,0].text(j,i,f'{arr[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if arr[i,j]<.6*np.nanmax(arr) else 'black')
    sign=np.array([[x['accuracy'] for x in methods[k]['sign_accuracy']] for k in keys],float)
    ims=ax[0,1].imshow(sign,cmap='viridis',vmin=0,vmax=1);fig.colorbar(ims,ax=ax[0,1],label='Sign accuracy',shrink=.8);ax[0,1].set(yticks=range(len(keys)),yticklabels=[pretty(k) for k in keys],xticks=range(6),xticklabels=SHORT,title='(b) Sign correctness (non-negligible truth)');ax[0,1].tick_params(axis='x',rotation=30)
    for ai,name in zip(ax[1],[ch['direct_fused'],ch['geo_fused']]):
        p=np.load(OUT/'results'/f'{name}_final_test.npz')['mean'];mask=np.array([m['class_id']==0 for m in meta]);ai.scatter(truth[mask,0],p[mask,0],s=5,alpha=.25,c='#4b6f9b');ai.plot([-50,50],[-50,50],'k--',lw=1);ai.set(xlabel='True vx (px/s)',ylabel='Estimated vx (px/s)',title=pretty(name)+' | all global cases')
    fig.suptitle(f'Frozen test: N={len(meta)} triples; one recording\nMAE includes rejected finite predictions; unavailable PIV fits have separate denominators',fontsize=13)
    save(fig,'Fig_F10V_3_six_parameter_recovery',['results/analysis.json'])
    fig,ax=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    for k,c in zip(keys,COLORS):
        r=methods[k]['risk_coverage'];ax[0].plot(r['coverage'],r['risk'],label=pretty(k),color=c,lw=1.6)
    ax[0].axhline(.05,color='k',ls='--',lw=.8);ax[0].set(xlabel='Accepted fraction of all test triples',ylabel='Any-parameter / class failure fraction',ylim=(0,1),xlim=(0,1),title='(a) Risk–coverage (same test denominator)');ax[0].legend(fontsize=7)
    xx=np.arange(len(keys));ax[1].bar(xx,[methods[k]['acceptance_rate'] for k in keys],color=COLORS);ax[1].set(xticks=xx,xticklabels=[pretty(k) for k in keys],ylim=(0,1),ylabel='Accepted fraction',title='(b) PIV quality / AI validation risk gate');ax[1].tick_params(axis='x',rotation=65)
    for i,k in enumerate(keys):ax[1].text(i,methods[k]['acceptance_rate']+.015,str(methods[k]['n_accepted']),ha='center',fontsize=8)
    modes=['','_drop','_false','_jitter','_severe']
    for base,c in [('direct_fused','#2986ac'),('direct_fused_robust','#3a9174')]:
        values=[np.mean(np.array(methods[ch[base]+m]['mae_all_finite_global'])/LIMITS) for m in modes];ax[2].plot(range(5),values,'o-',color=c,label=pretty(ch[base]))
    ax[2].set(xticks=range(5),xticklabels=['Original','50% drop','False points','Jitter','Severe'],ylabel='Mean normalized MAE',title='(c) Errors propagated from detector heat');ax[2].tick_params(axis='x',rotation=40);ax[2].legend(fontsize=8)
    save(fig,'Fig_F10V_4_heatmap_risk_coverage',['results/analysis.json'])
    fig,ax=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    for k,c in zip([ch['direct_gray'],ch['direct_fused'],ch['geo_fused']],COLORS[1:]):
        for level,style in [('0.9','--'),('0.95','-')]:
            interval=methods[k]['intervals'][level];ax[0].plot(range(6),interval['coverage'],style+'o',ms=3,c=c,label=pretty(k)+' '+level)
        ax[1].plot(range(6),np.array(methods[k]['intervals']['0.95']['mean_width'])/LIMITS,'o-',c=c,label=pretty(k))
    ax[0].axhline(.9,c='k',ls=':',lw=.7);ax[0].axhline(.95,c='k',ls=':',lw=.7);ax[0].set(xticks=range(6),xticklabels=SHORT,ylim=(0,1.03),ylabel='Empirical marginal coverage',title='(a) Independent synthetic calibration');ax[0].tick_params(axis='x',rotation=45);ax[0].legend(fontsize=7,loc='lower left')
    ax[1].set(xticks=range(6),xticklabels=SHORT,ylabel='Full 95% width / error tolerance',yscale='log',title='(b) Coverage may require wide intervals');ax[1].tick_params(axis='x',rotation=45);ax[1].legend(fontsize=8)
    conditions=['translation_x','rotation_cw','fixed_background','defocus','dropout','symmetric'];vals=np.array([methods[ch['geo_fused']]['intervals']['0.95']['by_case'][c]['coverage'] for c in conditions]);im=ax[2].imshow(vals,cmap='viridis',vmin=0,vmax=1);ax[2].set(yticks=range(6),yticklabels=conditions,xticks=range(6),xticklabels=SHORT,title='(c) GeoMotion 95% coverage by condition');ax[2].tick_params(axis='x',rotation=45);fig.colorbar(im,ax=ax[2],shrink=.6)
    fig.suptitle('Calibration groups share one recording: no real-domain or independent-experiment guarantee',fontsize=12)
    save(fig,'Fig_F10V_5_uncertainty',['results/analysis.json'])
    runtime={n:json.loads((OUT/'results'/f'real_replay_{n}.json').read_text()) for n in ['PIV_SINGLE','direct_gray','direct_fused','geo_fused']};fig,ax=plt.subplots(2,4,figsize=(18,8),layout='constrained');ns=list(runtime);labels=[pretty(n) for n in ns];xx=np.arange(4)
    for d,key in [(-.18,'p50'),(.18,'p95')]:ax[0,0].bar(xx+d,[runtime[n]['full_pipeline_ms'][key] for n in ns],.36,label=key)
    ax[0,0].set(ylabel='ms',title='(a) Full pipeline; ready causal windows');ax[0,0].legend()
    ax[0,1].bar(xx,[runtime[n]['fresh_six_hz'] for n in ns],color=COLORS[:4]);ax[0,1].set(ylabel='Fresh accepted six-parameter Hz',title='(b) Validity and age both required')
    ax[0,2].bar(xx,[runtime[n]['drop_fraction']*100 for n in ns],color=COLORS[:4]);ax[0,2].set(ylabel='Dropped frames (%)',title='(c) Single latest-frame slot')
    ax[0,3].bar(xx,[runtime[n]['processed_hz'] for n in ns],color=COLORS[:4]);ax[0,3].set(ylabel='Processed frames / active second',title='(d) Actual causal replay throughput')
    ax[1,0].bar(xx,[runtime[n]['parameter_count']/1000 for n in ns],color=COLORS[:4]);ax[1,0].set(ylabel='Learned parameters (thousand)',title='(e) Model size')
    ax[1,1].bar(xx,[runtime[n]['peak_cuda_allocated_mb'] for n in ns],color=COLORS[:4]);ax[1,1].set(ylabel='Peak allocated CUDA MiB',title='(f) Actual replay memory')
    ax[1,2].bar(xx,[runtime[n]['resource']['cpu_core_equivalents'] for n in ns],color=COLORS[:4]);ax[1,2].set(ylabel='Parent CPU core equivalents',title='(g) Software process CPU time')
    ax[1,3].bar(xx,[np.mean([s.get('gpu_util_percent',0) for s in runtime[n]['resource']['samples']]) for n in ns],color=COLORS[:4]);ax[1,3].set(ylabel='Device GPU utilization (%)',title='(h) 1 s samples during paced replay')
    for a in ax.flat:a.set_xticks(xx,labels,rotation=25,ha='right')
    fig.suptitle('Same RTX 4090 server and same real ROI | actual PTS; no camera end-to-end measurement',fontsize=12)
    save(fig,'Fig_F10V_6_runtime_resources',['results/real_replay_'+n+'.json' for n in ns])
    fig,ax=plt.subplots(2,4,figsize=(14,7),layout='constrained');name=ch['geo_fused'];pr=np.load(OUT/'results'/f'{name}_final_test.npz');conditions=['fixed_background','defocus','dropout','symmetric','non_affine','start','low_texture','collinear']
    for a,case in zip(ax.flat,conditions):
        ids=[i for i,m in enumerate(meta) if m['case']==case];i=min(ids,key=lambda j:pr['score'][j]);image_axes(a,raw[i,2],case.replace('_',' '));err=np.max(np.abs(pr['mean'][i]-truth[i])/LIMITS);status=('Non-affine motion' if meta[i]['class_id']==1 else 'Not observable; known generator') if meta[i]['class_id'] else f'max error / tol = {err:.1f}'
        a.set_xlabel(status+f'\nscore={pr["score"][i]:.1f} (lower = more confident)',fontsize=8)
    fig.suptitle('Failure audit | lowest-score example per preregistered difficult condition\nRaw estimates are diagnostics; rejected values are not valid observations',fontsize=12)
    save(fig,'Fig_F10V_7_difficult_cases',['results/analysis.json',f'results/{name}_final_test.npz'])
    dump(OUT/'figures/manifest.json',dict(identity=identity(),figures=FILES))
if __name__=='__main__':main()
