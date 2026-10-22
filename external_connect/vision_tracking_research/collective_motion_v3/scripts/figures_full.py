"""Full-file PIV figures from saved records. No new inference or smoothing."""
import sys,json,collections
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from full_records import full_records
import pyarrow.parquet as pq
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import cv2
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'font.family':'DejaVu Sans','pdf.fonttype':42,'figure.dpi':150})
COLORS=['#5064a7','#c07b38','#5a9684','#9663a5','#a0a0a0']
LABELS=['Translation x (px/s)','Translation y (px/s)','Divergence (1/s)','Vorticity (1/s)','Shear 1 (1/s)','Shear 2 (1/s)']
def save(fig,name,caption):
    fig.suptitle(caption,fontsize=13);fig.tight_layout(rect=[0,0,1,.96]);(OUT/'figures').mkdir(exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(OUT/'figures'/f'{name}.{ext}',bbox_inches='tight',dpi=200)
    plt.close(fig)
def main():
    p=pq.read_table(OUT/'results/full_video_piv_per_frame.parquet').to_pylist();dist=pq.read_table(OUT/'results/full_video_distribution.parquet').to_pylist();main=[r for r in p if r['method']=='candidate_3'];t=np.array([r['pts_s'] for r in main]);valid=np.array([r['valid_all'] for r in main]);counts=np.array([r['visible_count'] for r in dist]);centers=np.array([r['group_center_xy'] for r in dist]);density=np.array([r['density_summary']['mean_per_px2'] for r in dist]);boundary=np.array([r['ingress_egress_indicators']['boundary_candidate_fraction'] for r in dist]);values=np.array([[r['values'][n] if r['values'][n] is not None else np.nan for n in NAMES] for r in main]);reasons=collections.Counter(r['reason'] or 'VALID' for r in main)
    fig,ax=plt.subplots(2,2,figsize=(11,6));ax[0,0].plot(t,np.arange(len(t)),color=COLORS[0]);ax[0,0].set(xlabel='Source PTS (s)',ylabel='Decoded / analyzed frame ID',title='937 / 937 analyzed; no missing interval');ax[0,1].plot(t[1:],np.diff(t)*1000);ax[0,1].set(xlabel='Source PTS (s)',ylabel='Frame interval (ms)',ylim=(30,36),title='Monotonic PTS; no duplicate pictures');ax[1,0].scatter(t,np.ones(len(t)),s=3,c=np.where(valid,COLORS[0],'#d57c68'));ax[1,0].set(xlabel='Source PTS (s)',yticks=[],title='Blue: quality accepted; coral: rejected / warmup');short={'VALID':'Accepted','WARMUP':'Warmup','NON_AFFINE_OR_INCONSISTENT_VECTORS':'Inconsistent vectors','INSUFFICIENT_OR_DEGENERATE_SUPPORT':'Insufficient support'};labels=[short[k] for k in reasons];ax[1,1].barh(labels,list(reasons.values()),color=COLORS[0]);ax[1,1].set(xlabel='Frames',title='Acceptance is not ground-truth accuracy');save(fig,'Fig_F10W_1_full_video_coverage','Full source coverage and PIV observation status')
    fig,ax=plt.subplots(2,2,figsize=(11,7));ax[0,0].plot(t,counts,color=COLORS[0]);ax[0,0].set(ylabel='Visible candidates',xlabel='Source PTS (s)');ax[0,1].plot(t,centers[:,0],label='x');ax[0,1].plot(t,centers[:,1],label='y');ax[0,1].legend();ax[0,1].set(ylabel='Candidate center in ROI (px)',xlabel='Source PTS (s)');ax[1,0].plot(t,density*1e4);ax[1,0].set(xlabel='Source PTS (s)',ylabel='Candidates / 10,000 px²');sc=ax[1,1].scatter(centers[:,0],centers[:,1],c=t,s=5,cmap='viridis');ax[1,1].set(xlabel='Center x in ROI (px)',ylabel='Center y in ROI (px)',title='Distribution center; no material identity');ax[1,1].invert_yaxis();fig.colorbar(sc,ax=ax[1,1],label='Source PTS (s)');save(fig,'Fig_F10W_2_full_distribution','Full-duration visible candidate distribution')
    fig,ax=plt.subplots(3,2,figsize=(12,9),sharex=True)
    for j,a in enumerate(ax.flat):
        a.plot(t,values[:,j],lw=.7,color=COLORS[0]);a.set(ylabel=LABELS[j],xlabel='Source PTS (s)');a.scatter(t[~valid],np.full((~valid).sum(),.02),transform=a.get_xaxis_transform(),marker='|',color='#d57c68',s=20)
    save(fig,'Fig_F10W_3_full_six_parameters','Frozen PIV estimates: gaps remain null; real accuracy NOT_VERIFIED')
    methods=['candidate_3','candidate_5','candidate_10','raw_3','local_3'];comparison={};fig,ax=plt.subplots(2,2,figsize=(12,7))
    for k,m in enumerate(methods):
        rr=[r for r in p if r['method']==m];v=np.array([r['values']['translation_vx'] if r['values']['translation_vx'] is not None else np.nan for r in rr]);ax[0,0].plot(t,v,lw=.6,label=m,color=COLORS[k]);comparison[m]=dict(valid=sum(r['valid_all'] for r in rr),valid_fraction=np.mean([r['valid_all'] for r in rr]),residual=stats([r['residual'] if r['residual'] is not None else np.nan for r in rr]),algorithm_ms=stats([r['algorithm_ms'] for r in rr]),paired_vx_difference_from_main=stats(np.abs(v-values[:,0])))
    ax[0,0].legend(fontsize=8,ncol=2);ax[0,0].set(xlabel='Source PTS (s)',ylabel='Translation x (px/s)');ax[0,1].bar(methods,[comparison[m]['valid_fraction'] for m in methods],color=COLORS);ax[0,1].set(ylabel='Quality-accepted fraction',ylim=(0,1));ax[0,1].tick_params(axis='x',rotation=20)
    ax[1,0].plot(t,[r['valid_vectors'] for r in main]);ax[1,0].set(xlabel='Source PTS (s)',ylabel='Accepted PIV vectors');ax[1,1].plot(t,[r['coverage'] for r in main]);ax[1,1].axhline(.18,color='#a45b55',ls='--',label='Frozen support limit');ax[1,1].legend();ax[1,1].set(xlabel='Source PTS (s)',ylabel='Convex-hull coverage / ROI area');save(fig,'Fig_F10W_4_piv_comparators','Window / background comparison and spatial support (no truth ranking)')
    failures=[]
    for reason in ['NON_AFFINE_OR_INCONSISTENT_VECTORS','INSUFFICIENT_OR_DEGENERATE_SUPPORT']:
        ids=[r['frame_id'] for r in main if r['reason']==reason];groups=[]
        for _,g in itertools_group(ids):groups.append(g)
        failures.append(dict(reason=reason,frame_ids=ids,intervals=[dict(start_frame=g[0],end_frame=g[-1],first_pts_s=t[g[0]],last_pts_s=t[g[-1]]) for g in groups]))
    ids=[468]+[f['frame_ids'][0] for f in failures];source=old.read_frames(ids,VIDEO);records={r['frame_id']:r for r in full_records() if r['frame_id'] in ids};fig,ax=plt.subplots(1,3,figsize=(12,7))
    for a,fid in zip(ax,ids):
        bgr,_=source[fid];r=records[fid];x,y,w,h=r['roi'];a.imshow(cv2.cvtColor(bgr[y:y+h,x:x+w],cv2.COLOR_BGR2GRAY),cmap='gray',vmin=0,vmax=255);vv=np.array(r['methods']['candidate_3']['vectors'],float);good=vv[:,4]>0
        a.quiver(vv[good,0],vv[good,1],vv[good,2]*8,vv[good,3]*8,angles='xy',scale_units='xy',scale=1,color='#51ce9c',width=.005);a.scatter(vv[~good,0],vv[~good,1],s=5,marker='x',color='#ce8170');a.set(title=f'Frame {fid}, {t[fid]:.3f} s\n'+short[main[fid]['reason'] or 'VALID'],xlabel='ROI x (px)',ylabel='ROI y (px)')
    save(fig,'Fig_F10W_5_observed_failures','Representative saved local vectors; arrows ×8, rejected vectors marked')
    regions=pq.read_table(OUT/'results/full_video_regions.parquet').to_pylist();nr=9;nc=12;cs=np.zeros((nr,nc));supported=np.zeros_like(cs);nn=np.zeros_like(cs)
    for r in regions:
        i,j=r['y']//128,r['x']//128;cs[i,j]+=r['candidate_count'];nn[i,j]+=1;supported[i,j]+=r['status']=='CANDIDATE_TEMPORAL_SUPPORT_UNCALIBRATED'
    image=old.read_frames([0],VIDEO)[0][0];fig,ax=plt.subplots(1,3,figsize=(14,5));ax[0].imshow(cv2.cvtColor(image,cv2.COLOR_BGR2RGB));ax[0].add_patch(Rectangle((650,300),256,384,fill=False,color='#3c9bc3',lw=2));ax[0].axvline(450,color='#c0a96d',ls='--');ax[0].axvline(1330,color='#c0a96d',ls='--');ax[0].set(title='Original frame 0; ROI + visual channel guides');a=ax[1].pcolormesh(np.r_[np.arange(0,1476,128),1476],np.r_[np.arange(0,1120,128),1120],cs/nn,cmap='Purples');ax[1].invert_yaxis();fig.colorbar(a,ax=ax[1],label='Mean detected candidates / tile');ax[1].set(title='Fixed texture may also generate candidates');a=ax[2].pcolormesh(np.r_[np.arange(0,1476,128),1476],np.r_[np.arange(0,1120,128),1120],supported/nn,vmin=0,vmax=1,cmap='viridis');ax[2].invert_yaxis();fig.colorbar(a,ax=ax[2],label='Fraction with candidate + temporal evidence');ax[2].set(title='Heuristic evidence, not validated velocity');
    for a in ax:a.set(xlabel='Full-frame x (px)',ylabel='Full-frame y (px)')
    save(fig,'Fig_F10W_6_fullframe_scope','Full-frame candidate scan; channel guides are uncalibrated image annotations')
    fig,ax=plt.subplots(2,2,figsize=(11,7));ax[0,0].plot(t,boundary);ax[0,0].set(xlabel='Source PTS (s)',ylabel='Boundary-candidate fraction');ax[0,1].scatter(counts,[r['valid_vectors'] for r in main],s=7,c=np.where(valid,COLORS[0],'#d57c68'));ax[0,1].set(xlabel='Visible candidate count',ylabel='Accepted PIV vectors');ax[1,0].plot(t,[r['main_pipeline_component_ms'] for r in main]);ax[1,0].set(xlabel='Source PTS (s)',ylabel='Main offline components (ms)');ax[1,1].scatter(boundary,[r['residual'] for r in main],s=7,c=np.where(valid,COLORS[0],'#d57c68'));ax[1,1].axhline(1.5,ls='--',color='#a45b55');ax[1,1].set(xlabel='Boundary-candidate fraction',ylabel='Affine residual (px)');save(fig,'Fig_F10W_7_boundary_support_runtime','Boundary / support diagnostics; count change is not material flux')
    summary=dict(identity=identity(),n=len(main),valid=int(valid.sum()),failures=failures,count_range=[int(counts.min()),int(counts.max())],count_first_last=[int(counts[0]),int(counts[-1])],center_roi_first_last=[centers[0],centers[-1]],center_range=np.c_[centers.min(0),centers.max(0)],density_range=[density.min(),density.max()],parameters={n:stats(values[:,j]) for j,n in enumerate(NAMES)},comparisons=comparison,selected_failure_frames=ids,fullframe_scope='candidate/temporal heuristic only; no calibrated channel or material velocity mask',true_accuracy='NOT_VERIFIED')
    dump(OUT/'results/full_video_analysis_summary.json',summary)
def itertools_group(ids):
    import itertools
    for k,g in itertools.groupby(enumerate(ids),lambda a:a[1]-a[0]):yield k,[x[1] for x in g]
if __name__=='__main__':main()
