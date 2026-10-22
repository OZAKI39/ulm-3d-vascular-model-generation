"""Publication figures from saved computed results only; English, PNG and PDF."""
import sys,json,time,textwrap
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from channel import Channel
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Rectangle,Ellipse
from matplotlib.lines import Line2D
import cv2
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
COL=['#697687','#416ca0','#d38337','#5e9d89','#805b9f'];LABEL=['Old OpenPIV','SciPy','Rust','TorchPIV','CuPy'];c=Channel();manifest=[]
def read(p):return json.loads((ROOT/p).read_text())
def save(fig,name,note):
 fig.text(.01,.012,textwrap.fill(note,170),fontsize=8,color='#40464c');fig.subplots_adjust(bottom=.22 if fig.get_size_inches()[1]<5.5 else .16)
 for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'{name}.{ext}',dpi=210,bbox_inches='tight')
 manifest.append(dict(name=name,note=note,files={ext:sha(ROOT/'figures'/f'{name}.{ext}') for ext in ['png','pdf']}));plt.close(fig)
def scene(ax,fid=450):
 p=ROOT/'evidence'/f'source_{fid:04d}.png';g=cv2.imread(str(p),0);ax.imshow(g,cmap='gray',vmin=0,vmax=255);ax.set(xlim=(0,1476),ylim=(1120,0),xlabel='Image x (px)',ylabel='Image y (px)');return g
def boundary(ax):
 ax.add_patch(Polygon(c.config['channel_polygon'],fill=False,edgecolor='#d89522',linewidth=1.5))
 for p in c.config['excluded_polygons']:ax.add_patch(Polygon(p,facecolor='#bb4a4a',alpha=.22,edgecolor='#ac3030'))
def vectors(ax,s,stride=1):
 w=s['windows'];good=np.array(w['valid']);x=np.array(w['x']);y=np.array(w['y']);dt=s['time_window']['dt_s'] or 0;u=np.array(w['vx'],float)*dt;v=np.array(w['vy'],float)*dt;mask=good&((np.arange(len(x))%stride)==0)
 ax.scatter(x[~good&np.array(w['safe'])],y[~good&np.array(w['safe'])],s=1,c='#9c7d6d',alpha=.4)
 q=ax.quiver(x[mask],y[mask],6*u[mask],6*v[mask],np.hypot(u[mask],v[mask])/(dt or 1),cmap='viridis',angles='xy',scale_units='xy',scale=1,width=.0022);return q
def main():
 start=time.perf_counter();a=read('results/delivery_analysis.json');samples=read('results/figure_sample_states.json');audit=read('results/channel_audit.json');diag=read('results/diagnostics_v1/post_freeze_diagnostics.json');s=samples['450']
 fig,axs=plt.subplots(1,2,figsize=(12,5.9));scene(axs[0]);boundary(axs[0]);axs[0].add_patch(Rectangle((650,300),256,384,fill=False,edgecolor='#208fc0',lw=2));axs[0].set_title('(a) Full source view and historical ROI')
 axs[0].legend(handles=[Line2D([],[],color='#d89522',label='Editable channel boundary'),Line2D([],[],color='#bb4a4a',label='Excluded wall / fixed structure'),Line2D([],[],color='#208fc0',label='Historical 256 x 384 px ROI')],loc='lower left',fontsize=8)
 axs[1].imshow(c.valid,cmap='Greys',vmin=0,vmax=1);axs[1].scatter((c.x+c.crop[0])[c.safe],(c.y+c.crop[1])[c.safe],s=1,c='#3b92a4');axs[1].set(title='(b) Whole-window support mask',xlabel='Image x (px)',ylabel='Image y (px)');axs[1].text(.04,.04,'Observation mask: 805,495 px\n1,283 whole 48-px windows\n8.19 x historical ROI area',transform=axs[1].transAxes,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
 save(fig,'Fig_F10X_1_channel_geometry_mask','Image-defined geometry, not calibrated boundary truth. Windows require 100% valid pixels; central texture remains a background-risk region.')
 fig,axs=plt.subplots(1,3,figsize=(15,5.8));scene(axs[0]);boundary(axs[0]);q=np.array(s['candidates_global_xy']);axs[0].scatter(q[:,0],q[:,1],s=5,facecolors='none',edgecolors='#008bb4',linewidths=.6);axs[0].set_title('(a) Full-channel candidates, frame 450')
 w=s['windows'];density=np.array(w['candidate_count']).reshape(s['window_shape']);density[~c.safe]=np.nan;im=axs[1].imshow(density,cmap='magma',origin='upper');fig.colorbar(im,ax=axs[1],label='Candidates per 48-px window',shrink=.65);axs[1].set(title='(b) Candidate support, all safe windows',xlabel='Grid column',ylabel='Grid row')
 rr=audit['reviewed_rows'];n=sum(r['n_human_bubbles'] for r in rr);det=sum(sum(r['detected_within_8px']) for r in rr);covered=sum(sum(r['covered_by_selected_window']) for r in rr);axs[2].bar(['Detected\nwithin 8 px','Covered by\nactive window'],[det/n,covered/n],color=['#369ab4','#5e9d89']);axs[2].set(ylim=(0,1.08),ylabel='Fraction of reviewed bubble centers',title='(c) Human-label diagnostic');axs[2].text(0,.65,f'{det}/{n}',ha='center');axs[2].text(1,.65,f'{covered}/{n}',ha='center');axs[2].text(.02,.08,'Six partially reviewed frames\nNot full-channel ground truth\n2 confirmed centers missed',transform=axs[2].transAxes,fontsize=9)
 save(fig,'Fig_F10X_2_candidate_activity','Candidates are detections, not bubble identities. Activity selector is tested against preserved human labels; unreviewed pixels are not negative truth.')
 fig,axs=plt.subplots(1,3,figsize=(15,5));bt=a['backend'];regs=['historical','medium','full_visible_bbox']
 for j,reg in enumerate(regs):
  rows=[next(r for r in bt if r['backend']==f'B{i}' and r['region']==reg) for i in range(5)];x=np.arange(5);axs[j].bar(x,[r['timing']['p50'] for r in rows],color=COL);axs[j].scatter(x,[r['timing']['p95'] for r in rows],marker='_',s=180,c='black',label='p95');axs[j].set(xticks=x,xticklabels=LABEL,ylabel='Measured pair pipeline (ms)',title=f'({chr(97+j)}) {reg.replace("_"," ")}\n{rows[0]["windows"]} interrogation windows');axs[j].tick_params(axis='x',rotation=35);axs[j].legend()
 save(fig,'Fig_F10X_3_cpu_gpu_benchmark','Same source pairs and 48/24 windows; warm CUDA-synchronized runs. Includes preprocessing/PIV/validation/one six-fit, excludes decode/full-channel scan/multiple groups.')
 final=[r for r in a['accuracy'] if r['split']=='final_test' and r['mode']=='causal' and r['gap']==3 and not r['multipass']];final=sorted(final,key=lambda r:r['backend']);fig,axs=plt.subplots(1,3,figsize=(15,4.8));axs[0].bar(LABEL,[r['epe_valid_mean'] for r in final],color=COL);axs[0].set(ylabel='Mean valid-vector EPE (px)',title='(a) Held-out procedural motion');axs[0].tick_params(axis='x',rotation=35)
 mat=np.array([r['six_mae_accepted'] for r in final])/[5,5,.2,.2,.1,.1];im=axs[1].imshow(mat,cmap='magma',vmin=0);axs[1].set(xticks=range(6),xticklabels=['vx','vy','div','vort','s1','s2'],yticks=range(5),yticklabels=LABEL,title='(b) Accepted six MAE / tolerance');fig.colorbar(im,ax=axs[1],shrink=.7,label='Normalized error')
 axs[2].bar(np.arange(5)-.18,[r['accepted'] for r in final],.36,color=COL,label='Accepted');axs[2].bar(np.arange(5)+.18,[r['false_accept'] for r in final],.36,color='#bd5d53',label='False accepted');axs[2].set(xticks=range(5),xticklabels=LABEL,ylabel='Cases out of 44',title='(c) Quality gates do not prove truth');axs[2].tick_params(axis='x',rotation=35);axs[2].legend()
 save(fig,'Fig_F10X_4_numerics_ground_truth','New procedural seeds from one recording; no independent physical velocity truth. Backend supports differ. Rust/SciPy displacement difference <= 3.6e-15 on same-pair checks.')
 fig,axs=plt.subplots(1,3,figsize=(15,4.8));modes=['raw','local','fixed','causal','clahe'];matrix=np.array([[next(r['false_accept'] for r in a['accuracy'] if r['split']=='validation' and r['backend']=='B4' and r['mode']==m and r['gap']==gap and r['ws']==48 and not r['multipass']) for gap in [1,2,3]] for m in modes]);im=axs[0].imshow(matrix,cmap='YlOrRd',vmin=0);axs[0].set(xticks=range(3),xticklabels=['33.3','66.6','99.9'],yticks=range(5),yticklabels=modes,xlabel='Endpoint interval (ms)',title='(a) Validation false accepts / 22')
 for iy in range(5):
  for ix in range(3):axs[0].text(ix,iy,str(matrix[iy,ix]),ha='center',va='center')
 rr=[next(r for r in a['accuracy'] if r['split']=='final_test' and r['backend']=='B4' and r['mode']==mode and r['multipass']==multi) for mode,multi in [('local',False),('causal',False),('causal',True)]];names=['Local\nsingle gap 2','Causal\nsingle gap 3','Causal\n64 -> 32 gap 3'];axs[1].bar(names,[r['false_accept'] for r in rr],color=['#b79d71','#805b9f','#5e9d89']);axs[1].set(ylabel='False accepts / 44',title='(b) Frozen final comparisons');axs[2].bar(names,[r['mean_total_ms'] for r in rr],color=['#b79d71','#805b9f','#5e9d89']);axs[2].set(ylabel='Mean synthetic pair time (ms)',title='(c) Accuracy / cost tradeoff')
 save(fig,'Fig_F10X_5_preprocessing_multipass','Primary configuration frozen before final: causal EMA alpha=0.08, 48/24, gap 3. Final multipass result is a diagnostic, not used to retune. Midpoint lag = 49.95 ms.')
 fig,axs=plt.subplots(1,3,figsize=(15,5.7))
 for j,fid in enumerate([100,450,936]):
  scene(axs[j],fid if fid in [450,936] else 100);boundary(axs[j]);s0=samples[str(fid)];v=vectors(axs[j],s0);axs[j].set_title(f'({chr(97+j)}) Frame {fid}, t={s0["source_PTS"]:.2f} s');fig.colorbar(v,ax=axs[j],label='Apparent speed (px/s)',shrink=.6)
 save(fig,'Fig_F10X_6_full_channel_motion','Only quality-accepted supported windows have arrows; arrow lengths show interval displacement x6. Brown points are unobserved, never filled as zero velocity.')
 fig=plt.figure(figsize=(16,7));gs=fig.add_gridspec(2,4,width_ratios=[1.6,1,1,1]);ax=fig.add_subplot(gs[:,0]);scene(ax);boundary(ax);ax.set_xlim(600,1320);available=[r for r in a['uncertainty'] if r['estimate'].get('uncertainty')];selected=[available[i] for i in np.linspace(0,len(available)-1,4).astype(int)]
 for i,r in enumerate(selected):
  gp=next(g for g in s['groups'] if g['observation_key']==r['key']);eig,vec=np.linalg.eigh(gp['covariance']);ang=np.degrees(np.arctan2(vec[1,-1],vec[0,-1]));ax.add_patch(Ellipse(gp['center_global_xy'],4*np.sqrt(eig[-1]),4*np.sqrt(eig[0]),angle=ang,fill=False,edgecolor=COL[i],lw=2));ax.text(*gp['center_global_xy'],str(i),color='black',fontsize=12)
 ax.set_title('(a) Spatial sets: channel crop')
 for j,n in enumerate(NAMES):
  ax=fig.add_subplot(gs[j//3,1+j%3])
  for i,r in enumerate(selected):
   e=r['estimate'];lo,hi=e['uncertainty']['intervals'][n];val=e['values'][n];ax.vlines(i,lo,hi,color=COL[i],lw=2);ax.plot(i,val,'o',color=COL[i])
  ax.set(xticks=range(len(selected)),xlabel='Current-frame spatial region',ylabel=('px/s' if j<2 else '1/s'),title=n.replace('_',' '));ax.axhline(0,color='gray',lw=.6)
 fig.subplots_adjust(wspace=.5,hspace=.4);save(fig,'Fig_F10X_7_local_group_motion','Frame 450. Nominal quadrant-bootstrap 95% intervals, coverage NOT_CALIBRATED. Regions have no physical identity; local non-affine false acceptance remains unresolved.')
 perf=a['performance'];off=[r for r in perf if r['name'].startswith('original_offline')];on=[r for r in perf if 'replay' in r['name']];fig,axs=plt.subplots(2,2,figsize=(13,8));axs[0,0].bar([r['name'].replace('original_offline_','') for r in off],[r['p95'] for r in off],color='#805b9f');axs[0,0].set(ylabel='Complete offline p95 (ms)',title='(a) Same original fit: scheduling cost');axs[0,0].axhline(33.3,color='#bf594f',ls='--',label='30 fps budget');axs[0,0].legend()
 names=[r['name'].replace('verified_replay_','') for r in on];axs[0,1].bar(names,[r['p95'] for r in on],color='#4e8daa');axs[0,1].axhline(33.3,color='#bf594f',ls='--');axs[0,1].set(ylabel='Complete paced pipeline p95 (ms)',title='(b) Actual-PTS replay; latest-slot = 1')
 axs[1,0].bar(np.arange(len(on))-.18,[r['fresh_group_frame_hz'] for r in on],.36,label='Fresh valid group-frame Hz',color='#619c87');axs[1,0].bar(np.arange(len(on))+.18,[r['full_channel_scan_hz'] for r in on],.36,label='Processed full scans / s',color='#4e8daa');axs[1,0].set(xticks=range(len(on)),xticklabels=names,ylabel='Rate (Hz)',title='(c) Fresh output vs processed frames');axs[1,0].legend(fontsize=8)
 axs[1,1].plot([r['mean_valid_windows']/1283 for r in off],[r['p95'] for r in off],'o-',color='#805b9f')
 for r in off:axs[1,1].annotate(r['name'].replace('original_offline_',''),(r['mean_valid_windows']/1283,r['p95']),xytext=(4,3),textcoords='offset points',fontsize=8)
 axs[1,1].set(xlabel='Mean fresh valid-window fraction',ylabel='Complete offline p95 (ms)',title='(d) Speedup can reduce fresh coverage');fig.subplots_adjust(hspace=.4);save(fig,'Fig_F10X_8_runtime_validity','Offline and paced replay are separate. Drops, warmup, JSON transport and state ages retained. Low GPU utilization is observed; no 30/60 fps or camera latency claim.')
 fig,axs=plt.subplots(2,2,figsize=(13,8));sel=diag['selection'];cases=list(dict.fromkeys(r['case'] for r in sel));rates=[sum(r['covered_by_active'] for r in sel if r['case']==k)/max(1,sum(r['visible'] for r in sel if r['case']==k)) for k in cases];axs[0,0].bar(cases,rates,color='#5e9d89');axs[0,0].tick_params(axis='x',rotation=30);axs[0,0].set(ylim=(0,1.05),ylabel='Known visible centers covered',title='(a) New-location selector diagnostic')
 for gap in [1,2,3]:
  rr=[r for r in diag['step'] if r['gap']==gap];axs[0,1].plot([r['end'] for r in rr],[r['estimate']['values']['translation_vx'] for r in rr],'-o',label=f'gap {gap}');axs[0,1].plot([r['end'] for r in rr],[r['six_truth'][0] for r in rr],'--',alpha=.6)
 axs[0,1].set(xlabel='Frame after abrupt motion onset',ylabel='Interval vx (px/s)',title='(b) Longer interval mixes motion states');axs[0,1].legend();axs[0,1].text(.04,.10,'Missing points: rejected',transform=axs[0,1].transAxes,fontsize=9)
 rr=diag['wall_rules'];axs[1,0].bar([str(r['fraction']) for r in rr],[r['n_wall_crossing_allowed'] for r in rr],color='#bd5d53');axs[1,0].set(xlabel='Minimum valid-pixel fraction',ylabel='Admitted wall-crossing windows',title='(c) Whole-window geometry safety')
 rr=[r for r in read('results/research_v1/results/accuracy_final_test.json')['rows'] if r['backend']=='B4' and r['mode']=='causal' and not r['multipass'] and r['false_accept']];axs[1,1].axis('off');axs[1,1].text(.02,.95,'Retained failure evidence',va='top',fontsize=13,weight='bold');axs[1,1].text(.02,.8,'Final single-pass false accepts:\n'+ '\n'.join(f"  {r['case']} (seed {r['seed']})" for r in rr)+'\n\nReal wall-speed error: NOT_VERIFIED\nReal PZT-OFF segment: NOT_IDENTIFIED\nCross-recording generalization: NOT_VERIFIED',va='top',fontsize=10)
 fig.subplots_adjust(hspace=.5);save(fig,'Fig_F10X_9_failure_and_latency','Synthetic diagnostic uses already exposed development scenes; not another blind test. Dashed step curves are interval-equivalent truth. No physical wall or PZT ground truth inferred.')
 inputs=['results/delivery_analysis.json','results/figure_sample_states.json','results/channel_audit.json','results/diagnostics_v1/post_freeze_diagnostics.json','results/research_v1/results/accuracy_final_test.json','configs/channel_mask.json'];dump(ROOT/'figures/figure_manifest.json',dict(identity=identity(),figures=manifest,input_sha256={p:sha(ROOT/p) for p in inputs},script_sha256=sha(Path(__file__)),wall_s=time.perf_counter()-start,visual_review='PENDING'));print('figures',len(manifest),flush=True)
if __name__=='__main__':main()
