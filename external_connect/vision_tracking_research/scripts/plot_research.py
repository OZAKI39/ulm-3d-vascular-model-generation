from pathlib import Path
import json,collections
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def read(p):return json.loads((ROOT/p).read_text())
def msgs(tag):return [json.loads(s) for s in (ROOT/f'results/{tag}/states.jsonl').read_text().splitlines()]
def save(fig,name):
 fig.savefig(OUT/f'{name}.png',dpi=250,bbox_inches='tight');fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
full=msgs('full_cc_target_roi');agg=read('results/comparison_aggregate.json');accuracy=read('results/accuracy_review.json');review=read('results/track_review.json')
# 1: actual current-frame observations and real IDs, no invented biological labels.
fid=138;im=np.asarray(Image.open(ROOT/f'evidence/frame_{fid:05d}.png'));state=next(m for m in msgs('cc_norfair_120_target_roi') if m['frame_id']==fid)
fig,axes=plt.subplots(1,2,figsize=(11.5,6.5),gridspec_kw={'width_ratios':[1.35,1]},layout='constrained')
axes[0].imshow(im);axes[0].add_patch(Rectangle((650,350),180,350,fill=False,edgecolor='#008bb3',linewidth=2));axes[0].set_title('(a) Original sample and analysis ROI');axes[0].set_xlabel('Encoded image x (px)');axes[0].set_ylabel('Encoded image y (px)')
axes[1].imshow(im);axes[1].set_xlim(650,830);axes[1].set_ylim(700,350);axes[1].set_title('(b) Online Norfair IDs at frame 138');axes[1].set_xlabel('Encoded image x (px)');axes[1].set_ylabel('Encoded image y (px)')
for t in state['targets']:
 if t['status']!='OBSERVED':continue
 x,y=t['xy_px'];col='#008bb3' if t['visual_valid'] else '#c77818';axes[1].plot(x,y,'o',mfc='none',mec=col,ms=8,mew=1);axes[1].annotate(str(t['id']),(x,y),xytext=(4,-3),textcoords='offset points',fontsize=7,color=col)
axes[1].axhline(600,color='#963d73',ls='--',lw=1,label='Target line y = 600 px')
axes[1].plot([],[],'o',mfc='none',color='#008bb3',label='Quality-screened visual feature');axes[1].plot([],[],'o',mfc='none',color='#c77818',label='Untrusted candidate');axes[1].legend(loc='upper left',bbox_to_anchor=(0,-.1),fontsize=8)
fig.suptitle('Real microscopy sample: detection and online IDs',fontsize=15)
fig.supxlabel('CAUSAL_REPLAY_TEST | Image features are not certified bubble identities; no spatial calibration.',fontsize=9)
save(fig,'Fig1_detection_online_ids')
# 2: measured annotation recall and actual identity instability/loss.
fig,axes=plt.subplots(2,2,figsize=(11.8,8),layout='constrained')
a=axes[0,0];a.imshow(im);a.set_xlim(680,815);a.set_ylim(515,350);a.set_title('(a) Reviewed short visual landmarks');a.set_xlabel('x (px)');a.set_ylabel('y (px)')
ann=read('evidence/visual_annotations.json');colors=['#007c91','#c77818','#804e9c']
for color,(name,ps) in zip(colors,ann['short_tracks'].items()):
 arr=np.array(list(ps.values()));a.plot(arr[:,0],arr[:,1],'-o',color=color,ms=4,label=name.replace('landmark_',''))
a.legend(title='Visual label',fontsize=8)
a=axes[0,1]
for color,name in zip(colors,ann['short_tracks']):
 t=next(t for t in accuracy['short_track_review'] if t['method']=='cc_norfair' and t['label']==name);pp=t['points'];times=[p['frame_id']*.0333 for p in pp];a.plot(times,[p['gt_xy'][1] for p in pp],'-',color=color,alpha=.4)
 for x,p in zip(times,pp):
  a.plot(x,p['gt_xy'][1],'o' if p['matched'] else 'x',color=color,ms=7)
  if p['matched']:a.annotate(f"ID {p['id']}",(x,p['gt_xy'][1]),xytext=(3,5),textcoords='offset points',fontsize=8,color=color)
a.set_title('(b) Detections do not guarantee stable identity');a.set_xlabel('Replay PTS (s)');a.set_ylabel('Reviewed landmark y (px)');a.plot([],[],'o',color='k',label='Matched online observation');a.plot([],[],'x',color='k',label='Missing / not initialized');a.legend(fontsize=8)
a=axes[1,0];obs=[sum(t['status']=='OBSERVED' for t in m['targets']) for m in full];lost=[sum(t['status']=='LOST' for t in m['targets']) for m in full];tt=[m['source_pts_s'] for m in full];a.plot(tt,obs,lw=1,color='#007c91',label='Observed tracks');a.plot(tt,lost,lw=1,color='#c77818',label='Temporarily lost tracks');a.set_title('(c) Full replay: births, losses and reappearances');a.set_xlabel('Replay PTS (s)');a.set_ylabel('Track count in ROI');a.legend(fontsize=8)
a=axes[1,1];ks=accuracy['keyframe_detection'];methods=['cc','blob'];vals=[sum(k['matched'] for k in ks if k['detector']==m) for m in methods];a.bar(['Connected components','SimpleBlobDetector'],vals,color=['#007c91','#c77818']);a.set_ylim(0,28);a.set_ylabel('Matched selected visual landmarks / 21');a.set_title('(d) Positive-only visual review')
for i,v in enumerate(vals):a.text(i,v+.4,f'{v}/21',ha='center')
a.text(.5,.96,'Biological precision and group-shape error:\nNOT_EVALUATED',transform=a.transAxes,ha='center',va='top',fontsize=9)
fig.suptitle('Accuracy limits and real online track discontinuities',fontsize=15);fig.supxlabel('Assistant visual review; ±3 px coordinate uncertainty. Reappearance of an ID does not establish physical identity.',fontsize=9);save(fig,'Fig2_tracks_and_accuracy')
# 3: same methods/clips within each distinct ROI condition.
fig,axes=plt.subplots(1,3,figsize=(14,5.4),layout='constrained');labels=['CC + NF','Blob + NF','CC + NF + LK','CC + NF + FB'];cols=['#007c91','#c77818','#7762a8','#78914b']
for ax,roi,title in zip(axes[:2],['wide','target_roi'],['(a) Broad ROI: 620 × 960 px','(b) Target ROI: 180 × 350 px']):
 rows=[r for r in agg if r['roi']==roi]
 for j,(r,c) in enumerate(zip(rows,cols)):
  s=r['stages']['processing_ms'];ax.bar(j,s['mean'],color=c,alpha=.85);ax.plot([j,j],[s['mean'],s['p95']],color='#333333',lw=1.5);ax.plot(j,s['p95'],'_',color='black',ms=12)
 ax.axhline(16.67,color='#a83b45',ls='--',lw=1,label='60 fps budget (reference)');ax.axhline(33.3,color='#555555',ls=':',lw=1,label='Actual encoded interval');ax.set_xticks(range(4),labels,rotation=30,ha='right');ax.set_title(title);ax.set_ylabel('Processing time (ms): mean bar, p95 cap');ax.legend(fontsize=8)
a=axes[2];x=np.arange(4)
for off,roi,col in [(-.18,'wide','#aa7d52'),(.18,'target_roi','#007c91')]:
 a.bar(x+off,[r['valid_visual_hz'] for r in agg if r['roi']==roi],.34,label='Broad ROI' if roi=='wide' else 'Target ROI',color=col)
a.axhline(1/.0333,color='#555555',ls=':');a.set_xticks(x,labels,rotation=30,ha='right');a.set_ylabel('Quality-screened states received (Hz)');a.set_title('(c) Mock consumer output');a.legend(fontsize=8)
fig.suptitle('Measured causal replay throughput — i7-13700HX CPU',fontsize=15);fig.supxlabel('Three identical 90-frame clips per ROI/method; decode and transport reported separately. NF: Norfair; FB: Farnebäck.',fontsize=9);save(fig,'Fig3_latency_throughput')
# 4: measured frame ages + bounded failure behavior, not a fabricated camera chain.
fig,axes=plt.subplots(2,2,figsize=(12,7.7),layout='constrained')
for ax,tag,title in [(axes[0,0],'full_cc_target_roi','(a) Full replay: receive-to-mock age'),(axes[0,1],'fault_cc_target_roi','(b) Injected slow consumer / outage / PTS jump')]:
 data=read(f'results/{tag}/mock_received.json');start=data[0]['received_mono'];timev=np.array([v['received_mono']-start for v in data]);age=np.array([v['age_ms'] for v in data]);valid=np.array([v['valid_visual'] for v in data]);ax.plot(timev,age,lw=.8,color='#555');ax.scatter(timev[valid],age[valid],s=5,color='#007c91',label='Accepted visual state');ax.scatter(timev[~valid],age[~valid],s=10,color='#bf5d39',label='Rejected / no reliable target');ax.axhline(100,color='#963d73',ls='--',label='100 ms admission limit');ax.set_xlabel('Measured replay wall time (s)');ax.set_ylabel('Decoded-frame age (ms)');ax.set_title(title);ax.legend(fontsize=7)
a=axes[1,0];methods=['Broad\nCC + NF','Target\nCC + NF','Slow-consumer\ntest'];stats=[next(r for r in agg if r['roi']=='wide' and r['method']=='cc_norfair')['drop_fraction'],read('results/full_cc_target_roi/summary.json')['drop_fraction'],read('results/fault_cc_target_roi/summary.json')['drop_fraction']];a.bar(methods,np.array(stats)*100,color=['#c77818','#007c91','#bf5d39']);a.set_ylabel('Decoded frames skipped (%)');a.set_title('(c) Latest-frame overwrite is reported')
for i,v in enumerate(stats):a.text(i,v*100+1,f'{100*v:.1f}%',ha='center')
a.set_ylim(0,90)
a=axes[1,1];co=read('results/cotracker_online/windows.json');xx=np.arange(len(co));a.bar(xx,[v['latest_frame_age_ms'] for v in co],color='#7762a8');a.axhline(100,color='#963d73',ls='--');a.set_xticks(xx,['First window\n(cold inference)','After overrun\n(reinitialized)']);a.set_ylabel('Newest input frame age at output (ms)');a.set_title('(d) CoTracker3 Online: 2 bounded emissions')
for j,v in enumerate(co):a.text(j,v['latest_frame_age_ms']+12,f"{v['latest_frame_age_ms']:.0f} ms",ha='center')
fig.suptitle('Data freshness and admission — hardware control disabled',fontsize=15);fig.supxlabel('One Python latest-frame slot; all expired fault-test messages rejected. Camera-to-controller end-to-end latency: NOT_MEASURED.',fontsize=9);save(fig,'Fig4_frame_age_and_admission')
print('Saved four PNG/PDF figure pairs')
