"""Data-only publication plots and PTS-preserving videos. No image generation models."""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import cv2,av
from fractions import Fraction
from distribution import CandidateDetector,distribution
from synthetic import read_frames,build_sequence
from methods import piv
from affine import parameters

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelcolor':'#222222','text.color':'#222222','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','pdf.fonttype':42})
COLORS=['#2362a4','#c96726','#6f469c'];LABELS=['v0x (px/s)','v0y (px/s)','Divergence (1/s)','Vorticity (1/s)','Shear 1 (1/s)','Shear 2 (1/s)']

def savefig(fig,name):
    fig.savefig(ROOT/'figures'/f'{name}.png',dpi=220,bbox_inches='tight');fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)

def overlay(ax,im,objects,dist,vectors=None):
    ax.imshow(im,cmap='gray',vmin=0,vmax=255)
    if objects:
        xy=np.array([o['xy'] for o in objects]);ax.scatter(*xy.T,s=10,facecolor='none',edgecolor='#00bdd6',linewidth=.65,label='Candidates')
    if dist.get('group_center_xy') is not None:
        center=np.array(dist['group_center_xy']);ax.scatter(*center,marker='+',s=70,c='#d42a72',linewidth=1.7)
        if dist.get('covariance_2x2') is not None:
            cov=np.array(dist['covariance_2x2']);ev,vec=np.linalg.eigh(cov);angle=np.degrees(np.arctan2(vec[1,-1],vec[0,-1]));ax.add_patch(Ellipse(center,2*np.sqrt(ev[-1]),2*np.sqrt(ev[0]),angle=angle,fill=False,color='#d42a72',lw=1.5))
    if vectors:
        v=np.array(vectors);v=v[v[:,4]>0];ax.quiver(v[:,0],v[:,1],v[:,2],v[:,3],angles='xy',scale_units='xy',scale=1/3,color='#00e280',width=.004)
    ax.set_xlabel('Image x (px)');ax.set_ylabel('Image y (px)');ax.set_aspect('equal')

def plots():
    real=json.loads((ROOT/'results/real_observation.json').read_text());frames=np.load(ROOT/'results/real_frames.npz');rows=real['rows'];m=json.loads((ROOT/'results/metrics.json').read_text());r=rows[40];im=frames['images'][40]
    fig,axs=plt.subplots(1,3,figsize=(11,5),constrained_layout=True);overlay(axs[0],cv2.cvtColor(im,cv2.COLOR_BGR2GRAY),r['objects'],r['distribution']);axs[0].set_title('(a) Visible candidate distribution')
    _,density=distribution(r['objects'],im.shape);z=axs[1].imshow(density,cmap='magma');fig.colorbar(z,ax=axs[1],label='Candidates / px²',shrink=.6);axs[1].set_title('(b) Spatial density (bandwidth 12 px)');axs[1].set_xlabel('Image x (px)');axs[1].set_ylabel('Image y (px)')
    t=np.array([x['pts'] for x in rows]);cent=np.array([x['distribution']['group_center_xy'] for x in rows]);axs[2].plot(t,cent[:,0]-cent[0,0],label='Visible mean x');axs[2].plot(t,cent[:,1]-cent[0,1],label='Visible mean y');axs[2].set_xlabel('Source PTS (s)');axs[2].set_ylabel('Change of visible mean (px)');axs[2].legend();axs[2].set_title('(c) Distribution changes, not material velocity')
    fig.suptitle('Recorded Chronos | no persistent bubble IDs | uncalibrated 2D distribution',fontsize=12);savefig(fig,'Fig1_collective_distribution')
    fig,axs=plt.subplots(3,3,figsize=(12,10),constrained_layout=True);det=CandidateDetector()
    for ax,case in zip(axs[0],['rotation_cw','contraction','shear']):
        z=np.load(ROOT/f'results/synthetic/test/7201_{case}.npz');ims=z['frames'];o1,g1=det(ims[11]);o2,g2=det(ims[15]);d,_=distribution(o2,g2.shape);q=piv(g1,g2,.1332,'local',o1,o2,multipass=True)
        overlay(ax,g2,o2,d,q['vectors']);ax.set_title('Synthetic: '+case.replace('_',' ')+'\nDisplacement arrows ×3')
    for j,ax in enumerate(axs[1:].ravel()):
        for n,c in zip(P['windows'],COLORS):
            y=[x['methods'].get('PIV_'+str(n),{}).get('estimate',{}).get('values',{}).get(NAMES[j],np.nan) for x in rows];ax.plot(t,y,color=c,label=f'{n} frames',lw=1)
        ax.set_ylabel(LABELS[j]);ax.set_xlabel('Recorded source PTS (s)');ax.grid(alpha=.2)
    axs[1,0].legend(fontsize=8);fig.suptitle('Six-parameter interval estimates | real-video absolute motion accuracy NOT VERIFIED',fontsize=12);savefig(fig,'Fig2_motion_field_six_parameters')
    fig,axs=plt.subplots(2,3,figsize=(11,6),constrained_layout=True)
    methods=['PIV_multipass','CPD'];comparison=m['truth_tests']
    if 'matched_AI_comparison' in comparison:methods+=['SEA_RAFT_S'];comparison={k:{'overall':v} for k,v in comparison['matched_AI_comparison'].items()}
    for j,ax in enumerate(axs.ravel()):
        vals=[comparison[k]['overall']['per_parameter'][NAMES[j]]['mae'] for k in methods]
        ax.bar(np.arange(len(methods)),[v if v is not None else np.nan for v in vals],color=COLORS[:len(methods)]);ax.set_xticks(range(len(methods)),['PIV','CPD','SEA-RAFT'][:len(methods)]);ax.set_ylabel('MAE: '+LABELS[j]);ax.grid(axis='y',alpha=.2)
        for i,k in enumerate(methods):ax.text(i,0,f"n={comparison[k]['overall']['per_parameter'][NAMES[j]]['valid_n']}",ha='center',va='bottom',fontsize=8)
    fig.suptitle('Known synthetic truth | same endpoint pairs | accepted outputs only; unequal rejection rates',fontsize=12);savefig(fig,'Fig3_ground_truth_recovery')
    fig,axs=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
    for method,c in zip(['PIV_multipass','CPD'],COLORS):
        w=m['truth_tests'][method]['windows'];xs=P['windows'];axs[0].plot(xs,[w[str(n)]['all_six_valid_fraction'] for n in xs],'-o',label=method,color=c)
        errs=[np.mean([p['mae']/p['error_limit'] for p in w[str(n)]['per_parameter'].values() if p['mae'] is not None]) for n in xs];axs[1].plot(xs,errs,'-o',label=method,color=c)
    axs[0].set_ylabel('All-six output fraction');axs[1].set_ylabel('Mean normalized absolute error');axs[0].legend(fontsize=8)
    resp=json.loads((ROOT/'results/window_response.json').read_text())
    for n,c in zip(P['windows'],COLORS):
        rr=[r for r in resp['rows'] if r['window']==n];axs[2].plot([r['end_t'] for r in rr],[r['estimate']['values']['translation_vx'] for r in rr],'-o',ms=2,color=c,label=f'{n} frames')
    axs[2].axvline(resp['onset_t'],c='k',ls=':');axs[2].axhline(30,c='k',ls='--',lw=1);axs[2].set_xlabel('Synthetic time (s)');axs[2].set_ylabel('Step response v0x (px/s)');axs[2].legend(fontsize=8)
    for ax in axs[:2]:ax.set_xlabel('Number of frames');ax.set_xticks(P['windows']);ax.grid(alpha=.2)
    fig.suptitle('Causal windows: 66.6 / 133.2 / 299.7 ms | longer is not automatically more accurate');savefig(fig,'Fig4_window_tradeoff')
    fig,axs=plt.subplots(1,3,figsize=(13,5),constrained_layout=True);runtime=json.loads((ROOT/'results/causal_replay.json').read_text());fast=json.loads((ROOT/'results/causal_replay_fast.json').read_text());names=['PIV 2-pass\nfull pipeline','PIV 1-pass\nfull pipeline','CPD\nalgorithm'];timings=[runtime['software_processing_ms'],fast['software_processing_ms'],m['truth_tests']['CPD']['overall']['algorithm_ms']]
    if (ROOT/'results/ai_comparison.json').exists():
        ai=json.loads((ROOT/'results/ai_comparison.json').read_text())
        if ai['status']=='COMPLETE':names+=['SEA-RAFT\nforward + reverse'];timings+=[ai['timing_ms']]
    ix=np.arange(len(names));axs[0].bar(ix-.15,[s['p50'] for s in timings],.3,label='p50');axs[0].bar(ix+.15,[s['p95'] for s in timings],.3,label='p95');axs[0].set_xticks(ix,names,fontsize=8);axs[0].set_ylabel('Measured processing time (ms)');axs[0].legend();axs[0].axhline(100,c='r',ls=':',label='stale limit')
    cases=['translation','rotation_cw','fixed_background','dropout','defocus','sparse','symmetric','local_opposite'];mat=np.array([[m['truth_tests'][met]['cases'][c]['all_six_valid_fraction'] for c in cases] for met in ['PIV_multipass','CPD']]);z=axs[1].imshow(mat,aspect='auto',vmin=0,vmax=1,cmap='viridis');axs[1].set_xticks(range(len(cases)),[s.replace('_','\n') for s in cases],rotation=75,fontsize=8);axs[1].set_yticks([0,1],['PIV','CPD']);axs[1].set_title('Output fraction (not correctness)');fig.colorbar(z,ax=axs[1],shrink=.6)
    axs[2].axis('off');axs[2].text(0,1,'Admission and limits\n\n• Synthetic affine recovery: LIMITED PASS\n• Fixed texture / dropout: failures retained\n• Nominal confidence coverage: FAIL\n• Real absolute six-parameter truth: unavailable\n• Physical 3D aggregation: NOT INFERRED\n• Camera-to-controller latency: NOT MEASURED\n• PZT actions: NOT RUN\n• training_ready = false',va='top',fontsize=10,linespacing=1.6)
    fig.suptitle('Measured throughput and conservative scientific admission');savefig(fig,'Fig5_runtime_and_admission')

class Writer:
    def __init__(self,name):
        self.path=ROOT/'videos'/name;self.c=av.open(str(self.path),'w');self.s=self.c.add_stream('libx264',rate=Fraction(10000,333));self.s.width=1280;self.s.height=720;self.s.pix_fmt='yuv420p';self.s.time_base=Fraction(1,1000000);self.s.options={'crf':'18','preset':'fast'};self.pts=[]
    def add(self,bgr,t):
        f=av.VideoFrame.from_ndarray(bgr,format='bgr24');f.pts=round(t*1e6);f.time_base=Fraction(1,1000000);self.pts.append(t)
        for p in self.s.encode(f):self.c.mux(p)
    def close(self):
        for p in self.s.encode():self.c.mux(p)
        self.c.close()
        with av.open(str(self.path)) as c:ts=[float(f.time) for f in c.decode(video=0)]
        if len(ts)!=len(self.pts) or not np.all(np.diff(ts)>0):raise RuntimeError('Video complete decode / monotonic PTS failed')
        return {'file':str(self.path),'sha256':sha(self.path),'frames':len(ts),'first_pts':ts[0],'last_pts':ts[-1],'max_pts_error_s':float(np.max(np.abs(np.array(ts)-self.pts))),'complete_decode':True,'interpolation':False}

def text(canvas,s,xy,scale=.65,color=(225,225,225)):
    cv2.putText(canvas,s,xy,cv2.FONT_HERSHEY_SIMPLEX,scale,color,1,cv2.LINE_AA)

def panel(im,objects,dist,vec=None):
    if im.ndim==2:im=cv2.cvtColor(im,cv2.COLOR_GRAY2BGR)
    image=im.copy()
    for o in objects:cv2.circle(image,tuple(np.round(o['xy']).astype(int)),2,(255,200,0),1,cv2.LINE_AA)
    if dist.get('group_center_xy') is not None:
        c=tuple(np.round(dist['group_center_xy']).astype(int));cv2.drawMarker(image,c,(160,40,255),cv2.MARKER_CROSS,10,1)
        if dist.get('covariance_2x2') is not None:
            ev,evec=np.linalg.eigh(dist['covariance_2x2']);ax=tuple(np.maximum(1,np.round(np.sqrt(ev[::-1]))).astype(int));angle=np.degrees(np.arctan2(evec[1,-1],evec[0,-1]));cv2.ellipse(image,c,ax,angle,0,360,(160,40,255),1,cv2.LINE_AA)
    for row in vec or []:
        x,y,u,v,good=row
        if good and np.isfinite([x,y,u,v]).all():cv2.arrowedLine(image,(round(x),round(y)),(round(x+3*u),round(y+3*v)),(0,230,130),1,cv2.LINE_AA,tipLength=.25)
    scale=min(560/im.shape[1],560/im.shape[0]);return cv2.resize(image,None,fx=scale,fy=scale,interpolation=cv2.INTER_LINEAR)

def video_real():
    z=np.load(ROOT/'results/real_frames.npz');rows=json.loads((ROOT/'results/real_observation.json').read_text())['rows'];writer=Writer('Video1_real_collective_motion.mp4')
    for i,(im,r) in enumerate(zip(z['images'],rows)):
        can=np.full((720,1280,3),18,np.uint8);q=r['methods'].get('PIV_3',{});est=q.get('estimate',{});tile=panel(im,r['objects'],r['distribution'],q.get('vectors'));can[120:120+tile.shape[0],40:40+tile.shape[1]]=tile
        text(can,'F10U | Collective motion without persistent bubble IDs',(35,40),.85);text(can,'CAUSAL_REPLAY_TEST | recorded Chronos | offline rendering',(35,76),.7,(80,215,255));text(can,f"Source PTS {r['pts']:.4f} s   Frame {r['frame_id']}   Visible candidates {len(r['objects'])}",(480,132))
        text(can,'3-frame causal endpoint estimate; no future smoothing',(480,175),.62)
        for j,k in enumerate(NAMES):
            v=est.get('values',{}).get(k);s=f'{v:+.4f}' if v is not None else 'NOT_OBSERVABLE';text(can,f'{LABELS[j]:<20} {s}',(500,235+j*51),.71)
        text(can,'Cyan: candidate centers   Magenta: visible covariance',(480,570),.6);text(can,'Green: measured displacement arrows x3',(480,602),.6);text(can,'2D visible statistics; absolute motion accuracy unverified',(480,638),.61,(90,180,255));text(can,'Source playback cadence preserved; software latency measured separately',(35,700),.57)
        writer.add(can,float(z['pts'][i]-z['pts'][0]))
        if i in [0,40,119]:cv2.imwrite(str(ROOT/f'evidence/video1_frame_{i}.png'),can)
    return writer.close()

def video_synthetic():
    source=read_frames([840])[840];det=CandidateDetector();writer=Writer('Video2_synthetic_truth_vs_estimate.mp4');records=[];outdir=ROOT/'results/synthetic_visualization';outdir.mkdir(exist_ok=True);case_names=['translation','rotation_cw','contraction','shear','fixed_background','local_opposite'];counter=0
    for name in case_names:
        ims,meta=build_sequence(name,9001,source,length=48);H=meta.pop('H');np.savez_compressed(outdir/f'{name}.npz',frames=ims,H=H,points=meta.pop('points'));meta.pop('instantaneous_generator');dump(outdir/f'{name}.json',meta);past=[]
        for i,im in enumerate(ims):
            o,g=det(im);d,_=distribution(o,g.shape);past.append((g,o));truth=None;q={}
            if i>=2:
                q=piv(past[-3][0],g,.0666,'local',past[-3][1],o,multipass=True)
                if name!='local_opposite':HH=H[i]@np.linalg.inv(H[i-2]);truth=parameters(HH[:2,:2],HH[:2,2],.0666,np.array([128.,128.]))
            est=q.get('estimate',{});can=np.full((720,1280,3),18,np.uint8);tile=panel(im,o,d,q.get('vectors'));can[120:120+tile.shape[0],35:35+tile.shape[1]]=tile
            text(can,'SYNTHETIC_MOTION_TEST | real Chronos patches + analytic motion',(30,42),.82,(80,215,255));text(can,f'{name} | t = {i*.0333:.4f} s | 3-frame window',(30,82),.72)
            text(can,'Parameter                 Truth          Estimate',(620,165),.66)
            for j,k in enumerate(NAMES):
                tv=f'{truth[j]:+.3f}' if truth is not None else 'N/A';val=est.get('values',{}).get(k);v=f'{val:+.3f}' if val is not None else 'NOT_OBSERVABLE';text(can,f'{LABELS[j]}',(620,220+j*55),.61);text(can,tv,(900,220+j*55),.63);text(can,v,(1015,220+j*55),.58)
            text(can,'One global affine model is invalid here' if name=='local_opposite' else 'Translation: px/s; rates: 1/s',(620,610),.62,(90,180,255));text(can,'Errors and rejected estimates retained; no acoustic or PZT experiment',(620,650),.49);text(can,'Displacement arrows x3; display zoom preserves image aspect ratio',(30,704),.56)
            writer.add(can,counter*.0333);counter+=1
            q.pop('vectors',None);records.append({'case':name,'frame':i,'truth':truth,'result':q})
            if i==30:cv2.imwrite(str(ROOT/f'evidence/video2_{name}.png'),can)
        print('synthetic video',name,flush=True)
    dump(ROOT/'results/synthetic_video_estimates.json',{'scope':'post-benchmark visualization cases; no thresholds changed; not included in heldout benchmark metrics','identity':identity(),'rows':records})
    return writer.close()

if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('mode',choices=['plots','videos']);mode=a.parse_args().mode
    if mode=='plots':plots()
    else:dump(ROOT/'evidence/video_decode_verification.json',[video_real(),video_synthetic()])
