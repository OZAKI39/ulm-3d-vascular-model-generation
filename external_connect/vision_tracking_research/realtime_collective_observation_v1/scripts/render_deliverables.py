import sys,json,gzip,argparse,ast,functools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump,sha,clean
import numpy as np,av,cv2
from PIL import Image,ImageDraw,ImageFont
from fractions import Fraction
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
W,H=1600,1000
# Reuse previous scientific video writer, preserving source PTS.
writer_path=ROOT.parent/'collective_motion_v2/scripts/videos.py';node=next(n for n in ast.parse(writer_path.read_text()).body if isinstance(n,ast.ClassDef) and n.name=='Writer');scope=dict(av=av,np=np,Fraction=Fraction,sha=sha,W=W,H=H);exec(compile(ast.Module(body=[node],type_ignores=[]),str(writer_path),'exec'),scope);HistoricalWriter=scope['Writer']
class Writer(HistoricalWriter):
 def __init__(self,path):
  super().__init__(path);self.s.codec_context.time_base=Fraction(1,1000000)
font=functools.lru_cache()(lambda size:ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',size))
def label(d,xy,text,size=22,color='#e9eff4'):d.text(xy,text,font=font(size),fill=color)
def loadrows(p):
 with gzip.open(p,'rt') as f:return [json.loads(line) for line in f]
def background(title,subtitle):
 im=Image.new('RGB',(W,H),'#101820');d=ImageDraw.Draw(im);label(d,(30,20),title,32);label(d,(30,65),subtitle,19,'#b3c3d0');return im,d

def real_panel(gray,r):
 im,d=background('Chronos | Full-channel visible microbubble motion','RECORDED VIDEO ANALYSIS  |  No verified PZT events  |  Pixel units')
 a=cv2.cvtColor(gray,cv2.COLOR_GRAY2RGB);g=r['geometry'];x,y,_,_=g['crop_xywh'];cv2.polylines(a,[np.array(g['channel_polygon'],np.int32)],True,(80,210,190),2)
 for px,py in r['candidates']:cv2.circle(a,(round(x+px),round(y+py)),4,(90,220,250),1)
 for px,py,u,v in r['vectors'][::2]:
  delta=np.array([u,v])*.10;length=np.linalg.norm(delta);delta*=min(1.,35/max(length,1e-9));cv2.arrowedLine(a,(round(x+px),round(y+py)),(round(x+px+delta[0]),round(y+py+delta[1])),(245,217,125),1,tipLength=.3)
 scale=min(1090/a.shape[1],810/a.shape[0]);img=Image.fromarray(a).resize((round(a.shape[1]*scale),round(a.shape[0]*scale)));im.paste(img,(25,115));d=ImageDraw.Draw(im);r0=r['aggregation']['regions'][0]
 values=[('Frame / PTS',f"{r['frame_id']} / {r['source_pts']:.4f} s"),('Visible candidates',str(r0['count'])),('PIV valid / safe',f"{r['valid_windows']} / {r['safe_windows']}"),('NN distance',f"{r0['nn_px']:.2f} px" if r0['nn_px'] else 'null'),('Distribution Rg',f"{r0['rg_px']:.2f} px" if r0['rg_px'] else 'null'),('Local visible regions',str(len(r['aggregation']['clusters']))),('Aggregation trend',r['aggregation']['status'])]
 for i,(name,value) in enumerate(values):label(d,(1140,120+i*92),name,18,'#9eafbe');label(d,(1140,150+i*92),value.replace('_','\n') if len(value)>24 else value,20)
 label(d,(30,945),'Cyan: candidates | Gold: 100 ms PIV displacement (35 px cap) | Teal: channel',20);return im

def validate_video(path,expected=None):
 pts=[];size=None;shots=[]
 with av.open(str(path)) as c:
  for i,f in enumerate(c.decode(video=0)):
   pts.append(float(f.time));size=[f.width,f.height]
   if i in [0,20,300,600,936]:
    shot=ROOT/'evidence'/f'{path.stem}_frame{i}.png';Image.fromarray(f.to_ndarray(format='rgb24')).save(shot);shots.append(str(shot.relative_to(ROOT)))
 result=dict(path=str(path.relative_to(ROOT)),sha256=sha(path),frames=len(pts),first_pts=pts[0],last_pts=pts[-1],presentation_duration_s=pts[-1]-pts[0]+float(np.median(np.diff(pts))),resolution=size,pts_monotonic=bool(np.all(np.diff(pts)>0)),full_decode=True,screenshots=shots)
 if expected is not None:result.update(expected_frames=len(expected),max_pts_error_s=float(np.max(np.abs(np.array(pts)-expected))) if len(pts)==len(expected) else None)
 return result

def synthetic_video(syn,sr):
 clips=np.load(syn/'clips.npz');sequence=json.load(gzip.open(syn/'sequences.json.gz','rt'));lookup={r['case']:r['records'] for r in sequence if r['seed']==sr[0]['seed']};writer=Writer(ROOT/'videos/Video_F10Y_3_synthetic_aggregation_validation.mp4');count=0
 for case in clips.files:
  for i,frame in enumerate(clips[case]):
   record=lookup[case][i];im,d=background('Synthetic aggregation validation', 'SYNTHETIC_AGGREGATION_TEST  |  Real image patches + known 2D mathematical motion');im.paste(Image.fromarray(frame).convert('RGB').resize((800,800)),(30,120));label(d,(880,140),case,30);label(d,(880,210),f"Time: {record['time_s']:.2f} s",24);label(d,(880,275),'Image-estimated trend:',22,'#b3c3d0');label(d,(880,315),record['measured']['status'].replace('_','\n'),25);label(d,(880,465),'Mathematical truth-input trend:',22,'#b3c3d0');label(d,(880,510),record['gold']['status'].replace('_','\n'),25);label(d,(880,680),'No acoustic or physical fusion claim',19);label(d,(880,725),'Errors and NOT_OBSERVABLE retained',19);writer.add(im,count/30);count+=1
 writer.close();return validate_video(ROOT/'videos/Video_F10Y_3_synthetic_aggregation_validation.mp4',np.arange(count)/30)

def main():
 p=argparse.ArgumentParser();p.add_argument('--video',required=True);p.add_argument('--study',required=True);p.add_argument('--synthetic',required=True);p.add_argument('--stability');a=p.parse_args();study=Path(a.study);syn=Path(a.synthetic);rows=loadrows(study/'offline.jsonl.gz');t=np.array([r['source_pts'] for r in rows]);reg=[r['aggregation']['regions'][0] for r in rows];cfg=rows[0]['geometry'];outputs=[]
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
 def save(fig,name):
  for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'{name}.{ext}',bbox_inches='tight')
  plt.close(fig);outputs.append(name)
 fig,ax=plt.subplots(figsize=(12,6));ax.axis('off');boxes=[(.05,.75,'Chronos → LatestFrameStream\nWindows receiver: reused\nLIVE hardware NOT TESTED'),(.05,.30,'Real MP4 → PyAV / actual PTS\nREALTIME_REPLAY: tested'),(.40,.55,'FramePacket + latest slot\ncontinuity / session / age'),(.70,.75,'Whole-channel detection\nCuPy PIV + 2D trends'),(.70,.30,'TTL publisher → Qt UI\nasync recorder / local ZMQ')]
 for x,y,s in boxes:ax.text(x,y,s,transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.7',fc='#eef4f8',ec='#56788b'),va='center')
 for start,end in [((.30,.75),(.40,.58)),((.30,.30),(.40,.50)),((.60,.55),(.70,.72)),((.85,.65),(.85,.40))]:ax.annotate('',xy=end,xytext=start,xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#406a80',lw=1.8))
 ax.set_title('F10Y | Unified observation architecture',loc='left');ax.text(.05,.05,'Lab: Windows + WSL CPU (Rust). Remote GPU: server-local replay; no untested WAN raw-video path.',transform=ax.transAxes);save(fig,'Fig_F10Y_1_realtime_architecture')
 sample=None
 with av.open(a.video) as c:
  for i,f in enumerate(c.decode(video=0)):
   if i==450:sample=f.to_ndarray(format='gray');break
 r=rows[450];q=np.array(r['candidates']);vv=np.array(r['vectors']);x,y,_,_=cfg['crop_xywh'];fig,ax=plt.subplots(figsize=(10,8));ax.imshow(sample,cmap='gray');ax.plot(*np.array(cfg['channel_polygon']+cfg['channel_polygon'][:1]).T,color='teal',lw=1.5);ax.scatter(q[:,0]+x,q[:,1]+y,s=5,facecolors='none',edgecolors='#32a7d2',linewidth=.5)
 if len(vv):ax.quiver(vv[:,0]+x,vv[:,1]+y,vv[:,2],vv[:,3],color='#d29322',angles='xy',scale_units='xy',scale=12,width=.002)
 ax.set(xlabel='Image x / px',ylabel='Image y / px',title=f"Full observed channel | frame 450 | {len(q)} candidates, {len(vv)} valid PIV vectors");save(fig,'Fig_F10Y_2_full_channel_observation')
 fig,axs=plt.subplots(2,2,figsize=(12,7),layout='constrained');axs[0,0].plot(t,[v['count'] for v in reg]);axs[0,0].set_ylabel('Visible candidate count');cent=np.array([v['center'] for v in reg]);axs[0,1].plot(t,cent[:,0],label='x');axs[0,1].plot(t,cent[:,1],label='y');axs[0,1].legend();axs[0,1].set_ylabel('Visible center / px (crop)');
 for key in ['all','0:300','300:300','600:300']:
  vals=[next((v['mean_velocity'][1] for v in r['aggregation']['regions'] if v['key']==key and v['mean_velocity'] is not None),np.nan) for r in rows];axs[1,0].plot(t,vals,label=key,lw=.7)
 axs[1,0].legend(fontsize=8,ncol=2);axs[1,0].set_ylabel('Regional image-y velocity / px s⁻¹');axs[1,1].plot(t,[r['valid_windows'] for r in rows]);axs[1,1].set_ylabel('Valid PIV windows')
 for ax in axs.flat:ax.set_xlabel('Source PTS / s');ax.grid(alpha=.2)
 fig.suptitle('Population motion | 937 analyzed frames | no long-term bubble IDs');save(fig,'Fig_F10Y_3_population_motion')
 fig,axs=plt.subplots(2,3,figsize=(13,7),layout='constrained');keys=['density_per_10k_px','nn_px','rg_px','radial_px_s','divergence_s'];names=['Density / 10⁴ px²','Median nearest neighbor / px','Visible distribution Rg / px','Relative radial motion / px s⁻¹','Local divergence / s⁻¹']
 for ax,k,name in zip(axs.flat,keys,names):ax.plot(t,[v[k] for v in reg],lw=.8);ax.set(xlabel='Source PTS / s',ylabel=name);ax.grid(alpha=.2)
 statuses=[r['aggregation']['status'] for r in rows];labels=sorted(set(statuses));axs[1,2].barh(labels,[statuses.count(s) for s in labels]);axs[1,2].set_xlabel('Frames');fig.suptitle('Real recording: 2D visible trends; physical aggregation accuracy NOT VERIFIED');save(fig,'Fig_F10Y_4_aggregation_metrics')
 synsum=json.loads((syn/'summary.json').read_text());sr=synsum['rows'];fig,axs=plt.subplots(1,2,figsize=(13,6),layout='constrained');names=[f"{r['case']} / {r['seed']}" for r in sr];yy=np.arange(len(sr));axs[0].barh(yy,[1 if r['image_correct'] else -1 for r in sr],color=['#76b4a3' if r['image_correct'] else '#cb6e61' for r in sr]);axs[0].set(yticks=yy,yticklabels=names,xlabel='Correct direction / rejection (+1), error (−1)');axs[1].barh(yy,[r['image_response_s'] or 0 for r in sr],color='#8b9dc3');axs[1].set(yticks=yy,yticklabels=[],xlabel='First convergence indication / s (0 = none)');fig.suptitle('SYNTHETIC_AGGREGATION_TEST | Real-patch images, known mathematical motion');save(fig,'Fig_F10Y_5_synthetic_validation')
 summaries={n:json.loads((study/n/'summary.json').read_text()) for n in ['replay_1x','stress_2x','debug_half']};fig,axs=plt.subplots(2,2,figsize=(12,7),layout='constrained');keys=['scan_ms','piv_ms','aggregation_ms','analysis_ms','end_to_end_software_ms'];r=summaries['replay_1x'];xx=np.arange(len(keys));axs[0,0].bar(xx-.18,[r['stages'][k]['p50'] for k in keys],.36,label='p50');axs[0,0].bar(xx+.18,[r['stages'][k]['p95'] for k in keys],.36,label='p95');axs[0,0].set(xticks=xx,xticklabels=[k.replace('_ms','') for k in keys],ylabel='ms');axs[0,0].tick_params(axis='x',rotation=25);axs[0,0].legend()
 ns=list(summaries);axs[0,1].bar(ns,[summaries[n]['fresh_valid_hz'] for n in ns]);axs[0,1].set_ylabel('Fresh valid motion Hz');axs[1,0].bar(ns,[summaries[n]['analysis_skipped'] for n in ns]);axs[1,0].set_ylabel('Analysis skipped frames');rr=r['resources'];axs[1,1].plot([v['elapsed_s'] for v in rr],[v['rss_mib'] for v in rr],label='RSS MiB');axs[1,1].plot([v['elapsed_s'] for v in rr],[v.get('gpu_mib') for v in rr],label='GPU device MiB');axs[1,1].legend();axs[1,1].set_xlabel('Wall time / s');fig.suptitle('Measured real-time replay performance | 2× is synthetic rate stress');save(fig,'Fig_F10Y_6_realtime_performance')
 # Full real video: every source frame, original PTS.
 writer=Writer(ROOT/'videos/Video_F10Y_2_real_bubble_motion_and_distribution.mp4')
 with av.open(a.video) as c:
  for i,f in enumerate(c.decode(video=0)):writer.add(real_panel(f.to_ndarray(format='gray'),rows[i]),float(f.pts*f.time_base))
 writer.close();records=[validate_video(ROOT/'videos/Video_F10Y_2_real_bubble_motion_and_distribution.mp4',t)]
 records.append(synthetic_video(syn,sr))
 ui=ROOT/'videos/Video_F10Y_1_realtime_replay_observation.mp4'
 if ui.exists():records.append(validate_video(ui))
 dump(ROOT/'evidence/video_validation.json',records);dump(ROOT/'evidence/figures_validation.json',{'figures':outputs,'png_pdf_count':len(outputs)*2,'input_study':str(study),'input_synthetic':str(syn),'writer_reuse_sha':sha(writer_path),'visual_review':'PENDING'})
 print('rendered',len(outputs),len(records),flush=True)
if __name__=='__main__':main()
