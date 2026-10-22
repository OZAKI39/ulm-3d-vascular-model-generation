"""Render saved causal predictions; explicit PTS, no invented/interpolated frames."""
import sys,json,time
from pathlib import Path
from fractions import Fraction
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,cv2,av
from PIL import Image,ImageDraw,ImageFont
from vcommon import ROOT,OUT,P,NAMES,LIMITS,dump,sha,identity
from legacy import read_frames

W,H=1600,1000
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf';BOLD='/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
def font(n,b=False):return ImageFont.truetype(BOLD if b else FONT,n)
def text(draw,xy,s,size=24,color='#e2e7ef',bold=False):draw.text(xy,s,font=font(size,bold),fill=color)
def canvas(title,subtitle):
    im=Image.new('RGB',(W,H),'#080c13');d=ImageDraw.Draw(im);text(d,(36,22),title,31,bold=True);text(d,(36,65),subtitle,20,color='#aab8ca');return im,d

def microscopy(im,gray,objects,x=35,y=115,w=570,h=805):
    rgb=cv2.cvtColor(gray,cv2.COLOR_GRAY2RGB) if gray.ndim==2 else cv2.cvtColor(gray,cv2.COLOR_BGR2RGB)
    scale=min(w/rgb.shape[1],h/rgb.shape[0]);dw=round(rgb.shape[1]*scale);dh=round(rgb.shape[0]*scale)
    image=Image.fromarray(rgb).resize((dw,dh),Image.Resampling.LANCZOS);im.paste(image,(x,y));d=ImageDraw.Draw(im)
    for o in objects:
        xx,yy=o['xy'];cx=x+xx*scale;cy=y+yy*scale;d.ellipse((cx-3,cy-3,cx+3,cy+3),outline='#39c9ff',width=2)
    if objects:
        xy=np.array([o['xy'] for o in objects]);mx,my=xy.mean(0)*scale+[x,y];d.line((mx-9,my,mx+9,my),fill='#ffcf75',width=2);d.line((mx,my-9,mx,my+9),fill='#ffcf75',width=2)
    return scale,dw,dh

def motion_panel(draw,x,y,title,estimate,pipeline=None,truth=None):
    text(draw,(x,y),title,24,bold=True);y+=38
    if estimate is None:text(draw,(x,y),'NO CURRENT ESTIMATE',22,color='#a7afbd');return
    good=all(estimate.get('valid',{}).values()) and bool(estimate.get('valid'))
    status=estimate.get('observability_status','NOT_OBSERVABLE')
    text(draw,(x,y),'ACCEPTED' if good else 'REJECTED — diagnostic estimate only',20,color='#6bd7a2' if good else '#ffac80');y+=30
    raw=estimate.get('raw_estimate');raw=raw if raw is not None else [estimate.get('values',{}).get(n) for n in NAMES]
    names=['vx (px/s)','vy (px/s)','div (1/s)','vort (1/s)','s1 (1/s)','s2 (1/s)']
    for j,(label,value) in enumerate(zip(names,raw)):
        xx=x+(j%3)*285;yy=y+(j//3)*59;text(draw,(xx,yy),label,17,color='#94a9c3');text(draw,(xx,yy+22),'null' if value is None else f'{value:+.3f}',23)
    y+=125
    info=[]
    if pipeline is not None:info.append(f'processing {pipeline:.1f} ms')
    if truth is not None and all(v is not None for v in raw):info.append(f'max |error| / tolerance: {np.max(np.abs(np.array(raw)-truth)/LIMITS):.2f}')
    if estimate.get('confidence_score') is not None:info.append(f'risk score: {estimate["confidence_score"]:.2f}')
    if info:text(draw,(x,y),' | '.join(info),17,color='#a7b9cf')

class Writer:
    def __init__(self,path):
        self.path=path;self.c=av.open(str(path),'w');self.s=self.c.add_stream('libx264',rate=Fraction(10000,333));self.s.width=W;self.s.height=H;self.s.pix_fmt='yuv420p';self.s.options={'crf':'18','preset':'medium','threads':'4'};self.s.time_base=Fraction(1,1000000);self.pts=[]
    def add(self,image,pts):
        frame=av.VideoFrame.from_ndarray(np.asarray(image),format='rgb24');frame.pts=round(pts*1000000);frame.time_base=Fraction(1,1000000)
        for p in self.s.encode(frame):self.c.mux(p)
        self.pts.append(pts)
    def close(self):
        for p in self.s.encode():self.c.mux(p)
        self.c.close();return dict(file=self.path.name,sha256=sha(self.path),expected_frames=len(self.pts),expected_pts_s=self.pts,resolution=[W,H],interpolation=False)

def validate(record):
    path=OUT/'videos'/record['file'];pts=[];sample={0,record['expected_frames']//2,record['expected_frames']-1};shots=[]
    with av.open(str(path)) as c:
        for i,f in enumerate(c.decode(video=0)):
            pts.append(float(f.time));assert f.width==W and f.height==H
            if i in sample:
                target=OUT/'evidence'/f'{path.stem}_frame{i:04}.png';Image.fromarray(f.to_ndarray(format='rgb24')).save(target);shots.append(str(target.relative_to(OUT)))
    assert len(pts)==record['expected_frames'];assert np.all(np.diff(pts)>0);err=float(np.max(np.abs(np.array(pts)-record['expected_pts_s'])));assert err<.0001
    return dict(**record,fully_decoded=True,decoded_frames=len(pts),pts_monotonic=True,max_pts_error_s=err,first_pts_s=pts[0],last_pts_s=pts[-1],screenshots=shots,visual_review='PENDING')

def real_video():
    names=['PIV_SINGLE','direct_fused','geo_fused'];rr={n:json.loads((OUT/'results'/f'real_replay_{n}.json').read_text()) for n in names};by={n:{r['frame_id']:r for r in rr[n]['rows']} for n in names};begin,end=P['real_replay']['frame_range'];frames=read_frames(range(begin,end),P['source_video']);x,y,w,h=P['real_replay']['roi_xywh'];first=frames[begin][1];writer=Writer(OUT/'videos/Video_F10V_1_real_causal_replay.mp4')
    for fid,(bgr,pts) in sorted(frames.items()):
        im,d=canvas('Real Chronos | PIV and learned collective motion','REAL_RECORDING_CAUSAL_REPLAY  •  pixel units  •  real six-parameter accuracy NOT_VERIFIED')
        available=next((by[n].get(fid) for n in names if fid in by[n]),None);objects=available['objects'] if available else []
        gray=cv2.cvtColor(bgr[y:y+h,x:x+w],cv2.COLOR_BGR2GRAY);scale,dw,dh=microscopy(im,gray,objects)
        r=by['PIV_SINGLE'].get(fid)
        if r:
            cv=np.array(im)
            for xx,yy,u,v,good in r['vectors']:
                if not good:continue
                pt=(round(35+xx*scale),round(115+yy*scale));endpt=(round(pt[0]+u*scale*8),round(pt[1]+v*scale*8));cv2.arrowedLine(cv,pt,endpt,(255,207,117),1,tipLength=.3)
            im=Image.fromarray(cv)
        d=ImageDraw.Draw(im);text(d,(35,940),f'PTS {pts:.4f} s  |  frame {fid}  |  candidates {len(objects)}',21);text(d,(35,967),'Cyan: detections; gold: visible center / PIV arrows x8',17,color='#aab8ca')
        for j,n in enumerate(names):
            r=by[n].get(fid);motion_panel(d,655,120+j*270,{'PIV_SINGLE':'PIV single pass','direct_fused':'Direct6 + heatmap','geo_fused':'GeoMotion'}[n],r['state']['motion'] if r else None,r['pipeline_ms'] if r else None)
        text(d,(655,942),'No current estimate is carried across a skipped frame.',18,color='#aab8ca');writer.add(im,pts-first)
    return validate(writer.close())

def synthetic_videos():
    groups=[('Video_F10V_2_synthetic_truth.mp4',['translation','rotation_cw','contraction']),('Video_F10V_3_difficult_ablation.mp4',['fixed_background','defocus','dropout','step'])];outputs=[]
    for filename,cases in groups:
        writer=Writer(OUT/'videos'/filename);offset=0
        for case in cases:
            record=json.loads((OUT/'results/synthetic_demo'/f'{case}.json').read_text());data=np.load(OUT/'results/synthetic_demo'/f'{case}.npz');dt=record['dt_s']
            for row,frame in zip(record['records'],data['frames']):
                difficult='3_' in filename;im,d=canvas(('Difficult microscopy | heatmap ablation' if difficult else 'Synthetic motion | truth and frozen estimates'),f'SYNTHETIC_MOTION_TEST  •  {case.replace("_"," ")}  •  independent scene cuts; no interpolated frames')
                microscopy(im,frame,row['objects'],w=570,h=570);d=ImageDraw.Draw(im);text(d,(35,710),f'Scene time {row["time_s"]:.4f} s  |  frame {row["frame"]}',23)
                text(d,(35,748),'Known interval-equivalent truth',21,bold=True)
                if row['truth'] is not None:
                    for j,(lab,v) in enumerate(zip(['vx','vy','div','vort','s1','s2'],row['truth'])):text(d,(35+(j%3)*190,788+(j//3)*40),f'{lab}: {v:+.3f}',21)
                else:text(d,(35,790),'Causal window warming up',21,color='#aab8ca')
                text(d,(35,904),'Translation: px/s; deformation: 1/s',18,color='#aab8ca');text(d,(35,940),'Real-patch generator; not physical bubble simulation',17,color='#aab8ca')
                names=['direct_gray','direct_fused','direct_fused_corrupt'] if difficult else ['PIV_SINGLE','direct_fused','geo_fused']
                labels={'PIV_SINGLE':'PIV single pass','direct_gray':'Direct6 gray','direct_fused':'Direct6 + detected heatmap','geo_fused':'GeoMotion','direct_fused_corrupt':'Direct6 + corrupted detector heat'}
                for j,n in enumerate(names):motion_panel(d,655,120+j*270,labels[n],row['predictions'].get(n),truth=row['truth'])
                if difficult:text(d,(655,942),'Corruption: 60% detected-center dropout + 20 false points.',18,color='#aab8ca')
                writer.add(im,offset+row['time_s'])
            offset+=len(record['records'])*dt
        outputs.append(validate(writer.close()))
    return outputs

if __name__=='__main__':
    start=time.perf_counter();(OUT/'videos').mkdir(exist_ok=True);rows=[real_video()]+synthetic_videos();dump(OUT/'videos/manifest.json',dict(identity=identity(),wall_s=time.perf_counter()-start,videos=rows,scope='saved actual predictions; exact mapped PTS; no LIVE or hardware control claim'))
