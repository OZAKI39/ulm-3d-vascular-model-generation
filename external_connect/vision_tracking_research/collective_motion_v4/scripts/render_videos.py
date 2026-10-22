"""Full source PTS, saved actual results; reuse F10V PyAV codec and validator."""
import sys,json,gzip,time,functools,importlib.util
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
import cv2,av
from PIL import Image,ImageDraw
spec=importlib.util.spec_from_file_location('f10v_video_helpers',ROOT.parent/'collective_motion_v2/scripts/videos.py');vv=importlib.util.module_from_spec(spec);spec.loader.exec_module(vv);vv.OUT=ROOT;vv.font=functools.lru_cache(64)(vv.font)
text=vv.text;cfg=json.loads((ROOT/'configs/channel_mask.json').read_text());NAMES2=['vx','vy','div','vort','s1','s2']
def overlay(g,s,groups=False):
 rgb=cv2.cvtColor(g,cv2.COLOR_GRAY2RGB);cv2.polylines(rgb,[np.array(cfg['channel_polygon'],np.int32)],True,(245,204,91),2)
 for p in cfg['excluded_polygons']:cv2.polylines(rgb,[np.array(p,np.int32)],True,(230,124,112),2)
 w=s['windows'];dt=s['time_window']['dt_s']
 for i,(x,y) in enumerate(zip(w['x'],w['y'])):
  p=(round(x),round(y))
  if not w['safe'][i]:continue
  if w['valid'][i]:
   u,v=w['vx'][i]*dt,w['vy'][i]*dt;cv2.arrowedLine(rgb,p,(round(x+u*6),round(y+v*6)),(45,196,98),1,tipLength=.25)
  else:cv2.circle(rgb,p,1,(155,115,75) if w['active'][i] else (95,95,95),-1)
 for x,y in s['candidates_global_xy']:cv2.circle(rgb,(round(x),round(y)),2,(35,185,241),1)
 if groups:
  for i,gp in enumerate(s['groups']):
   if not all(gp['six']['valid'].values()):continue
   eig,vec=np.linalg.eigh(np.array(gp['covariance']));axes=tuple(max(2,round(2*np.sqrt(z))) for z in eig[::-1]);ang=np.degrees(np.arctan2(vec[1,-1],vec[0,-1]));center=tuple(round(v) for v in gp['center_global_xy']);cv2.ellipse(rgb,center,axes,ang,0,360,(242,175,65),2);cv2.putText(rgb,str(i),center,cv2.FONT_HERSHEY_SIMPLEX,.5,(24,40,70),2,cv2.LINE_AA)
 return Image.fromarray(rgb)
def curves(d,fid,compact,x=1040,y=740):
 t=np.array([r['pts'] for r in compact]);v=np.array([r['cell_1_1_six'] for r in compact],float)
 for j in range(6):
  xx=x+(j%2)*260;yy=y+(j//2)*77;w=230;h=42;finite=v[:,j][np.isfinite(v[:,j])];lo,hi=(finite.min(),finite.max()) if len(finite) else (-1,1);pad=max((hi-lo)*.1,.001);lo-=pad;hi+=pad
  text(d,(xx,yy),NAMES2[j]+f' [{lo:.2g}, {hi:.2g}]',14,color='#aec0d1');yy+=21;d.line((xx,yy+h,xx+w,yy+h),fill='#536579');points=[]
  for k in range(fid+1):
   if np.isfinite(v[k,j]):points.append((xx+w*t[k]/t[-1],yy+h*(hi-v[k,j])/(hi-lo)))
   else:
    if len(points)>1:d.line(points,fill='#52c9ef',width=1)
    points=[]
  if len(points)>1:d.line(points,fill='#52c9ef',width=1)
  cursor=xx+w*t[fid]/t[-1];d.line((cursor,yy,cursor,yy+h),fill='#edbd6a')
def main():
 start=time.perf_counter();analysis=json.loads((ROOT/'results/delivery_analysis.json').read_text());compact=analysis['compact'];path=ROOT/'results/full_video_verified_v5/results/offline_regular/records.jsonl.gz';writers=[vv.Writer(ROOT/'videos'/n) for n in ['Video_F10X_1_full_channel_PIV.mp4','Video_F10X_2_channel_and_collective_motion.mp4']];mapping=[]
 with av.open(str(VIDEO)) as container,gzip.open(path,'rt') as f:
  for fid,(frame,line) in enumerate(zip(container.decode(video=0),f,strict=True)):
   s=json.loads(line);pts=float(frame.time);assert s['frame_id']==fid and abs(pts-s['source_PTS'])<1e-8;g=frame.to_ndarray(format='gray');w=s['windows'];r=compact[fid];mapping.append(dict(frame=fid,source_pts_s=pts,output_pts_s=pts))
   im,d=vv.canvas('F10X | Full-channel PIV','OFFLINE 937/937 source frames | causal background | apparent pattern motion, pixel units')
   im.paste(overlay(g,s).resize((1080,820)),(24,114));d=ImageDraw.Draw(im);x=1130
   for yy,msg,size,col in [(119,f'Frame {fid} / 936',24,'#ffffff'),(156,f't = {pts:.4f} s',24,'#ffffff'),(210,f'Candidates: {r["candidates"]}',22,'#5bc8ef'),(250,f'Valid windows: {r["valid_windows"]}',22,'#69d799'),(288,'Safe windows: 1283',20,'#d0dbea'),(326,f'Support footprint: {100*r["supported_footprint_fraction"]:.1f}%',20,'#d0dbea'),(364,'No pixelwise velocity filling',18,'#a7b5c8'),(418,'Backend: CuPy / RTX 4090',20,'#d0dbea'),(454,'Window 48 px; overlap 24 px',19,'#d0dbea'),(490,'Interval: 3 frames / 99.9 ms',19,'#d0dbea'),(526,'Midpoint lag: 49.95 ms',19,'#d0dbea'),(575,'Unobserved window reasons',20,'#ffffff')]:text(d,(x,yy),msg,size,color=col)
   labels={1:'Outside safe geometry',2:'Insufficient support',3:'Window warmup',4:'Update deferred',5:'PIV quality rejected',6:'Registration failed',7:'Backend failed'}
   for j,k in enumerate([1,2,3,4,5,6,7]):text(d,(x,613+j*30),f'{labels[k]}: {w["reason_code"].count(k)}',17,color='#b9c4d4')
   text(d,(x,851),'Invalid velocity = null',20,color='#f1bf77');text(d,(x,887),'Physical accuracy: NOT_VERIFIED',16,color='#f1bf77');text(d,(24,955),'Cyan: candidates | Green: measured interval displacement x6 | Brown/gray: unobserved | Gold: editable channel boundary',18)
   writers[0].add(im,pts)
   im,d=vv.canvas('F10X | Channel and local collective motion','SPATIAL OBSERVATION REGIONS, NOT TRACKED PHYSICAL IDs | no inferred fusion or water velocity')
   im.paste(overlay(g,s,True).resize((980,744)),(24,120));d=ImageDraw.Draw(im);text(d,(24,893),f'Frame {fid} | t = {pts:.4f} s | valid spatial groups: {s["fresh_valid_group_count"]}',22);text(d,(24,932),'Ellipses: 2 x covariance axes | Labels are current-frame row indices',18);text(d,(24,963),'Arrows: measured displacement x6 | No estimates carried across missing observations',17)
   text(d,(1035,114),'Current local six parameters',23,bold=True);text(d,(1035,148),'#',16,color='#a7bed3')
   for k,n in enumerate(NAMES2):text(d,(1070+k*78,148),n,16,color='#a7bed3')
   text(d,(1035,172),'Translation: px/s; deformation: 1/s',15,color='#a7bed3')
   for j,gp in enumerate(s['groups'][:23]):
    vals=gp['six']['values'];text(d,(1035,199+j*21),str(j),15,color='#f5c477')
    for k,n in enumerate(NAMES):text(d,(1070+k*78,199+j*21),'null' if vals[n] is None else f'{vals[n]:+.2f}',14,color='#c6e4d6' if vals[n] is not None else '#ae9090')
   if len(s['groups'])>23:text(d,(1035,685),f'{len(s["groups"])-23} further regions saved in machine-readable state',14)
   text(d,(1035,713),'Fixed spatial cell (1,1); curves stop at current PTS',15,color='#e7c278');curves(d,fid,compact);writers[1].add(im,pts)
   if fid%200==0:print('full video render',fid,flush=True)
 outputs=[vv.validate(w.close()) for w in writers];dump(ROOT/'videos/source_pts_mapping.json',mapping);dump(ROOT/'videos/full_video_manifest.json',dict(identity=identity(),source_sha256=sha(VIDEO),videos=outputs,wall_s=time.perf_counter()-start,record_sha256=sha(path),encoder=dict(codec='libx264',crf=18,preset='medium',threads=4,rate='10000/333',time_base='1/1000000'),scope='Full offline source coverage; no live or real-time claim',visual_review='PENDING'))
 print('full videos fully decoded',flush=True)
if __name__=='__main__':main()
