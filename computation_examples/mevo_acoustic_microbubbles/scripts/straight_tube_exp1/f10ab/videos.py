"""Actual saved physical frames; unchanged true-radius orthographic sphere renderer."""
import json,time,numpy as np,h5py,av
from PIL import Image,ImageDraw
from scripts.straight_tube_exp1.f10aa.videos import scene,font
from .metrics import measure
from .figures import all_rows
from .common import *

def raw(row):
 with h5py.File(BASE/row['trajectory_path']) as h:return h['time_s'][:],h['positions_m'][:]
def main():
 complete=all_rows();rows=[r for r in complete if r['stage']=='discovery' and r['status']=='PASS'];out=BASE/'videos';out.mkdir(exist_ok=True);manifest=[]
 def pick(state,action):return next((r for r in rows if r['state_id']==state and r['action']==action),None)
 compare=[pick('S2_r0',a) for a in ['A0_OFF','A7_X210','A5_X180']];compare=[r for r in compare if r]
 chosen=pick('S2_r0','A7_X210') or compare[-1]
 limits=[pick(s,'A7_X210') for s in ['S0_r0','S2_r0','S6_r0']];limits=[r for r in limits if r]
 target_path=BASE/'analysis/TARGET.json';target=read(target_path) if target_path.exists() else None
 if target:
  limits=[next(r for r in complete if r['run_id']==c['run_id']) for c in sorted(target['cases'],key=lambda c:c['state_id'])]
 jobs=[('Video_F10AB_1_collective_shape_responses',compare,'Same initial shape; different acoustic actions',False),
       ('Video_F10AB_2_rigid_motion_vs_deformation',[chosen]*3,'Fixed IDs: raw / centered / Kabsch-aligned',True),
       ('Video_F10AB_3_shape_control_or_limitations',limits,'Target: +3% aspect in 2 ms | drift <20 nm | connected endpoint' if target else 'Actual limitations: initial-shape dependence; no target success implied',False)]
 for name,group,title,align in jobs:
  data=[raw(r) for r in group];end=min(t[-1] for t,x in data);physical=np.arange(0,end+1e-13,1e-5)
  path=out/(name+'.mp4');writer=av.open(str(path),'w');stream=writer.add_stream('libx264',rate=30);stream.width=1800;stream.height=880;stream.pix_fmt='yuv420p';stream.options={'crf':'19','preset':'medium'};mapping=[];begin=time.perf_counter()
  for frame,time_s in enumerate(physical):
   im=Image.new('RGB',(1800,880),'#080d16');d=ImageDraw.Draw(im)
   d.text((24,14),title,font=font(25),fill='#edf4fb');d.text((24,52),f'Physical time {time_s*1e3:.3f} ms | 2.5 MHz | R = 1.78 um | zero background flow',font=font(21),fill='#bacadb')
   sources=[]
   for k,(row,(tt,xx)) in enumerate(zip(group,data)):
    index=int(np.searchsorted(tt,time_s-1e-13));assert abs(tt[index]-time_s)<1e-13
    x=xx[index];m=measure(xx[0],x);center=x.mean(0);display=x
    label=row['state_id']+' | '+row['action']
    if target and name.endswith('3_shape_control_or_limitations'):
     admitted=next(c['admitted'] for c in target['cases'] if c['run_id']==row['run_id']);label+=' | final '+('PASS' if admitted else 'FAIL')
    if align:
     label=['Original positions','Translation removed','Rigid-aligned; orange: proximity'][k]
     if k==1:display=x-center
     if k==2:display=m['aligned_positions_m']
    ii,jj=np.triu_indices(len(x),1);neighbors=(np.linalg.norm(x[ii]-x[jj],axis=1)<=2.1*R).astype(float)
    im.paste(scene(display,size=(560,560),contact=neighbors if align and k==2 else None),(20+600*k,130));d=ImageDraw.Draw(im)
    d.text((24+600*k,95),label,font=font(21),fill='#8cc6f1')
    d.text((24+600*k,692),f'COM (um): {center[0]*1e6:.3f}, {center[1]*1e6:.3f}, {center[2]*1e6:.3f}',font=font(17),fill='#d9e3ed')
    d.text((24+600*k,719),f'Rg = {m["Rg_3d_m"]*1e6:.3f} um | aspect = {m["aspect_3d"]:.3f}',font=font(18),fill='#d9e3ed')
    axis=m['principal_axes_3d'][:,0];orientation='ORIENTATION_UNIDENTIFIABLE' if not np.isfinite(axis).all() else 'Major axis: '+', '.join(f'{v:.2f}' for v in axis)
    if align:orientation='Q eigenvalues (um^2): '+', '.join(f'{v*1e12:.2f}' for v in m['eigenvalues_3d_m2'])
    d.text((24+600*k,746),orientation,font=font(16),fill='#d9e3ed')
    d.text((24+600*k,771),f'Nonrigid RMS = {m["similarity_rmsd_3d_m"]*1e9:.1f} nm | components = {m["proximity_components"]}',font=font(16),fill='#d9e3ed')
    schedule=read((BASE/row['trajectory_path']).parent/'INPUT.json')['schedule'];s=schedule['segments'][0]
    amps='/'.join(f'{v:.0f}' for v in s['amplitude_Pa_equivalent']);phases='/'.join(f'{np.degrees(v):.0f}' for v in s['phase_rad'])
    d.text((24+600*k,797),'Y-/Y+/X-/X+ A: '+amps+' Pa-eq',font=font(15),fill='#b6c8da')
    d.text((24+600*k,820),'Phase: '+phases+' deg',font=font(15),fill='#b6c8da')
    sources.append(dict(run_id=row['run_id'],trajectory=row['trajectory_path'],saved_index=index,sha256=row['trajectory_sha256']))
   d.text((24,851),'EXPLORATORY | Uncalibrated homogeneous-water proxy | Fixed physical frames; slow playback; geometric proximity is not fusion',font=font(16),fill='#b1bbc7')
   for packet in stream.encode(av.VideoFrame.from_image(im)):writer.mux(packet)
   if frame in [0,len(physical)//2,len(physical)-1]:im.save(out/(name+f'_{frame:04d}.png'))
   mapping.append(dict(frame=frame,time_s=float(time_s),sources=sources))
  for packet in stream.encode():writer.mux(packet)
  writer.close();pts=[]
  with av.open(str(path)) as decoder:
   for f in decoder.decode(video=0):pts.append(float(f.pts*f.time_base))
  assert len(pts)==len(physical) and np.allclose(pts,np.arange(len(pts))/30,atol=1e-6)
  record=dict(name=name,sha256=sha(path),fps=30,frames=len(pts),resolution=[1800,880],physical_duration_s=float(end),playback_s=len(pts)/30,fully_decoded=True,wall_s=time.perf_counter()-begin,mapping=mapping)
  write(out/(name+'.json'),record);manifest.append({k:v for k,v in record.items() if k!='mapping'})
 write(out/'MANIFEST.json',manifest)
if __name__=='__main__':main()
