"""True coordinates in a fixed oblique 3D view; no synthetic physical motion."""
import time,sys
import numpy as np,h5py,av
from PIL import Image,ImageDraw,ImageFont
from .common import *
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
def font(n):return ImageFont.truetype(FONT,n)
def scene(x,roots,size=(900,680),capture=False):
 w,h=size;ss=2;scale=min(w/.00135,h/.0011);im=Image.new('RGB',(w*ss,h*ss),'#080d16');d=ImageDraw.Draw(im)
 # Orthographic camera, orthonormal projected axes; no anisotropic particle transform.
 az=np.radians(-58);el=np.radians(32);right=np.array([-np.sin(az),np.cos(az),0]);up=np.array([-np.sin(el)*np.cos(az),-np.sin(el)*np.sin(az),np.cos(el)]);depth=np.cross(right,up)
 def proj(q):return ((w/2+scale*np.dot(q,right))*ss,(h/2-scale*np.dot(q,up))*ss)
 for gx in np.arange(-.0005,.00051,.0001):d.line([proj([gx,-.0004,-.00008]),proj([gx,.0004,-.00008])],fill='#263346',width=1)
 for gy in np.arange(-.0004,.00041,.0001):d.line([proj([-.0005,gy,-.00008]),proj([.0005,gy,-.00008])],fill='#263346',width=1)
 for root in roots:
  corners=[[root,-.0004,-.00008],[root,.0004,-.00008],[root,.0004,.00008],[root,-.0004,.00008],[root,-.0004,-.00008]]
  d.line([proj(q) for q in corners],fill='#756783',width=ss)
 radius=R*scale*ss
 for q in x[np.argsort(x@depth)]:
  px,py=proj(q);col='#66cdb4' if capture and len(roots) and abs(q[0]-roots).min()<20e-6 else '#55a9ec';d.ellipse((px-radius,py-radius,px+radius,py+radius),fill=col);rr=.45*radius;d.ellipse((px-rr-.2*radius,py-rr-.2*radius,px+rr-.2*radius,py+rr-.2*radius),fill='#d5e8fa')
 im=im.resize(size,Image.Resampling.LANCZOS);d=ImageDraw.Draw(im)
 origin=np.array([-.00052,-.0004,-.00008]);o=tuple(np.array(proj(origin))/ss)
 for j,(name,color) in enumerate([('X','#dcac9e'),('Y','#a8cdb8'),('Z','#adbedb')]):
  e=tuple(np.array(proj(origin+np.eye(3)[j]*100e-6))/ss);d.line([o,e],fill=color,width=2);d.text(e,name,font=font(17),fill=color)
 d.text((25,h-40),'3D orthographic | sphere radius = 1.78 µm | axes = 100 µm',font=font(16),fill='#c5d1e0');return im

def main():
 begin=time.perf_counter();out=BASE/'videos';out.mkdir(exist_ok=True);s=read(BASE/'SCIENTIFIC_RESULTS.json');good=[r for r in s['rows'] if r['kind']=='HF' and r['status']=='PASS'];st=read(BASE/'STATIC_RESULTS.json');roots=np.array([r['position_m'][0] for a in st['rows'] if a['frequency_Hz']==2500000 and a['phase_deg']==0 for r in a['roots']]);rows=[]
 def hf(id,label):
  r=next(r for r in good if r['id']==id);p=BASE/'dataset'/(id+'.h5')
  with h5py.File(p) as h:tt=h['time_s'][:];xx=h['positions_m'][:]
  z=np.load(BASE/'series'/('HF_'+id+'.npz'));return dict(id=id,label=label,time=tt,x=xx,science=r,series=z,path=p,scope='FULL_F10AA_PHYSICS')
 with h5py.File(BASE/'single/atlas.h5') as h:tt=h['time_s'][:];xx=h['positions_m'][:]
 # 501 saved uniform dense samples -> exact 501-frame video, no time rescaling per segment.
 single=dict(id='single_atlas',label='Independent single bubbles | F1 + Stokes drag',time=tt,x=xx,path=BASE/'single/atlas.h5',scope='SINGLE_BUBBLE_DIAGNOSTIC')
 # Camera is fixed, displayed coordinates and radii have the same scale.
 jobs=[('Video_F10AG_1_3d_acoustic_capture','3D capture dynamics',[single]),('Video_F10AG_2_long_time_stripe_dynamics','Full multibubble evolution',[hf('long150_random_0','Random initial cloud A'),hf('long150_random_1','Independent initial cloud B')]),('Video_F10AG_3_rotation_mechanism','Opposite phase orders: actual particle motion',[hf('ring50_ROT_PLUS','ROT_PLUS | prescribed ring'),hf('ring50_ROT_MINUS','ROT_MINUS | prescribed ring')])]
 selected=set(sys.argv[1:]);prior=read(out/'MANIFEST.json')['videos'] if selected and (out/'MANIFEST.json').exists() else [];rows=[r for r in prior if Path(r['file']).stem not in selected]
 for name,title,panels in jobs:
  if selected and name not in selected:continue
  start=time.perf_counter();times=panels[0]['time'];assert np.allclose(np.diff(times),np.diff(times)[0],rtol=0,atol=1e-12)
  for p in panels[1:]:assert np.allclose(p['time'],times,rtol=0,atol=1e-12)
  path=out/(name+'.mp4');writer=av.open(str(path),'w');stream=writer.add_stream('libx264',rate=60);stream.width=1920;stream.height=1080;stream.pix_fmt='yuv420p';stream.options={'crf':'18','preset':'medium'};mapping=[]
  for i,t in enumerate(times):
   im=Image.new('RGB',(1920,1080),'#080d16');d=ImageDraw.Draw(im);d.text((35,20),title,font=font(31),fill='#f1f5fb');d.text((35,68),f'Physical time: {t*1000:.3f} ms | 2.5 MHz | zero background flow',font=font(23),fill='#cad9e9')
   for j,p in enumerate(panels):
    left=35+j*950 if len(panels)>1 else 200;size=(900,680) if len(panels)>1 else (1500,680);x=p['x'][i];im.paste(scene(x,roots if 'rotation' not in name else [],size,capture=True),(left,145));d=ImageDraw.Draw(im);d.text((left,112),p['label'],font=font(24),fill='#87c6f6')
    if p['scope']=='SINGLE_BUBBLE_DIAGNOSTIC':lines=['Y± OFF; X− / X+ = 1000 / 1000 Pa-eq; phase = 0 / 0 deg',f'Independent starts: {len(x)} | within 20 µm of X band: {int((abs(x[:,0,None]-roots).min(1)<=20e-6).sum())}','Green: instantaneous band proximity; held capture evaluated separately','Numerical status: PASS | F1-only diagnostic, not multibubble aggregation']
    else:
     z=p['series'];r=p['science'];task=r['task'];a=' / '.join(f'{v:g}' for v in task['amplitudes']);phase=' / '.join(f'{np.degrees(v):g}' for v in task['phases']);lines=[f'Y− / Y+ / X− / X+ : {a} Pa-eq',f'Phase: {phase} deg | N = {r["N"]} | numerical gates: PASS']
     if 'rotation' in name:
      k=int(np.argmin(abs(z['omega_time_s']-t)));angle=float(np.interp(t,z['angle_time_s'],z['angle_deg']));lines += [f'Omega_z = {z["omega_s"][k]:+.5f} rad/s | angle = {angle:+.4f} deg',f'Rigid rotation fraction = {z["rigid_fraction"][k]:.5f} (criterion 0.8)']
     else:lines += [f'Peak S = {z["S"][i]:.3f} | Rg = {z["Rg_m"][i]*1e6:.2f} µm',f'Force-root bands occupied (<20 µm) = {int((abs(x[:,0,None]-roots).min(axis=0)<=20e-6).sum())}']
    for k,line in enumerate(lines):d.text((left,835+32*k),line,font=font(18),fill='#d0ddeb')
   status='SINGLE_BUBBLE_DIAGNOSTIC | 3D local stability does not imply full-cloud capture' if len(panels)==1 else ('COLLECTIVE_ROTATION_UNCONFIRMED | no injected streaming or vortex flow' if 'rotation' in name else 'STRIPE_NOT_CONFIRMED | actual full-physics trajectories')
   d.text((35,987),status,font=font(21),fill='#e4bd8d');d.text((35,1038),'EXPLORATORY | HOMOGENEOUS_WATER_CONTROL_PROXY | PHYSICS_UNCALIBRATED | training_ready=false',font=font(18),fill='#b8c6d7')
   for packet in stream.encode(av.VideoFrame.from_image(im)):writer.mux(packet)
   if i in [0,len(times)//2,len(times)-1]:im.save(out/f'{name}_{i:04d}.png')
   mapping.append(dict(frame=i,physical_time_s=float(t),sources=[dict(id=p['id'],frame=i,scope=p['scope']) for p in panels]))
  for packet in stream.encode():writer.mux(packet)
  writer.close();pts=[]
  with av.open(str(path)) as reader:
   for frame in reader.decode(video=0):assert (frame.width,frame.height)==(1920,1080);pts.append(float(frame.pts*frame.time_base))
  assert len(pts)==len(times) and np.allclose(pts,np.arange(len(times))/60,atol=1e-7)
  r=dict(file=path.name,frames=len(times),fps=60,resolution=[1920,1080],physical_duration_s=float(times[-1]),playback_duration_s=len(times)/60,uniform_physical_dt_s=float(np.diff(times)[0]),sha256=sha(path),sources=[dict(id=p['id'],path=str(p['path'].relative_to(BASE)),sha256=sha(p['path']),scope=p['scope']) for p in panels],mapping=mapping,fully_decoded=True,wall_s=time.perf_counter()-start,command='PYTHONPATH=src:. .venv/bin/python -m scripts.straight_tube_exp1.f10ag.videos',identity=identity());write(out/(name+'.json'),r);rows.append({k:v for k,v in r.items() if k!='mapping'});print('VIDEO',name,flush=True)
 write(out/'MANIFEST.json',dict(videos=rows,wall_s=sum(r['wall_s'] for r in rows),last_render_wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
