"""Exact saved-state rendering. Playback slows physical time; never interpolates it."""
import argparse,json,math
from fractions import Fraction
from pathlib import Path
import av,numpy as np,h5py
from PIL import Image,ImageDraw,ImageFont
from .common import *
W,H=1800,1060
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
def text(d,xy,s,size=23,fill='#e2e8ed'):d.text(xy,s,font=ImageFont.truetype(FONT,size),fill=fill)
def load(name):
 with h5py.File(BASE/'dataset'/(name+'.h5')) as h:
  r={k:h[k][:] for k in h if isinstance(h[k],h5py.Dataset)};r['force_samples']={k:h['force_samples'][k][:] for k in h['force_samples']};r['schedule']=json.loads(h.attrs['action_schedule_json'])
 return r

def panel(d,q,box,limit,title,centered=False,force=None):
 left,top,width,height=box;scale=min(width,height)/(2*limit);center=np.array([left+width/2,top+height/2]);q=q.copy()
 if centered:q-=q.mean(0)
 for k in [-limit,-limit/2,0,limit/2,limit]:
  a=center+np.array([k*scale,0]);b=center+np.array([0,k*scale]);d.line([(a[0],top),(a[0],top+height)],fill='#182431',width=1);d.line([(left,b[1]),(left+width,b[1])],fill='#182431',width=1)
  text(d,(a[0]-15,top+height+4),f'{k:g}',15,'#a3b2c1')
 coords=q[:,:2]*1e6*np.array([1,-1])*scale+center
 for j in np.argsort(q[:,2]):
  a,b=coords[j];radius=R*1e6*scale
  d.ellipse([a-radius,b-radius,a+radius,b+radius],fill='#164e88',outline='#77b7eb',width=1)
  for f in [.9,.7,.5,.3]:
   z=radius*f;offset=(1-f)*radius*.25;color=(int(30+45*(1-f)),int(102+72*(1-f)),int(171+60*(1-f)))
   d.ellipse([a-z-offset,b-z-offset,a+z-offset,b+z-offset],fill=color)
 if force is not None:
  v=force[:,:2];vmax=np.linalg.norm(v,axis=1).max();factor=55/max(vmax,1e-30)
  for (x,y),(vx,vy) in zip(coords,v):
   dx,dy=vx*factor,-vy*factor;d.line([(x,y),(x+dx,y+dy)],fill='#ecc681',width=2)
   length=math.hypot(dx,dy)
   if length>3:
    ux,uy=dx/length,dy/length;d.polygon([(x+dx,y+dy),(x+dx-6*ux+3*uy,y+dy-6*uy-3*ux),(x+dx-6*ux-3*uy,y+dy-6*uy+3*ux)],fill='#ecc681')
 text(d,(left,top-37),title,23);text(d,(left+width/2-32,top+height+30),'X (µm)',17);text(d,(left+6,top+8),'Y (µm)',16,'#a3b2c1')

def plot(d,t,y,idx,box,title,color='#65b0ef'):
 l,b,w,h=box;y=np.array(y,float);finite=np.isfinite(y)
 if not finite.any():text(d,(l,b),title+' NOT OBSERVABLE',18);return
 low,high=np.nanmin(y),np.nanmax(y)
 if high-low<1e-9:low-=.5;high+=.5
 else:pad=.1*(high-low);low-=pad;high+=pad
 points=np.c_[l+(t-t[0])/max(t[-1]-t[0],1e-15)*w,b+h-(y-low)/(high-low)*h];d.line([tuple(v) for v in points[finite]],fill=color,width=2);j=min(idx,len(y)-1);x,yy=points[j];d.ellipse([x-4,yy-4,x+4,yy+4],fill='white');text(d,(l,b-28),f'{title}: {y[j]:.4g}',18);text(d,(l,b+h+3),f'{t[0]*1e3:g} ms',15,'#adbac6');text(d,(l+w-70,b+h+3),f'{t[-1]*1e3:g} ms',15,'#adbac6')

def base(title,sub,time_s,slow):
 im=Image.new('RGB',(W,H),'#05090e');d=ImageDraw.Draw(im);text(d,(35,22),title,32);text(d,(35,68),sub,21,'#acbdce');text(d,(1380,25),f't = {time_s*1e3:.3f} ms',29);text(d,(1380,70),f'Playback: {slow:g}× slower',18,'#acbdce');return im,d

def encode(path,frames,fps,mapping):
 path.parent.mkdir(exist_ok=True);out=av.open(str(path),'w');s=out.add_stream('libx264',rate=fps);s.width=W;s.height=H;s.pix_fmt='yuv420p';s.options={'crf':'19','preset':'fast'};s.time_base=Fraction(1,fps);s.codec_context.time_base=Fraction(1,fps)
 for k,im in enumerate(frames):
  frame=av.VideoFrame.from_ndarray(np.asarray(im),format='rgb24');frame.pts=k;frame.time_base=Fraction(1,fps)
  for packet in s.encode(frame):out.mux(packet)
 for packet in s.encode():out.mux(packet)
 out.close();pts=[];shots=[]
 with av.open(str(path)) as c:
  for k,f in enumerate(c.decode(video=0)):
   pts.append(float(f.pts*f.time_base))
   if k in [0,len(mapping)//2,len(mapping)-1]:
    shot=BASE/'evidence'/(path.stem+f'_frame{k}.png');f.to_image().save(shot);shots.append(str(shot.relative_to(BASE)))
 assert len(pts)==len(mapping) and np.all(np.diff(pts)>0)
 expected=np.arange(len(pts))/fps;error=float(np.max(abs(np.array(pts)-expected)));assert error<1e-5
 record=dict(path=str(path.relative_to(BASE)),frame_count=len(pts),fps=fps,first_pts_s=pts[0],last_pts_s=pts[-1],duration_s=len(pts)/fps,resolution=[W,H],sha256=sha(path),fully_decoded=True,PTS_monotonic=True,max_PTS_error_s=error,source_mapping=mapping,screenshots=shots,no_physical_interpolation=True)
 sources={s['run'] for row in mapping for s in (row.get('sources') or [row])};runs=read(BASE/'analysis/RUN_ANALYSIS.json');record['source_dataset_sha256']={n:sha(BASE/'dataset'/(n+'.h5')) for n in sources};record['source_raw_sha256']={n:runs[n]['trajectory_sha256'] for n in sources};record['render_identity']=identity()
 write(path.with_suffix('.json'),record);return record

def main():
 videos=BASE/'videos';videos.mkdir(exist_ok=True);records=[];all_runs=read(BASE/'analysis/RUN_ANALYSIS.json')
 names=[f'medium_train__{a}__2000us__DT1' for a in ['OFF','XM','PAIR']];ds=[load(n) for n in names];indices=np.arange(0,len(ds[0]['time_s']),10);fps=30;slow=1/fps/(ds[0]['time_s'][10]-ds[0]['time_s'][0]);mapping=[dict(frame=k,time_s=float(ds[0]['time_s'][j]),sources=[dict(run=n,saved_index=int(j)) for n in names]) for k,j in enumerate(indices)]
 def first():
  for j in indices:
   im,d=base('F10Z | Dispersed microbubble dynamics','Local XY projection of 32 true 3D bubbles | R = 1.78 µm | No imposed flow',ds[0]['time_s'][j],slow)
   for i,(r,action) in enumerate(zip(ds,['OFF','X− only: 400 Pa-eq','X pair: 1036.20 Pa-eq, 210°'])):
    l=40+i*590;panel(d,r['positions_m'][j],(l,165,530,480),25,action)
    text(d,(l,700),f"ΔCOM X: {(r['centroid_m'][j,0]-r['centroid_m'][0,0])*1e6:+.5f} µm",21);text(d,(l,735),f"Rg: {r['rg_m'][j]*1e6:.4f} µm | NN: {r['nn_m'][j]*1e6:.4f} µm",19);text(d,(l,768),f"3D proximity components: {r['components'][j]} | 32 bubbles",19)
    meta=all_runs[names[i]];text(d,(l,795),'Step gates: '+meta['status']+' | dt/cache: '+('PASS' if meta.get('refinement_admitted') else 'NOT ADMITTED'),16,'#a6b7c7')
    plot(d,r['time_s'],(r['rg_m']-r['rg_m'][0])*1e9,j,(l,842,500,125),'Actual ΔRg (nm)')
   text(d,(35,1025),'HOMOGENEOUS_WATER_CONTROL_PROXY | Numerically computed states; physical calibration remains unverified',18,'#8194a7');yield im
 records.append(encode(videos/'Video_F10Z_1_dispersed_bubble_dynamics.mp4',first(),fps,mapping))
 ds2=[load('original__PAIR__2000us__DT1'),ds[2]];names2=['original__PAIR__2000us__DT1',names[2]];mapping2=[dict(frame=k,time_s=float(ds2[0]['time_s'][j]),sources=[dict(run=n,saved_index=int(j)) for n in names2]) for k,j in enumerate(indices)]
 def second():
  for j in indices:
   im,d=base('F10Z | Translation versus deformation','X pair: 1036.20 Pa-equivalent each, phase difference 210° | 2.5 MHz',ds2[0]['time_s'][j],slow)
   for i,(r,label,limit) in enumerate(zip(ds2,['Prescribed near-contact group','Initially dispersed group'],[9,25])):
    l=40+i*900;panel(d,r['positions_m'][j],(l,175,380,380),limit,label)
    panel(d,r['positions_m'][j],(l+430,175,380,380),limit,'Centroid removed',centered=True)
    text(d,(l,625),f"32 true bubbles | {r['components'][j]} proximity components",22)
    meta=all_runs[names2[i]];text(d,(l,651),'Step gates: '+meta['status']+' | dt/cache: '+('PASS' if meta.get('refinement_admitted') else 'NOT ADMITTED'),16,'#a6b7c7')
    plot(d,r['time_s'],(r['centroid_m'][:,0]-r['centroid_m'][0,0])*1e6,j,(l,700,365,210),'COM X displacement (µm)')
    plot(d,r['time_s'],(r['rg_m']-r['rg_m'][0])*1e9,j,(l+430,700,365,210),'Internal ΔRg (nm)')
   text(d,(35,1000),'Actual radius and positions; XY projections of different Z layers can overlap. No deformation magnification.',19,'#a6b7c7');text(d,(35,1030),'Geometric proximity is not fusion. Prescribed initial contact is not an acoustically formed aggregate.',17,'#8194a7');yield im
 records.append(encode(videos/'Video_F10Z_2_translation_vs_deformation.mp4',second(),fps,mapping2))
 runs=read(BASE/'analysis/RUN_ANALYSIS.json');candidates=[r for r in runs.values() if r['initial']=='medium_train' and r['action']=='PAIR' and r['name'].endswith('DT1') and r.get('refinement_admitted')];r=max(candidates,key=lambda a:a['completed_s']);d3=load(r['name']);fs=d3['force_samples'];force_indices=np.arange(len(fs['time_s']));indices3=[int(np.argmin(abs(d3['time_s']-t))) for t in fs['time_s']];fps3=10;slow3=(1/fps3)/(fs['time_s'][1]-fs['time_s'][0]);mapping3=[dict(frame=k,time_s=float(fs['time_s'][k]),run=r['name'],saved_index=int(j),force_sample_index=k) for k,j in enumerate(indices3)]
 def third():
  for k,j in enumerate(indices3):
   im,d=base('F10Z | Forces and internal convergence','X pair | Saved force and position samples at the same physical time | Arrows normalized within each panel',d3['time_s'][j],slow3)
   for i,(field,title) in enumerate([('primary_force_N','External primary force'),('secondary_force_N','Internal F2 pair force'),('total_force_N','Total acoustic force')]):
    l=40+i*590;f=fs[field][k];panel(d,d3['positions_m'][j],(l,180,530,470),25,title,force=f);text(d,(l,720),f'Mean |F| = {np.linalg.norm(f,axis=1).mean()*1e12:.4g} pN',23);text(d,(l,755),f'Net |ΣF| = {np.linalg.norm(f.sum(0))*1e12:.4g} pN',23)
   radial=d3['radial_mean_m_s'][min(j,len(d3['radial_mean_m_s'])-1)]*1e6;quant=d3['pair_quantiles_m'][j]*1e6
   text(d,(35,790),f'Relative radial speed: {radial:+.3f} µm/s | Pair-distance P10 / P50 / P90: {quant[0]:.2f} / {quant[1]:.2f} / {quant[2]:.2f} µm | Step and dt/cache gates: PASS',18)
   plot(d,d3['time_s'],(d3['rg_m']-d3['rg_m'][0])*1e9,j,(40,850,490,120),'ΔRg (nm)');plot(d,d3['time_s'],d3['nn_m']*1e6,j,(640,850,490,120),'Median 3D NN (µm)');plot(d,d3['time_s'],d3['components'],j,(1240,850,490,120),'Proximity components')
   text(d,(35,1025),'Model trajectory only | Internal force magnitude is not net COM propulsion | No unsimulated time is displayed',18,'#8194a7');yield im
 records.append(encode(videos/'Video_F10Z_3_force_and_aggregation_trends.mp4',third(),fps3,mapping3));write(BASE/'evidence/VIDEO_VALIDATION.json',records)
 print([(v['path'],v['frame_count'],v['duration_s']) for v in records])
if __name__=='__main__':main()
