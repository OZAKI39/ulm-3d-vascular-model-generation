"""Exact saved frames, fixed camera/true spheres, uniform playback, full decode."""
import time,numpy as np,h5py,av
from PIL import Image,ImageDraw
from scripts.straight_tube_exp1.f10aa.videos import scene,font
from scripts.straight_tube_exp1.f10ab.metrics import measure
from .common import *

def main(only_name=None):
 out=BASE/'videos';out.mkdir(exist_ok=True);selected=read(BASE/'ELITE_VALIDATION_PREREGISTRATION.json')['selected'];comparison=read(BASE/'COMPARISON_V2.json');best=[]
 allrows={r['id']:r for r in read(BASE/'ALL_RUNS.json')};elites=read(BASE/'ELITE_MANIFEST.json')
 for method in comparison['methods']:
  matches=[e for e in elites if e['archive']==f'search_v2/{method}/RESULT.json' and e['task_type']=='archive']
  pick=max(matches,key=lambda e:abs(e['measures'][1])) if matches else None
  best.append((method,allrows[pick['run_id']] if pick else None))
 long=read(AD/'LONG.json');checks=read(AD/'LONG_CONVERGENCE.json');freq=[]
 for f in [2400000.,2500000.,2600000.]:
  r=next((r for r in long if r['task']['frequency']==f and r['id'].endswith('Y90') and r['status']=='PASS'),None)
  cc=[c for c in checks if r and c['reference']==r['id']]
  freq.append((f'{f/1e6:g} MHz',r if len(cc)==2 and all(c['admitted'] for c in cc) else None))
 discovery=[(f'Elite {i}',r) for i,r in enumerate(selected)]
 if (BASE/'OUTSIDE_PREREGISTRATION.json').exists():discovery[1]=('Outside archive; seed checks FAIL',read(BASE/'OUTSIDE_PREREGISTRATION.json')['chosen'])
 jobs=[('Video_F10AC_1_shape_diversity_discovery',discovery,'Archived and retained outside-archive responses | fixed IDs, no remapping'),('Video_F10AC_2_algorithm_best_elites',best,'Archived elite with largest |delta log aspect| per method | same initial state and window'),('Video_F10AD_3_frequency_dependent_dynamics',freq,'Y pair: 700 Pa-eq each, 90 deg phase difference | admitted fixed frequencies')]
 manifest=[r for r in read(out/'MANIFEST.json') if r['name']!=only_name] if only_name and (out/'MANIFEST.json').exists() else []
 for name,panels,title in jobs:
  if only_name and name!=only_name:continue
  data=[]
  for label,r in panels:
   if r:
    path=BASE/r['trajectory_path']
    if not path.exists():path=BASE/'dataset'/(r['id']+'.h5')
    with h5py.File(path) as h:
     if 'raw_sha256' in h.attrs:assert h.attrs['raw_sha256']==r['trajectory_sha256']
     t=h['time_s'][:];index=h['raw_frame_index'][:] if 'raw_frame_index' in h else np.arange(len(t));data.append((t,h['positions_m'][:],index))
   else:data.append(None)
  end=min([v[0][-1] for v in data if v is not None],default=.002);times=np.arange(0,end+1e-13,1e-5);path=out/(name+'.mp4');writer=av.open(str(path),'w');stream=writer.add_stream('libx264',rate=30);stream.width=1800;stream.height=880;stream.pix_fmt='yuv420p';stream.options={'crf':'19','preset':'medium'};mapping=[];begin=time.perf_counter()
  for iframe,tt in enumerate(times):
   im=Image.new('RGB',(1800,880),'#080d16');d=ImageDraw.Draw(im);d.text((24,15),title,font=font(24),fill='#ecf2fa');d.text((24,51),f'Physical time {tt*1e3:.3f} ms | Equilibrium R = 1.78 um | zero background flow | slow playback',font=font(21),fill='#bdcddd');sources=[]
   for k,((label,r),arr) in enumerate(zip(panels,data)):
    xbase=24+k*600;d.text((xbase,95),label,font=font(22),fill='#91c8f3')
    if arr is None:
     d.text((xbase,325),'DYNAMIC_NOT_ADMITTED' if name.startswith('Video_F10AD') else 'NO_VALID_ELITE',font=font(21),fill='#e2b177')
     if name.startswith('Video_F10AD'):
      f=int(float(label.split()[0])*1e6);fd=np.load(AD/str(f)/'FIELD_DIAGNOSTICS.npz');from matplotlib import colormaps
      p=abs(fd['xy_pressure'][0])*700;img=(colormaps['magma'](p/max(p.max(),1e-10))[...,:3]*255).astype('uint8');im.paste(Image.fromarray(img).resize((480,350)),(xbase,390))
     continue
    t,x,raw_index=arr;ix=int(np.searchsorted(t,tt-1e-13));assert abs(t[ix]-tt)<1e-13;m=measure(x[0],x[ix]);im.paste(scene(x[ix],size=(560,560)),(20+600*k,126));d=ImageDraw.Draw(im)
    lines=[f'f = {r["task"].get("frequency",2500000.)/1e6:g} MHz | N = {len(x[0])}',f'Rg = {m["Rg_3d_m"]*1e6:.3f} um | 3D aspect = {m["aspect_3d"]:.3f}',f'Geometric components: {m["proximity_components"]} | largest: {m["largest_component"]}',f'Nonrigid RMS: {m["similarity_rmsd_3d_m"]*1e9:.1f} nm','A [Y-/Y+/X-/X+]: '+' / '.join(f'{v:.0f}' for v in r['task']['amplitudes'])+' Pa-eq','Phase: '+' / '.join(f'{np.degrees(v):.1f}' for v in r['task']['phases'])+' deg']
    for j,line in enumerate(lines):d.text((xbase,691+26*j),line,font=font(16),fill='#d1dfec')
    sources.append(dict(id=r['id'],trajectory=r['trajectory_path'],sha256=r['trajectory_sha256'],saved_index=int(raw_index[ix])))
   d.text((24,853),'EXPLORATORY | Uncalibrated homogeneous water | Actual sphere size and separation | Geometric components do not mean fusion',font=font(16),fill='#aab7c7')
   for packet in stream.encode(av.VideoFrame.from_image(im)):writer.mux(packet)
   if iframe in [0,len(times)//2,len(times)-1]:im.save(out/f'{name}_{iframe:04d}.png')
   mapping.append(dict(frame=iframe,physical_time_s=float(tt),sources=sources))
  for packet in stream.encode():writer.mux(packet)
  writer.close();pts=[]
  with av.open(str(path)) as dec:
   for frame in dec.decode(video=0):pts.append(float(frame.pts*frame.time_base))
  assert len(pts)==len(times) and np.allclose(pts,np.arange(len(times))/30,atol=1e-6)
  record=dict(name=name,sha256=sha(path),fps=30,frames=len(pts),physical_duration_s=float(end),playback_s=len(pts)/30,fully_decoded=True,resolution=[1800,880],wall_s=time.perf_counter()-begin,mapping=mapping,identity=identity());write(out/(name+'.json'),record);manifest.append({k:v for k,v in record.items() if k!='mapping'})
 write(out/'MANIFEST.json',manifest)
if __name__=='__main__':
 import sys
 main(sys.argv[1] if len(sys.argv)>1 else None)
