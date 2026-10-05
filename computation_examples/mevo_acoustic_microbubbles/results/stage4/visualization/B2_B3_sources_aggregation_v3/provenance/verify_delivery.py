from pathlib import Path
import csv,hashlib,json,subprocess
from html.parser import HTMLParser
import imageio.v2 as io
import imageio_ffmpeg
import numpy as np
from PIL import Image
root=Path('/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles')
out=root/'results/stage4/visualization/B2_B3_sources_aggregation_v3'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
m=json.loads((out/'RENDER_MANIFEST.json').read_text())
assert m['status']=='PASS' and sha(root/'scripts/stage4/render_aggregation.py')==m['renderer_sha256']
assert sha(out/'provenance/stage4_render_aggregation.py')==m['renderer_sha256']
assert all(sha(p)==s for p,s in m['protected_inputs'].items())
holds=np.array(m['microscope']['video_frames_per_state'])
starts=np.r_[0,np.cumsum(holds)[:-1]];ends=np.cumsum(holds)-1
n=int(holds.sum());assert n==1812
readers={};checked={}
for name,info in m['videos'].items():
 path=out/info['path'];assert sha(path)==info['sha256']
 r=io.get_reader(path);readers[name]=r;meta=r.get_meta_data()
 assert r.count_frames()==n and meta['fps']==60 and meta['size']==tuple(info['resolution'])
 assert meta['codec']=='h264'
 decoded=subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-threads','2','-i',str(path),'-f','null','-'],capture_output=True,text=True,check=True)
 assert not decoded.stderr,decoded.stderr
 checked[name]=dict(full_decode='PASS',frames=n,fps=60,resolution=info['resolution'],duration_s=n/60)
 print(name,'decode PASS',flush=True)
for name in ['B2','B3']:
 r=readers[name]
 raw=[io.imread(out/'keyframes'/f'{name}_{k:02d}.png') for k in range(len(holds))]
 timers=np.stack([a[20:80,1640:1880].astype(np.float32) for a in raw])
 errors=[]
 for k in range(len(holds)):
  for frame in [int(starts[k]),int(ends[k])]:
   im=r.get_data(frame)
   score=np.abs(timers-im[20:80,1640:1880]).mean(axis=(1,2,3))
   assert int(score.argmin())==k,(name,frame,k,score.tolist())
   mae=float(np.abs(im.astype(float)-raw[k]).mean());assert mae<2.,mae
   errors.append(dict(frame=frame,state=k,mae=mae))
  if not m['states'][name][k]['anchor_active']:
   assert m['states'][name][k]['visible']==0 and m['states'][name][k]['visible_particle_ids']==[]
   crop=raw[k][325:695,1270:1710].astype(int)
   assert not np.any((crop[:,:,2]>120)&(crop[:,:,2]-crop[:,:,0]>50)), 'Ghost bubble after crossing'
 checked[name]['snapshot_boundaries']=errors
 Image.fromarray(r.get_data(int(starts[24 if name=='B2' else 17]))).resize((1536,864),Image.Resampling.LANCZOS).save(f'/tmp/mevo_aggregation_full_encoded_{name}.png')
for frame in [0,120,1320,n-1]:
 actual=readers['B2_B3'].get_data(frame)
 expected=np.concatenate([readers[c].get_data(frame) for c in ['B2','B3']],axis=1)
 assert np.abs(actual.astype(float)-expected).mean()<2.
for r in readers.values():r.close()
rows=list(csv.DictReader((out/'frame_map.csv').open()));assert len(rows)==2*n
state_index=np.repeat(np.arange(len(holds)),holds)
for row in rows:
 k=int(state_index[int(row['frame'])]);s=m['states'][row['case']][k]
 assert int(row['saved_state'])==k and float(row['physical_time_s'])==s['t']
 assert int(row['original_saved_state'])==s['source_state']
 assert (row['anchor_active']=='True')==s['anchor_active']
 for key,col in [('anchor_cluster_size','local_group'),('active_particles','active_particles'),('A','A'),('B','B')]:assert int(row[col])==s[key]
 assert s['active_particles']+s['A']+s['B']==256
 assert int(row['all_cluster_count'])==s['stats']['cluster_count']
 assert float(row['near_contact_fraction'])==s['stats']['near_contact_fraction']
 assert float(row['median_nearest_neighbor_m'])==s['stats']['median_nearest_neighbor_m']
assert [m['states']['B2'][-1][k] for k in ['A','B','active_particles']]==[76,0,180]
assert [m['states']['B3'][-1][k] for k in ['A','B','active_particles']]==[0,88,168]
assert m['states']['B2'][-1]['anchor_cluster_size']==75
assert m['states']['B3'][17]['anchor_cluster_size']==22
class Links(HTMLParser):
 def handle_starttag(self,tag,attrs):
  for key,val in attrs:
   if key in ['src','href','poster']:assert (out/val).is_file(),val
Links().feed((out/'OPEN_RESULTS.html').read_text())
report=dict(status='PASS',videos=checked,source_hashes_unchanged=True,frame_map_rows=2*n,
 original_metric_states_checked=50,code_commit='f481c86',radius_um=1.78,field_color_scale_kPa=[0,39],
 original_saved_state_indices=m['microscope']['original_saved_indices'],scope='Display/data identity and temporal checks, not new physical validation',
 anchor_exit_handling='B3 microscope empties after recorded crossing, no frozen inactive particle rendered as live')
(out/'DELIVERY_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
files={str(p.relative_to(out)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file() and p.name not in ['DELIVERY_SHA256.json','render.log']}
(out/'DELIVERY_SHA256.json').write_text(json.dumps(files,indent=2)+'\n')
print(json.dumps(dict(status='PASS',files=len(files),bytes=sum(v['bytes'] for v in files.values()))))
