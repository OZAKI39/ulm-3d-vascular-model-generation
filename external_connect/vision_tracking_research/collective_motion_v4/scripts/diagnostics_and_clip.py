"""Bounded post-freeze diagnostics; never retune frozen scientific configuration."""
import sys,time,json,gzip
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from backends import Backend
from channel import Channel
from studies import score
import cv2
from scipy.spatial import cKDTree

def main():
 cv2.setNumThreads(2);start=time.perf_counter();mon=ResourceMonitor();engines={n:Backend(n) for n in ['B0','B1','B2','B3','B4']};folder=OUT/'results';folder.mkdir(exist_ok=True)
 # A contiguous real clip, independently executed for all five backends.
 images={i:(t,g) for i,t,g in frames(range(358,420))};clip=[];x,y,w,h=650,300,256,384
 for n,e in engines.items():
  for _ in range(3):e.run(images[358][1][y:y+h,x:x+w],images[360][1][y:y+h,x:x+w],.0666,mode='local',actual_coordinates=True)
  for fid in range(360,420):
   pt,g=images[fid];prev=images[fid-2];r=e.run(prev[1][y:y+h,x:x+w],g[y:y+h,x:x+w],pt-prev[0],mode='local',actual_coordinates=True)
   clip.append(dict(backend=n,frame_id=fid,source_pts=pt,dt_s=pt-prev[0],roi=[x,y,w,h],**{k:r[k] for k in ['x','y','u','v','valid','snr','estimate','total_ms','core_ms','validation_ms']}))
  print('clip',n,flush=True)
 dump(folder/'backend_clip.json',dict(identity=identity(),rows=clip,scope='SELECTED_CLIP frames 360..419, historical ROI, local enhancement, gap 2; no physical truth'))
 # Actual data-transfer sizes, synchronized standalone probe (not subtracted from total).
 cp=engines['B4'].cp if hasattr(engines['B4'],'cp') else __import__('cupy');arr=np.zeros((2,1120,912),np.float32);transfer=[]
 for i in range(25):
  t=time.perf_counter();z=cp.asarray(arr);cp.cuda.Stream.null.synchronize();up=(time.perf_counter()-t)*1000;t=time.perf_counter();b=z.get();cp.cuda.Stream.null.synchronize();down=(time.perf_counter()-t)*1000
  if i>=5:transfer.append(dict(upload_ms=up,download_ms=down))
 # Truth-based selector probe embeds the existing generator at a new spatial location.
 index=json.loads((SERVER/'runs/research_v1/data/development/index.json').read_text());base=images[360][1].copy();selection=[];steps=[]
 for meta in index:
  data=np.load(SERVER/'runs/research_v1'/meta['file']);f=data['frames'];c=Channel();offset=np.array([1020,700]);c.scan(base)
  if meta['case'] in ['sparse','low_contrast','static','dropout','fixed_background','translation_x','boundary']:
   for fid in [0,1,15]:
    full=base.copy();full[700:956,1020:1276]=f[fid];s=c.scan(full);ids=meta['meta']['visible_ids_generator_only'][fid];truth=data['points'][fid][ids]+offset;q=s['candidates']+np.array(c.crop[:2]);dist=cKDTree(q).query(truth)[0] if len(q) else np.full(len(truth),np.inf);coverage=[]
    for xx,yy in truth:
     inside=(np.abs(c.x+c.crop[0]-xx)<=24)&(np.abs(c.y+c.crop[1]-yy)<=24);coverage.append(bool(np.any(inside&s['selected'])))
    selection.append(dict(case=meta['case'],fid=fid,visible=len(truth),detected_within8=int((dist<=8).sum()),covered_by_active=int(sum(coverage)),new_location=[1020,700],scan_ms=s['scan_ms'],registration=s['registration'],scope='Known rendered visible particles, synthetic patch insertion; no PZT state inferred'))
  if meta['case']=='step':
   backgrounds=[];bg=f[0].astype(np.float32)
   for a in f:backgrounds.append(bg.copy());cv2.accumulateWeighted(a,bg,.08)
   for end in [7,8,9,10]:
    for gap in [1,2,3]:
     beg=end-gap;r=engines['B4'].run(f[beg],f[end],float(data['times'][end]-data['times'][beg]),mode='causal',background_a=backgrounds[beg],background_b=backgrounds[end],actual_coordinates=True)
     steps.append(dict(end=end,gap=gap,**score(r,data,meta,beg,end)))
 # Pixel-window safety diagnostic: weaker rules retain wall-crossing windows.
 c=Channel();g0=images[360][1];g1=images[362][1];r=engines['B4'].run(c.cut(g0),c.cut(g1),.0666,mode='local',do_fit=False,actual_coordinates=True);wall=[]
 for fraction in [.5,.9,1.]:
  t=time.perf_counter();safe=c.coverage>=fraction;valid=safe&r['valid'];wall.append(dict(fraction=fraction,n_allowed=int(safe.sum()),n_wall_crossing_allowed=int((safe&~c.safe).sum()),n_valid=int(valid.sum()),median_apparent_speed_px_s=float(np.median(np.hypot(r['u'][valid],r['v'][valid])/.0666)),rule_ms=(time.perf_counter()-t)*1000,physical_wall_velocity_error=None))
 # Pure PDMS/background diagnostic has no bubble truth; do not label candidates as bubbles.
 c=Channel();s=c.scan(base);excluded_before=len(s['excluded_candidates']);static=c.scan(base)
 dump(folder/'post_freeze_diagnostics.json',dict(identity=identity(),selection=selection,step=steps,wall_rules=wall,static_background=dict(excluded_candidates_before_geometry=excluded_before,selected_windows=int(static['selected'].sum()),mean_activity=float(static['activity'].mean()),note='Static candidates may remain active; zero motion does not mean no bubble. PZT-OFF real segment NOT_IDENTIFIED.'),transfer=dict(bytes=arr.nbytes,rows=transfer,upload=stats([r['upload_ms'] for r in transfer]),download=stats([r['download_ms'] for r in transfer]),scope='Synchronized separate transfer microbenchmark, not total solver transfer attribution'),resources=mon.close(),wall_s=time.perf_counter()-start,configuration_changed=False))
 print('diagnostics complete',flush=True)
if __name__=='__main__':main()
