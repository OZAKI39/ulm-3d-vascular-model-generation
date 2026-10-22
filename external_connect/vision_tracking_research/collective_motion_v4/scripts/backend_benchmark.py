"""First-stage same-pair backend experiment; no channel system is needed."""
import sys,time,argparse,traceback,json,importlib.metadata as md
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from backends import Backend
import cv2
def main():
 p=argparse.ArgumentParser();p.add_argument('--names',default='B0,B1,B2,B3,B4');p.add_argument('--repeats',type=int,default=12);a=p.parse_args();cv2.setNumThreads(2)
 assert sha(VIDEO)==VIDEO_SHA
 ids=[90,92,360,362,700,702];images={i:(t,g) for i,t,g in frames(ids)}
 rois={'historical':[650,300,256,384],'medium':[600,240,600,624],'full_visible_bbox':[480,0,840,1120]}
 mon=ResourceMonitor();rows=[];failures=[];fields=OUT/'results/backend_fields';fields.mkdir(exist_ok=True)
 for name in a.names.split(','):
  try:
   backend=Backend(name)
   for region,(xx,yy,w,h) in rois.items():
    for first in [90,360,700]:
     ta,aa=images[first];tb,bb=images[first+2];aa=aa[yy:yy+h,xx:xx+w];bb=bb[yy:yy+h,xx:xx+w]
     cold=time.perf_counter();warm=backend.run(aa,bb,tb-ta);cold=(time.perf_counter()-cold)*1000
     for _ in range(2):backend.run(aa,bb,tb-ta)
     for rep in range(a.repeats):
      r=backend.run(aa,bb,tb-ta);ts=time.perf_counter();wire=json.dumps(clean(r['estimate']),allow_nan=False);output_ms=(time.perf_counter()-ts)*1000
      row={k:r[k] for k in ['core_ms','preprocess_ms','validation_ms','fit_ms','total_ms','dtype','upload_ms','download_ms','quality_scope']}
      row.update(name=name,region=region,frame_a=first,frame_b=first+2,dt_s=tb-ta,rep=rep,n_windows=r['valid'].size,n_valid=int(r['valid'].sum()),estimate=r['estimate'],state_output_ms=output_ms,end_to_end_ms=r['total_ms']+output_ms,cold_ms=cold,roi=[xx,yy,w,h],ws=48,overlap=24)
      rows.append(row)
     np.savez_compressed(fields/f'{name}_{region}_{first}.npz',**{k:r[k] for k in ['x','y','u','v','valid','texture']},snr=r['snr'] if r['snr'] is not None else np.full_like(r['u'],np.nan))
     backend.profile=True;pr=backend.run(aa,bb,tb-ta);backend.profile=False
     dump(fields/f'{name}_{region}_{first}_profile.json',dict(components=pr['components_ms'],core_ms=pr['core_ms'],note='synchronized instrumentation; overhead excluded from main benchmark'))
     print(name,region,first,'p95',stats([x['end_to_end_ms'] for x in rows if x['name']==name and x['region']==region])['p95'],flush=True)
  except Exception as e:
   err=dict(name=name,error=repr(e),traceback=traceback.format_exc());failures.append(err);print(err['traceback'],flush=True)
 dump(OUT/'results/backend_benchmark.json',dict(identity=identity(),rows=rows,failures=failures,resources=mon.close(),dependencies={k:md.version(k) for k in ['numpy','scipy','openpiv','torch','cupy-cuda12x']},scope='Same preprocessed raw pairs; different algorithms not claimed equivalent. Historical coordinates retained for first-stage comparison. No decoder in per-pair timer.'))
if __name__=='__main__':main()
