"""Finite full-channel integration profile, selected/regular and update rates."""
import sys,time,json,traceback
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from observer import Observer
import cv2
def main():
 cv2.setNumThreads(2);mon=ResourceMonitor();items=list(frames(list(range(360,390))));rows=[];failures=[]
 for name,selected,period in [('B2',False,1),('B2',True,1),('B3',False,1),('B3',True,1),('B4',False,1),('B4',True,1),('B4',True,3),('B4',True,6)]:
  obs=Observer(name,selected=selected,schedule=period)
  for i,pts,g in items:
   t=time.perf_counter()
   try:
    state=obs.process(i,pts,g);wire=json.dumps(state,allow_nan=False,separators=(',',':'));elapsed=(time.perf_counter()-t)*1000
    rows.append(dict(backend=name,selected=selected,period=period,frame_id=i,elapsed_ms=elapsed,algorithm_ms=state['algorithm_ms'],scan_ms=state['scan_ms'],group_ms=state['group_ms'],components=state['component_times'],n_candidates=len(state['candidates_global_xy']),valid_windows=sum(state['windows']['valid']),computed_windows=sum(state['windows']['fresh']),selected_windows=sum(state['windows']['active']),n_groups=len(state['groups']),valid_groups=state['fresh_valid_group_count'],failure=state['backend_failure'],registration=state['registration']))
    if i==389:dump(OUT/f'results/pilot_{name}_{selected}_{period}.json',state)
   except Exception as e:failures.append(dict(backend=name,selected=selected,period=period,frame=i,error=repr(e),traceback=traceback.format_exc()))
  print(name,selected,period,'p95',stats([r['elapsed_ms'] for r in rows if r['backend']==name and r['selected']==selected and r['period']==period and r['frame_id']>=365]),flush=True)
 dump(OUT/'results/channel_pilot.json',dict(identity=identity(),rows=rows,failures=failures,resources=mon.close(),scope='30 contiguous frames; warmup excluded only from warmed summary, all records retained. No decode in timed profile.'))
if __name__=='__main__':main()
