"""Local read-only reuse of human reviewed ROIs; never treat partial labels as background."""
import sys,json,time,shutil
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from channel import Channel
import cv2
from scipy.spatial import cKDTree
ANNOT=Path('/home/lzy/projects/external_connect/vision_tracking_research/smart_online_tracking_v2')
sys.path.insert(0,str(ANNOT/'src'))
from annotation_regions import reviewed_boxes
def main():
 cv2.setNumThreads(2);start=time.perf_counter();c=Channel();annotations={};manifest=[]
 for pointer in (ANNOT/'data/annotations').glob('video_*/latest.json'):
  meta=json.loads(pointer.read_text());mfile=pointer.parent/meta['mask_file'];assert sha(mfile)==meta['mask_sha256']
  masks=dict(np.load(mfile));boxes=reviewed_boxes(meta,masks)
  if not boxes:continue
  fid=int(pointer.parent.name.split('_')[-1]);annotations[fid]=(meta,masks,boxes)
  target=OUT/'evidence/human_annotation_snapshot'/pointer.parent.name;target.mkdir(parents=True,exist_ok=True);shutil.copy2(pointer,target/'latest.json');shutil.copy2(mfile,target/mfile.name)
  manifest.append(dict(frame_id=fid,meta_sha256=sha(pointer),mask_sha256=sha(mfile),review_scope=meta['review_scope'],boxes=boxes,source_sha256=meta['source_sha256']))
 rows=[];scans=[];first=None
 for fid,pts,g in frames():
  if first is None:first=g.copy()
  s=c.scan(g);scans.append(dict(frame_id=fid,pts=pts,scan_ms=s['scan_ms'],candidates=len(s['candidates']),selected=int(s['selected'].sum()),registration=s['registration']))
  if fid in [0,450,936]:
   canvas=cv2.cvtColor(g,cv2.COLOR_GRAY2BGR);cv2.polylines(canvas,[np.asarray(c.config['channel_polygon'],np.int32)],True,(20,220,240),3)
   for p in c.config['excluded_polygons']:cv2.polylines(canvas,[np.asarray(p,np.int32)],True,(80,80,240),2)
   cv2.rectangle(canvas,(650,300),(906,684),(255,170,40),3)
   cv2.imwrite(str(OUT/f'evidence/channel_overlay_{fid:04d}.png'),canvas)
  if fid not in annotations:continue
  meta,masks,boxes=annotations[fid];review=np.zeros(g.shape,bool)
  for y0,x0,y1,x1 in boxes:review[y0:y1,x0:x1]=True
  truths=[]
  labels=masks['SINGLE_BUBBLE']
  for label in np.unique(labels[review]):
   if not label:continue
   yy,xx=np.nonzero((labels==label)&review);truths.append([float(xx.mean()),float(yy.mean())])
  truths=np.asarray(truths);q=s['candidates']+np.array(c.crop[:2]);distance=cKDTree(q).query(truths)[0] if len(q) else np.full(len(truths),np.inf)
  safe=[];active=[]
  for x,y in truths:
   lx,ly=x-c.crop[0],y-c.crop[1];inside=(np.abs(c.x-lx)<=c.ws/2)&(np.abs(c.y-ly)<=c.ws/2)
   safe.append(bool(np.any(inside&c.safe)));active.append(bool(np.any(inside&s['selected'])))
  included=c.valid[np.rint(truths[:,1]).astype(int),np.rint(truths[:,0]).astype(int)].astype(bool)
  rows.append(dict(frame_id=fid,reviewed_boxes=boxes,n_human_bubbles=len(truths),human_centers=truths,nearest_candidate_distance=distance,detected_within_8px=distance<=8,inside_observation_mask=included,has_safe_window=safe,covered_by_selected_window=active,scope='Previously exposed human partial ROI labels; diagnostic candidate recall only, no velocity truth or cross-recording holdout'))
 drift=[]
 for shift in [0,3,9,20]:
  ch=Channel();ch.scan(first);moved=cv2.warpAffine(first,np.float32([[1,0,shift],[0,1,0]]),(first.shape[1],first.shape[0]),borderMode=cv2.BORDER_REFLECT);rr=ch.scan(moved)['registration'];drift.append(dict(applied_x_px=shift,expected_valid=shift<=6,observed=rr))
 coverages={str(f):int((c.coverage>=f).sum()) for f in [.5,.9,1.]}
 cv2.imwrite(str(OUT/'evidence/channel_mask.png'),c.valid*255)
 dump(OUT/'results/channel_audit.json',dict(identity=identity(),wall_s=time.perf_counter()-start,source_sha256=sha(VIDEO),annotation_manifest=manifest,reviewed_rows=rows,registration_tests=drift,full_scan_rows=scans,geometry=dict(channel_pixels=int(c.channel.sum()),observation_pixels=int(c.valid.sum()),historical_roi_pixels=256*384,excluded_pixels=int((c.channel*c.excluded).sum()),grid_shape=c.x.shape,window_coverage_rule_counts=coverages,physical_boundary_truth=False),limitations=['human ROIs do not cover whole channel','central-ridge exclusion may hide real bubbles; count reported','PZT state not known from MP4; cannot identify a real PZT-OFF segment']))
 print('Completed scans',len(scans),'human labelled frames',len(rows),'geometry',coverages)
if __name__=='__main__':main()
