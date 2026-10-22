from pathlib import Path
import json,sys,collections,csv
import numpy as np,cv2
from scipy.optimize import linear_sum_assignment
ROOT=Path(__file__).resolve().parents[1];sys.dont_write_bytecode=True
from replay_benchmark import Detector,CFG,diststats,dump
summaries=[json.loads(p.read_text()) for p in sorted((ROOT/'results').glob('*/summary.json')) if 'cotracker_online' not in str(p)]
agg=[]
for roi in ['wide','target_roi']:
 for method in ['cc_norfair','blob_norfair','cc_norfair_lk','cc_norfair_farneback']:
  tags=[f'{method}_{s}'+('_target_roi' if roi=='target_roi' else '') for s in [120,360,720]];ss=[s for s in summaries if s['tag'] in tags];frames=[];mock=[]
  for tag in tags:
   frames+=json.loads((ROOT/'results'/tag/'frames.json').read_text());mock+=json.loads((ROOT/'results'/tag/'mock_received.json').read_text())
  wall=sum(s['wall_s'] for s in ss); n=sum(s['decoded_frames'] for s in ss)
  a={'roi':roi,'method':method,'n_source':n,'n_processed':len(frames),'processed_hz':len(frames)/wall,'valid_visual_hz':sum(m['valid_visual'] for m in mock)/wall,'drop_fraction':1-len(frames)/n,'wall_s':wall,'mean_candidates':float(np.mean([f['n_detections'] for f in frames])),'max_candidates':max(f['n_detections'] for f in frames),'cpu_pct_weighted':sum(s['cpu_one_core_pct']*s['wall_s'] for s in ss)/wall,'max_rss_mb':max(s['max_rss_mb'] for s in ss),'stages':{k:diststats([f[k] for f in frames]) for k in ['decode_ms','detect_ms','track_ms','shape_quality_ms','flow_ms','serialize_send_ms','processing_ms']},'receive_to_mock_ms':diststats([m['age_ms'] for m in mock]),'transfer_ms':diststats([m['transfer_ms'] for m in mock]),'mean_receive_queue_ms':float(np.mean([f['frame_age_at_start_ms'] for f in frames]))};agg.append(a)
dump(ROOT/'results/comparison_aggregate.json',agg)
annotations=json.loads((ROOT/'evidence/visual_annotations.json').read_text());CFG['roi_xywh']=[650,350,180,350];acc=[]
for kind in ['cc','blob']:
 det=Detector(kind)
 for fid,gt in annotations['positive_only_keyframes'].items():
  image=cv2.imread(str(ROOT/f'evidence/frame_{int(fid):05d}.png'));objects,*_=det(image);pts=np.array([d['xy'] for d in objects]);truth=np.array(gt);matches=[]
  if len(pts):
   distances=np.linalg.norm(truth[:,None]-pts[None,:],axis=2);rr,cc=linear_sum_assignment(distances)
   for i,j in zip(rr,cc):
    if distances[i,j]<=6:matches.append({'annotation_index':int(i),'det_index':int(j),'error_px':float(distances[i,j])})
  acc.append({'detector':kind,'frame_id':int(fid),'positive_landmarks':len(gt),'matched':len(matches),'missed_selected_landmarks':len(gt)-len(matches),'matches':matches,'n_candidate_detections':len(pts),'positions':pts.tolist(),'ground_truth_type':annotations['annotation_type'],'false_positive_rate':'NOT_EVALUATED_POSITIVE_ONLY_LABELS','biological_accuracy':'NOT_EVALUATED'})
trackacc=[]
for method in ['cc_norfair','blob_norfair','cc_norfair_lk']:
 msgs=[json.loads(s) for s in (ROOT/f'results/{method}_120_target_roi/states.jsonl').read_text().splitlines()];byid={m['frame_id']:m for m in msgs}
 for label,gt in annotations['short_tracks'].items():
  points=[]
  for fid,xy in gt.items():
   fid=int(fid);m=byid.get(fid);cands=[o for o in m['targets'] if o['status']=='OBSERVED'] if m else []
   if cands:
    d=np.array([np.linalg.norm(np.array(o['xy_px'])-xy) for o in cands]);j=int(d.argmin());o=cands[j]
    points.append({'frame_id':fid,'gt_xy':xy,'matched':bool(d[j]<=6),'id':o['id'] if d[j]<=6 else None,'xy':o['xy_px'] if d[j]<=6 else None,'error_px':float(d[j]) if d[j]<=6 else None,'visual_valid':o['visual_valid'] if d[j]<=6 else False})
   else:points.append({'frame_id':fid,'gt_xy':xy,'matched':False,'id':None,'xy':None,'error_px':None,'reason':'not initialized/missed'})
  valid=[p for p in points if p['matched']]; ids=[p['id'] for p in valid];velocity_errors=[]
  for a,b in zip(valid,valid[1:]):
   if a['id']==b['id']:
    dt=(b['frame_id']-a['frame_id'])*.0333
    velocity_errors.append(float(np.linalg.norm((np.array(b['xy'])-a['xy'])/dt-(np.array(b['gt_xy'])-a['gt_xy'])/dt)))
  trackacc.append({'method':method,'label':label,'points':points,'matched_observations':len(valid),'reviewed_id_changes':sum(a!=b for a,b in zip(ids,ids[1:])),'reviewed_span_s':(valid[-1]['frame_id']-valid[0]['frame_id'])*.0333 if len(valid)>1 else 0,'centroid_error_px':diststats([p['error_px'] for p in valid]),'velocity_error_px_per_replay_s':diststats(velocity_errors),'velocity_reference_uncertainty_note':'Coordinate annotation uncertainty ±3 px; finite-difference velocity error is descriptive, not expert accuracy.'})
accout={'keyframe_detection':acc,'short_track_review':trackacc,'NOT_EVALUATED':['biological bubble precision/recall','global ID switches','false physical identity reconnections','group segmentation shape error','expert crossing-event false positives/negatives','3D velocity','true liquid velocity','um-scale calibration'],'review_scope':'21 selected positive visual landmarks and 3 short visual tracks; assistant visual review, not expert labels'}
dump(ROOT/'results/accuracy_review.json',accout)
# Real observed temporal losses/reappearances in tracker output, not physical reidentification truth.
full=[json.loads(s) for s in (ROOT/'results/full_cc_target_roi/states.jsonl').read_text().splitlines()];tracks=collections.defaultdict(list)
for m in full:
 for t in m['targets']:tracks[t['id']].append({'frame_id':m['frame_id'],'pts_s':m['source_pts_s'],**t})
reappeared=[]
for tid,seq in tracks.items():
 states=[s['status'] for s in seq]
 for i in range(1,len(states)):
  if states[i-1]=='LOST' and states[i]=='OBSERVED':reappeared.append({'id':tid,'frame_id':seq[i]['frame_id'],'identity_reconfirmed':False})
dump(ROOT/'results/track_review.json',{'tracked_ids':len(tracks),'software_reappearance_events':reappeared,'longest_tracks':[{'id':tid,'n_observations':sum(p['status']=='OBSERVED' for p in seq),'span_s':seq[-1]['pts_s']-seq[0]['pts_s']} for tid,seq in sorted(tracks.items(),key=lambda p:len(p[1]),reverse=True)[:10]],'crossing_events':[{'frame_id':m['frame_id'],'id':t['id'],'xy':t['xy_px']} for m in full for t in m['targets'] if t['target_crossing']=='OBSERVED']})
print(json.dumps([{k:a[k] for k in ['roi','method','processed_hz','valid_visual_hz','drop_fraction']} for a in agg],indent=2));print('matches',[(a['detector'],a['frame_id'],a['matched'],a['positive_landmarks']) for a in acc]);print('software reappearances',len(reappeared))
