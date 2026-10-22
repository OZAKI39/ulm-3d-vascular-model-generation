"""Summarize actual artifacts, preserve nulls, and audit the exact-Jacobian change."""
import sys,json,gzip,time,collections
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from channel import Channel
from fit_adapter import fit_affine
import cv2
def records(p):
 with gzip.open(p,'rt') as f:
  for line in f:yield json.loads(line)
def resource(d):
 s=d['samples'];return dict(cpu_core_equivalents=d['cpu_core_equivalents'],rss_peak_mb=max(v['rss_mb'] for v in s),gpu_memory_peak_mb=max(v.get('gpu_device_memory_mb',0) or 0 for v in s),gpu_util_mean=float(np.mean([v.get('gpu_device_utilization_percent',0) or 0 for v in s])),gpu_util_max=max(v.get('gpu_device_utilization_percent',0) or 0 for v in s),scope=d['gpu_scope'])
def main():
 start=time.perf_counter();mainpath=ROOT/'results/full_video_verified_v5/results/offline_regular';c=Channel();compact=[];samples={};equiv=[];old=records(ROOT/'results/full_video_v1/offline_regular/records.jsonl.gz');changes=0;failures=[]
 for s,o in zip(records(mainpath/'records.jsonl.gz'),old,strict=True):
  assert s['frame_id']==o['frame_id'];w=s['windows'];good=np.array(w['valid']);assert w['valid']==o['windows']['valid'];assert w['vx']==o['windows']['vx'] and w['vy']==o['windows']['vy']
  if s['backend_failure']:failures.append([s['frame_id'],s['backend_failure']])
  assert len(s['groups'])==len(o['groups'])
  for g,h in zip(s['groups'],o['groups']):
   changes+=g['six']['valid']!=h['six']['valid'];aa=g['six'].get('raw_estimate');bb=h['six'].get('raw_estimate')
   if aa is not None and bb is not None:equiv.append(float(np.max(np.abs(np.array(aa)-bb))))
  footprint=np.zeros(c.local.shape,np.uint8)
  for x,y in zip(c.x0.ravel()[good],c.y0.ravel()[good]):footprint[y:y+48,x:x+48]=1
  selected=[g for g in s['groups'] if g['observation_key'].startswith('cell_1_1_')];g=max(selected,key=lambda z:z['visible_count']) if selected else None
  vals=[g['six']['values'][n] for n in NAMES] if g else [None]*6
  compact.append(dict(frame_id=s['frame_id'],pts=s['source_PTS'],valid_windows=int(good.sum()),valid_fraction=float(good.sum()/c.safe.sum()),supported_footprint_fraction=float((footprint*c.local).sum()/c.local.sum()),candidates=len(s['candidates_global_xy']),valid_groups=s['fresh_valid_group_count'],groups=len(s['groups']),cell_1_1_six=vals,cell_1_1_count=g['visible_count'] if g else 0,cell_1_1_key=g['observation_key'] if g else None,mask_valid=s['channel_mask_valid'],scan_ms=s['scan_ms'],group_ms=s['group_ms']))
  if s['frame_id'] in [0,3,100,450,936]:samples[str(s['frame_id'])]=s
 print('equivalence audit',len(compact),changes,max(equiv),flush=True)
 assert len(compact)==937 and changes==0 and max(equiv)<1e-4
 # Bootstrap on existing actual vectors only. Nominal intervals, not calibrated confidence.
 s=samples['450'];w=s['windows'];x=np.array(w['x']);y=np.array(w['y']);good=np.array(w['valid']);dt=s['time_window']['dt_s'];uu=np.array(w['vx'],float)*dt;vv=np.array(w['vy'],float)*dt;unc=[]
 for g in s['groups']:
  rx,ry,rw,rh=g['reference_roi_global'];inside=good&(x>=rx)&(x<=rx+rw)&(y>=ry)&(y<=ry+rh)
  e=fit_affine(np.c_[x[inside]-rx,y[inside]-ry],np.c_[uu[inside],vv[inside]],dt,(rh,rw),'POSTHOC_BOOTSTRAP',bootstrap=24)
  unc.append(dict(key=g['observation_key'],center=g['center_global_xy'],estimate=e,scope='Nominal spatial resampling; uncertainty calibration NOT_VERIFIED'))
 bench=[]
 for p in ['stage_b_v2','stage_b_checks']:
  d=json.loads((ROOT/f'results/{p}/backend_benchmark.json').read_text());bench+=d['rows']
 bt=[]
 for n in ['B0','B1','B2','B3','B4']:
  for reg in ['historical','medium','full_visible_bbox']:
   rr=[r for r in bench if r['name']==n and r['region']==reg]
   if not rr:continue
   bt.append(dict(backend=n,region=reg,windows=rr[0]['n_windows'],timing=stats([r['end_to_end_ms'] for r in rr]),components={k:stats([r[k] for r in rr if r.get(k) is not None]) for k in ['core_ms','preprocess_ms','validation_ms','fit_ms','upload_ms','state_output_ms']},scope='Pair + PIV + six + serialization; excludes decode/channel scan/multiple groups'))
 perf=[]
 paths=sorted((ROOT/'results/full_video_v1').glob('*/summary.json'))+[ROOT/'results/full_video_verified_v5/results/offline_regular/summary.json']+sorted((ROOT/'results/replay_final_v6/results').glob('*/summary.json'))
 for p in paths:
  d=json.loads(p.read_text());rows=d['rows'];perf.append(dict(name=('original_' if 'full_video_v1/' in str(p) else 'verified_')+p.parent.name,path=str(p.relative_to(ROOT)),p50=d['full_pipeline_ms']['p50'],p95=d['full_pipeline_ms']['p95'],max=d['full_pipeline_ms']['max'],frames=len(rows),mean_valid_windows=float(np.mean([r['valid_windows'] for r in rows])),mean_computed_windows=float(np.mean([r['computed_windows'] for r in rows])),mean_groups=float(np.mean([r['valid_groups'] for r in rows])),full_channel_scan_hz=d.get('full_channel_scan_hz'),scan_kernel_hz=1000/np.mean([r['scan_ms'] for r in rows]),fresh_group_frame_hz=d.get('fresh_valid_group_motion_frame_hz'),fresh_group_states_hz=d.get('fresh_valid_group_states_hz'),local_piv_frame_hz=d.get('local_piv_update_frame_hz'),dropped=len(d.get('dropped_frames',[])),resources=resource(d['resources']),components={k:stats([r[k] for r in rows]) for k in ['decode_ms','scan_ms','group_ms']},age=stats([r['age_ms'] for r in rows]) if 'age_ms' in rows[0] else None))
 accuracy=[]
 for split in ['development','validation','final_test']:
  d=json.loads((ROOT/f'results/research_v1/results/accuracy_{split}.json').read_text());keys=sorted(set((r['backend'],r['mode'],r['gap'],r['ws'],r['multipass']) for r in d['rows']))
  for key in keys:
   rr=[r for r in d['rows'] if (r['backend'],r['mode'],r['gap'],r['ws'],r['multipass'])==key];accepted=[r for r in rr if r['accepted']];errs=[r['six_absolute_error'] for r in accepted if r['six_absolute_error'] is not None];sign=[]
   for r in accepted:
    if r['six_truth'] is not None:
     mask=np.abs(r['six_truth'])>1e-6;sign.extend(np.array(r['six_sign_correct'])[mask].tolist())
   accuracy.append(dict(split=split,backend=key[0],mode=key[1],gap=key[2],ws=key[3],multipass=key[4],n=len(rr),accepted=len(accepted),false_accept=sum(r['false_accept'] for r in rr),correct_reject=sum(r['correct_reject'] for r in rr),false_cases=[r['case'] for r in rr if r['false_accept']],epe_valid_mean=float(np.mean([r['epe_valid'] for r in rr if r['epe_valid'] is not None])),six_mae_accepted=np.mean(errs,axis=0).tolist() if errs else None,nonzero_truth_sign_accuracy=float(np.mean(sign)) if sign else None,over_tolerance_among_accepted=sum(bool(r['over_tolerance']) for r in accepted),mean_valid_fraction=float(np.mean([r['valid_vectors']/r['n_vectors'] for r in rr])),correct_rejection_cases=[r['case'] for r in rr if r['correct_reject']],mean_total_ms=float(np.mean([r['total_ms'] for r in rr]))))
 dump(ROOT/'results/delivery_analysis.json',dict(identity=identity(),wall_s=time.perf_counter()-start,compact=compact,backend=bt,performance=perf,accuracy=accuracy,coverage_stats=stats([r['supported_footprint_fraction'] for r in compact]),coverage_mean=float(np.mean([r['supported_footprint_fraction'] for r in compact])),valid_window_mean=float(np.mean([r['valid_fraction'] for r in compact])),exact_jacobian_equivalence=dict(frames=937,windows_exact=True,group_gates_changed=changes,max_six_difference=max(equiv),checked_groups=len(equiv)),backend_failures=failures,uncertainty=unc))
 dump(ROOT/'results/figure_sample_states.json',samples);print('analysis done; coverage',np.mean([r['supported_footprint_fraction'] for r in compact]),'max six difference',max(equiv),flush=True)
if __name__=='__main__':main()
