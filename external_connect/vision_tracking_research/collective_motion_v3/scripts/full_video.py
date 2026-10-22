"""Every decoded frame, frozen OpenPIV variants; verified atomic resume chunks."""
import sys,time,json,collections,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
import av,cv2,pyarrow as pa,pyarrow.parquet as pq

def frame_regions(gray,previous,objects):
    # Full image scan is only candidate evidence, not a segmented physical channel.
    h,w=gray.shape;xy=np.array([o['xy'] for o in objects]).reshape(-1,2);tiles=[]
    for y in range(0,h,128):
        for x in range(0,w,128):
            crop=gray[y:y+128,x:x+128];mask=(xy[:,0]>=x)&(xy[:,0]<x+crop.shape[1])&(xy[:,1]>=y)&(xy[:,1]<y+crop.shape[0]);n=int(mask.sum())
            delta=float(np.mean(np.abs(crop.astype(float)-previous[y:y+128,x:x+128]))) if previous is not None else None
            status='WARMUP' if delta is None else ('LOW_CANDIDATE_SUPPORT' if n<8 else ('STATIC_OR_BACKGROUND_AMBIGUITY' if delta<1 else 'CANDIDATE_TEMPORAL_SUPPORT_UNCALIBRATED'))
            tiles.append(dict(x=x,y=y,w=crop.shape[1],h=crop.shape[0],candidate_count=n,mean_abs_frame_change=delta,status=status))
    return tiles

def main():
    audit=json.loads((OUT/'results/full_video_source_audit.json').read_text());assert sha(VIDEO)==audit['sha256']==P['video_sha256']
    chunks=OUT/'results/full_video_chunks';chunks.mkdir(exist_ok=True)
    # Resume requires same numerical source and config, independent of later documentation commits.
    dependency={str(p.relative_to(ROOT.parent)):sha(p) for d in [ROOT/'src',ROOT.parent/'collective_motion_v1/src',ROOT.parent/'scripts',ROOT.parent/'configs'] for p in d.glob('*.py')}
    key=dict(video_sha=sha(VIDEO),protocol_sha=sha(ROOT/'configs/protocol.json'),runner_sha=sha(__file__),dependencies=dependency)
    keypath=chunks/'RUN_IDENTITY.json'
    if keypath.exists():assert json.loads(keypath.read_text())==key,'INCOMPATIBLE_RESUME'
    else:dump(keypath,key)
    loaded={}
    for p in sorted(chunks.glob('chunk_*.json')):
        info=json.loads(p.with_suffix('.sha.json').read_text()) if not p.name.endswith('.sha.json') else None
        if info is None:continue
        assert sha(p)==info['sha256']
        for row in json.loads(p.read_text()):loaded[row['frame_id']]=row
    start=time.perf_counter();monitor=ResourceMonitor();det=old.CandidateDetector();history=collections.deque(maxlen=10);pending=[];allrows=[];previous=None;x,y,w,h=P['roi_xywh']
    variants=[('candidate',3),('candidate',5),('candidate',10),('raw',3),('local',3)]
    def save():
        if not pending:return
        path=chunks/f"chunk_{pending[0]['frame_id']:06d}_{pending[-1]['frame_id']:06d}.json"
        if path.exists():raise RuntimeError('REFUSE_OVERWRITE_CHUNK')
        dump(path,pending);dump(path.with_suffix('.sha.json'),dict(sha256=sha(path),frames=len(pending)));pending.clear()
    with av.open(str(VIDEO)) as c:
        for fid,f in enumerate(c.decode(video=0)):
            tick=time.perf_counter();im=f.to_ndarray(format='bgr24');pts=float(f.time);roi=im[y:y+h,x:x+w]
            objects,gray=det(roi);history.append((fid,pts,gray,objects));fullgray=cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
            if fid in loaded:allrows.append(loaded[fid]);previous=fullgray;continue
            detection_ms=(time.perf_counter()-tick)*1000;dist,_=old.distribution(objects,gray.shape);methods={}
            for mode,n in variants:
                method=f'{mode}_{n}'
                if len(history)<n:result=dict(estimate=old.empty('WARMUP',method),vectors=[],valid_vector_count=0,algorithm_ms=0)
                elif history[-1][1]-history[-n][1]<=0 or any(history[j+1][1]-history[j][1]>.15 for j in range(len(history)-n,len(history)-1)):
                    result=dict(estimate=old.empty('PTS_DISCONTINUITY',method),vectors=[],valid_vector_count=0,algorithm_ms=0)
                else:
                    a=history[-n];b=history[-1]
                    try:result=old.piv(a[2],b[2],b[1]-a[1],mode=mode,objects_a=a[3],objects_b=b[3],bootstrap=0)
                    except Exception as e:result=dict(estimate=old.empty('PROCESSING_EXCEPTION',method),exception=repr(e),vectors=[],valid_vector_count=0,algorithm_ms=0)
                methods[method]=result
            main=methods['candidate_3'];main_ms=detection_ms+main.get('algorithm_ms',0)+main['estimate'].get('fit_ms',0)
            t=time.perf_counter();fullobjects,_=det(im);tiles=frame_regions(fullgray,previous,fullobjects)
            row=dict(frame_id=fid,pts_s=pts,source_pts=f.pts,roi=P['roi_xywh'],processing_status='WARMUP' if fid<2 else ('VALID' if all(main['estimate']['valid'].values()) else 'NOT_OBSERVABLE'),objects=objects,distribution=dist,methods=methods,full_frame_candidates=fullobjects,full_frame_tiles=tiles,focus_laplacian_var=float(cv2.Laplacian(gray,cv2.CV_32F).var()),detection_decode_ms=detection_ms,main_pipeline_component_ms=main_ms,fullframe_scan_ms=(time.perf_counter()-t)*1000,total_offline_ms=(time.perf_counter()-tick)*1000)
            pending.append(row);allrows.append(row);previous=fullgray
            if len(pending)>=P['chunk_frames']:save();print('saved through',fid,flush=True)
    save();resource=monitor.close();assert [r['frame_id'] for r in allrows]==list(range(audit['decoded_frames']))
    summary=[];distributions=[];vectors=[];regions=[]
    for row in allrows:
        dist=row['distribution'];distributions.append(dict(frame_id=row['frame_id'],pts_s=row['pts_s'],**dist))
        for name,result in row['methods'].items():
            e=result['estimate'];summary.append(dict(frame_id=row['frame_id'],pts_s=row['pts_s'],method=name,roi=row['roi'],candidate_count=len(row['objects']),status=e['observability_status'],reason=next((v for v in e['reasons'].values() if v is not None),None),valid_all=all(e['valid'].values()),valid_vectors=result.get('valid_vector_count',0),coverage=e.get('coverage_fraction'),condition=e.get('condition_number'),residual=e.get('residual'),values=e['values'],valid=e['valid'],reasons=e['reasons'],algorithm_ms=result.get('algorithm_ms'),fit_ms=e.get('fit_ms'),main_pipeline_component_ms=row['main_pipeline_component_ms']))
            for v in result.get('vectors',[]):vectors.append(dict(frame_id=row['frame_id'],method=name,x=v[0],y=v[1],dx=v[2],dy=v[3],valid=bool(v[4])))
        for tile in row['full_frame_tiles']:regions.append(dict(frame_id=row['frame_id'],pts_s=row['pts_s'],**tile))
    for name,rows in [('full_video_piv_per_frame',summary),('full_video_distribution',distributions),('full_video_vectors',vectors),('full_video_regions',regions)]:pq.write_table(pa.Table.from_pylist(clean(rows)),OUT/'results'/f'{name}.parquet',compression='zstd')
    valid=sum(all(r['methods']['candidate_3']['estimate']['valid'].values()) for r in allrows);reasons=collections.Counter(next((v for v in r['methods']['candidate_3']['estimate']['reasons'].values() if v), 'VALID') for r in allrows)
    validity=dict(decoded_frames=audit['decoded_frames'],analyzed_frames=len(allrows),frame_coverage=len(allrows)/audit['decoded_frames'],valid_motion_frames=valid,invalid_motion_frames=len(allrows)-valid,warmup_frames=2,valid_motion_fraction=valid/len(allrows),reasons=dict(reasons),missing_frame_ids=[],source_first_pts_s=allrows[0]['pts_s'],source_last_pts_s=allrows[-1]['pts_s'],real_accuracy='NOT_VERIFIED',identity=identity())
    dump(OUT/'results/full_video_validity.json',validity);dump(OUT/'results/full_video_runtime.json',dict(identity=identity(),resource=resource,invocation_wall_s=time.perf_counter()-start,resumed_frames=len(loaded),main_pipeline_components_ms=stats([r['main_pipeline_component_ms'] for r in allrows]),all_variants_fullscan_offline_ms=stats([r['total_offline_ms'] for r in allrows]),scope='offline timing is not paced real-time throughput'))
    print(json.dumps(validity),flush=True)
if __name__=='__main__':main()
