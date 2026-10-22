"""Full decode, integer PTS and duplicate-picture audit; never substitute sources."""
import sys,time,json,hashlib,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
import av,cv2,pyarrow as pa,pyarrow.parquet as pq
def main():
    start=time.perf_counter();h=sha(VIDEO);assert h==VIDEO_SHA,'VIDEO_IDENTITY_MISMATCH'
    rows=[];seen={}
    with av.open(str(VIDEO)) as c:
        s=c.streams.video[0];meta=dict(codec=s.codec_context.name,width=s.width,height=s.height,declared_frames=s.frames,nominal_fps=float(s.average_rate),stream_time_base=str(s.time_base),container_duration_s=c.duration/1e6 if c.duration else None)
        for fid,f in enumerate(c.decode(video=0)):
            gray=f.to_ndarray(format='gray');fh=hashlib.sha256(gray.tobytes()).hexdigest()
            pts=float(f.pts*f.time_base) if f.pts is not None else None
            rows.append(dict(frame_id=fid,pts=f.pts,time_base_num=f.time_base.numerator,time_base_den=f.time_base.denominator,pts_s=pts,duration_pts=f.duration,gray_sha256=fh,duplicate_of_frame=seen.get(fh),decode_status='OK'))
            seen.setdefault(fh,fid)
            if fid==0:cv2.imwrite(str(OUT/'evidence/source_first.png'),f.to_ndarray(format='bgr24'))
        cv2.imwrite(str(OUT/'evidence/source_last.png'),f.to_ndarray(format='bgr24'))
    ts=np.array([r['pts_s'] for r in rows]);dt=np.diff(ts);median=float(np.median(dt));monotonic=bool(np.isfinite(ts).all() and (dt>0).all())
    pq.write_table(pa.Table.from_pylist(rows),OUT/'results/full_video_frame_timestamps.parquet',compression='zstd')
    record=dict(identity=identity(),path=str(VIDEO),sha256=h,bytes=VIDEO.stat().st_size,**meta,decoded_frames=len(rows),first_pts_s=float(ts[0]),last_pts_s=float(ts[-1]),presentation_duration_s=float(ts[-1]-ts[0]+median),strictly_monotonic_pts=monotonic,frame_intervals_s=stats(dt),abnormal_interval_frame_ids=(np.flatnonzero((dt<median*.5)|(dt>median*1.5))+1).tolist(),identical_decoded_picture_count=sum(r['duplicate_of_frame'] is not None for r in rows),dropped_capture_frames='NOT_DETERMINABLE_FROM_ENCODED_FILE_WITHOUT_CAPTURE_COUNTER',full_decode_pass=True,elapsed_s=time.perf_counter()-start)
    dump(OUT/'results/full_video_source_audit.json',record);print(json.dumps(record,default=str),flush=True)
if __name__=='__main__':main()
