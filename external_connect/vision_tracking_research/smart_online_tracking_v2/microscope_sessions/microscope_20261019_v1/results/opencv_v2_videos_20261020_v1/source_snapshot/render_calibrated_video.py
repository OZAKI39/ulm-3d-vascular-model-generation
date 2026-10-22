"""Frame-complete offline OpenCV overlays; exact source presentation timing.

The camera rectangle excludes already baked text/PZT graphics. It is not a
verified vessel mask. Experimental profiles require explicit offline selection.
"""
from pathlib import Path
from fractions import Fraction
import argparse, gzip, hashlib, json, os, subprocess, sys, time

BASE=Path(__file__).resolve().parents[1]
LIVE=BASE.parent/'realtime_collective_observation_v1'
sys.path.insert(0,str(LIVE/'src'))
import av
import cv2
import numpy as np
from microscope_detector import detect,validate_config
from video_text import draw_text,text_width,font_record


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def dump(path,value):
    path=Path(path);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False));tmp.replace(path)


def overlay(image,objects,box):
    x,y,w,h=box;result=image.copy()
    # Same candidate-center symbol and palette as the existing observation UI.
    for o in objects:
        px,py=o['xy_local']
        center=(round((px+x)*16),round((py+y)*16))
        cv2.circle(result,center,4*16,(250,233,115),1,cv2.LINE_AA,shift=4)
    text=f'Candidates: {len(objects)}';size=20;width=text_width(text,size)
    left=round(x+w-width-14)
    cv2.rectangle(result,(left-7,0),(x+w-1,29),(20,22,20),-1)
    draw_text(result,text,left,23,size,(230,235,235))
    return result


def run(source,out,profile_path,box,allow_experimental):
    started=time.monotonic();cv2.setNumThreads(2)
    profile=json.loads(profile_path.read_text());params=validate_config(profile['parameters'])
    if not profile.get('admitted') and not allow_experimental:
        raise ValueError('Experimental profile needs explicit offline authorization')
    out.mkdir(parents=True,exist_ok=False)
    movie=out/f'{source.stem}_opencv_calibrated.mp4';partial=movie.with_suffix('.partial.mp4')
    frame_log=out/'detections.jsonl.gz'
    identity=dict(source=str(source),source_sha256=sha(source),profile_path=str(profile_path),
        profile_sha256=sha(profile_path),parameters=params,profile_globally_admitted=profile.get('admitted',False),
        use='USER_REQUESTED_OFFLINE_VISUALIZATION',channel_mask='NOT_DRAWN',
        detection_scope='Camera rectangle, including unmasked wall/outside-lumen areas; not vessel bubble counts',
        camera_crop_xywh=box,objects='Unreviewed image candidates, not physical bubble identities or boundaries',
        display='Original camera and baked PZT panel; cyan radius-4 centers; Arial candidate count',
        temporal_processing='Every decoded frame; no temporal smoothing, skipped frames or invented trajectories',
        pzt='Preserve baked source pixels; UNKNOWN readback remains UNKNOWN',
        inference='OpenCV calibrated dark components; no new PIV or inferred channel mask',
        source_hashes={str(p):sha(p) for p in [Path(__file__),LIVE/'src/microscope_detector.py',LIVE/'src/video_text.py']},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),font=font_record(),
        status='RUNNING')
    dump(out/'manifest.json',identity)
    expected_times=[];expected_durations={};counts=[];det_ms=[];first_pts=None;last_time=None;unchanged=True
    with av.open(str(source)) as src:
        stream=src.streams.video[0];tb=stream.time_base
        x,y,w,h=box
        if (stream.width,stream.height)!=(1964,1024) or box!=[0,40,1280,984]:
            raise ValueError('This overlay layout must be visually inspected before processing')
        if src.streams.audio:raise ValueError('Audio-bearing video needs an explicit preservation path')
        with av.open(str(partial),'w',options={'movflags':'+faststart'}) as dst, gzip.open(frame_log,'wt') as log:
            enc=dst.add_stream('libx264',rate=Fraction(30,1));enc.width=stream.width;enc.height=stream.height
            enc.pix_fmt='yuv420p';enc.time_base=tb;enc.codec_context.time_base=tb
            enc.options={'crf':'18','preset':'fast','threads':'2','bf':'0','g':'60'}
            def mux(packet):
                pts=Fraction(packet.pts)*packet.time_base
                packet.duration=round(expected_durations[pts]/packet.time_base)
                dst.mux(packet)
            for index,frame in enumerate(src.decode(stream)):
                if frame.pts is None:raise ValueError('Missing source presentation timestamp')
                absolute=Fraction(frame.pts)*frame.time_base
                if first_pts is None:first_pts=absolute
                presentation=absolute-first_pts
                if last_time is not None and presentation<=last_time:raise ValueError('Nonmonotonic source PTS')
                last_time=presentation;image=frame.to_ndarray(format='bgr24')
                tick=time.perf_counter();gray=cv2.cvtColor(image[y:y+h,x:x+w],cv2.COLOR_BGR2GRAY)
                _,_,objects=detect(gray,params);det_ms.append((time.perf_counter()-tick)*1000)
                rendered=overlay(image,objects,box)
                unchanged &= bool(np.array_equal(rendered[:,1280:],image[:,1280:]))
                output=av.VideoFrame.from_ndarray(rendered,format='bgr24')
                output.pts=round(presentation/tb);output.time_base=tb
                if not frame.duration:raise ValueError('Missing source frame duration')
                expected_durations[presentation]=Fraction(frame.duration)*frame.time_base
                output.duration=round(expected_durations[presentation]/tb)
                for packet in enc.encode(output):mux(packet)
                record=dict(frame=index,source_pts=frame.pts,source_time_base=str(frame.time_base),
                    output_pts=output.pts,output_time_base=str(tb),source_time_s=float(absolute),
                    relative_time_s=float(presentation),count=len(objects),
                    centers_xy=[[o['xy_local'][0]+x,o['xy_local'][1]+y] for o in objects],areas_px=[o['area'] for o in objects])
                log.write(json.dumps(record,separators=(',',':'))+'\n')
                counts.append(len(objects));expected_times.append(presentation)
                if index in (0,stream.frames//2,max(0,stream.frames-1)):
                    cv2.imwrite(str(out/f'preview_{index:05d}.jpg'),rendered,[cv2.IMWRITE_JPEG_QUALITY,94])
                if index%300==0:
                    dump(out/'progress.json',dict(frames=index+1,expected_source_frames=stream.frames,
                        source_time_s=float(presentation),wall_s=time.monotonic()-started))
                    print(json.dumps(dict(video=source.name,frames=index+1,time_s=float(presentation))),flush=True)
            for packet in enc.encode():mux(packet)
    # Full decode validates playability and per-frame presentation identity.
    deviations=[];decoded=0
    with av.open(str(partial)) as c:
        s=c.streams.video[0]
        for index,f in enumerate(c.decode(s)):
            if index>=len(expected_times):raise ValueError('Extra encoded frame')
            deviations.append(float(abs(Fraction(f.pts)*f.time_base-expected_times[index])))
            decoded+=1
        output_duration=float(s.duration*s.time_base) if s.duration is not None else None
    expected_duration=float(expected_times[-1]+expected_durations[expected_times[-1]])
    if decoded!=len(expected_times) or max(deviations,default=0)>1e-6 or not unchanged or abs(output_duration-expected_duration)>1e-6:
        raise RuntimeError('Frame/timestamp/PZT preservation validation failed')
    partial.replace(movie)
    identity.update(status='COMPLETED',output=str(movie),output_sha256=sha(movie),
        detections_sha256=sha(frame_log),decoded_input_frames=len(expected_times),decoded_output_frames=decoded,
        skipped_frames=0,maximum_pts_error_s=max(deviations,default=0),output_duration_s=output_duration,expected_duration_s=expected_duration,
        pzt_pixels_unchanged_before_lossy_encode=unchanged,
        candidate_counts=dict(min=min(counts),max=max(counts),median=float(np.median(counts))),
        detector_ms=dict(p50=float(np.median(det_ms)),p95=float(np.percentile(det_ms,95))),
        wall_s=time.monotonic()-started,training_ready=False,control_permission=False)
    dump(out/'manifest.json',identity);print(json.dumps(identity),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--profile',type=Path,required=True);p.add_argument('--camera-box',nargs=4,type=int,default=[0,40,1280,984])
    p.add_argument('--allow-experimental-profile',action='store_true');a=p.parse_args()
    run(a.source.resolve(),a.output.resolve(),a.profile.resolve(),a.camera_box,a.allow_experimental_profile)
