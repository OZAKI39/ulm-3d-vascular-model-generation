"""Bounded, timestamped experiment recording; never captures the Qt window.

Two H.264 files share the same received images and presentation timestamps.
Annotations use an exact analyzed frame when available, otherwise an explicitly
labelled, recent *past* result. They never borrow a future frame or ROI.
"""
from collections import OrderedDict, deque
from dataclasses import replace
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import json
import queue
import threading
import time
import traceback

import av
import cv2
import numpy as np

from ycommon import clean, dump
from piv_display import piv_arrows,arrow_shape,draw_cv_arrows,DISPLAY_SPEC,REFERENCE_SPEED,GAIN_SECONDS,MAX_LENGTH_PX,display_vectors,display_is_smoothed
from overlay_style import DEFAULT_LAYERS,DISPLAY_POLICY,overlay_opacity,overlay_caption
from pzt_video import PZTVideoRenderer, sample_id, EXPERIMENT_MODEL

TIME_BASE = Fraction(1, 1_000_000)


def render_overlay(packet, state, geometry, layers, match, lag_ms, utc, elapsed,ttl_ms=250.):
    """Draw in original image pixels, with the same overlay colours as Canvas."""
    a = packet.image
    image = cv2.cvtColor(a, cv2.COLOR_GRAY2BGR) if a.ndim == 2 else a.copy()
    origin = np.asarray((packet.roi or [0, 0])[:2], float)
    g = geometry or (state.get('geometry') if state else None)
    if g and tuple(g['image_wh']) == packet.original_resolution and layers.get('Mask', True):
        polygon = np.rint(np.asarray(g['channel_polygon'])-origin).astype(np.int32)
        cv2.polylines(image, [polygon], True, (208,223,85), 2, cv2.LINE_AA)
        if g.get('excluded_polygons'):
            shade = image.copy()
            for poly in g['excluded_polygons']:
                cv2.fillPoly(shade, [np.rint(np.asarray(poly)-origin).astype(np.int32)], (40,30,30))
            image = cv2.addWeighted(shade, .3, image, .7, 0)
    opacity=overlay_opacity((lag_ms or 0.)+state.get('age_upper_bound_at_arrival_ms',0.),ttl_ms) if state else 0.
    if state and state.get('geometry') and opacity>0:
        overlay = image if opacity==1. else image.copy()
        offset = np.asarray(state['geometry']['crop_xywh'][:2], float)-origin
        def point(v):return tuple(np.rint(np.asarray(v)+offset).astype(int))
        if layers.get('Candidates', True):
            for q in state.get('candidates', []):cv2.circle(overlay, point(q), 4, (250,233,115), 1, cv2.LINE_AA)
        if layers.get('Validity', False):
            for q, ok in zip(state.get('window_grid_xy', []),state.get('window_valid', [])):
                if not ok:cv2.circle(overlay,point(q),1,(90,195,245),-1)
        if (layers.get('Regions', False) or layers.get('Shape', False)) and state.get('aggregation'):
            for cluster in state['aggregation']['clusters']:
                x,y,x1,y1 = cluster['extent']
                if layers.get('Regions', False):cv2.rectangle(overlay,point([x,y]),point([x1,y1]),(216,166,227),1,cv2.LINE_AA)
                if layers.get('Shape', False):
                    cv2.circle(overlay,point(cluster['center']),max(0,round(cluster['rg_px'])),(216,166,227),1,cv2.LINE_AA)
        if state.get('valid') and layers.get('PIV', True):
            smoothed=display_is_smoothed(state,layers)
            arrows=piv_arrows(display_vectors(state,smoothed),offset,sample=not smoothed)
            draw_cv_arrows(overlay,arrows)
        image = overlay if opacity==1. else cv2.addWeighted(overlay,opacity,image,1.-opacity,0)
        if state.get('valid') and layers.get('PIV', True):
            if arrows:
                y=image.shape[0]-18
                cv2.rectangle(image,(6,y-16),(260,y+12),(25,25,25),-1)
                draw_cv_arrows(image,[arrow_shape([16,y],[REFERENCE_SPEED,0])])
                cv2.putText(image,f'{REFERENCE_SPEED:g} px/s | '+('smoothed' if smoothed else 'raw'),(49,y+5),cv2.FONT_HERSHEY_SIMPLEX,.45,(235,235,235),1,cv2.LINE_AA)
                v=np.asarray(display_vectors(state,smoothed)).reshape(-1,4)
                if len(v) and np.any(np.linalg.norm(v[:,2:],axis=1)*GAIN_SECONDS>MAX_LENGTH_PX):
                    cv2.putText(image,'Long arrows capped',(image.shape[1]-170,y+5),cv2.FONT_HERSHEY_SIMPLEX,.4,(255,226,59),1,cv2.LINE_AA)
    # Camera caption; the PZT sidebar is composed separately, without UI controls.
    lines = [f"UTC {utc.replace('T',' ')[:23]} | t={elapsed:.3f} s"]
    if packet.source_mode=='REALTIME_REPLAY':lines[0]+=' | Replay'
    if state:
        count = len(state.get('candidates', []))
        lines.append(f'Candidates: {count} | '+overlay_caption(lag_ms or 0.))
        if not state.get('valid'):lines.append('Motion unavailable')
        elif state.get('motion_reference')=='IMAGE_SPACE_UNREGISTERED':lines.append('Image motion; camera drift not corrected')
    else:
        lines.append('Waiting for current analysis')
    font_scale = max(.28, min(.56, image.shape[1]/1200))
    line_h = max(13, round(font_scale*29))
    cv2.rectangle(image,(0,0),(image.shape[1]-1,min(image.shape[0]-1,line_h*3+8)),(25,25,25),-1)
    for i,line in enumerate(lines):
        # Fit narrow camera crops without changing image size.
        width = cv2.getTextSize(line,cv2.FONT_HERSHEY_SIMPLEX,font_scale,1)[0][0]
        scale = font_scale*min(1.,(image.shape[1]-12)/max(width,1))
        cv2.putText(image,line,(6,(i+1)*line_h),cv2.FONT_HERSHEY_SIMPLEX,scale,(225,225,225),1,cv2.LINE_AA)
    return image


class VideoPair:
    def __init__(self, directory, index, shape, overlay_shape=None):
        self.shape = shape[:2]
        self.outputs = []
        suffix = '' if index == 1 else f'_part{index:03d}'
        self.names = [f'experiment_raw{suffix}.mp4', f'experiment_overlay{suffix}.mp4']
        h,w = self.shape
        self.encoded_wh = [w+w%2,h+h%2]  # Pad odd dimensions; never crop pixels.
        oh,ow = (overlay_shape or shape)[:2]
        self.overlay_encoded_wh = [ow+ow%2,oh+oh%2]
        try:
            for name,wh in zip(self.names,(self.encoded_wh,self.overlay_encoded_wh)):
                container = av.open(str(directory/name),'w',options={'movflags':'+frag_keyframe+empty_moov+default_base_moof'})
                stream = container.add_stream('libx264',rate=30)
                self.outputs.append((container,stream))
                stream.width,stream.height = wh
                stream.pix_fmt = 'yuv420p'
                stream.time_base = TIME_BASE
                stream.codec_context.time_base = TIME_BASE
                stream.options = {'crf':'18','preset':'ultrafast','tune':'zerolatency','threads':'2','bf':'0','g':'30'}
        except Exception:
            for container,_ in self.outputs:container.close()
            raise
        self.last_pts = -1
        self.frames = 0

    def write(self, raw, overlay, pts):
        pts = max(self.last_pts+1, pts)
        duration = pts-self.last_pts if self.last_pts >= 0 else 33333
        for a,(container,stream) in zip((raw,overlay),self.outputs):
            h,w = a.shape[:2]
            if h%2 or w%2:a = cv2.copyMakeBorder(a,0,h%2,0,w%2,cv2.BORDER_REPLICATE)
            frame = av.VideoFrame.from_ndarray(a,format='gray' if a.ndim==2 else 'bgr24')
            frame.pts,frame.time_base,frame.duration = pts,TIME_BASE,duration
            for packet in stream.encode(frame):container.mux(packet)
        self.last_pts = pts
        self.frames += 1
        return pts

    def close(self):
        errors = []
        for container,stream in self.outputs:
            try:
                for packet in stream.encode():container.mux(packet)
            except Exception as exc:errors.append(str(exc))
            finally:
                try:container.close()
                except Exception as exc:errors.append(str(exc))
        if errors:raise RuntimeError('; '.join(errors))


class ExperimentVideoRecorder:
    def __init__(self, directory, config, start_monotonic):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True,exist_ok=True)
        self.start_monotonic = start_monotonic
        self.start_utc = datetime.now(timezone.utc).isoformat(timespec='microseconds')
        self.ttl_ms = float(config.get('ttl_ms',250))
        self.diagnostics = bool(config.get('diagnostics',False))
        self.wait_s = min(float(config.get('video_annotation_wait_ms',250))/1000,self.ttl_ms/1000)
        # Bounded startup cushion: CUDA warmup can briefly stall Python and
        # release a burst of received frames. This does not extend analysis TTL.
        self.queue = queue.Queue(max(1,int(config.get('video_queue',48))))
        self.condition = threading.Condition()
        self.states = OrderedDict()
        self.done = threading.Event()
        self.closed = False
        self.error = None
        self.submitted = self.written = self.dropped = 0
        self.drop_examples = deque(maxlen=100)
        self.parts = []
        self.first_arrival = None
        self.last_elapsed_s = 0.
        self.match_counts = dict(EXACT_FRAME=0,RECENT_PAST_FRAME=0,NO_CURRENT_ANALYSIS=0)
        self.pzt_renderer = PZTVideoRenderer()
        self.pzt_sample_count = 0
        self.thread = threading.Thread(target=self.run,name='experiment-video-encoder',daemon=True)
        self.thread.start()

    def submit(self, packet, geometry_revision, geometry, layers, pzt=None):
        if self.closed or self.done.is_set():return False
        self.submitted += 1
        if self.error or self.queue.full():
            self.dropped += 1
            self.drop_examples.append(dict(session_id=packet.session_id,frame_id=packet.frame_id,
                arrival_monotonic=packet.arrival_monotonic,reason='ENCODER_ERROR' if self.error else 'QUEUE_FULL'))
            return False
        item = dict(packet=replace(packet,image=packet.image.copy()),revision=geometry_revision,
            geometry=geometry,layers=dict(layers),enqueued=time.monotonic(),
            receive_utc=datetime.now(timezone.utc).isoformat(timespec='microseconds'),
            pzt=clean(pzt or dict(frame_unix_s=time.time(),sample=dict(scope='NOT_RECORDED'))))
        try:self.queue.put_nowait(item)
        except queue.Full:
            self.dropped += 1
            self.drop_examples.append(dict(session_id=packet.session_id,frame_id=packet.frame_id,reason='QUEUE_FULL'))
            return False
        return True

    def add_analysis(self,state):
        if self.closed or self.error:return
        with self.condition:
            self.states[(state['session_id'],state['frame_id'])] = state
            while len(self.states)>256:self.states.popitem(last=False)
            self.condition.notify_all()

    def choose_state(self,item):
        p = item['packet']
        with self.condition:
            deadline = item['enqueued']+self.wait_s
            while (p.session_id,p.frame_id) not in self.states and not self.done.is_set():
                left = deadline-time.monotonic()
                if left <= 0:break
                self.condition.wait(left)
            for state in reversed(list(self.states.values())):
                if state['session_id']!=p.session_id or state['frame_id']>p.frame_id:continue
                if state.get('geometry_revision')!=item['revision']:continue
                if not state.get('geometry'):continue
                lag = (p.arrival_monotonic-state['arrival_monotonic'])*1000
                if lag < -1e-5 or lag+state.get('age_upper_bound_at_arrival_ms',0)>self.ttl_ms:continue
                match = 'EXACT_FRAME' if state['frame_id']==p.frame_id else 'RECENT_PAST_FRAME'
                return state,match,max(lag,0.)
        return None,'NO_CURRENT_ANALYSIS',None

    def status(self):
        return dict(state='ERROR' if self.error else 'SAVED' if self.closed and self.written else
                    'NO_FRAMES' if self.closed else 'RECORDING' if self.written else 'WAITING_FOR_FRAMES',
            written_frames=self.written,submitted_frames=self.submitted,dropped_frames=self.dropped,
            queue_size=self.queue.qsize(),error=self.error,directory=str(self.directory))

    def report(self):
        status=self.status()
        if not self.diagnostics:status.pop('queue_size',None)
        return dict(**status,schema='EXPERIMENT_VIDEO_1.1',start_utc=self.start_utc,
            start_monotonic=self.start_monotonic,first_received_monotonic=self.first_arrival,
            last_frame_elapsed_s=self.last_elapsed_s,time_base='1/1000000',
            timeline='WSL/local monotonic receive intervals; variable frame rate; no synthetic frames',
            utc_scope='local host receive wall clock; not camera exposure time',
            recording_scope='Received preview frames only; not the camera internal high-speed recording',
            raw_encoding='H.264 CRF18 yuv420p; unannotated, not lossless sensor data',
            annotations='Exact frame preferred; recent past overlay delay labelled; future/stale/other-ROI results excluded',
            piv_display=DISPLAY_SPEC,overlay_display=DISPLAY_POLICY,
            pzt_recording=dict(samples='pzt_readbacks.jsonl',geometry='pzt_experiment_model.json',
                sample_count=self.pzt_sample_count,frame_reference='video_frames.jsonl: pzt_sample_id, pzt_frame_unix_s',
                association='Latest monitor observation available at local frame arrival; no future readback; not exposure synchronized',
                replay_policy='No live AFG values assigned to historical replay',
                overlay='Native camera pixels plus black PZT status/layout sidebar; raw video dimensions unchanged'),
            source_modes=sorted({part['source_mode'] for part in self.parts}),
            match_counts=self.match_counts,parts=self.parts,dropped_frame_examples=list(self.drop_examples),
            metadata='video_frames.jsonl',video_complete=bool(self.closed and not self.error and self.written and not self.dropped))

    def run(self):
        pair = None
        metadata = None
        pzt_metadata = None
        last_pzt_id = None
        try:
            dump(self.directory/'video_manifest.json',self.report())
            dump(self.directory/'pzt_experiment_model.json',EXPERIMENT_MODEL)
            metadata = (self.directory/'video_frames.jsonl').open('x',encoding='utf-8',buffering=1)
            pzt_metadata = (self.directory/'pzt_readbacks.jsonl').open('x',encoding='utf-8',buffering=1)
            while not self.done.is_set() or not self.queue.empty():
                try:item = self.queue.get(timeout=.05)
                except queue.Empty:continue
                p = item['packet']
                state,match,lag = self.choose_state(item)
                if self.first_arrival is None:self.first_arrival = p.arrival_monotonic
                elapsed = p.arrival_monotonic-self.first_arrival
                self.last_elapsed_s = elapsed
                overlay = render_overlay(p,state,item['geometry'],item['layers'],match,lag,
                    item['receive_utc'],elapsed,self.ttl_ms)
                overlay = self.pzt_renderer.compose(overlay,item['pzt'])
                pzt_id = sample_id(item['pzt']['sample'])
                if pzt_id != last_pzt_id:
                    pzt_metadata.write(json.dumps(dict(sample_id=pzt_id,**item['pzt']['sample']),separators=(',',':'),allow_nan=False)+'\n')
                    self.pzt_sample_count += 1
                    last_pzt_id = pzt_id
                if pair is None or pair.shape!=p.image.shape[:2]:
                    if pair:pair.close();pair=None
                    pair = VideoPair(self.directory,len(self.parts)+1,p.image.shape,overlay.shape)
                    self.parts.append(dict(part=len(self.parts)+1,raw=pair.names[0],overlay=pair.names[1],
                        source_mode=p.source_mode,image_wh=[p.image.shape[1],p.image.shape[0]],encoded_wh=pair.encoded_wh,
                        overlay_encoded_wh=pair.overlay_encoded_wh,
                        overlay_camera_xy=[0,(overlay.shape[0]-p.image.shape[0])//2],
                        start_elapsed_s=elapsed,frames=0))
                part = self.parts[-1]
                pts = round((elapsed-part['start_elapsed_s'])*1e6)
                encoded_pts = pair.write(p.image,overlay,pts)
                row = dict(record_index=self.written,part=part['part'],part_frame_index=pair.frames-1,
                    video_pts_us=encoded_pts,requested_pts_us=pts,pts_adjusted=encoded_pts!=pts,
                    elapsed_since_first_frame_s=elapsed,elapsed_since_start_s=p.arrival_monotonic-self.start_monotonic,
                    receive_utc=item['receive_utc'],arrival_monotonic=p.arrival_monotonic,
                    source_mode=p.source_mode,source_name=p.source_name,source_session_id=p.session_id,source_frame_id=p.frame_id,
                    source_pts_s=p.source_pts,source_timestamp_quality=p.timestamp_quality,
                    capture_timestamp=p.capture_timestamp,capture_timestamp_status=p.capture_timestamp_status,
                    host_decode_timestamp=p.host_decode_timestamp,host_receive_timestamp=p.host_receive_timestamp,
                    age_upper_bound_at_arrival_ms=p.age_upper_bound_at_arrival_ms,original_resolution=p.original_resolution,
                    source_roi=p.roi,geometry_revision=item['revision'],layers=item['layers'],
                    annotation_match=match,analysis_frame_id=state['frame_id'] if state else None,
                    analysis_session_id=state['session_id'] if state else None,annotation_lag_ms=lag,
                    analysis_completed_monotonic=state.get('completed_monotonic') if state else None,
                    motion_valid=bool(state and state.get('valid')),motion_reference=state.get('motion_reference') if state else None,
                    piv_display_smoothed=bool(state and display_is_smoothed(state,item['layers'])),
                    candidate_count=len(state.get('candidates',[])) if state else None,
                    pzt_sample_id=pzt_id,pzt_frame_unix_s=item['pzt']['frame_unix_s'])
                metadata.write(json.dumps(clean(row),separators=(',',':'),allow_nan=False)+'\n')
                self.written += 1;part['frames'] += 1;self.match_counts[match] += 1
                if self.written%30==1:dump(self.directory/'video_manifest.json',self.report())
        except Exception:self.error = traceback.format_exc()
        finally:
            if pair:
                try:pair.close()
                except Exception:self.error = (self.error or '')+'\n'+traceback.format_exc()
            if metadata:metadata.close()
            if pzt_metadata:pzt_metadata.close()
            self.closed = True
            try:dump(self.directory/'video_manifest.json',self.report())
            except Exception:self.error = (self.error or '')+'\n'+traceback.format_exc()

    def close(self):
        self.done.set()
        with self.condition:self.condition.notify_all()
        self.thread.join(25)
        if self.thread.is_alive():
            self.error = 'VIDEO_ENCODER_STOP_TIMEOUT; output is not finalized'
            raise RuntimeError(self.error)
        return self.report()
