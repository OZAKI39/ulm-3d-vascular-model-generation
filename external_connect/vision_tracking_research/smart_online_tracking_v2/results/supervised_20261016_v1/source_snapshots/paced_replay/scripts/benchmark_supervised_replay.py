"""Native-rate, latest-slot saved-video benchmark. Never a camera latency test."""
from pathlib import Path
import sys, json, time, threading, copy
import av, cv2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from runtime import MockConsumer, ResourceMonitor, enforce_freshness, invalid_state, identity
from vision_engine import VisionEngine, stats
from annotation_store import sha
from prepare_candidates import VIDEO


class PacedReplay:
    def __init__(self, roi):
        self.roi = roi
        self.lock = threading.Lock()
        self.latest = None
        self.done = threading.Event()
        self.headers = []
        self.error = None
        self.thread = threading.Thread(target=self.run, daemon=True)

    def run(self):
        try:
            x, y, w, h = self.roi
            with av.open(str(VIDEO)) as container:
                stream = container.streams.video[0]
                stream.codec_context.thread_count = 1
                start = first_pts = None
                for index, frame in enumerate(container.decode(stream)):
                    if index < 360:
                        continue
                    if index > 659:
                        break
                    pts = float(frame.time)
                    if start is None:
                        start, first_pts = time.monotonic(), pts
                    release = start + pts - first_pts
                    time.sleep(max(0, release - time.monotonic()))
                    image = frame.to_ndarray(format='bgr24')[y:y+h, x:x+w].copy()
                    arrived = time.monotonic()
                    header = dict(seq=index, source_time_s=pts, session='PACED_CAUSAL_REPLAY',
                                  source='CAUSAL_REPLAY', roi_xywh=self.roi,
                                  source_time_semantics='recorded PTS; paced on host, not exposure time',
                                  arrival_mono=arrived,
                                  age_upper_bound_at_arrival_ms=1000*(arrived-release),
                                  scheduled_release_mono=release)
                    with self.lock:
                        self.latest = (header, image)
                    self.headers.append(header)
        except Exception as exc:
            self.error = repr(exc)
        finally:
            self.done.set()

    def snapshot(self):
        with self.lock:
            return self.latest


def main():
    study = ROOT / 'results/supervised_20261016_v1'
    out = study / 'paced_replay'
    out.mkdir(exist_ok=False)
    # Inherited original saved-video ROI, registered before either timed route.
    roi = [650, 350, 180, 350]
    configs = {r: json.loads((study/f'selected_{r}.json').read_text()) for r in ('opencv','cellpose')}
    (out/'preregistered.json').write_text(json.dumps(dict(
        source='CAUSAL_REPLAY', frames=[360,659], roi_xywh=roi, configs=configs,
        source_sha256=sha(VIDEO), source_rate='native PTS', queue='latest slot, no FIFO backlog',
        freshness_ms=100, no_hardware=True, execution=identity(ROOT)), indent=2))
    summaries = {}
    for route, config in configs.items():
        engine = VisionEngine(route, origin=roi[:2], detector_config=config)
        warm = cv2.imread(str(study/'replay/images/000360.png'))[170:520, :180]
        # Warmup detection only; tracker reset and warmup excluded from timing.
        warm_t = time.monotonic()
        engine.detect(warm)
        warm_s = time.monotonic()-warm_t
        engine.reset('warmup_complete')
        monitor, mock = ResourceMonitor(), MockConsumer()
        source = PacedReplay(roi)
        rows, last = [], -1
        source.thread.start()
        while True:
            item = source.snapshot()
            if item is None or item[0]['seq'] == last:
                if source.done.is_set():
                    break
                time.sleep(.001)
                continue
            header, frame = item
            last = header['seq']
            begin = time.monotonic()
            msg, timing, _ = engine.process(frame, header)
            count = len(msg['candidates'])
            fresh = enforce_freshness(msg)
            if not fresh:
                engine.reset('stale_result')
            timing['serialize_send_ms'] = mock.send(msg)
            timing['software_ms'] = 1000*(time.monotonic()-begin)
            rows.append(dict(frame_id=last, source_pts_s=header['source_time_s'],
                             header=header, candidates=count, fresh_at_send=fresh,
                             age_at_send_ms=msg['age_upper_bound_ms'], timing=timing))
        end = time.monotonic()
        mock.send(invalid_state('PACED_CAUSAL_REPLAY',last,'REPLAY_ENDED'))
        mock.close(); monitor.close(); source.thread.join(2)
        if source.error:
            raise RuntimeError(source.error)
        span = end-source.headers[0]['scheduled_release_mono']
        valid = [r for r in mock.rows if r['candidate_valid']]
        summary = dict(source='CAUSAL_REPLAY', route=route, roi_xywh=roi, warmup_s=warm_s,
                       native_source_hz=30, source_frames=len(source.headers), processed_frames=len(rows),
                       latest_slot_overwrites=len(source.headers)-len(rows), measured_wall_s=span,
                       processed_hz=len(rows)/span, fresh_candidate_hz=len(valid)/span,
                       fresh_candidate_states=len(valid), validated_bubble_hz=0, validated_group_hz=0,
                       stale_result_rejected=sum(not r['fresh_at_send'] for r in rows),
                       stage_ms={k:stats([r['timing'][k] for r in rows]) for k in rows[0]['timing']},
                       mock_transport_ms=stats([r['transport_ms'] for r in mock.rows]),
                       mock_consume_ms=stats([r['consume_ms'] for r in mock.rows]),
                       mock_age_ms=stats([r['age_upper_bound_ms'] for r in mock.rows if r['new_frame']]),
                       mock_send_drops=mock.send_drops, mock_received=len(mock.rows),
                       resources={k:stats([r[k] for r in monitor.rows if k in r]) for k in
                                  ['cpu_pct_one_core_basis','rss_mb','gpu_device_util_pct','gpu_device_memory_mb']},
                       resource_scope='WSL process CPU/RAM; device-wide GPU includes other workloads',
                       torch_peak_allocated_mb=engine.torch.cuda.max_memory_allocated()/1e6 if engine.torch else 0,
                       CAMERA_TO_CONTROLLER_END_TO_END_LATENCY='NOT_MEASURED', HARDWARE_CONTROL_DISABLED=True,
                       tracking_accuracy='NOT_EVALUATED_NO_HUMAN_TRACK_IDS', reset_events=engine.reset_events)
        summaries[route] = summary
        for name, data in [('summary',summary),('resources',monitor.rows),('mock',mock.rows),('source',source.headers)]:
            (out/f'{route}_{name}.json').write_text(json.dumps(data,indent=2))
        (out/f'{route}_frames.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
        print(json.dumps({k:summary[k] for k in ['route','processed_frames','processed_hz','fresh_candidate_hz','stale_result_rejected']}),flush=True)
    (out/'summary.json').write_text(json.dumps(summaries,indent=2))


if __name__ == '__main__':
    main()
