"""Decode every video frame; audit frozen inputs, device files, artifacts and metadata."""
from pathlib import Path
import sys,json,hashlib,subprocess,ast
import av,cv2,numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from annotation_store import sha
from runtime import identity


def main():
    tag='supervised_20261016_v1';out=ROOT/'results'/tag
    checks=out/'visual_check';checks.mkdir(exist_ok=True)
    source=[json.loads(s)['source_pts_s'] for s in (out/'replay/records.jsonl').read_text().splitlines()
            if json.loads(s)['route']=='opencv']
    videos={};samples=[]
    for video in sorted((ROOT/'videos'/tag).glob('*.mp4')):
        times=[]
        with av.open(str(video)) as container:
            stream=container.streams.video[0]
            metadata=dict(codec=stream.codec_context.name,width=stream.width,height=stream.height,
                          average_rate=str(stream.average_rate),time_base=str(stream.time_base),
                          duration_s=float(container.duration/av.time_base),library='PyAV '+av.__version__)
            for i,frame in enumerate(container.decode(stream)):
                times.append(float(frame.time))
                if i in (0,30,59):
                    image=frame.to_ndarray(format='bgr24');path=checks/f'{video.stem}_{i:03d}.png'
                    cv2.imwrite(str(path),image);samples.append((video.stem,i,image))
        assert metadata['codec']=='h264' and (metadata['width'],metadata['height'])==(1280,960)
        assert len(times)==60 and all(b>a for a,b in zip(times,times[1:]))
        error=max(abs((a-source[0])-b) for a,b in zip(source,times))
        assert error<1/90000
        videos[video.name]=dict(metadata=metadata,decoded_frames=len(times),maximum_pts_error_s=error,
                                source_pts_preserved=True,sha256=sha(video))
    (out/'video_validation.json').write_text(json.dumps(videos,indent=2))
    # Inspection sheet is derived from actual encoded video, not render inputs.
    thumbs=[cv2.resize(a,(640,480)) for _,_,a in samples]
    sheet=np.vstack([np.hstack(thumbs[:3]),np.hstack(thumbs[3:])]);cv2.imwrite(str(checks/'encoded_first_middle_last.png'),sheet)
    freeze=ROOT/'results/expert_freezes/vision2_supervised_20261016_v1'
    frozen=json.loads((freeze/'frozen.json').read_text())
    frozen_check={p:sha(freeze/p)==v for p,v in frozen['snapshot_sha256'].items()}
    assert all(frozen_check.values())
    original=json.loads((out/'protected_sources_before.json').read_text())
    protection={p:sha(p)==v for p,v in original.items()}
    assert all(protection.values())
    (out/'integrity.json').write_text(json.dumps(dict(frozen_files=frozen_check,protected_device_files=protection,
        protected_scope='36 original top-level Python/config files recorded before this study; not a whole-disk audit'),indent=2))
    # Preserve the executed source of the finite study and original paced benchmark.
    snapshots=out/'source_snapshots';snapshots.mkdir(exist_ok=True)
    studies={'static':json.loads((out/'preregistered.json').read_text()),
             'ordered_replay':json.loads((out/'replay/preregistered.json').read_text()),
             'paced_replay':json.loads((out/'paced_replay/preregistered.json').read_text())}
    needed={'static':['scripts/run_supervised_study.py','src/vision_engine.py','src/vision_evaluation.py','src/runtime.py',
                      'src/annotation_store.py','src/annotation_regions.py'],
            'ordered_replay':['scripts/replay_supervised.py','scripts/prepare_candidates.py','src/vision_engine.py','src/runtime.py'],
            'paced_replay':['scripts/benchmark_supervised_replay.py','scripts/prepare_candidates.py','src/vision_engine.py','src/runtime.py']}
    source_check={}
    for run,prereg in studies.items():
        for name in needed[run]:
            content=(ROOT/name).read_bytes();expected=prereg['execution']['source_sha256'][name]
            if name=='scripts/benchmark_supervised_replay.py':
                # Preserve exactly the version that ran; its nominal 30 Hz metadata
                # was later corrected from logged PTS without altering measurements.
                s=content.decode().replace('import sys, json, time, threading\n','import sys, json, time, threading, copy\n')
                s=s.replace("native_source_hz=(len(source.headers)-1)/(source.headers[-1]['source_time_s']-source.headers[0]['source_time_s']),\n                       source_frames=len(source.headers), processed_frames=len(rows),",
                            'native_source_hz=30, source_frames=len(source.headers), processed_frames=len(rows),')
                content=s.encode()
            digest=hashlib.sha256(content).hexdigest();assert digest==expected,(run,name)
            path=snapshots/run/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(content)
            source_check[f'{run}/{name}']=digest
    # Older shared detector dependency retained explicitly in addition to Git.
    legacy=ROOT.parent/'scripts/replay_benchmark.py'
    (snapshots/'legacy_replay_benchmark.py').write_bytes(legacy.read_bytes())
    source_check['legacy_replay_benchmark.py']=sha(legacy)
    (snapshots/'verification.json').write_text(json.dumps(source_check,indent=2))
    # Verify the new deliverable sources parse without touching active GUI bytecode.
    scripts=['run_supervised_study.py','replay_supervised.py','benchmark_supervised_replay.py',
             'render_supervised_delivery.py','write_supervised_reports.py','verify_supervised_delivery.py']
    for script in scripts:ast.parse((ROOT/'scripts'/script).read_text())
    (out/'final_execution_identity.json').write_text(json.dumps(identity(ROOT),indent=2))
    print(json.dumps(dict(videos={k:v['metadata'] for k,v in videos.items()},frozen_files_unchanged=len(frozen_check),
                         original_device_files_unchanged=len(protection),executed_source_hashes_match=len(source_check)),indent=2))


if __name__=='__main__':main()
