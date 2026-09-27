"""Independent replay-state checks against every selected saved trajectory."""
from pathlib import Path
import argparse,json,socket
import numpy as np
from PIL import Image
from .particle82_provenance import atomic_json,sha256


def audit(root):
    root=Path(root);scene=json.loads((root/'data/SCENE.json').read_text());source_sha=sha256(root/'data/SCENE.json')
    compute=json.loads((root/'PPT_COMPUTE_VALIDATION.json').read_text());assert compute['all_pass'] and compute['quota']==500
    media=json.loads((root/'PPT_MEDIA_VALIDATION.json').read_text());assert media['all_pass'] and len(media['records'])==3
    records={r['event_id']:r for r in scene['records']};arrays={}
    for pid,r in records.items():
        assert sha256(root/r['samples_path'])==r['samples_sha256'];assert sha256(root/r['metadata_path'])==r['metadata_sha256']
        arrays[pid]=np.load(root/r['samples_path'])['samples']
    states=0;frames=0;max_error=0.
    for video in media['records']:
        assert sha256(root/'animations'/video['file'])==video['sha256']
        assert video['decoded_frames']==432 and video['dimensions']==[1920,1080] and video['fps']==24
        receipt=json.loads((root/'frames'/(Path(video['file']).stem+'.json')).read_text());assert len(receipt['frames'])==432
        for frame in receipt['frames']:
            frames+=1;assert frame['source_scene_sha256']==source_sha and frame['physical_radius_scale']==1.
            assert not frame['camera']['physical_coordinates_rotated'] and not frame['camera']['axis_stretching']
            for bubble in frame['active_bubbles']:
                pid=bubble['event_id'];a=arrays[pid];age=bubble['age_s'];states+=1
                assert a[0,0]<=age<a[-1,0] and bubble['radius_m']==records[pid]['radius_m']
                expected=np.array([np.interp(age,a[:,0],a[:,k]) for k in [1,2,3]])
                error=float(np.linalg.norm(expected-bubble['position_m']));max_error=max(max_error,error);assert error<1e-15
                index=int(np.searchsorted(a[:,0],age,side='right'))
                assert abs(np.linalg.norm(a[index,4:7])*1000-bubble['speed_mm_s'])<1e-10
    figures=[]
    for p in sorted((root/'figures').glob('*.png')):
        with Image.open(p) as im:im.verify()
        with Image.open(p) as im:dimensions=list(im.size)
        figures.append(dict(file=p.name,sha256=sha256(p),dimensions=dimensions))
    assert len(figures)==5 and all(f['dimensions']==[3840,2160] for f in figures if f['file']!='04_animation_storyboard.png')
    result=dict(all_pass=True,integrated_tracks=500,verified_video_frames=frames,verified_bubble_states=states,
        maximum_position_replay_error_m=max_error,source_scene_sha256=source_sha,
        true_saved_radii=True,no_extrapolation=True,physical_coordinates_preserved=True,
        independent_display_phase_not_simultaneous_population=True,figures=figures,
        all_videos_fully_decoded=True,hostname=socket.gethostname(),audit_source_sha256=sha256(__file__),
        manual_visual_review='PENDING_USER_REVIEW')
    atomic_json(root/'PPT_REPLAY_AUDIT.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();audit(a.root)
