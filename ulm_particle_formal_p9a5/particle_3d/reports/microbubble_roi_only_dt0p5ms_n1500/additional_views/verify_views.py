"""Check delivered videos, unchanged source results and camera/style consistency."""
from pathlib import Path
import hashlib,json
import numpy as np
import imageio_ffmpeg
from PIL import Image

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def main():
    manifest=json.loads((HERE/'RENDER_MANIFEST.json').read_text())
    lock=json.loads((HERE/'audit/SOURCE_LOCK.json').read_text())
    original=json.loads((ROOT/'data/render_manifest.json').read_text())
    for path,digest in lock['protected_files'].items():assert sha(ROOT/path)==digest,path
    assert manifest['particle_integrations']==0 and manifest['CFD_calls']==0
    assert manifest['source_tracks_unchanged'] and manifest['original_movie_unchanged']
    assert manifest['source_track_files_verified']==3000
    assert manifest['original_renderer_sha256']==sha(ROOT/'scripts/render_results.py')
    assert len(manifest['views'])==3
    ages=np.linspace(0,original['age_end_s'],original['video_frames'])
    expected_age_sha=hashlib.sha256(ages.tobytes()).hexdigest()
    original_direction=np.array([-.72,.69,.11]);original_direction/=np.linalg.norm(original_direction)
    directions=[original_direction]
    result=[]
    for view in manifest['views']:
        path=HERE/view['video'];assert sha(path)==view['sha256']
        assert view['text_language']=='English' and view['all_tracks_used']
        assert view['zoom']==original['animation_zoom_factor']==1.3
        assert view['marker_size_px']==9.1
        assert view['age_grid_sha256']==expected_age_sha
        assert sum(view['count_by_group'].values())==1500
        assert view['layout']['parallel_scale']==original['animation_layout']['parallel_scale']
        assert np.allclose(view['layout']['projected_size_ratios'],1.3,atol=1e-10,rtol=0)
        bounds=np.asarray(view['layout']['enlarged_bounds_px'])
        assert np.all(bounds[0]>=[40,108]) and np.all(bounds[1]<=[1880,1008])
        direction=np.asarray(view['camera_direction']);assert np.isclose(np.linalg.norm(direction),1)
        for prior in directions:assert not np.allclose(direction,prior)
        directions.append(direction)
        assert any('NVIDIA GeForce RTX 4090' in line for line in view['render_backend'])
        decoder=imageio_ffmpeg.read_frames(str(path),pix_fmt='rgb24')
        metadata=next(decoder);metadata['decoded_frames']=sum(1 for _ in decoder)
        assert metadata['decoded_frames']==view['frames']==288
        assert tuple(metadata['size'])==(1920,1080) and metadata['fps']==24
        assert metadata['duration']==12
        with Image.open(HERE/'figures'/f"{view['name']}_preview.png") as im:
            assert im.size==(1920,1080) and min(im.info['dpi'])>299
        result.append(dict(name=view['name'],path=view['video'],sha256=sha(path),video=metadata))
    report=dict(PASS=True,videos=result,original_delivery_files_unchanged=len(lock['protected_files']),
                same_time_axis=True,same_scale_and_styles=True,distinct_camera_directions=True,
                projected_surface_within_content_rectangle=True,source_integrity=manifest['source_tracks_unchanged'],
                GPU=manifest['GPU'],encoder=manifest['encoder'])
    (HERE/'VALIDATION.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)

if __name__=='__main__':main()
