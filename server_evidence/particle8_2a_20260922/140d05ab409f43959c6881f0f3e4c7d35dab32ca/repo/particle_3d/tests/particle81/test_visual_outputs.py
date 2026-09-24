import hashlib
import numpy as np
import pytest
from PIL import Image
import imageio_ffmpeg
from particle_3d.particle81_visuals import ANIMATIONS,FOOTERS,TITLE,VesselView,compose_frame
from particle_3d.particle81_replay import frame_schedule
from particle_3d.particle81_figures import FIGURES
from particle_3d.particle8_replay import read,canonical_hash


@pytest.mark.parametrize('kind',ANIMATIONS)
def test_every_mp4_gif_frame_and_saved_receipt(root,scene,kind):
    name='particle8_1_anim_'+kind;record=read(root/'frames'/(name+'_frames.json'));schedule=frame_schedule(scene,kind)
    assert record['source_scene_sha256']==scene.sha256 and record['frame_count']==len(schedule)
    reader=imageio_ffmpeg.read_frames(str(root/'animations'/(name+'.mp4')),pix_fmt='rgb24');meta=next(reader)
    assert meta['size']==(1600,1000) and meta['fps']==15
    hashes=set();count=0
    for pixels in reader:hashes.add(hashlib.sha256(pixels).digest());count+=1
    assert count==len(schedule) and len(hashes)>40
    gif=Image.open(root/'animations'/(name+'.gif'))
    for k in range(gif.n_frames):gif.seek(k);gif.convert('RGB').load()
    assert gif.n_frames>30
    for f,spec in zip(record['frames'],schedule):
        expected=scene.snapshot(spec['time_s'],spec['only_ids'])
        assert f['state']==expected and f['state_sha256']==canonical_hash(expected)
        assert f['count_scope']==('SELECTED MB ACCOUNTING' if spec['only_ids'] is not None else 'ACQUISITION ACCOUNTING')
        assert f['rendered_completed_ids']==expected['completed_ids']
        assert f['rendered_active_ids']==[p['particle_id'] for p in expected['active']]
        assert not f['physical_coordinates_rotated'] and not f['physical_axis_stretching']
        assert f['display_interpolation_only'] and f['time_compression_or_cuts_explicit']
        assert TITLE in f['labels'] and set(FOOTERS).issubset(f['labels'])
    for index in record['keyframe_indices']:
        with Image.open(root/'inspection'/f'{name}_decoded_{index:04d}.png') as image:assert image.size==(1600,1000)


@pytest.mark.parametrize('index',range(10))
def test_all_static_figures_bound_to_same_saved_scene(root,scene,index):
    name='particle8_1_fig_'+FIGURES[index]
    record=read(root/'data'/(name+'_source.json'))
    assert record['source_scene_sha256']==scene.sha256 and TITLE in record['scientific_labels']
    assert set(FOOTERS).issubset(record['scientific_labels'])
    with Image.open(root/'figures'/(name+'.png')) as image:
        assert min(image.size)>=1000;image.verify()


def test_repeat_actual_render_pixels_and_physical_state_are_identical(scene,root):
    kind='03_single_journeys';schedule=frame_schedule(scene,kind);index=min(45,len(schedule)-1)
    images=[];records=[]
    for _ in range(2):
        view=VesselView(scene);im,record=compose_frame(scene,view,kind,index,len(schedule),schedule[index]);view.close()
        images.append(np.asarray(im));records.append(record)
    np.testing.assert_array_equal(images[0],images[1]);assert records[0]==records[1]
    from particle_3d.particle81_simulation import dump
    dump(root/'data/render_determinism.json',dict(repeated_independent_render_count=2,frame_kind=kind,frame_index=index,
        raw_rgb_sha256=hashlib.sha256(images[0].tobytes()).hexdigest(),exact_pixels_equal=True,
        exact_state_equal=True,source_scene_sha256=scene.sha256,renderer_rng_used=False))
