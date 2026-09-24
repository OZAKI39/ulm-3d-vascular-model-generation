import json
import numpy as np
import pytest
from PIL import Image
from particle_3d.particle8_replay import read,digest,snapshot,canonical_hash
from particle_3d.particle8_visuals import FIGURES,ANIMATIONS,ReplayView


@pytest.mark.parametrize('index',range(1,6))
def test_encoded_mp4_full_decode_and_motion(output,index):
    import imageio_ffmpeg
    name=ANIMATIONS[index-1];manifest=read(output/'data'/(name+'_frames.json'))
    reader=imageio_ffmpeg.read_frames(str(output/'animations'/(name+'.mp4')),pix_fmt='rgb24')
    meta=next(reader);assert meta['size']==(1440,900) and meta['fps']==15
    assert 'h264' in meta['codec'];n=0;hashes=set()
    for frame in reader:
        import hashlib
        hashes.add(hashlib.sha256(frame).hexdigest());n+=1
    assert n==manifest['frame_count'] and len(hashes)>30
    gif=Image.open(output/'animations'/(name+'.gif'));assert gif.n_frames>=25


@pytest.mark.parametrize('index',range(1,6))
def test_every_frame_is_replayable_without_simulation(output,scenes,index):
    name=ANIMATIONS[index-1];d=read(output/'data'/(name+'_frames.json'));prev=-1
    key={1:'real_single_mb',2:'real_mixed',3:'synthetic',4:'synthetic'}
    for i,f in enumerate(d['frames']):
        assert f['frame_index']==i and f['video_time_s']==i/15
        t=f['physical_or_display_time_s'];assert t>=prev;prev=t
        if index<=4:
            expected=snapshot(scenes[key[index]],t)
            for k,v in expected.items(): assert f['state'][k]==v
        else:
            assert f['state']['diagnostic_proposals_not_admitted_births']
            assert f['clock_role']=='DISPLAY_CLOCK_NOT_PHYSICAL_TIME'
        assert canonical_hash(f['state'])==f['state_sha256']
    if index==3:
        assert any(s['state']['counts']['MB']['active'] and s['state']['counts']['RBC']['active'] for s in d['frames'])
        assert len(set(s['window_label'] for s in d['frames']))==4
        # Moving centers must be resolved across consecutive frames, not teleporting dots.
        movements=0
        for a,b in zip(d['frames'],d['frames'][1:]):
            old={p['particle_id']:p for p in a['state']['active']}
            for p in b['state']['active']:
                if p['particle_id'] in old and p['position_m']!=old[p['particle_id']]['position_m']:movements+=1
        assert movements>300
    if index==4: assert all(f['state']['mismatch_count']==0 and f['state']['continuous_state_sha256']==f['state']['restarted_state_sha256'] for f in d['frames'])


def test_static_pngs_and_source_maps(output):
    for name in FIGURES:
        im=Image.open(output/'figures'/(name+'.png'));assert im.width>=1800;im.verify()
        p=read(output/'data'/(name+'_plot.json'))
        assert 'RBC PASSAGE NOT ESTABLISHED' in ' '.join(p['labels'])
        assert 'MODEL ASSUMPTION' in ' '.join(p['labels'])
        for path,h in p['source_sha256'].items():assert digest(output/path)==h


def test_keyframe_reproduces_exact_pixels(output):
    import matplotlib.pyplot as plt
    idx=3;name=ANIMATIONS[idx-1];m=read(output/'data'/(name+'_frames.json'));k=m['keyframe_indices'][2]
    f=m['frames'][k];v=ReplayView(output,idx);v.draw(f['physical_or_display_time_s'],f['window_label']);v.fig.canvas.draw()
    actual=np.asarray(v.fig.canvas.buffer_rgba())[:,:,:3]
    saved=np.asarray(Image.open(output/'frames'/f'{name}_frame_{k:04d}.png'))
    np.testing.assert_array_equal(actual,saved);plt.close(v.fig)


def test_density_uses_true_triangle_area_without_numpy_cross_2d():
    from particle_3d.particle8_visuals import triangle_area_2d
    triangles=np.array([[[0,0],[2,0],[0,3]],[[0,3],[2,0],[0,0]]],float)
    np.testing.assert_array_equal(triangle_area_2d(triangles),[3,3])
