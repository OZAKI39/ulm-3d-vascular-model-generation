import hashlib
import numpy as np
import pytest
import imageio_ffmpeg
from PIL import Image
from particle_3d.particle8_replay import read,canonical_hash
from particle_3d.particle8_full3d_data import CASES,replay_frame,camera_angles


@pytest.mark.parametrize('case',list(CASES))
def test_all_mp4_and_gif_frames_decode(root,case):
    name='particle8_full3d_'+case;m=read(root/'frames'/(name+'_manifest.json'))
    reader=imageio_ffmpeg.read_frames(str(root/'mp4'/(name+'.mp4')),pix_fmt='rgb24');meta=next(reader)
    assert meta['size']==(1600,1000) and meta['fps']==15 and 'h264' in meta['codec']
    hashes=set();n=0
    for frame in reader:hashes.add(hashlib.sha256(frame).digest());n+=1
    assert n==m['frame_count'] and len(hashes)>40
    im=Image.open(root/'gif'/(name+'.gif'))
    for i in range(im.n_frames):im.seek(i);im.convert('RGB').load()
    assert im.n_frames>=30
    for i in [0,n//2,n-1]:
        with Image.open(root/'keyframes'/f'{name}_frame_{i:04d}.png') as image:assert image.size==(1600,1000)


@pytest.mark.parametrize('case',list(CASES))
def test_every_rendered_frame_and_tail_matches_p8(root,scenes,case):
    name='particle8_full3d_'+case;m=read(root/'frames'/(name+'_manifest.json'));scene=scenes[CASES[case]['scene']]
    for i,f in enumerate(m['frames']):
        expected=replay_frame(scene,f['physical_time_s'],CASES[case]['tail_window_s'])
        assert f['state']==expected and f['state_sha256']==canonical_hash(expected)
        ids=[p['particle_id'] for p in expected['particles']]
        assert f['rendered_ids_main']==f['rendered_ids_secondary']==ids
        assert not set(ids)&set(expected['pending_ids']) and not expected['pending_drawn_ids']
        assert f['coordinates_rotated'] is False
        assert f['camera_main']['azimuth_deg']==camera_angles(case,i,m['frame_count'])['azimuth_deg']
        assert f['classification']==CASES[case]['classification']
        if scene['source_classification']=='real_frozen':
            assert 'open section' not in f['context_note'] and 'straddle' not in f['context_note']
        labels=' '.join(f['labels']);assert 'DISPLAY-ONLY INTERPOLATION' in labels and 'Real RBC passage NOT established' in labels
        assert 'neighbor settings NOT FROZEN' in labels and 'No CFD' in labels
        if case.startswith('D_'):assert f['restart_mismatches']==0 and f['restarted_state_sha256']==f['state_sha256']


def test_fixed_and_rotating_views_have_identical_physical_replay(root):
    a=read(root/'frames/particle8_full3d_A_real_single_mb_orbit_manifest.json')
    b=read(root/'frames/particle8_full3d_A_real_single_mb_fixed_manifest.json')
    assert [f['state_sha256'] for f in a['frames']]==[f['state_sha256'] for f in b['frames']]
    assert a['frames'][-1]['camera_main']!=b['frames'][-1]['camera_main']


def test_synthetic_has_simultaneous_moving_rbc_mb_and_deletion(root):
    m=read(root/'frames/particle8_full3d_C_synthetic_lifecycle_orbit_manifest.json');frames=m['frames']
    assert any({p['species'] for p in f['state']['particles']}=={'MB','RBC'} for f in frames)
    c=frames[-1]['state']['counts'];assert c['MB']['deleted']==2 and c['RBC']['deleted']==2352
    for sp,r in c.items():assert r['scheduled']==r['pending']+r['active']+r['deleted']
    moved=0
    for a,b in zip(frames,frames[1:]):
        old={p['particle_id']:p for p in a['state']['particles']}
        moved+=sum(p['particle_id'] in old and p['position_m']!=old[p['particle_id']]['position_m'] for p in b['state']['particles'])
    assert moved>250


def test_storyboard_and_decode_receipts(root):
    m=read(root/'data/media_audit.json')
    assert m['all_pass'] and len(m['animations'])==5
    for p in ['storyboard/full3d_storyboard.png','storyboard/scientific_scope_summary.png']:
        with Image.open(root/p) as im:im.verify()
