import json
import imageio_ffmpeg
from PIL import Image
from particle_3d.particle8_replay import digest

def test_actual_exports_decode_and_bind_to_reviewed_scene(saved_scene):
    s=saved_scene;audit=json.loads((s.root/'data/media_audit.json').read_text())
    assert audit['all_pass'] and len(audit['animations'])==4
    for r in audit['animations']:
        stem='particle8_2_anim_'+r['kind'];mp4=s.root/'animations'/(stem+'.mp4');gif=s.root/'animations'/(stem+'.gif')
        assert digest(mp4)==r['mp4_sha256'] and digest(gif)==r['gif_sha256']
        assert r['source_scene_sha256']==s.sha256 and r['all_pass']
        reader=imageio_ffmpeg.read_frames(str(mp4));meta=next(reader)
        assert sum(1 for _ in reader)==r['decoded_mp4_frames']
        image=Image.open(gif)
        for i in range(image.n_frames):image.seek(i);image.load()
        assert image.n_frames==r['decoded_gif_frames']
