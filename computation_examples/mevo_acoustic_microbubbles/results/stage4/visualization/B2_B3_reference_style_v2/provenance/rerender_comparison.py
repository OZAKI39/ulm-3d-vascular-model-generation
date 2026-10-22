"""Correct comparison label margins; retain verified, pixel-identical singles."""
import json
from pathlib import Path
import shutil
import time
import numpy as np
import imageio.v2 as imageio
from scripts.stage4 import render_motion as r

source = r.ROOT / 'results/stage4/visualization/B2_B3_reference_style_v1'
out = r.ROOT / 'results/stage4/visualization/B2_B3_reference_style_v2'
shutil.copytree(source, out, ignore=shutil.ignore_patterns('stage4_B2_B3_motion_3d_60fps.mp4', 'render.log', 'DELIVERY_SHA256.json', 'DELIVERY_VERIFICATION.json'))
shutil.move(out / 'provenance/stage4_render_motion.py', out / 'provenance/stage4_render_motion_v1.py')
shutil.copy2(r.__file__, out / 'provenance/stage4_render_motion.py')
shutil.copy2(__file__, out / 'provenance/rerender_comparison.py')
font = Path('/mnt/c/Windows/Fonts/arial.ttf')
r.font_manager.fontManager.addfont(str(font))
r.plt.rcParams.update({'font.family': 'Arial', 'font.weight': 'normal', 'font.size': 14,
                      'axes.labelweight': 'normal', 'axes.titleweight': 'normal'})
times = np.r_[np.zeros(30), np.arange(1201) * .005, np.full(90, 6.)]
data, protected = r.load_data(times)
protected[str(font)] = r.sha(font)
wall = r.wall_sections(data['B2']['surface'])
manifest = json.loads((source / 'RENDER_MANIFEST.json').read_text())
old_renderer = manifest['renderer_sha256']
checks = {}
started = time.perf_counter()
for names in [['B2'], ['B3'], ['B2', 'B3']]:
    stem = 'stage4_' + '_'.join(names) + '_motion_3d_60fps'
    fig, render, checks[stem] = r.make_scene(data, names, wall)
    for frame in [30, 230, 630, 1230]:
        raw = render(frame, times[frame])
        key = out / 'keyframes' / f'{stem}_{frame:04d}.png'
        if len(names) == 1:
            assert np.array_equal(raw, imageio.imread(key)), 'Reused individual pixels changed'
        else:
            imageio.imwrite(key, raw)
    if len(names) == 1:
        assert r.sha(out / f'{stem}.mp4') == manifest['videos'][stem]['sha256']
        manifest['videos'][stem]['renderer_sha256'] = old_renderer
        manifest['videos'][stem]['reuse'] = 'Byte-identical v1 video; all four regenerated keyframes are pixel-identical with the current renderer'
        print(json.dumps(dict(video=stem, verified_reuse=True)), flush=True)
        r.plt.close(fig)
        continue
    with imageio.get_writer(out / f'{stem}.mp4', fps=r.FPS, codec='libx264', quality=8,
                            macro_block_size=1, pixelformat='yuv420p',
                            ffmpeg_params=['-movflags', '+faststart']) as writer:
        for frame, now in enumerate(times):
            writer.append_data(render(frame, now))
            if frame % 180 == 0 or frame == len(times) - 1:
                print(json.dumps(dict(video=stem, frame=frame, total=len(times), elapsed_s=time.perf_counter()-started)), flush=True)
    r.plt.close(fig)
    reader = imageio.get_reader(out / f'{stem}.mp4')
    assert reader.count_frames() == 1321 and reader.get_meta_data()['fps'] == 60
    reader.close()
    manifest['videos'][stem].update(sha256=r.sha(out / f'{stem}.mp4'), renderer_sha256=r.sha(r.__file__))
assert all(r.sha(path) == digest for path, digest in protected.items())
manifest.update(renderer_sha256=r.sha(r.__file__), screen_checks=checks,
                reused_individual_videos_source_manifest_sha256=r.sha(source / 'RENDER_MANIFEST.json'),
                layout_correction='Comparison Z labels moved within the frame; camera, data, glyphs, trajectory timing and all individual-video pixels unchanged',
                elapsed_s=time.perf_counter()-started)
(out / 'RENDER_MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
(out / 'PREVIEW_CHECKS.json').write_text(json.dumps(checks, indent=2)+'\n')
r.write_viewer(out)
reproduce = (out / 'REPRODUCE.txt').read_text().replace('B2_B3_reference_style_v2', 'B2_B3_reference_style_v3').replace('e7e3d2f', 'adf0c63')
reproduce += '\nV2：修正对照版 Z 轴标签边界。单组视频沿用 V1 原文件，4 个关键帧重新绘制后逐像素一致；仅对照版重新编码。复现命令完整重绘三段。\n'
(out / 'REPRODUCE.txt').write_text(reproduce)
verifier = (out / 'provenance/verify_delivery.py').read_text().replace('B2_B3_reference_style_v1', 'B2_B3_reference_style_v2').replace("renderer_commit='e7e3d2f'", "renderer_commit='adf0c63'")
(out / 'provenance/verify_delivery.py').write_text(verifier)
print(json.dumps(dict(status='PASS', output=str(out), elapsed_s=manifest['elapsed_s'])), flush=True)
