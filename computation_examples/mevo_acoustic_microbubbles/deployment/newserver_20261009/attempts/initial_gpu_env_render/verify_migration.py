"""Bounded migration checks; preserve delivered source and scientific results."""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
import subprocess
import sys
import time
import traceback

ROOT = Path('/workspace/mevo_acoustic_microbubbles')
PROJECT = ROOT / 'repo/computation_examples/mevo_acoustic_microbubbles'
F5 = PROJECT / 'results/fast_pzt/f5'
OUT = ROOT / 'migration/validation'


def write(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def file_integrity():
    rows = json.loads((ROOT / 'migration/inventory/SOURCE_FILES.json').read_text())
    for row in rows:
        p = ROOT / 'repo' / row['path']
        assert p.stat().st_size == row['size'], p
        assert digest(p) == row['sha256'], p
    links = []
    for p in PROJECT.rglob('*'):
        if p.is_symlink():
            assert p.exists(), p
            assert p.resolve().is_relative_to(ROOT), p
            links.append(dict(path=str(p.relative_to(PROJECT)), target=str(p.readlink())))
    return dict(files=len(rows), bytes=sum(r['size'] for r in rows), links=links,
                all_delivered_files_unchanged=True)


def environment():
    import numpy, scipy, h5py, pyvista, vtk, warp, yaml, PIL, imageio_ffmpeg
    assert sys.prefix == str(ROOT / '.venv')
    assert sys.prefix != sys.base_prefix
    assert 'include-system-site-packages = false' in (ROOT / '.venv/pyvenv.cfg').read_text()
    packages = {d.metadata['Name']: d.version for d in metadata.distributions()}
    for name, version in dict(numpy='2.2.6', scipy='1.15.3', h5py='3.16.0',
                              pyvista='0.46.5', vtk='9.5.2', **{'warp-lang':'1.18.0+cu12'}).items():
        assert metadata.version(name) == version, (name, metadata.version(name))
    warp.init()
    assert warp.is_cuda_available(), 'GPU required; CPU fallback is not a pass'
    gpu = warp.get_device('cuda:0')
    a = warp.array(numpy.arange(8, dtype=numpy.float64), dtype=warp.float64, device=gpu)
    numpy.testing.assert_array_equal(a.numpy(), numpy.arange(8, dtype=numpy.float64))
    return dict(python=sys.version, prefix=sys.prefix, base_prefix=sys.base_prefix,
                packages=packages, gpu=str(gpu), gpu_name=gpu.name,
                cuda_available=True, default_environment_modified=False)


def saved_data():
    import h5py
    import numpy as np
    from mevo_mb.flow.fields import SimVascularFlowField
    manifests = json.loads((F5 / 'INPUT_MANIFEST.json').read_text())
    for rel, row in manifests.items():
        assert digest(F5 / rel) == row['sha256'], rel
    datasets = []
    for p in sorted((PROJECT / 'results/fast_pzt').rglob('*.h5')):
        # Read every dataset, including compressed chunks. Inputs and actual runs only;
        # historical absolute provenance strings inside files remain unchanged.
        with h5py.File(p, 'r') as h:
            def visit(name, obj):
                if isinstance(obj, h5py.Dataset):
                    value = obj[()]
                    datasets.append(dict(file=str(p.relative_to(PROJECT)), dataset=name,
                                         shape=list(obj.shape), dtype=str(obj.dtype)))
            h.visititems(visit)
    with h5py.File(F5 / 'science_v1/main/trajectory.h5') as h, h5py.File(F5 / 'science_v1/branch_B/trajectory.h5') as b:
        x = h['positions_m'][:, 0]
        t = h['time_s'][:]
        np.testing.assert_array_equal(h['particle_id'][:], np.arange(32))
        np.testing.assert_allclose(t, np.arange(1201)*0.0005, atol=1e-15, rtol=0)
        np.testing.assert_array_equal(b['positions_m'][0, 0], x[600])
        with np.load(F5 / 'science_v1/branch_A_actual_suffix.npz') as a:
            np.testing.assert_array_equal(a['positions_m'], h['positions_m'][600:])
        assert np.isfinite(x).all() and h['active'][:].all() and b['active'][:].all()
        flow = SimVascularFlowField(F5 / 'inputs/result_121.vtu', spatial_method='tetra_p1')
        ids = np.arange(0, 1201, 60)
        u = flow.velocity(x[ids].reshape(-1, 3), 0).reshape(len(ids), 32, 3)
        error = float(np.linalg.norm(u-h['flow/velocity_m_s'][ids, 0], axis=-1).max())
        assert error < 2e-12
    return dict(input_hashes=len(manifests), hdf5_datasets=len(datasets),
                hdf5_files=len({r['file'] for r in datasets}), datasets=datasets,
                F5_frames=1201, F5_particles=32, archived_branch_state_exact=True,
                CPU_actual_saved_flow_points=672, CPU_flow_error_m_s=error)


def f5_checks():
    from scripts.fast_pzt_f5 import checks
    work = OUT / 'f5'
    work.mkdir(exist_ok=True)
    for name in ['science_v1', 'ACOUSTIC_RELEASE_PROTOCOL.json']:
        p = work / name
        if not p.exists():
            p.symlink_to(os.path.relpath(F5 / name, work))
    # Only output/input path bindings change. Original functions, kernels, tolerances,
    # assertions and historical source identity are left intact.
    checks.BASE = work
    checks.INPUT = F5 / 'inputs'
    label = 'checks_newserver_20261009'
    checks.main(label)
    result = json.loads((work / label / 'RESULT.json').read_text())
    assert result['status'] == 'PASS'
    return dict(label=label, checks=list(result['rows']), status=result['status'],
                source_checks_unmodified=True, original_assertions_unchanged=True,
                archived_ON_OFF_prefixes_bitwise=True)


def cpu_regression():
    done = subprocess.run([sys.executable, '-m', 'pytest', '-q',
                           'tests/test_stage4_flow.py',
                           '-o', 'cache_dir=' + str(OUT/'pytest-cache')],
                          cwd=PROJECT, text=True, capture_output=True)
    (OUT/'cpu-pytest.log').write_text(done.stdout+done.stderr)
    assert done.returncode == 0, done.stdout+done.stderr
    return dict(returncode=done.returncode, log='cpu-pytest.log')


def render_and_video():
    import h5py
    import numpy as np
    import pyvista as pv
    from PIL import Image
    from scripts.fast_pzt_f5 import style
    style.FONT = ROOT/'references/fonts/arial.ttf'
    assert style.FONT.exists()
    protocol = json.loads((F5/'ACOUSTIC_RELEASE_PROTOCOL.json').read_text())
    with h5py.File(F5/'science_v1/main/trajectory.h5') as h:
        states = h['positions_m'][[0, 600, 1200], 0]
    global_image = style.global_image(protocol['events'][0], states[0].mean(0), white=True, size=(900,740))
    global_image.save(OUT/'PZT_global_render_smoke.png')
    scene = style.local_scene(states, True, False, size=(800,650))
    try:
        style.image_state(scene, states[1]).save(OUT/'bubble_cluster_render_smoke.png')
        gl = scene[0].render_window.ReportCapabilities()
    finally:
        scene[0].close()
    for name in ['PZT_global_render_smoke.png', 'bubble_cluster_render_smoke.png']:
        with Image.open(OUT/name) as im:
            assert np.asarray(im).std() > 10, 'Empty render: '+name
    (OUT/'OPENGL.txt').write_text(gl)
    movie = F5/'visualization/Video_F5_acoustic_release_flow_60fps.mp4'
    probe = json.loads(subprocess.check_output(['ffprobe','-v','error','-select_streams','v:0',
                    '-count_frames','-show_entries',
                    'stream=width,height,r_frame_rate,nb_read_frames','-of','json',str(movie)],text=True))
    s=probe['streams'][0]
    assert s['width']==1920 and s['height']==1080
    assert s['r_frame_rate']=='60/1' and int(s['nb_read_frames'])==1201
    done=subprocess.run(['ffmpeg','-v','error','-i',str(movie),'-f','null','-'],capture_output=True,text=True)
    assert done.returncode==0 and not done.stderr, done.stderr
    # Encode a tiny test from actual rendered states, not a new scientific video.
    frames=OUT/'encode-smoke'
    frames.mkdir(exist_ok=True)
    for i in range(3):
        global_image.save(frames/f'{i:03d}.png')
    subprocess.run(['ffmpeg','-v','error','-y','-framerate','60','-i',str(frames/'%03d.png'),
                    '-c:v','libx264','-pix_fmt','yuv420p',str(OUT/'encode_smoke.mp4')],check=True)
    return dict(offscreen_render=True, reference_style_unmodified=True, font=str(style.FONT),
                saved_video=probe, full_decode=True, h264_encode=True,
                scope='Environment smoke only; original trajectories and delivered movie unchanged')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    phases={}
    for name, fn in [('integrity_before',file_integrity),('environment',environment),
                     ('saved_data',saved_data),('cpu_regression',cpu_regression),
                     ('f5_gpu_checks',f5_checks),('render_video',render_and_video),
                     ('integrity_after',file_integrity)]:
        started=time.perf_counter()
        print('BEGIN', name, flush=True)
        try:
            detail=fn()
            phases[name]=dict(status='PASS', elapsed_s=time.perf_counter()-started, detail=detail)
        except Exception:
            phases[name]=dict(status='FAIL', elapsed_s=time.perf_counter()-started,
                              traceback=traceback.format_exc())
            print(phases[name]['traceback'],flush=True)
        write(name+'.json',phases[name])
        print('END',name,phases[name]['status'],flush=True)
    status='PASS' if all(x['status']=='PASS' for x in phases.values()) else 'FAIL'
    write('VALIDATION.json',dict(status=status, phases={k:{x:v[x] for x in ('status','elapsed_s')} for k,v in phases.items()},
            source_commit='66d4cebec18555a26bc9a8035cc3082756d0e76f',
            original_science_commit='ee6b35dce83bd3583457a7e85c163b6886e59ce5'))
    return 0 if status=='PASS' else 1


if __name__=='__main__':
    raise SystemExit(main())
