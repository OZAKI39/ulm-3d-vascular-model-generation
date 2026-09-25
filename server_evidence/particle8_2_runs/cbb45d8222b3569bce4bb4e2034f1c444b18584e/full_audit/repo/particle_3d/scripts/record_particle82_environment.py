#!/usr/bin/env python3
from pathlib import Path
import argparse,importlib.metadata,json,platform,subprocess,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_provenance import atomic_json,collect_host,require_remote

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--host-provenance',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    host=json.loads(Path(a.host_provenance).read_text());require_remote(host,host['hostname'])
    packages={}
    for name in ['numpy','scipy','vtk','pyvista','matplotlib','lammps','Pillow','pytest','psutil','imageio-ffmpeg']:
        try:packages[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:packages[name]='NOT_INSTALLED'
    import imageio_ffmpeg
    ffmpeg=subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-version'],capture_output=True,text=True)
    freeze=subprocess.run(['uv','pip','freeze','--python',sys.executable],capture_output=True,text=True)
    compiler=subprocess.run(['g++','--version'],capture_output=True,text=True)
    atomic_json(a.output,dict(hostname=host['hostname'],source_git_commit=host['source_git_commit'],
        python=sys.version,executable=sys.executable,platform=platform.platform(),packages=packages,
        pip_freeze_command=['uv','pip','freeze','--python',sys.executable],pip_freeze=freeze.stdout,
        ffmpeg_executable=imageio_ffmpeg.get_ffmpeg_exe(),ffmpeg_version=ffmpeg.stdout,
        native_point_compiler=compiler.stdout,native_point_flags=['-O3','-fno-fast-math','-ffp-contract=off'],
        solver_gpu_usage='GPU_AVAILABLE_BUT_NOT_USED_BY_CURRENT_SOLVER',
        source_and_frozen_inputs_recorded_in='host_provenance.json'))
