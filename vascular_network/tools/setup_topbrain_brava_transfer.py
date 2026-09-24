#!/usr/bin/env python3
"""Install only the pinned binary transfer tools after a resolver safety check."""
import argparse
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import venv
import ctypes
import hashlib
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
PINS = {'open3d':'0.20.0', 'POT':'0.9.7.post1'}
CORE = ['numpy','scipy','vtk','pyvista','nibabel','scikit-image','networkx','matplotlib','pandas']


def run(command, **kw):
    print(' '.join(map(str,command)),flush=True)
    kw.setdefault('cwd',ROOT)
    result = subprocess.run(list(map(str,command)),**kw)
    if result.returncode and kw.get('capture_output'): print(result.stderr, file=sys.stderr)
    result.check_returncode()
    return result


def installed(python):
    code = 'import importlib.metadata as m,json;print(json.dumps({d.metadata["Name"].lower():d.version for d in m.distributions()}))'
    return json.loads(run([python,'-c',code],capture_output=True,text=True).stdout)


def ensure_usb_runtime():
    try:
        ctypes.CDLL('libusb-1.0.so.0')
        return
    except OSError:
        pass
    runtime=ROOT/'outputs/topbrain_brava_transfer/runtime'
    if (runtime/'lib/libusb-1.0.so.0').exists():return
    archive=runtime.parent/'libusb.conda'
    url='https://api.anaconda.org/download/conda-forge/libusb/1.0.29/linux-64/libusb-1.0.29-h73b1eb8_0.conda'
    digest='89c84f5b26028a9d0f5c4014330703e7dff73ba0c98f90103e9cef6b43a5323c'
    data=urllib.request.urlopen(url,timeout=60).read()
    if hashlib.sha256(data).hexdigest()!=digest:raise RuntimeError('libusb binary checksum mismatch')
    archive.write_bytes(data)
    extractor=Path(sys.base_prefix)/'bin/python'
    run([extractor,'-c','from conda_package_handling.api import extract;import sys;extract(sys.argv[1],dest_dir=sys.argv[2])',archive,runtime])
    (runtime.parent/'libusb_runtime_source.json').write_text(json.dumps(dict(url=url,sha256=digest,version='1.0.29',scope='project-local binary only'),indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,default=ROOT/'outputs/topbrain_brava_transfer')
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    python=Path(sys.executable);before=installed(python);isolated=ROOT/'.topbrain_brava_transfer_env'
    baseline_path=a.output_dir/'main_environment_before.json'
    if not baseline_path.exists():baseline_path.write_text(json.dumps(dict(python=sys.version,executable=sys.executable,packages=before),indent=2))
    baseline={k.lower():v for k,v in json.loads(baseline_path.read_text())['packages'].items()}
    if (isolated/'bin/python').exists():
        available=installed(isolated/'bin/python')
        if all(available.get(k.lower())==v for k,v in PINS.items()):python=isolated/'bin/python'
    current=installed(python)
    wanted=[k+'=='+v for k,v in PINS.items() if current.get(k.lower())!=v]
    if wanted:
        constraint=a.output_dir/'environment_constraints.txt'
        constraint.write_text(''.join(k+'=='+before[k.lower()]+'\n' for k in CORE if k.lower() in before))
        report=a.output_dir/'pip_dry_run.json'
        command=[python,'-m','pip','install','--index-url','https://pypi.org/simple','--only-binary=:all:',*wanted,'--dry-run','--report',report,'-c',constraint]
        try:
            run(command)
            plan=json.loads(report.read_text())['install']
            changed=[x['metadata']['name'] for x in plan if x['metadata']['name'].lower() in {n.lower() for n in CORE} and before.get(x['metadata']['name'].lower())!=x['metadata']['version']]
            if changed:raise RuntimeError('Core dependency replacement: '+str(changed))
        except (subprocess.CalledProcessError,RuntimeError):
            if not isolated.exists():venv.EnvBuilder(with_pip=True,system_site_packages=True).create(isolated)
            python=isolated/'bin/python'
            command=[python,'-m','pip','install','--index-url','https://pypi.org/simple','--only-binary=:all:',*wanted,'--dry-run','--report',report]
            run(command)
        command=[x for x in command if x!='--dry-run']
        run(command)
    ensure_usb_runtime()
    code='''import sys,json,inspect,importlib.metadata as m,hashlib
from vascular_processing.transfer_runtime import prepare_runtime
prepare_runtime()
import open3d as o3d,ot
f=ot.gromov.partial_fused_gromov_wasserstein
print(json.dumps(dict(python=sys.version,executable=sys.executable,open3d=m.version('open3d'),POT=m.version('POT'),pot_signature=str(inspect.signature(f)),pot_source_sha256=hashlib.sha256(inspect.getsource(f).encode()).hexdigest(),icp_scaling=o3d.pipelines.registration.TransformationEstimationPointToPoint(with_scaling=True).with_scaling),indent=2))'''
    verification=json.loads(run([python,'-c',code],capture_output=True,text=True).stdout)
    assert all(verification[k]==v for k,v in PINS.items())
    after=installed(Path(sys.executable))
    verification.update(isolated_env=python==isolated/'bin/python',main_core_unchanged=all(baseline.get(k.lower())==after.get(k.lower()) for k in CORE),
                        main_existing_packages_changed={k:[v,after.get(k)] for k,v in baseline.items() if after.get(k)!=v},
                        main_added_packages={k:v for k,v in after.items() if k not in baseline})
    library=ROOT/'outputs/topbrain_brava_transfer/runtime/lib/libusb-1.0.so.0'
    if library.exists():verification.update(runtime_library=str(library.relative_to(ROOT)),runtime_library_sha256=hashlib.sha256(library.read_bytes()).hexdigest())
    (a.output_dir/'environment_summary.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(verification,indent=2))
    return 0

if __name__=='__main__':raise SystemExit(main())
