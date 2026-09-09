"""Idempotent pinned source installation; never delete a preexisting tree."""
from pathlib import Path
import datetime
import hashlib
import json
import subprocess
import tarfile
import time
import shutil
import os
import fcntl

ROOT=Path(__file__).resolve().parents[1]
COMMIT='5a410848bd5c57d5ae1c171112e78eab4a82e650'
PAL='05712164d940a42e06afdd705249912fa0c49f14'
ARCHIVE_SHA='5a9c4f22c169ad16de259a72ef2f128b233afaf81c251e6c3baa8eb34d14838f'
URL=f'https://gitlab.com/unigespc/palabos/-/archive/{PAL}/palabos-{PAL}.tar.gz'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    for folder in ['vendor','build','downloads','metadata','logs']:(ROOT/folder).mkdir(exist_ok=True)
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    record={'commands':[],'counted_as_solver':False,'created_utc':stamp}
    log_path=ROOT/'logs'/('setup_'+stamp+'.log')
    def run(cmd,cwd=ROOT,timeout=1800):
        start=time.monotonic()
        with log_path.open('ab') as log:r=subprocess.run(cmd,cwd=cwd,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        record['commands'].append({'argv':cmd,'cwd':str(cwd),'elapsed_monotonic_s':time.monotonic()-start,'exit_code':r.returncode})
        if r.returncode:raise RuntimeError('Command failed; see '+str(log_path))
    missing=[x for x in ['git','curl','patch','make','cmake','gcc','g++','mpicxx.openmpi','h5pcc'] if not shutil.which(x)]
    if missing:
        # No apt upgrade; user can supply sudo locally, never paste a password.
        packages=set()
        mapping={'git':'git','curl':'curl','patch':'patch','make':'build-essential','cmake':'cmake','gcc':'build-essential','g++':'build-essential','mpicxx.openmpi':'libopenmpi-dev openmpi-bin','h5pcc':'libhdf5-openmpi-dev'}
        for x in missing:packages.update(mapping[x].split())
        raise RuntimeError('Missing '+str(missing)+'; run exactly: sudo apt-get install --no-install-recommends '+' '.join(sorted(packages)))
    v=ROOT/'vendor/HemoCell'
    if not v.exists():
        run(['git','clone','--no-checkout','https://github.com/UvaCsl/HemoCell.git',str(v)],timeout=180)
        run(['git','checkout','--detach',COMMIT],cwd=v)
    head=subprocess.check_output(['git','-C',str(v),'rev-parse','HEAD'],text=True).strip()
    if head!=COMMIT:raise RuntimeError('Existing HemoCell version differs; no checkout/reset/delete performed')
    dirty=subprocess.check_output(['git','-C',str(v),'diff','--name-only','HEAD'],text=True).strip()
    if dirty:raise RuntimeError('Tracked upstream files changed; inspect before reuse: '+dirty)
    archive=ROOT/'downloads'/f'palabos-{PAL}.tar.gz'
    if not archive.exists():run(['curl','--fail','--location','--retry','2','--connect-timeout','20','--max-time','180','--output',str(archive),URL],timeout=200)
    if sha(archive)!=ARCHIVE_SHA:raise RuntimeError('PALABOS_ARCHIVE_HASH_MISMATCH')
    dependency=v/'palabos';sources=ROOT/'metadata/sources.json'
    if dependency.exists():
        if not sources.exists():raise RuntimeError('Existing Palabos has no task provenance; refusing overwrite or repeated patch')
        old=json.loads(sources.read_text())
        if old['palabos_commit']!=PAL or old['archive_sha256']!=ARCHIVE_SHA:raise RuntimeError('DEPENDENCY_VERSION_MISMATCH')
        run(['patch','--batch','--dry-run','--reverse','-d',str(dependency),'-p1','-i',str(v/'patch/palabos.patch')])
    else:
        with tarfile.open(archive) as tar:tar.extractall(v,filter='data')
        (v/f'palabos-{PAL}').rename(dependency)
        run(['patch','--batch','--forward','-d',str(dependency),'-p1','-i',str(v/'patch/palabos.patch')])
        with sources.open('x') as f:json.dump({'hemocell_commit':COMMIT,'palabos_commit':PAL,'url':URL,'archive_sha256':ARCHIVE_SHA,'patch_sha256':sha(v/'patch/palabos.patch'),'setup_sha256':sha(v/'setup.sh'),'patch_applied_once':True,'upstream_setup_executed':False,'hemocell_license':'AGPL-3.0-or-later','local_core_changes':False},f,indent=2)
    binary=ROOT/'build/benchmark/pure_fluid_benchmark';library=ROOT/'build/upstream/libhemocell.a'
    manifest=ROOT/'metadata/installed.json'
    def identity():
        return {'hemocell_commit':COMMIT,'palabos_commit':PAL,'binary':str(binary),'binary_sha256':sha(binary),'library_sha256':sha(library),
                'case_sha256':{p.name:sha(p) for p in (ROOT/'cases/pure_fluid_benchmark').iterdir() if p.suffix in ['.cpp','.txt']},
                'precision_bits':64,'MPI':'OpenMPI 4.1.6','HDF5':'parallel C and HL, system 1.10.10','core_modified':False}
    if manifest.exists() and binary.exists() and library.exists() and json.loads(manifest.read_text())==identity():
        record['status']='REUSED_VERIFIED_NO_COMPILATION'
    else:
        if manifest.exists():raise RuntimeError('Installed identity changed. Inspect mismatch before explicit rebuild; no automatic overwrite.')
        run(['/usr/bin/cmake','-S',str(v),'-B',str(ROOT/'build/upstream'),'-DCMAKE_BUILD_TYPE=Release','-DBUILD_TESTING=OFF','-DCMAKE_C_COMPILER=/usr/bin/gcc','-DCMAKE_CXX_COMPILER=/usr/bin/g++','-DMPI_C_COMPILER=/usr/bin/mpicc.openmpi','-DMPI_CXX_COMPILER=/usr/bin/mpicxx.openmpi','-DHDF5_PREFER_PARALLEL=ON'])
        run(['/usr/bin/cmake','--build',str(ROOT/'build/upstream'),'--target','hemocell','--parallel','2'])
        run(['/usr/bin/cmake','-S',str(ROOT/'cases/pure_fluid_benchmark'),'-B',str(ROOT/'build/benchmark'),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_C_COMPILER=/usr/bin/gcc','-DCMAKE_CXX_COMPILER=/usr/bin/g++'])
        run(['/usr/bin/cmake','--build',str(ROOT/'build/benchmark'),'--target','pure_fluid_benchmark','--parallel','2'])
        with manifest.open('x') as f:json.dump(identity(),f,indent=2)
        record['status']='INSTALLED'
    with (ROOT/'metadata'/('setup_'+stamp+'.json')).open('x') as f:json.dump(record,f,indent=2)
    print(record['status']);print(binary)


if __name__=='__main__':
    lock=Path('/home/lzy/projects/mirheo_starter/runs/fluid_calibration/.gpu_exclusive.lock')
    lock.parent.mkdir(parents=True,exist_ok=True)
    with lock.open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        main()
