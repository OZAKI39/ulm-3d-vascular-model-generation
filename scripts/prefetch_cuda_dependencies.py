"""Download pinned Torch CUDA wheels from PyPI with four bounded workers.

This optional installer helper verifies PyPI's SHA256 for each wheel. It does
not use another project's environment and performs no scientific computation.
"""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time
import urllib.request
from packaging.requirements import Requirement
from packaging.tags import sys_tags
from packaging.utils import parse_wheel_filename

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'validation/cache/cuda_wheels'
TAG_ORDER={tag:i for i,tag in enumerate(sys_tags())}
TAGS=set(TAG_ORDER)

def fetch(requirement):
    version=next(iter(requirement.specifier)).version
    with urllib.request.urlopen('https://pypi.org/pypi/'+requirement.name+'/'+version+'/json',timeout=60) as response:
        metadata=json.load(response)
    files=[f for f in metadata['urls'] if f['filename'].endswith('.whl')
           and parse_wheel_filename(f['filename'])[3]&TAGS]
    if not files:raise ValueError('No compatible wheel for '+requirement.name)
    # Some pinned releases publish both legacy and dual manylinux tags.
    # Follow this interpreter's compatibility ordering, then prefer dual tags.
    files.sort(key=lambda f:(min(TAG_ORDER[t] for t in parse_wheel_filename(f['filename'])[3]&TAGS),
                             -len(parse_wheel_filename(f['filename'])[3]),f['filename']))
    record=files[0];path=CACHE/record['filename'];part=path.with_suffix('.part')
    expected=record['digests']['sha256']
    for attempt in range(5):
        try:
            if path.exists() and hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==expected:break
            print('DOWNLOAD',path.name,record['size'],flush=True)
            offset=part.stat().st_size if part.exists() else 0
            if offset<record['size']:
                request=urllib.request.Request(record['url'],headers={'Range':f'bytes={offset}-'} if offset else {})
                with urllib.request.urlopen(request,timeout=120) as response:
                    mode='ab' if offset and response.status==206 else 'wb'
                    with part.open(mode) as stream:
                        while block:=response.read(1024**2):stream.write(block)
            with part.open('rb') as stream:actual=hashlib.file_digest(stream,'sha256').hexdigest()
            if actual!=expected:
                part.unlink()
                raise ValueError('Wheel hash mismatch: '+path.name)
            part.replace(path);break
        except Exception as exc:
            print('RETRY',path.name,repr(exc),flush=True)
            if attempt==4:raise
            time.sleep(2)
    print('VERIFIED',path.name,flush=True)
    return dict(file=path.name,source=record['url'],size=record['size'],sha256=expected,status='PASS')

def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    requirements=[Requirement(r) for r in importlib.metadata.requires('torch')]
    pinned=[r for r in requirements if (r.name.startswith('nvidia-') or r.name=='triton')
            and (not r.marker or r.marker.evaluate())]
    with ThreadPoolExecutor(max_workers=4) as workers:records=list(workers.map(fetch,pinned))
    (ROOT/'validation/CUDA_WHEEL_DOWNLOADS.json').write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':main()
