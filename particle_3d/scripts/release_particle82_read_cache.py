#!/usr/bin/env python3
"""Advise eviction of this stage's clean output cache, never drop system caches."""
from pathlib import Path
import argparse,json,os,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_provenance import require_remote,atomic_json

p=argparse.ArgumentParser();p.add_argument('--directory',action='append',required=True)
p.add_argument('--host-provenance',required=True);p.add_argument('--output',required=True);a=p.parse_args()
h=json.loads(Path(a.host_provenance).read_text());require_remote(h,h['hostname'])
def memory():
    return {name:(Path('/sys/fs/cgroup')/name).read_text() for name in ['memory.current','memory.stat','memory.events']}
before=memory();count=total=0;start=time.time()
for directory in a.directory:
    folder=Path(directory).resolve()
    if not folder.is_relative_to('/root/particle8_2_runs'):raise ValueError('Only this stage output directories may be advised')
    if folder.name not in ['shards','point_paths']:raise ValueError('Only saved point-path arrays may be advised')
    for file in folder.glob('*.npz'):
        prior=file.stat();fd=os.open(file,os.O_RDONLY)
        try:os.posix_fadvise(fd,0,0,os.POSIX_FADV_DONTNEED)
        finally:os.close(fd)
        after=file.stat()
        if (prior.st_size,prior.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise ValueError('Input file metadata unexpectedly changed')
        count+=1;total+=prior.st_size
atomic_json(a.output,dict(hostname=h['hostname'],count=count,bytes_advised=total,seconds=time.time()-start,
    before=before,after=memory(),method='READ_ONLY_POSIX_FADV_DONTNEED_ON_OWN_CLEAN_NPZ_OUTPUTS',
    no_system_drop_caches=True,no_file_content_changes=True,no_scientific_data_deletion=True))
print(json.dumps(dict(count=count,bytes_advised=total,before=before['memory.current'],after=memory()['memory.current'])))
