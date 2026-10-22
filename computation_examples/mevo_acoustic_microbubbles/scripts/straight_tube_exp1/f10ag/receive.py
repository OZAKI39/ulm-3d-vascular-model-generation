"""Read-only receipt: raw full-step data and recomputable caches stay on server."""
import subprocess,shlex,time
from .common import *
from .remote import SSH,SERVER

def main():
    start=time.perf_counter();BASE.mkdir(parents=True,exist_ok=True)
    cmd=['rsync','-az','--exclude=hf/*/runs/*/trajectory.h5','--exclude=hf/*/runs/*/LAST_STEP.json','--exclude=fields','--exclude=fields/*/coarse.h5','--exclude=fields/*/fine.h5','--exclude=evidence/*.tar.gz','-e',shlex.join(SSH[:-1]),SSH[-1]+':'+SERVER+'/results/straight_tube_exp1/f10ag/',str(BASE)+'/'];subprocess.run(cmd,check=True)
    subprocess.run(['rsync','-az','-e',shlex.join(SSH[:-1]),SSH[-1]+':'+SERVER+'/logs/',str(BASE/'logs')+'/'],check=True)
    files={str(p.relative_to(BASE)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BASE.rglob('*') if p.is_file() and p.name!='LOCAL_RECEIPT.json' and '.tar.gz' not in p.name}
    write(BASE/'LOCAL_RECEIPT.json',dict(files=files,wall_s=time.perf_counter()-start,retained_remotely=['Full per-step HF trajectory.h5','LAST_STEP.json','coarse/fine caches','execution source archives'],server_root=SERVER,identity=identity()))
if __name__=='__main__':main()
