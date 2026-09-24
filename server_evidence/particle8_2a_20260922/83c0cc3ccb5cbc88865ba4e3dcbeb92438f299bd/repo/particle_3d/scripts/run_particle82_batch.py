#!/usr/bin/env python3
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_batch import run_batch

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input-ledger',required=True);p.add_argument('--output-dir',required=True)
    p.add_argument('--workers',type=int,required=True);p.add_argument('--id-start',type=int,required=True)
    p.add_argument('--id-end',type=int,required=True);p.add_argument('--resume',action='store_true')
    p.add_argument('--shard-size',type=int,default=1);p.add_argument('--config',required=True)
    p.add_argument('--host-provenance',required=True)
    a=p.parse_args();read=lambda f:json.loads(Path(f).read_text())
    catalog,perf=run_batch(read(a.input_ledger),a.output_dir,a.workers,a.id_start,a.id_end,a.resume,a.shard_size,read(a.config),read(a.host_provenance))
    print(json.dumps(dict(scheduled=catalog['scheduled'],completed=catalog['completed'],outlets=catalog['outlet_counts'],wall_seconds=perf['wall_seconds'])))
