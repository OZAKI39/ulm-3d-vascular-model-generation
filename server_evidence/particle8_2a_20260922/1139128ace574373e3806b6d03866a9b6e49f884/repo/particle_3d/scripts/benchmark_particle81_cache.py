#!/usr/bin/env python3
"""Recompute four original cohort IDs with exact BVH/gap memoization."""
from pathlib import Path
import sys,multiprocessing
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle81_simulation import OUTPUT,environment,integrate_one,dump
from particle_3d.particle8_replay import read


def task(event):return integrate_one(event,output=OUTPUT/'data/cache_v2_benchmark',force=True)


def main():
    root=OUTPUT/'data/cache_v2_benchmark';births=read(OUTPUT/'data/birth_ledger.json')['events'];ids=[1,4,7,13]
    environment()
    with multiprocessing.get_context('fork').Pool(2) as pool:rows=list(pool.imap_unordered(task,[births[i-1] for i in ids]))
    result=[]
    for meta in rows:
        pid=meta['particle_id'];original=read(OUTPUT/f'trajectories/mb_{pid:06d}.json')
        a=np.load(OUTPUT/original['samples_path'])['samples'];b=np.load(root/meta['samples_path'])['samples']
        result.append(dict(particle_id=pid,raw_arrays_bitwise_identical=a.dtype==b.dtype and a.shape==b.shape and a.tobytes()==b.tobytes(),terminal_same=meta['end_reason']==original['end_reason'],
            wall_seconds=meta['wall_seconds'],old_wall_seconds=original['wall_seconds'],provider_calls_same=meta['provider_calls']==original['provider_calls']))
    dump(OUTPUT/'data/cache_v2_parity.json',dict(cases=result,
        all_pass=all(r['raw_arrays_bitwise_identical'] and r['terminal_same'] and r['provider_calls_same'] for r in result),
        comparison='SHAPE_DTYPE_AND_C_ORDER_RAW_BYTES_INCLUDING_SIGNED_ZERO',role='EXACT_BVH_AND_GEOMETRY_QUERY_MEMOIZATION_ONLY'))
    print(result)


if __name__=='__main__':main()
