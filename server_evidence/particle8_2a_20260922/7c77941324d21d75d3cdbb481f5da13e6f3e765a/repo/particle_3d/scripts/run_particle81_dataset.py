#!/usr/bin/env python3
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle81_simulation import run
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--count',type=int,default=2200);p.add_argument('--workers',type=int,default=8)
    p.add_argument('--first-id',type=int,default=1,help='Optional disjoint integration partition; ledger retains full prefix')
    args=p.parse_args();run(args.count,args.workers,args.first_id)
