#!/usr/bin/env python3
"""Read-only, resumable SSH synchronization; never deletes remote/local evidence."""
from pathlib import Path
import argparse,json,subprocess
REPO=Path(__file__).resolve().parents[2]
SSH='ssh -p 4159 -i /home/lzy/.ssh/vast_step3b_ed25519 -o IdentitiesOnly=yes -o BatchMode=yes -o ConnectTimeout=15 -o ServerAliveInterval=20'
FULL='/root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit'
POINT='/root/particle8_2_runs/a4820609c9ed71f09bcea968b5819701632a2a1e/diagnostic_pipeline/point_basin_100000'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['summaries','review','all'],default='summaries');a=p.parse_args()
    out=REPO/'particle_3d/reports/particle8_2';out.mkdir(exist_ok=True)
    # Full paths are discoverable in REMOTE_RAW_INVENTORY.json; large diagnostic
    # arrays can be fetched later using --mode all without recomputation.
    cmd=['rsync','-a','-z','--partial','--info=stats1','-e',SSH,'--exclude=repo/']
    if a.mode=='summaries':
        cmd += ['--exclude=trajectories/','--exclude=receipts/','--exclude=shards/','--exclude=point_paths/','--exclude=animations/','--exclude=frames/','--exclude=keyframes/','--exclude=inspection/','--exclude=PROPOSAL_INDEX.json','--exclude=*.npz','--exclude=*.csv.gz']
    elif a.mode=='review':
        cmd += ['--exclude=point_paths/','--exclude=point_refinement_100000/shards/','--exclude=scaling/workers_*/trajectories/','--exclude=continuation/guard_x*/trajectories/','--exclude=timestep/dt_div_*/trajectories/']
    subprocess.run(cmd+['root@50.115.148.16:'+FULL+'/',str(out)+'/'],check=True)
    point_args=['rsync','-a','-z','--partial','-e',SSH]
    if a.mode!='all':point_args+=['--exclude=shards/','--exclude=*.npz']
    subprocess.run(point_args+['root@50.115.148.16:'+POINT+'/',str(out/'point_basin_100000')+'/'],check=True)
