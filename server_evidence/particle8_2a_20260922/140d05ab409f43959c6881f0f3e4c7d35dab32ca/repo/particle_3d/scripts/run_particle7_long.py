#!/usr/bin/env python3
from pathlib import Path
import sys,argparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle7_long import run_long
p=argparse.ArgumentParser(); p.add_argument('--mb-target',type=int,default=1000); p.add_argument('--output',default='particle_3d/reports/particle7/data')
a=p.parse_args(); print(run_long(a.output,a.mb_target),flush=True)
