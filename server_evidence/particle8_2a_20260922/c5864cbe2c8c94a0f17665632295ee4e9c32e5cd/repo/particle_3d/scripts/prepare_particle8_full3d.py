#!/usr/bin/env python3
from pathlib import Path
import sys,argparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle8_full3d_data import prepare,OUTPUT
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=OUTPUT);a=p.parse_args()
print(prepare(a.output))
