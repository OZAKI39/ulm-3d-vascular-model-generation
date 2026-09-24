#!/usr/bin/env python3
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle8_full3d_data import OUTPUT,CASES
from particle_3d.particle8_full3d_visuals import render
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--cases',nargs='+',choices=CASES)
p.add_argument('--tail-mode',choices=['recent','full'],default='recent');p.add_argument('--preview',action='store_true');a=p.parse_args()
render(a.output,a.cases,a.tail_mode,a.preview)
