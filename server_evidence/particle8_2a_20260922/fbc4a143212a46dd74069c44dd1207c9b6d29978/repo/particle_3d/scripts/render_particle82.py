#!/usr/bin/env python3
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_visuals import render
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);render(p.parse_args().root)
