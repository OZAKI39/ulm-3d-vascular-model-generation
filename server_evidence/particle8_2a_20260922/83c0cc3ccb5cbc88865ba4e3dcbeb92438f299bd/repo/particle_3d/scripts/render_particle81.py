#!/usr/bin/env python3
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle81_visuals import render_animations,ANIMATIONS
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preview',action='store_true');p.add_argument('--kinds',nargs='+',choices=ANIMATIONS)
    a=p.parse_args();render_animations(kinds=a.kinds,preview=a.preview)
