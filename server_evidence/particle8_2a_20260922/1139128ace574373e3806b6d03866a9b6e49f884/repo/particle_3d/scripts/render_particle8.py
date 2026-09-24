#!/usr/bin/env python3
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle8_visuals import export_animations,export_figures
from particle_3d.particle8_replay import DEFAULT_OUTPUT
p=argparse.ArgumentParser(description='Render only from P8 saved replay caches; no simulation.')
p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
p.add_argument('--animations',nargs='*',type=int,choices=range(1,6))
p.add_argument('--figures-only',action='store_true');p.add_argument('--animations-only',action='store_true')
a=p.parse_args()
if not a.figures_only:export_animations(a.output,a.animations)
if not a.animations_only:export_figures(a.output)
