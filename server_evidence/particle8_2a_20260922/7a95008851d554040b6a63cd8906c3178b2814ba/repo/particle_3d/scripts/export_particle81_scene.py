#!/usr/bin/env python3
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle81_replay import assemble
if __name__=='__main__':
    c=assemble();print({k:c[k] for k in ['scheduled','admitted','completed','outlet_counts','physical_samples']})
