#!/usr/bin/env python3
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_results import export_natural

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--host-provenance',required=True);a=p.parse_args()
    c=export_natural(a.root,json.loads(Path(a.host_provenance).read_text()));print(c['outlet_counts'])
