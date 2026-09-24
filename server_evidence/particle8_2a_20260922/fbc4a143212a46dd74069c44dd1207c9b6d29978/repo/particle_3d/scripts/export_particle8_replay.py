#!/usr/bin/env python3
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle8_replay import build_bundle,DEFAULT_OUTPUT
p=argparse.ArgumentParser(description='Export read-only replay from accepted P7 evidence; no solver runs.')
p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
a=p.parse_args(); scenes=build_bundle(a.output)
for k,s in scenes.items(): print(k,len(s['records']),'particles',len(s['event_ledger']),'events',flush=True)
