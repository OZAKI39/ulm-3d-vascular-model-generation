#!/usr/bin/env python3
"""Flat CSV views of original diagnostic records, without changing classifications."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_diagnostics import csv_write
p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();root=Path(a.root)
read=lambda x:json.loads((root/x).read_text())
csv_write(root/'data/timestep_audit.csv',read('timestep/TIMESTEP_AUDIT.json')['rows'])
rows=[]
for factor in read('continuation/EXTENDED_GUARD_AUDIT.json')['factors']:
 rows.extend(dict(guard_factor=factor['guard_factor'],**r) for r in factor['rows'])
csv_write(root/'data/extended_guard_outcomes.csv',rows)
csv_write(root/'data/worker_scaling_results.csv',read('scaling/SERVER_SCALING_BENCHMARK.json')['worker_scaling_results'])
