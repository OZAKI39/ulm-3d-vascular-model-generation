#!/usr/bin/env python3
"""Freeze existing NN Pilot artifacts and calibrate four support gates offline."""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import argparse
import json

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,default=ROOT/'outputs/topbrain_brava_transfer')
    p.add_argument('--freeze-nn',action='store_true')
    p.add_argument('--calibrate-unknown',action='store_true')
    p.add_argument('--evaluate-existing',action='store_true')
    p.add_argument('--mevo-binary-comparison',action='store_true',help='Print preserved historical comparison; no recomputation')
    a=p.parse_args()
    if a.mevo_binary_comparison:
        report=json.loads((a.output_dir/'mevo_binary_comparison.json').read_text())
        print(json.dumps({k:report[k] for k in ['status','baseline','FGW','delta','paired_wins']},indent=2))
        return 0
    from vascular_processing.nn_support_gate import freeze_pilot,calibrate
    if a.freeze_nn:freeze_pilot(a.output_dir)
    if a.calibrate_unknown:calibrate(a.output_dir)
    if not (a.freeze_nn or a.calibrate_unknown):p.error('Choose --freeze-nn or --calibrate-unknown')
    return 0

if __name__=='__main__':raise SystemExit(main())
