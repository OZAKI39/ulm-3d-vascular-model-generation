"""Pure 0D bounded design. No CFD execution, calibration, or feedback."""
from pathlib import Path
import argparse
import json
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from network_1d0d.balance_design import run_design,freeze_design


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--freeze',action='store_true',help='Explicitly freeze only after all design gates and sensitivity checks')
    a=p.parse_args();summary,cache=run_design(a.config,a.output)
    if a.freeze:freeze_design(summary,cache,a.output/'frozen_balance_design.yaml')
    print(json.dumps({k:summary[k] for k in ('status','optimized','improvement','crosscheck','deterministic_reproducibility','bound_sensitivity')},indent=2))


if __name__=='__main__':main()
