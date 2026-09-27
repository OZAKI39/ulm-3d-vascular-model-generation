"""CLI: geometry-only forward regression or pure-0D fit. Never launches CFD."""
from pathlib import Path
import argparse
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from network_1d0d.parameter_config import run_configuration

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True,help='New output directory; never overwrites an earlier run')
    args=p.parse_args();r=run_configuration(args.config,args.output)
    print(r['status'],args.output/'fit_summary.json')
    if r['status']=='FIT_NOT_CONVERGED':raise SystemExit(2)

if __name__=='__main__':main()
