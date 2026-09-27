"""CLI: freeze an audited pure-0D result without overwriting any artifact."""
from pathlib import Path
import argparse
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from network_1d0d.parameter_freeze import freeze_parameters

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--fit-summary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();r=freeze_parameters(args.config,args.fit_summary,args.output)
    print(r['status'],r['artifact_sha256'],args.output)

if __name__=='__main__':main()
