import sys,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from dataset import make_banks,generate
from vcommon import dump,DATA,identity
p=argparse.ArgumentParser();p.add_argument('--banks',action='store_true');p.add_argument('--splits',nargs='+',default=[]);p.add_argument('--workers',type=int,default=8);p.add_argument('--limit',type=int);a=p.parse_args()
if a.banks:print(make_banks())
for s in a.splits:print(generate(s,a.workers,a.limit))
