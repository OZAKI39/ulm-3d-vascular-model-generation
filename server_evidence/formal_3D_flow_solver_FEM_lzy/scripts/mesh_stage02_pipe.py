#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fem3d.benchmark_geometry import build_pipe
p=argparse.ArgumentParser()
p.add_argument("--profile",required=True,choices=("pipe_coarse","pipe_medium","pipe_fine"))
build_pipe(ROOT,p.parse_args().profile)
