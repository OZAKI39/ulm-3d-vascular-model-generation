"""Reproduce the current provenance-gated audit, without modifying source cases."""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from network_1d0d.case_audit import run_audit

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'reports/a_network_1d0d_boundary_v1')
    parser.add_argument('--flow-root', type=Path, default=ROOT.parent/'ulm_flow_mean_2p0_mmps')
    parser.add_argument('--fem-root', type=Path, default=ROOT.parent/'formal_3D_flow_solver/FEM')
    args = parser.parse_args()
    run_audit(ROOT, args.output, args.flow_root, args.fem_root)
