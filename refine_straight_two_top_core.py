#!/usr/bin/env python3
"""Export the fine, upright straight ABS core; never invoke a slicer."""
import argparse
from pathlib import Path
from vascular_processing.straight_core_print_quality import ROOT, run

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path)
    parser.add_argument('--orientation-study', action='store_true',
                        help='Compare orientations and export low-support plus horizontal-clearance projects, without slicing')
    args = parser.parse_args()
    if args.orientation_study:
        from vascular_processing.straight_core_orientation import run_study
        run_study(args.config or ROOT/'config/straight_two_top_core_orientation.yaml')
    else:
        run(args.config or ROOT/'config/straight_two_top_core_print_quality.yaml')
