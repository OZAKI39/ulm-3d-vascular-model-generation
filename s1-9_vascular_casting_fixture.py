#!/usr/bin/env python3
"""V2: internal wall-grown PLA supports, closed walls and vertical trim guides."""
import argparse
from vascular_processing.casting_fixture_v2_pipeline import run
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config');p.add_argument('--output');p.add_argument('--geometry-only',action='store_true',help='Skip actual local slicing; status remains NEEDS_ADJUSTMENT');p.add_argument('--no-step',action='store_true',help='Diagnostic run: omit large faceted assembly STEP');p.add_argument('--no-render',action='store_true');a=p.parse_args()
    run(a.config,a.output,do_slice=not a.geometry_only,do_step=not a.no_step,do_render=not a.no_render)
