#!/usr/bin/env python3
"""Straight 1.6 mm ABS core and PLA fixtures; unsliced projects only."""
import argparse
from pathlib import Path
from vascular_processing.straight_channel_casting import ROOT, run

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT/'config/straight_channel_1p6_p1s_abs.yaml')
    parser.add_argument('--stage', choices=['geometry', 'project', 'all'], default='all',
                        help='geometry: models; project: package saved models; all: both, without slicing')
    args = parser.parse_args()
    run(args.config, args.stage)
