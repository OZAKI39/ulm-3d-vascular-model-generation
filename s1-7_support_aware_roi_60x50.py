#!/usr/bin/env python3
"""Independent frozen-source support-aware manufacturing ROI stage; local slicing only."""
import argparse
from pathlib import Path
from vascular_processing.support_aware_roi import DEFAULT_OUTPUT,prepare,model


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-root',type=Path,default=DEFAULT_OUTPUT)
    p.add_argument('--stage',choices=['all','prepare','model','baseline','baseline-upright','evaluate','report','box-stl'],default='all')
    p.add_argument('--box-thickness-mm',type=float,default=2.,help='Side-wall and bottom thickness for the box-stl stage')
    args=p.parse_args()
    if args.stage=='box-stl':
        from vascular_processing.support_aware_box_stl import export_box_stl
        export_box_stl(args.output_root,args.box_thickness_mm)
        return
    if args.stage in ('all','prepare'):prepare(args.output_root)
    if args.stage in ('all','model'):model(args.output_root)
    if args.stage == 'baseline':
        from vascular_processing.support_aware_evaluation import baseline
        baseline(args.output_root)
    if args.stage == 'baseline-upright':
        from vascular_processing.support_aware_evaluation import upright_baseline
        upright_baseline(args.output_root)
    if args.stage in ('all','evaluate'):
        from vascular_processing.support_aware_evaluation import evaluate
        evaluate(args.output_root)
        from vascular_processing.support_aware_evaluation import upright_baseline
        upright_baseline(args.output_root)
    if args.stage in ('all','report'):
        from vascular_processing.support_aware_review import report
        report(args.output_root)


if __name__=='__main__':main()
