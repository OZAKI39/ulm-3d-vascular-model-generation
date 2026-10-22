#!/usr/bin/env python3
"""BG001 frozen-model support planning; local slicing only, never printer control."""
import argparse

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',help='support_aware_roi_60x50 root; discovered automatically by default')
    parser.add_argument('--output',help='Separate support output directory')
    parser.add_argument('--config',help='Planner YAML configuration')
    parser.add_argument('--mode',choices=['rfus','v3.1','v3','v2','custom','benchmark'],default='rfus',help='RFUS full-underbody curtains in fixed O000. V2/V3/V3.1 are read-only research history; custom aliases RFUS.')
    parser.add_argument('--stage',choices=['all','orientation','plan','resume','refine','geometry','research','removal','slice','coupon','report'],default='all')
    args=parser.parse_args()
    if args.mode=='benchmark':
        from vascular_processing.support_inputs import discover
        from vascular_processing.vascular_support_pipeline import run
        output=args.output or str(discover(source=args.source)['root'].parent/'vascular_support_planner/robust_full_underbody/native_benchmark')
        from pathlib import Path
        if any(p in ('v1','v2','v3','v3_1','inverted_study','legacy_explicit_support_v1') for p in Path(output).parts):
            parser.error('Historical support outputs are read-only; use a new RFUS benchmark directory.')
        run(args.source,output,args.stage,args.config,'benchmark')
    elif args.mode in ('v2','v3','v3.1'):
        from vascular_processing.support_inputs import discover
        relative='v3_1/orientation_then_v3_1_report.md' if args.mode=='v3.1' else f'{args.mode}/vascular_support_{args.mode}_report.md'
        path=discover(source=args.source)['root'].parent/'vascular_support_planner'/relative
        print(args.mode.upper()+'_READ_ONLY_HISTORY: '+str(path))
        print(path.read_text())
    else:
        from vascular_processing.support_rfus_pipeline import run
        result=run(args.source,args.output,args.stage,args.config)
        print('RFUS status: '+result['status'])
        if result.get('output'):print('RFUS output: '+result['output'])

if __name__=='__main__':main()
