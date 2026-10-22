#!/usr/bin/env python3
"""Frozen vascular core: continuous roof + segmented removable ABS curtains."""
import argparse
from vascular_processing.core_support_pipeline import run

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output')
    p.add_argument('--printer',choices=['p1s-0.4','x1c-0.2'],default='p1s-0.4',help='Current production target is P1S / 0.4 mm / ABS; historical X1C outputs remain separate')
    p.add_argument('--stage',choices=['all','closure','prepare','geometry','slice','review','legacy-review','legacy-closure'],default='all')
    p.add_argument('--max-iterations',type=int,default=3,help='At most 3 practical repair rounds; P1S uses current topology revision 3 once; completed native slices are reused')
    args=p.parse_args()
    if not 1<=args.max_iterations<=3:p.error('--max-iterations must be between 1 and 3')
    if args.printer=='p1s-0.4':
        if args.stage not in ['all','closure','prepare','geometry','slice','review']:
            p.error('Legacy stages require --printer x1c-0.2')
        from vascular_processing.core_support_p1s import run as p1s_support
        s=p1s_support(args.output,args.stage,args.max_iterations)
    elif args.stage in ['all','closure','review']:
        from vascular_processing.core_support_practical_pipeline import run as practical_support
        s=practical_support(args.output,args.max_iterations,review_only=args.stage=='review')
    elif args.stage=='legacy-closure':
        from vascular_processing.core_support_closure import run as close_support
        s=close_support(args.output,args.max_iterations,review_only=True)
    else:s=run(args.output,'review' if args.stage=='legacy-review' else args.stage)
    raise SystemExit(2 if s.get('status')=='NEEDS_ADJUSTMENT' else 0)
