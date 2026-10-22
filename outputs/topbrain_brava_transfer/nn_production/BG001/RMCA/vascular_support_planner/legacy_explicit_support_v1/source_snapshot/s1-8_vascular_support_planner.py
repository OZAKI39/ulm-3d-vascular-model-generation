#!/usr/bin/env python3
"""BG001 frozen-model support planning; local slicing only, never printer control."""
import argparse
from vascular_processing.vascular_support_pipeline import run

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',help='support_aware_roi_60x50 root; discovered automatically by default')
    parser.add_argument('--output',help='Separate support output directory')
    parser.add_argument('--config',help='Planner YAML configuration')
    parser.add_argument('--mode',choices=['custom','benchmark'],default='custom',help='Explicit custom support is production; slicer trees are benchmark only')
    parser.add_argument('--stage',choices=['all','plan','geometry','research','removal','slice','coupon','report'],default='all')
    args=parser.parse_args();run(args.source,args.output,args.stage,args.config,args.mode)

if __name__=='__main__':main()
