"""Read in order: geometry -> network -> boundary -> FEM -> flow_field -> particles.

The default CLI is read-only. New scientific computations require a named compute
subcommand, explicit input/output, and cannot overwrite frozen campaigns.
"""
import argparse,json
from mouse_vessel_flow.paths import CASE_NAMES,case_config

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    sub.add_parser('list',help='List retained case identities; no computation')
    for name in ('inspect','analyze'):
        p=sub.add_parser(name);p.add_argument('--case',choices=CASE_NAMES,default='roi_only')
    p=sub.add_parser('validate');p.add_argument('--case',choices=(*CASE_NAMES,'all'),default='all')
    p.add_argument('--full-integrity',action='store_true',help='Hash all migrated scientific files')
    p.add_argument('--output',required=True,help='New report directory, never an old result directory')
    p=sub.add_parser('visualize');p.add_argument('--output',required=True)
    p=sub.add_parser('prepare-cfd');p.add_argument('--case',choices=CASE_NAMES,default='roi_only');p.add_argument('--output',required=True)
    p=sub.add_parser('compute-particles');p.add_argument('--case',choices=CASE_NAMES[:3],required=True)
    p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--output',required=True)
    p=sub.add_parser('compute-cfd');p.add_argument('--directory',required=True);p.add_argument('--build-json',required=True)
    args=parser.parse_args(argv)
    if args.command=='list':value=[case_config(n) for n in CASE_NAMES]
    elif args.command=='inspect':
        from mouse_vessel_flow.fem import inspect_case
        value=inspect_case(args.case)
    elif args.command=='analyze':
        from mouse_vessel_flow.particles.analysis import summarize
        value=summarize(args.case)
    elif args.command=='validate':
        from mouse_vessel_flow.validation import validate
        value=validate(args.case,args.output,full_integrity=args.full_integrity)
    elif args.command=='visualize':
        from mouse_vessel_flow.visualization import report_figures
        value=report_figures(args.output)
    elif args.command=='prepare-cfd':
        from mouse_vessel_flow.fem import prepare_case
        value={'prepared':str(prepare_case(args.case,args.output)),'solver_executed':False}
    elif args.command=='compute-particles':
        from mouse_vessel_flow.particles.trajectories import run_saved_cohort
        value=run_saved_cohort(args.case,args.output,args.ids)
    elif args.command=='compute-cfd':
        from mouse_vessel_flow.fem import run_case
        value=run_case(args.directory,args.build_json)
    print(json.dumps(value,ensure_ascii=False,indent=2))
    return 1 if isinstance(value,dict) and value.get('status')=='FAIL' else 0

if __name__=='__main__':raise SystemExit(main())
