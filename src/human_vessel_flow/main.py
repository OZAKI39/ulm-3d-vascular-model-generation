"""Read this entry first: geometry -> mesh -> boundary -> native FEM -> particles.

Default actions inspect frozen evidence. Computation has separate explicit commands.
"""
from pathlib import Path
import argparse, json
from human_vessel_flow.paths import ROOT

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    commands.add_parser('status',help='Read current scientific/result identities')
    commands.add_parser('native-status',help='Check configured native binary identity; do not solve')
    for name in ['verify','plot','boundary-design','prepare-case','prepare-geometry','prepare-mesh','prepare-population']:
        sub=commands.add_parser(name);sub.add_argument('--output',type=Path,required=True)
        if name=='prepare-population':sub.add_argument('--backend',choices=['cpu','cuda'],default='cuda')
    sub=commands.add_parser('solve');sub.add_argument('--case',type=Path,required=True)
    sub=commands.add_parser('postprocess');sub.add_argument('--case',type=Path,required=True)
    sub=commands.add_parser('run-population');sub.add_argument('--output',type=Path,required=True);sub.add_argument('--backend',choices=['cpu','cuda'],default='cuda')
    sub=commands.add_parser('audit-generated-mesh');sub.add_argument('--raw',type=Path,required=True);sub.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.command=='status':
        from human_vessel_flow.geometry import describe
        from human_vessel_flow.particles import historical_summary
        from human_vessel_flow.validation import scientific_gates
        result=dict(project=str(ROOT),geometry=describe(),historical_population=historical_summary(),scientific_gates=scientific_gates())
    elif args.command=='native-status':
        from human_vessel_flow.fem import native_status
        result=native_status()
    elif args.command=='verify':
        from human_vessel_flow.validation import run_checks
        result=run_checks(args.output)
    elif args.command=='plot':
        from human_vessel_flow.visualization import build_review
        result=build_review(args.output)
    elif args.command=='boundary-design':
        from human_vessel_flow.boundary import write_design
        result=write_design(args.output)
    elif args.command=='prepare-case':
        from human_vessel_flow.fem import prepare_case
        result=dict(case=str(prepare_case(args.output)),executed=False)
    elif args.command=='prepare-mesh':
        from human_vessel_flow.meshing import prepare
        result=dict(command=prepare(args.output),executed=False)
    elif args.command=='prepare-geometry':
        from human_vessel_flow.geometry import prepare_surface
        result=prepare_surface(args.output)
    elif args.command=='audit-generated-mesh':
        from human_vessel_flow.mesh_audit import audit_generated
        result=audit_generated(args.raw,args.output)
    elif args.command=='solve':
        from human_vessel_flow.fem import solve
        result=solve(args.case)
    elif args.command=='postprocess':
        from human_vessel_flow.postprocess import postprocess
        result=postprocess(args.case)
    elif args.command in {'prepare-population','run-population'}:
        from human_vessel_flow.particles import prepare_run,run
        command=prepare_run(args.output,args.backend)
        result=dict(command=command,executed=False)
        if args.command=='run-population':result.update(exit_code=run(command,args.output),executed=True)
    print(json.dumps(result,ensure_ascii=False,indent=2,default=str))
    return 1 if isinstance(result,dict) and any(v=='FAIL' for v in result.values() if isinstance(v,str)) else 0

if __name__=='__main__':raise SystemExit(main())
