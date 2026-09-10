"""CPU preflight/review or explicitly authorized short native repair controls."""
import argparse
import json
from py_scripts.single_rbc_repair.workflow import load_config, paths, verify_protection, short_diagnostics, require_authorization, DEFAULT_CONFIG


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default=str(DEFAULT_CONFIG))
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--preflight-only',action='store_true',help='Check CPU evidence, hashes and gates; no solver or compiler')
    mode.add_argument('--execute',action='store_true',help='Run the frozen short diagnostic queue only after explicit repair authorization')
    mode.add_argument('--review',action='store_true',help='Generate the offline page from stored records only')
    args = p.parse_args(); c = load_config(args.config); b,_,out = paths(c)
    if args.execute:
        result = short_diagnostics(c)
    elif args.review:
        from py_scripts.single_rbc_repair.reporting import export
        result = export(c)
    else:
        required = ['failure_timeline.json','root_cause_evidence.json','material_matching.json','candidate_comparability.json','frozen_benchmark_plan.json','authorization_request.json']
        missing = [n for n in required if not (b/n).exists()]
        try:require_authorization(c);authorized=True;authorization_reason=None
        except (RuntimeError,FileNotFoundError) as e:authorized=False;authorization_reason=str(e)
        result = dict(status='CPU_PREFLIGHT_READY' if not missing else 'INCOMPLETE_CPU_EVIDENCE',missing=missing,
                      protection=verify_protection(c),gpu_authorized=authorized,authorization_reason=authorization_reason,
                      solver_started=False,compiler_started=False,qualified_speedup=None,data=str(b),review=str(out/'comparison_review.html'))
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
