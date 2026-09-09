"""Compare native DPD/SDPD using the verified remainder of the original GPU budget."""
import argparse


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',default='py_scripts/fluid_model_comparison.yaml')
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--preflight-only',action='store_true');modes.add_argument('--execute',action='store_true')
    parser.add_argument('--only',nargs='+');parser.add_argument('--retry-failed',action='store_true')
    args=parser.parse_args()
    from py_scripts.fluid_comparison.experiments import load_config,prepare,execute
    c=load_config(args.config);p=prepare(c);print('PREPARED '+str(p),flush=True)
    if args.execute:execute(c,p,only=args.only,retry_failed=args.retry_failed)


if __name__=='__main__':main()
