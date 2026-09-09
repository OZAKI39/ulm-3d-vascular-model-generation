"""Mirheo vs HemoCell/Palabos pure fluid deployment audit; default is CPU preflight."""
import argparse
import json
from py_scripts.solver_benchmark.physics import load_config


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',default='py_scripts/solver_benchmark.yaml')
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--preflight-only',action='store_true')
    mode.add_argument('--execute',action='store_true')
    mode.add_argument('--analyze-only',action='store_true')
    mode.add_argument('--register-authorization',metavar='CONFIRMATION_JSON')
    args=parser.parse_args()
    from py_scripts.solver_benchmark.workflow import freeze,execute_cpu,execute_gpu,register_authorization
    c=load_config(args.config)
    if args.register_authorization:
        print(register_authorization(c,args.register_authorization));return
    frozen,budget=freeze(c)
    if args.execute:
        execute_cpu(c,frozen);print(json.dumps(execute_gpu(c,frozen),ensure_ascii=False))
    from py_scripts.solver_benchmark.reporting import export
    print(export(c,frozen))


if __name__=='__main__':main()
