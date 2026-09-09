"""Diagnose recorded SDPD evidence. GPU execution always requires an explicit flag."""
import argparse


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default='py_scripts/sdpd_diagnostics.yaml')
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--analyze-only',action='store_true')
    mode.add_argument('--execute-approved-probes',action='store_true')
    args=p.parse_args()
    from py_scripts.sdpd_diagnostics.workflow import load_config,export
    c=load_config(args.config);package=export(c)
    print('DIAGNOSIS '+str(package),flush=True)
    if args.execute_approved_probes:
        from py_scripts.sdpd_diagnostics.probes import execute
        execute(c,package)
        package=export(c)
    print('REPORT '+str(package/'report_zh.md'),flush=True)


if __name__=='__main__':main()
