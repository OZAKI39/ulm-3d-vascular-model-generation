"""One fixed unforced SDPD experiment; CPU preparation/viewing never starts CUDA."""
import argparse


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default='py_scripts/sdpd_equilibration.yaml')
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--preflight-only',action='store_true')
    mode.add_argument('--execute',action='store_true')
    mode.add_argument('--analyze-only',action='store_true')
    mode.add_argument('--register-authorization',metavar='USER_CONFIRMATION_JSON',
                      help='Only after a real user message explicitly authorizes additional seconds and this exact scope.')
    args=p.parse_args()
    from py_scripts.sdpd_diagnostics.equilibration import load_config,export,execute,register_authorization
    c=load_config(args.config)
    try:
        if args.register_authorization:
            print('AUTHORIZATION_RECORD '+str(register_authorization(c,args.register_authorization)))
            return
        package=execute(c) if args.execute else export(c,analyze=args.analyze_only)
        from py_scripts.fluid_physics.common import read_json
        from test_code.review_sdpd_equilibration import write_page
        page=write_page(c,package)
        print('STATUS '+read_json(package/'equilibration_summary.json')['status'])
        print('PACKAGE '+str(package))
        print('REPORT '+str(package/'report_zh.md'))
        print('HTML '+str(page))
    except PermissionError as exc:
        p.exit(3,str(exc)+'\n')


if __name__=='__main__':main()
