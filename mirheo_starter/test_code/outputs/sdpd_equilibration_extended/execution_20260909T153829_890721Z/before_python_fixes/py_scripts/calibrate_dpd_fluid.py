"""Execute only explicitly requested small tests under the persistent 3600 s ledger."""
import argparse
from .fluid_physics.common import config_read
from .fluid_physics.calibration import prepare,run_campaign


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default='py_scripts/fluid_physics.yaml');p.add_argument('--execute',action='store_true');p.add_argument('--retry-failed',action='store_true',help='explicit new attempts; preserves previous results and all charges');p.add_argument('--only',nargs='+',help='explicit configured task IDs; does not create a new campaign')
    a=p.parse_args();c=config_read(a.config)
    if not a.execute:p.error('use --execute to authorize GPU execution; --help does not initialize CUDA')
    if a.only and set(a.only)-{t['id'] for t in c['calibration']['tasks']}:p.error('unknown task id')
    prepared=prepare(c);run_campaign(c,prepared,retry_failed=a.retry_failed,only=a.only)


if __name__=='__main__':main()
