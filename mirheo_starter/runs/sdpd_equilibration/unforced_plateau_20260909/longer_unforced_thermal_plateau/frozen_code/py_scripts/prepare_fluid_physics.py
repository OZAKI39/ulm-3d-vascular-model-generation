"""CPU-only freeze of accepted physics, units and resource estimates."""
import argparse
from .fluid_physics.common import config_read
from .fluid_physics.calibration import prepare


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default='py_scripts/fluid_physics.yaml');p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();directory=prepare(config_read(a.config));print('PREPARED',directory);print('GPU_NOT_STARTED')


if __name__=='__main__':main()
