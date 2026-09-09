"""Read/open the existing offline report; never imports or launches a solver."""
import argparse
import subprocess
from pathlib import Path
from py_scripts.fluid_physics.common import PROJECT_ROOT,read_json,sha256_file


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default='py_scripts/solver_benchmark.yaml');p.add_argument('--open',action='store_true');a=p.parse_args()
    latest=read_json(PROJECT_ROOT/'data/solver_benchmark/LATEST.json');page=Path(latest['html'])
    if sha256_file(page)!=latest['html_sha256']:raise ValueError('REPORT_HASH_MISMATCH')
    config=Path(a.config);config=config if config.is_absolute() else PROJECT_ROOT/config
    if sha256_file(config)!=read_json(latest['results'])['config_sha256']:raise ValueError('CONFIG_REPORT_MISMATCH')
    print(page)
    if a.open:
        target=subprocess.check_output(['wslpath','-w',str(page)],text=True).strip()
        subprocess.Popen(['/mnt/c/Windows/explorer.exe',target],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)


if __name__=='__main__':main()
