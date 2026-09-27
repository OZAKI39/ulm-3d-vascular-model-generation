"""Run retained workflow checks; opt into complete trajectory re-integration."""
from pathlib import Path
import argparse, json, os, subprocess, sys, time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
TESTS = ['particle9a5_formal_trajectories', 'particle9a4_population_inlet',
         'rbc_mb_flow_rotation', 'rbc_mb_coflow', 'rbc_mb_normal', 'rbc_mb_interaction']
SLOW = 'worker_independent_fixed_subset or same_host_adapter_preserves_frozen_p9a1_samples'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--full-simulation', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', VTK_SMP_MAX_THREADS='1',
               PYTHONPATH=str(ROOT/'particle_3d/src'), SONOVUE_ROOT=str(ROOT/'sonovue_size_distribution_v0'))
    command = [sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider',
               '--junitxml='+str(out/'tests.xml')]
    command += ['particle_3d/tests/'+name for name in TESTS]
    if not args.full_simulation:
        command += ['-k', 'not ('+SLOW+')']
    start = time.monotonic()
    with (out/'tests.log').open('w') as log:
        result = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    tree = ET.parse(out/'tests.xml').getroot()
    suites = [tree] if tree.tag == 'testsuite' else list(tree)
    counts = {k:sum(int(s.get(k,0)) for s in suites) for k in ('tests','failures','errors','skipped')}
    receipt = dict(**counts, returncode=result.returncode, groups=TESTS,
                   mode='full' if args.full_simulation else 'current_inputs_and_short_integration',
                   excluded_full_trajectory_tests=[] if args.full_simulation else SLOW.split(' or '),
                   runtime_s=time.monotonic()-start, command=command)
    (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
    if result.returncode or not counts['tests'] or any(counts[k] for k in ('failures','errors','skipped')):
        raise SystemExit(1)

if __name__ == '__main__':
    main()
