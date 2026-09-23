"""Verify the curated GitHub snapshot without running CFD or integrating trajectories."""
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024**2), b''):
            result.update(chunk)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scientific-inputs', action='store_true',
                        help='Also use the original FEM and SonoVue contract loaders; requires scientific dependencies.')
    args = parser.parse_args()
    failed, checked = [], 0
    for line in (HERE / 'workflow_files.sha256').read_text().splitlines():
        expected, relative = line.split('  ', 1)
        path = (ROOT / relative).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file() or digest(path) != expected:
            failed.append(relative)
        checked += 1
    result = {'status': 'PASS' if not failed else 'FAIL',
              'checked_files': checked, 'failures': failed,
              'scope': 'Retained workflow files; sync metadata is excluded to avoid recursive hashes.'}
    if args.scientific_inputs and not failed:
        import sys
        sys.path.insert(0, str(ROOT / 'particle_3d/src'))
        from particle_3d.audit import read_frozen
        from particle_3d.sonovue_adapter import read_sonovue
        fem = read_frozen(ROOT / 'formal_3D_flow_solver/FEM_SimVascular')[0]
        sonovue = read_sonovue(ROOT / 'sonovue_size_distribution_v0')
        result['frozen_FEM_contract'] = fem['status']
        result['sonovue_contract'] = sonovue[0]['status']
        result['sonovue_hashed_files'] = len(sonovue[2])
    print(json.dumps(result, indent=2))
    raise SystemExit(bool(failed))


if __name__ == '__main__':
    main()
