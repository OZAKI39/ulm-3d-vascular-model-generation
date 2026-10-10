"""Read-only validation of this curated Git export; never rerun science/render jobs."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ast
import hashlib
import json
import os
import re

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
REPO = PROJECT.parents[1]


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check():
    manifest = read(HERE / 'PACKAGE_MANIFEST.json')
    copies = read(HERE / 'COPY_PROVENANCE.json')
    policy = read(HERE / 'EXPORT_POLICY.json')
    paths = []
    links = []
    for row in manifest['files']:
        p = REPO / row['path']
        assert not Path(row['path']).is_absolute() and '..' not in Path(row['path']).parts
        if row['mode'] == '120000':
            assert p.is_symlink() and str(p.readlink()) == row['target'], p
            assert p.exists() and p.resolve().is_relative_to(PROJECT), p
            links.append(row['path'])
        else:
            assert p.is_file() and not p.is_symlink(), p
            assert p.stat().st_size == row['bytes'] and sha(p) == row['sha256'], p
            assert p.stat().st_size <= policy['new_file_limit_bytes'], p
        paths.append(p)

    for row in copies:
        p = REPO / row['path']
        if 'symlink' in row:
            assert str(p.readlink()) == row['symlink'], p
            assert sha(p) == row['content_sha256'], p
        else:
            assert sha(p) == row['sha256'], p

    original = read(PROJECT / 'deployment/newserver_20261009/inventory/SOURCE_FILES.json')
    expected_change = 'computation_examples/mevo_acoustic_microbubbles/.gitignore'
    verified = 0
    for row in original:
        if row['path'] == expected_change:
            continue
        p = REPO / row['path']
        assert p.is_file() and sha(p) == row['sha256'], p
        verified += 1

    secrets = {
        'private_key': re.compile(r'-----BEGIN (?:OPENSSH|RSA|EC|DSA|ENCRYPTED)? ?PRIVATE KEY-----'),
        'github_token': re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b'),
        'aws_key_id': re.compile(r'\bAKIA[A-Z0-9]{16}\b'),
        'live_bearer_value': re.compile(r'(?i)Bearer\s+[A-Za-z0-9._~+/-]{32,}'),
    }
    py = js = scanned = 0
    forbidden = {'.venv','.venv-render','.python','.cache','__pycache__','node_modules','wheels'}
    for p in paths:
        assert not forbidden.intersection(p.relative_to(REPO).parts), p
        assert p.suffix not in {'.whl','.ttf','.pem','.key'}, p
        if p.is_symlink():
            continue
        if p.suffix == '.py':
            ast.parse(p.read_text(), filename=str(p))
            py += 1
        if p.suffix == '.json':
            read(p)
            js += 1
        if p.suffix in {'.py','.json','.md','.txt','.log','.sh','.conf','.toml','.yml','.yaml'} or p.name == '.gitignore':
            data=p.read_text()
            for kind, pattern in secrets.items():
                assert not pattern.search(data), f'{kind}: {p}'
            scanned += 1

    import numpy as np
    import h5py
    datasets=0
    trajectories=[]
    for p in paths:
        if p.suffix != '.h5' or p.is_symlink():
            continue
        with h5py.File(p, 'r') as h:
            def visit(name, obj):
                nonlocal datasets
                if isinstance(obj,h5py.Dataset):
                    obj[()]
                    datasets += 1
            h.visititems(visit)
            x=h['positions_m'][:]
            t=h['time_s'][:]
            np.testing.assert_array_equal(h['particle_id'][:],np.arange(32))
            np.testing.assert_array_equal(h['radii_m'][:],np.full(32,1.78e-6))
            assert x.shape[-2:]==(32,3) and np.isfinite(x).all()
            assert np.isfinite(t).all() and (np.diff(t)>0).all()
            assert h['active'][:].all()
            trajectories.append(dict(path=str(p.relative_to(REPO)),frames=len(t),particles=32,
                                     start_s=float(t[0]),end_s=float(t[-1]),sha256=sha(p)))

    validation=read(PROJECT/'deployment/newserver_20261009/validation/VALIDATION.json')
    assert validation['status']=='PASS'
    assert read(PROJECT/'deployment/newserver_20261009/LOCAL_BACKUP_DELETION.json')['status']=='DELETED_VERIFIED'
    assert read(PROJECT/'deployment/newserver_20261009/CURRENT_SERVER_CHECK_20261010.json')['status']=='PASS'
    return dict(status='PASS', checked_utc=datetime.now(timezone.utc).isoformat(),
                package_files=len(paths), copied_artifacts_checked=len(copies),
                inherited_baseline_files_unchanged=verified,
                inherited_baseline_exception={'path':expected_change,'reason':'Ignore new environment/cache/wheel files'},
                frozen_numerical_sources_configs_and_F1_F5_results_unchanged=True,
                relative_input_links=links, parsed_python_files=py, parsed_JSON_files=js,
                credential_pattern_scan_files=scanned, credential_pattern_matches=0,
                maximum_regular_file_bytes=max(p.stat().st_size for p in paths if not p.is_symlink()),
                full_HDF5_dataset_reads=datasets, original_visualization_trajectories=trajectories,
                migration_validation_date='2026-10-09', migration_validation_status='PASS',
                historical_log_bytes_preserved=True, new_scientific_simulations=0,
                note='Packaging, integrity and structured-data checks only; no new physics validation.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--write-report',action='store_true')
    args=ap.parse_args()
    result=check()
    if args.write_report:
        (HERE/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
