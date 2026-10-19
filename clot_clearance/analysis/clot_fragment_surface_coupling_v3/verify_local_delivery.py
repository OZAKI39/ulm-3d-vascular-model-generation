"""Read back the sealed V3 archive with Python's standard library only."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import time

TASK = 'clot_fragment_surface_coupling_v3'
RUN = '20261018T170400Z'


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 ** 2), b''):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / 'verification' / TASK / RUN
    manifest_path = output / 'DELIVERY_MANIFEST.json'
    before = digest(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    assert (manifest['task'], manifest['run']) == (TASK, RUN)
    assert manifest['file_count'] == len(manifest['files'])
    assert manifest['total_bytes'] == sum(v['bytes'] for v in manifest['files'].values())
    failures = []
    checked_bytes = 0
    begin = time.monotonic()
    for index, (relative, expected) in enumerate(manifest['files'].items(), 1):
        name = Path(relative)
        assert not name.is_absolute() and '..' not in name.parts
        path = root / name
        if path.is_symlink() or not path.is_file():
            failures.append(dict(path=relative, reason='missing or wrong file type'))
        elif path.stat().st_size != expected['bytes']:
            failures.append(dict(path=relative, reason='byte count mismatch'))
        elif digest(path) != expected['sha256']:
            failures.append(dict(path=relative, reason='SHA256 mismatch'))
        else:
            checked_bytes += expected['bytes']
        if index % 50 == 0:
            print('READBACK', index, '/', manifest['file_count'], 'files', checked_bytes, 'bytes', flush=True)
    for relative, expected in manifest['symlinks'].items():
        name = Path(relative)
        assert not name.is_absolute() and '..' not in name.parts
        path = root / name
        if not path.is_symlink() or os.readlink(path) != expected:
            failures.append(dict(path=relative, reason='symlink text mismatch'))
        elif not path.exists():
            failures.append(dict(path=relative, reason='broken symlink'))
    assert digest(manifest_path) == before, 'Manifest changed during readback'
    result = dict(status='PASS' if not failures else 'FAIL', task=TASK, run=RUN,
                  utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  project_root=str(root), manifest_sha256=before,
                  checked_files=manifest['file_count'], checked_symlinks=len(manifest['symlinks']),
                  verified_bytes=checked_bytes, wall_s=time.monotonic() - begin,
                  failures=failures, scope='Every sealed regular file byte count and SHA256; exact symlink text and resolvability. Historical V2 dependencies retain their separate archive manifest.')
    temporary = output / 'LOCAL_DELIVERY_READBACK.json.tmp'
    temporary.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(output / 'LOCAL_DELIVERY_READBACK.json')
    print(json.dumps(result, ensure_ascii=False), flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
