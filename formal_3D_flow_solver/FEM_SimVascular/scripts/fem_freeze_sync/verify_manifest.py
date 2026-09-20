"""Verify the committed handoff inventory, without source WSL access or CFD."""
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXCLUDED = {'external', '__pycache__', '.pytest_cache', '.venv', 'local_runtime'}
MANIFEST = 'sync_metadata/source_manifest.csv'


def payload_files():
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and
                  not EXCLUDED.intersection(p.relative_to(ROOT).parts) and p.suffix != '.pyc')


def verify():
    with (ROOT/MANIFEST).open(newline='') as f:
        rows = list(csv.DictReader(f))
    expected = {r['relative_path'] for r in rows}
    actual = {p.relative_to(ROOT).as_posix() for p in payload_files()} - {MANIFEST}
    assert len(rows) == len(expected) and expected == actual, (expected-actual, actual-expected)
    for r in rows:
        p = ROOT/r['relative_path']
        assert p.resolve().is_relative_to(ROOT)
        assert p.stat().st_size == int(r['size_bytes'])
        assert hashlib.sha256(p.read_bytes()).hexdigest() == r['sha256'], r['relative_path']
        assert p.stat().st_size < 95*2**20
    for line in (ROOT/'frozen_reference/SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split('  ', 1)
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    return len(rows)


if __name__ == '__main__':
    print(f'PASS: {verify()} handoff files match size/SHA256; canonical scientific hashes verified; no CFD.')
