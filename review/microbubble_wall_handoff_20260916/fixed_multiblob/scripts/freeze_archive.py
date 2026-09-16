#!/usr/bin/env python3
"""Freeze all archive files; keep verification receipts outside the archive."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt-dir', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    receipt_dir = args.receipt_dir.resolve()
    if receipt_dir == root or root in receipt_dir.parents:
        raise SystemExit('Verification receipts must be outside the archive.')
    manifest = root / 'SHA256SUMS'
    if manifest.exists():
        raise SystemExit('Refusing to replace an existing frozen manifest.')
    paths = sorted(root.rglob('*'))
    if any(path.is_symlink() for path in paths):
        raise SystemExit('Symlinks require an explicit archive policy.')
    files = [path for path in paths if path.is_file()]
    lines = []
    for path in files:
        relative = path.relative_to(root).as_posix()
        if '\n' in relative or '\\' in relative:
            raise SystemExit('Filename requires GNU checksum escaping.')
        lines.append(f'{digest(path)}  {relative}\n')
    with manifest.open('x', encoding='utf-8', newline='\n') as stream:
        stream.writelines(lines)
    receipt_dir.mkdir(parents=True, exist_ok=True)
    log = receipt_dir / 'FINAL_LOCAL_SHA_VERIFY.log'
    with log.open('w', encoding='utf-8') as stream:
        result = subprocess.run(
            ['sha256sum', '-c', 'SHA256SUMS'], cwd=root,
            stdout=stream, stderr=subprocess.STDOUT, check=False)
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*')
              if p.is_file() and p != manifest}
    expected = {p.relative_to(root).as_posix() for p in files}
    receipt = {
        'status': 'PASS' if result.returncode == 0 and actual == expected else 'FAIL',
        'archive': str(root),
        'manifest_sha256': digest(manifest),
        'verified_file_count': len(files),
        'exact_file_set_match': actual == expected,
        'sha256sum_return_code': result.returncode,
        'verification_log': str(log),
        'manifest_policy': 'All regular files, including nested manifests; only root SHA256SUMS excluded.'
    }
    (receipt_dir / 'FINAL_LOCAL_SHA_VERIFY.json').write_text(
        json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(receipt, indent=2), flush=True)
    if receipt['status'] != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
