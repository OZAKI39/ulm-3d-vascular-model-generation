"""Read-only candidate search; restore only a staged, exact SHA256 match."""
import concurrent.futures
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time
import zipfile

A = Path(__file__).resolve().parent
ROOT = A.parents[2]
OLD = A.parent / 'repository_cleanup_v1'
ROWS = json.loads((OLD / 'RECOVERY_UNAVAILABLE.json').read_text())
WANTED = {}
for row in ROWS:
    WANTED.setdefault(row['hash'], []).append(row)
SIZES = {r['size_bytes'] for r in ROWS}
FREEZE = Path(json.loads((A / 'BACKUP_RESULT.json').read_text())['destination'])
ERRORS, MATCHES = [], []

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def staged_restore(data, sha, source):
    assert hashlib.sha256(data).hexdigest() == sha
    staged = A / 'recovery_staging' / sha
    staged.parent.mkdir(exist_ok=True)
    if not staged.exists():
        staged.write_bytes(data)
    assert digest(staged) == sha
    for row in WANTED[sha]:
        target = ROOT / row['path']
        if target.exists():
            assert digest(target) == sha, str(target)
            continue
        assert len(data) == row['size_bytes']
        target.parent.mkdir(parents=True, exist_ok=True)
        # Exclusive creation avoids clobbering anything created by another process.
        with target.open('xb') as out:
            out.write(staged.read_bytes())
        assert digest(target) == sha
        record = dict(path=row['path'], sha256=sha, source=source,
                      staged=str(staged.relative_to(A)), status='RECOVERED_FROM_DUPLICATE')
        MATCHES.append(record)
        with (A / 'EXACT_RECOVERY_LOG.jsonl').open('a') as log:
            log.write(json.dumps(record) + '\n')

def candidates():
    count = 0
    # Full index for real directories under projects/temp_storage and current repo;
    # rest of home and /tmp are SHA candidates filtered by exact expected size.
    for scanroot in ['/home/lzy', '/tmp']:
        for parent, dirs, files in os.walk(scanroot, followlinks=False,
                                         onerror=lambda e: ERRORS.append(str(e))):
            dirs[:] = [d for d in dirs if d not in ('.ssh', '.gnupg')
                       and Path(parent) / d != FREEZE
                       and Path(parent) / d != A]
            for name in files:
                p = Path(parent) / name
                if p.is_symlink():
                    continue
                try:
                    st = p.stat()
                    if not p.is_file():
                        continue
                    count += 1
                    full = str(p).startswith('/home/lzy/projects/temp_storage/')
                    if full or st.st_size in SIZES:
                        yield p, st.st_size, 'FULL_SIBLING_INDEX' if full else 'EXPECTED_SIZE_FILTER'
                except OSError as e:
                    ERRORS.append(str(e))
    (A / 'HASH_SCAN_METADATA.json').write_text(json.dumps(dict(
        files_seen=count, scan_roots=['/home/lzy', '/tmp'],
        full_index_roots=['/home/lzy/projects/temp_storage'],
        excluded=['.ssh', '.gnupg', str(FREEZE), str(A)],
        symlinks='not followed; underlying real paths scanned when inside roots',
        hash_selection='all sibling files; other files with exact expected size', errors=ERRORS), indent=2))

def hash_one(item):
    p, size, scope = item
    try:
        return str(p), size, digest(p), scope, ''
    except OSError as e:
        return str(p), size, '', scope, str(e)

def scan_files():
    started = time.monotonic()
    count = total = 0
    iterator = iter(candidates())
    with (A / 'FILE_SHA256_INDEX.csv').open('w', newline='') as f, concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        writer = csv.writer(f)
        writer.writerow(['absolute_path', 'size', 'sha256', 'scope', 'error'])
        while True:
            batch = []
            for _ in range(256):
                try:
                    batch.append(next(iterator))
                except StopIteration:
                    break
            if not batch:
                break
            for result in pool.map(hash_one, batch):
                writer.writerow(result)
                path, size, sha, scope, error = result
                count += 1
                total += size
                if error:
                    ERRORS.append(error)
                if sha in WANTED and any(not (ROOT / r['path']).exists() for r in WANTED[sha]):
                    staged_restore(Path(path).read_bytes(), sha, path)
            if count % 4096 == 0:
                print(json.dumps(dict(stage='files', hashed=count, bytes=total,
                                      recovered=len(MATCHES), elapsed=time.monotonic()-started)), flush=True)
    return dict(files_hashed=count, bytes_hashed=total, elapsed_seconds=time.monotonic()-started)

def scan_git():
    discovery = json.loads((A / 'FILESYSTEM_DISCOVERY.json').read_text())
    records = []
    for root in discovery['git_roots']:
        env = dict(os.environ, GIT_OPTIONAL_LOCKS='0', GIT_NO_LAZY_FETCH='1')
        result = subprocess.run(['git', '-C', root, 'cat-file', '--batch-all-objects',
                                 '--batch-check=%(objectname) %(objecttype) %(objectsize)'],
                                env=env, capture_output=True, text=True)
        matched_size = 0
        for line in result.stdout.splitlines():
            oid, kind, size = line.split()
            if kind != 'blob' or int(size) not in SIZES:
                continue
            matched_size += 1
            blob = subprocess.run(['git', '-C', root, 'cat-file', 'blob', oid],
                                  env=env, capture_output=True)
            if blob.returncode:
                continue
            sha = hashlib.sha256(blob.stdout).hexdigest()
            if sha in WANTED:
                staged_restore(blob.stdout, sha, 'git:' + root + ':' + oid)
        records.append(dict(repository=root, exit_code=result.returncode,
                            object_count=len(result.stdout.splitlines()), size_candidate_blobs=matched_size,
                            error=result.stderr.strip()))
    (A / 'GIT_OBJECT_SEARCH.json').write_text(json.dumps(records, indent=2))
    return records

def scan_archives():
    paths = json.loads((A / 'FILESYSTEM_DISCOVERY.json').read_text())['candidate_archives']
    records = []
    for item in paths:
        path = item if isinstance(item, str) else item['path']
        record = dict(path=path, members_checked=0, size_candidates=0)
        try:
            if zipfile.is_zipfile(path):
                with zipfile.ZipFile(path) as z:
                    for member in z.infolist():
                        record['members_checked'] += 1
                        if not member.is_dir() and member.file_size in SIZES:
                            record['size_candidates'] += 1
                            data = z.read(member)
                            sha = hashlib.sha256(data).hexdigest()
                            if sha in WANTED:
                                staged_restore(data, sha, path + ':' + member.filename)
            elif tarfile.is_tarfile(path):
                with tarfile.open(path, 'r|*') as t:
                    for member in t:
                        record['members_checked'] += 1
                        if member.isfile() and member.size in SIZES:
                            record['size_candidates'] += 1
                            data = t.extractfile(member).read()
                            sha = hashlib.sha256(data).hexdigest()
                            if sha in WANTED:
                                staged_restore(data, sha, path + ':' + member.name)
            else:
                record['unsupported_format'] = True
        except Exception as e:
            record['error'] = repr(e)
        records.append(record)
    (A / 'ARCHIVE_SEARCH.json').write_text(json.dumps(records, indent=2))
    return records

if __name__ == '__main__':
    info = dict(files=scan_files())
    print('file search complete; git object search starting', flush=True)
    info['git_repositories'] = len(scan_git())
    print('git search complete; archive member search starting', flush=True)
    info['archives'] = len(scan_archives())
    info['new_exact_recoveries'] = MATCHES
    info['errors'] = ERRORS
    (A / 'HASH_RECOVERY_RESULT.json').write_text(json.dumps(info, indent=2))
    print(json.dumps(info), flush=True)
