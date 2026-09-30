"""Verify copied scientific files and chunks against the captured source hashes."""
import hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    m=json.loads((ROOT/'github_sync/SNAPSHOT_MANIFEST.json').read_text());errors=[];checked=0
    for e in m['files']:
        if e['disposition']=='omitted':continue
        if e['disposition']=='chunked':
            total=hashlib.sha256()
            for part in e['parts']:
                p=ROOT/part['path']
                if not p.is_file() or p.stat().st_size!=part['bytes'] or sha(p)!=part['sha256']:errors.append(part['path']);continue
                with p.open('rb') as f:
                    for b in iter(lambda:f.read(4*1024*1024),b''):total.update(b)
            if total.hexdigest()!=e['sha256']:errors.append(e['path'])
        else:
            p=ROOT/e['path']
            if e['kind']=='symlink':
                if not p.is_symlink() or str(p.readlink())!=e['target'] or not p.exists():errors.append(e['path'])
            elif not p.is_file() or p.stat().st_size!=e['bytes'] or sha(p)!=e['sha256']:errors.append(e['path'])
        checked+=1
    print(json.dumps(dict(status='FAIL' if errors else 'PASS',checked_source_entries=checked,errors=errors),indent=2))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
