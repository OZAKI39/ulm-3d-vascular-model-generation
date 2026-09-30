"""Restore authoritative NPZ files from verified byte chunks without overwrites."""
import argparse,hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def blocks(path):
    with path.open('rb') as f:
        yield from iter(lambda:f.read(4*1024*1024),b'')


def sha(path):
    h=hashlib.sha256()
    for b in blocks(path):h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify-only',action='store_true');a=p.parse_args()
    m=json.loads((ROOT/'github_sync/SNAPSHOT_MANIFEST.json').read_text());done=[]
    for e in m['files']:
        if e['disposition']!='chunked':continue
        output=ROOT/e['path'];total=hashlib.sha256();size=0
        for part in e['parts']:
            path=ROOT/part['path']
            if path.stat().st_size!=part['bytes'] or sha(path)!=part['sha256']:raise ValueError(f'Invalid chunk: {path}')
            for b in blocks(path):total.update(b);size+=len(b)
        if size!=e['bytes'] or total.hexdigest()!=e['sha256']:raise ValueError(f'Invalid concatenation: {output}')
        if not a.verify_only:
            if output.exists():
                if output.stat().st_size!=size or sha(output)!=e['sha256']:raise FileExistsError(f'Refusing to overwrite different data: {output}')
            else:
                temporary=output.with_name(output.name+'.restoring')
                with temporary.open('xb') as f:
                    for part in e['parts']:
                        for b in blocks(ROOT/part['path']):f.write(b)
                if sha(temporary)!=e['sha256']:raise ValueError('Restored file failed checksum')
                temporary.rename(output)
        done.append(dict(path=e['path'],bytes=size,sha256=e['sha256'],restored=not a.verify_only))
    print(json.dumps(dict(status='PASS',files=done),indent=2))


if __name__=='__main__':main()
