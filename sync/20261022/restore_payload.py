#!/usr/bin/env python3
"""Verify the snapshot; optionally restore exact-byte gzip exports in this checkout."""
import argparse,gzip,hashlib,json,os,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def digest(p,compressed=False):
    h=hashlib.sha256()
    with (gzip.open(p,'rb') if compressed else p.open('rb')) as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--restore',action='store_true');a=parser.parse_args()
    rows=json.loads((ROOT/'sync/20261022/COPIED_FILES.json').read_text());n=0
    for r in rows:
        p=ROOT/r['path'];assert digest(p)==r['sha256'],str(p)
        if r['encoding']=='gzip_exact_original_bytes':
            assert digest(p,True)==r['source_sha256'],str(p)
            if a.restore:
                target=ROOT/r.get('original_path',str(p.relative_to(ROOT))[:-3])
                if target.exists():assert digest(target)==r['source_sha256'],f'Refusing to overwrite changed file: {target}'
                else:
                    with tempfile.NamedTemporaryFile(dir=target.parent,delete=False) as out:
                        tmp=Path(out.name)
                        with gzip.open(p,'rb') as src:
                            for b in iter(lambda:src.read(4*1024*1024),b''):out.write(b)
                    os.chmod(tmp,int(r['source_mode'],8));tmp.replace(target)
                n+=1
    print(json.dumps(dict(verified_files=len(rows),restored_files=n)))
if __name__=='__main__':main()
