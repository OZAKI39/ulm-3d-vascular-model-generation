"""Copy an immutable local clot project; inventory omissions and split large data.

Run on a fresh export destination. Source files are never modified or followed
through symlinks. Chunking is byte-for-byte concatenation, not scientific data
compression, decimation or re-encoding.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

CHUNK_SIZE=40*1024*1024


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def omission(rel,size):
    parts=rel.parts
    if parts[0].startswith('build') or any(p in {'.tools','.pytest_cache','__pycache__','.git'} for p in parts):
        return 'Build products, installed tools, caches or temporary test files'
    if rel.suffix in {'.pyc','.nbc','.nbi'}:
        return 'Regenerable Python/Numba cache'
    if parts[0]=='references' and rel.name in {'1-s2.0-S2666496826000567-main.pdf','paper_extracted.txt'}:
        return 'Third-party paper full text; retain DOI and implementation mapping'
    if parts[:2]==('visualization','fragmentation_layout_draft'):
        return 'Superseded visualization layout draft'
    if 'vtk' in parts and parts[0] in {'results','runs','verification'}:
        return 'Regenerable per-state VTK export; canonical NPZ states retained'
    if rel.suffix=='.pvd' and parts[0] in {'results','runs','verification'}:
        return 'Index of omitted regenerable VTK exports'
    if parts[0]=='visualization' and rel.name=='display_samples.npz':
        return 'Derived display interpolation; raw solver states and frame maps retained'
    if parts[0]=='visualization' and rel.suffix=='.gif' and size>10*1024*1024:
        return 'Large duplicate legacy GIF; corresponding MP4 and renderer retained'
    return None


def export(source,destination):
    source=source.resolve();destination=destination.resolve()
    if source==destination or source in destination.parents:raise ValueError('Export outside the source project')
    destination.mkdir(parents=True,exist_ok=True);entries=[]
    for p in sorted(source.rglob('*')):
        rel=p.relative_to(source)
        if p.is_symlink():
            target=str(p.readlink());size=len(target.encode());digest=hashlib.sha256(target.encode()).hexdigest();kind='symlink'
        elif p.is_file():size=p.stat().st_size;digest=sha(p);kind='file'
        else:continue
        row=dict(path=rel.as_posix(),kind=kind,bytes=size,sha256=digest)
        reason=omission(rel,size)
        if reason:row.update(disposition='omitted',reason=reason)
        else:
            out=destination/rel
            if out.exists() or out.is_symlink():raise FileExistsError(out)
            out.parent.mkdir(parents=True,exist_ok=True)
            if kind=='symlink':
                row.update(disposition='copied',target=target);out.symlink_to(target)
            elif size>48*1024*1024:
                chunks=[]
                with p.open('rb') as stream:
                    k=0
                    for block in iter(lambda:stream.read(CHUNK_SIZE),b''):
                        part=out.with_name(out.name+f'.part{k:03d}')
                        if part.exists():raise FileExistsError(part)
                        part.write_bytes(block);chunks.append(dict(path=part.relative_to(destination).as_posix(),bytes=len(block),sha256=hashlib.sha256(block).hexdigest()));k+=1
                row.update(disposition='chunked',parts=chunks,reason='Lossless byte chunks below GitHub file limits; restore with github_sync/restore_large_files.py')
            else:
                shutil.copy2(p,out);assert sha(out)==digest
                row['disposition']='copied'
        entries.append(row)
    counts=Counter(e['disposition'] for e in entries);totals=defaultdict(int);omissions=defaultdict(lambda:dict(files=0,bytes=0))
    for e in entries:
        totals[e['disposition']]+=e['bytes']
        if e['disposition']=='omitted':omissions[e['reason']]['files']+=1;omissions[e['reason']]['bytes']+=e['bytes']
    manifest=dict(format_version=1,created_utc=datetime.now(timezone.utc).isoformat(),source_root=str(source),
                  canonical_data_policy='All scientific NPZ states retained exactly, with chunking when needed; no downsampling or changed values',
                  counts=dict(counts),bytes=dict(totals),omission_groups=dict(omissions),files=entries)
    support=destination/'github_sync';support.mkdir(exist_ok=True)
    (support/'SNAPSHOT_MANIFEST.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    ignored=['# Scope-local rules; do not change repository-wide exclusions.','build*/','.tools/','.pytest_cache/','__pycache__/','*.py[cod]','*.nbc','*.nbi','.venv/','*.restoring']
    ignored+=['!/'+e['path'] for e in entries if e['disposition']=='copied' and Path(e['path']).suffix in {'.npz','.vtp','.vtu','.stl','.swc'}]
    ignored+=['/'+e['path'] for e in entries if e['disposition']=='chunked']
    (destination/'.gitignore').write_text('\n'.join(ignored)+'\n')
    (destination/'.gitattributes').write_text('# Preserve archived bytes, including original CRLF prompts.\n* -text\n*.npz -diff\n*.part* -diff\n')
    print(json.dumps({k:manifest[k] for k in ['counts','bytes','omission_groups']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--destination',type=Path,required=True)
    a=p.parse_args();export(a.source,a.destination)
