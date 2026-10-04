"""Recover serialized copies only when the original SHA256 matches exactly.

No model/solver runs. Candidate changes are limited to gzip's stored timestamp
and literal fixture-directory prefixes in otherwise existing JSON files.
"""
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
import struct
from audit_inventory import AUDIT, PROJECT, ROOT, dump_csv, dump_json


def main():
    unavailable=json.loads((AUDIT/'RECOVERY_UNAVAILABLE.json').read_text())
    rows=json.loads((AUDIT/'CLASSIFIED_INVENTORY.json').read_text());mapping={r['path']:r for r in rows}
    temp=Path(json.loads((AUDIT/'VALIDATION_TEMP.json').read_text())['path'])
    by_size=defaultdict(list)
    for base in [PROJECT,temp]:
        for p in base.rglob('*.gz'):
            if p.is_file():by_size[p.stat().st_size].append(p)
    recovered=[]
    for row in unavailable:
        target=ROOT/row['path'];desired=row['hash'];candidates=[]
        if target.suffix=='.gz':
            mtime=mapping[row['path']]['mtime_ns']//1_000_000_000
            for p in by_size[row['size_bytes']]:
                raw=p.read_bytes()
                if raw[:3]!=b'\x1f\x8b\x08':continue
                for second in range(mtime-3,mtime+2):
                    candidates.append((p,raw[:4]+struct.pack('<I',second)+raw[8:],'Existing compressed payload; restored gzip mtime from inventory timestamp window.'))
        elif target.suffix=='.json':
            parts=target.parts
            index=next((i for i,part in enumerate(parts) if part.startswith('test_') or part.startswith('erosion_coupled') or part.startswith('equilibrated_erosion')),None)
            if index is not None:
                suffix=Path(*parts[index:]);old_prefix=str(Path(*parts[:index]))
                for prefix in [temp/'baseline_pytest',temp/'final_pytest',
                               PROJECT/'verification/damage_energy_audit_v1/tests/full',
                               PROJECT/'verification/mechanical_equilibrium_v1/adr_completion/test_artifacts']:
                    p=prefix/suffix
                    if p.is_file():
                        candidates.append((p,p.read_bytes().replace(str(prefix).encode(),old_prefix.encode()),'Existing JSON with literal fixture-root prefix restored.'))
        for source,data,method in candidates:
            if len(data)==row['size_bytes'] and hashlib.sha256(data).hexdigest()==desired:
                target.parent.mkdir(parents=True,exist_ok=True)
                with target.open('xb') as f:f.write(data)
                ns=mapping[row['path']]['mtime_ns'];os.utime(target,ns=(ns,ns))
                recovered.append(dict(path=row['path'],source=str(source),hash=desired,size_bytes=len(data),method=method,status='RESTORED_EXACT_HASH'))
                break
    remaining=[r for r in unavailable if not (ROOT/r['path']).exists()]
    dump_csv('SERIALIZED_ARTIFACT_RESTORATION_LOG.csv',recovered,['path','source','hash','size_bytes','method','status'])
    dump_json('RECOVERY_STILL_UNAVAILABLE.json',remaining)
    print('Exact serialized restorations:',len(recovered),'remaining:',len(remaining),flush=True)


if __name__=='__main__':main()
