#!/usr/bin/env python3
"""Check preserved files against the pre-migration inventory, without writes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from vascular_processing.project_paths import digest


def verify(full=False):
    audit=ROOT/'migration/20260926'
    plan=json.loads((audit/'migration_plan.json').read_text());old=Path(plan['source'])
    edits=json.loads((audit/'path_compatibility_edits.json').read_text())
    retained={r['path'] for r in json.loads((audit/'scope_corrections.json').read_text())}
    rows=[json.loads(line) for line in (audit/'inventory_before.jsonl').read_text().splitlines()]
    if not full:
        rows=[r for r in rows if r['path'].startswith(('vascular_processing/','config/','configs/','tools/','tests/','s1-',
            'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/final_abs_casting_mold/'))]
    def check(row):
        rel=row['path'];path=(old if rel in retained else ROOT)/rel
        if row['kind']=='symlink':
            same_target=path.is_symlink() and path.resolve()==(Path(row['target']) if Path(row['target']).is_absolute() else path.parent/row['target']).resolve()
            return rel,'UNCHANGED_TARGET' if same_target else 'MISMATCH'
        if not path.is_file():return rel,'MISSING'
        value=digest(path)
        if value==row['sha256']:return rel,'UNCHANGED'
        edit=edits.get(rel,{})
        if (edit.get('before_sha256')==row['sha256'] and edit.get('after_sha256')==value
            and digest(ROOT/edit['original_copy'])==row['sha256']):return rel,'DOCUMENTED_PATH_COMPATIBILITY_EDIT'
        return rel,'MISMATCH'
    with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(check,rows))
    changed=[dict(path=p,status=s) for p,s in results if s in ('MISSING','MISMATCH')]
    return dict(status='PASS' if not changed else 'FAIL',full_inventory=full,checked_files=len(rows),
        unchanged=sum(s.startswith('UNCHANGED') for _,s in results),
        documented_path_compatibility_edits=[p for p,s in results if s=='DOCUMENTED_PATH_COMPATIBILITY_EDIT'],
        unexpected_changes=changed)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--full',action='store_true',help='Also hash every original dataset and preserved output')
    args=parser.parse_args();result=verify(args.full);print(json.dumps(result,ensure_ascii=False,indent=2));raise SystemExit(0 if result['status']=='PASS' else 2)
