"""Qualify archived q-sign groups by actual bounds; no state/solver evaluation."""
from scripts import c7_early_branch as a
from pathlib import Path
import numpy as np
modified=[];cache={}
def visit(obj):
    if isinstance(obj,dict):
        label=obj.get('label','')
        if label.startswith(('branch/','c5_regression/')) and 'q_counts' in obj:
            if label not in cache:
                data=a.arrays(a.OUT/label/'STATE.npz');groups=a.bound_groups(data['s'],data['lower'],data['q_scaled'],.95)
                cache[label]=dict(q_counts={k:int(v.sum()) for k,v in groups.items()},q_IDs={k:np.flatnonzero(v).tolist() for k,v in groups.items()})
            obj.update(cache[label]);assert obj['q_counts']['WEAK_LOWER']==obj['weak_count'];assert obj['q_counts']['INTERIOR']==obj['interior_count']
        for value in obj.values():visit(value)
    elif isinstance(obj,list):
        for value in obj:visit(value)
for root in [a.OUT/'branch',a.OUT/'c5_regression']:
    for p in root.rglob('*.json'):
        r=a.read(p);before=a.json.dumps(r,sort_keys=True);visit(r)
        if a.json.dumps(r,sort_keys=True)!=before:a.save(p,r);modified.append(str(p.relative_to(a.OUT)))
a.save(a.OUT/'provenance/CONSTRAINT_COUNTS_RECONCILIATION.json',dict(files=modified,states=len(cache),scope='q_counts/q_IDs initially described derivative signs only. Positive-state near-zero q is now qualified by actual phase bounds, with INTERIOR explicit. Existing interior_count and weak_count agree exactly; no state vectors, residuals, acceptance, solver settings or physical values changed.'))
print('Constraint count metadata reconciled:',len(cache),'states')
