#!/usr/bin/env python3
"""Apply immutable eligibility/order; uncomputed quantities stay null."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import write_json,timestamp,sha256
from fem3d.cap_selection import select_candidate
R=ROOT/'reports/stage01_5';O=ROOT/'outputs/stage01_5'
policy=json.loads((R/'acceptance_policy.json').read_text());lock=json.loads((R/'acceptance_policy_lock.json').read_text())
assert sha256(R/'acceptance_policy.json')==lock['sha256']
rows=[]
for candidate in policy['candidates']:
    record=json.loads((O/candidate/'qc/surface_invariants.json').read_text())
    # This execution must not invent a missing volume result after geometry rejection.
    if record['status']=='PASS':
        raise RuntimeError('Geometry-admissible candidate requires the separately authorized volume/QC stage before selection')
    assert not (O/candidate/'mesh/fluid.msh').exists()
    failures={name:[k for k,v in p['geometry']['checks'].items() if not v] for name,p in record['ports'].items()}
    rows.append({'candidate':candidate,'geometry_status':record['status'],'surface_qc':record,'volume_qc':None,
        'volume_status':'NOT RUN: geometry gate failed','rejection_reasons_by_port':failures})
selection=select_candidate(rows,policy)
selection.update(timestamp=timestamp(),acceptance_policy_sha256=lock['sha256'],eligible_candidates=[],
    candidate_geometry={r['candidate']:{'status':r['geometry_status'],'reasons':r['rejection_reasons_by_port']} for r in rows},
    selection_rule_evaluation=[{'priority':1,'rule':'All geometry hard gates must pass','outcome':'A/B/C rejected by cap area invariant'},
        {'priority':2,'rule':'Minimum cap-adjacent q<0.1','outcome':'NOT EVALUATED: no volume meshing permitted'},
        {'priority':3,'rule':'Minimum total q<0.1','outcome':'NOT EVALUATED'},
        {'priority':4,'rule':'Maximum P1 minSICN','outcome':'NOT EVALUATED'},
        {'priority':5,'rule':'Maximum P5 minSICN','outcome':'NOT EVALUATED'},
        {'priority':6,'rule':'Smaller tetra count','outcome':'NOT EVALUATED'}],
    unperformed=['3D tetra meshing','new tetra QC','winner export','DOLFINx 1-rank reload','DOLFINx 2-rank reload'],
    parameter_search_expanded=False,fem_solved=False,stage3_started=False)
write_json(R/'candidate_selection.json',selection)
write_json(R/'quality_comparison.json',{'timestamp':timestamp(),'status':'FAIL','baseline':policy['baseline'],
    'candidates':{r['candidate']:r for r in rows},'selected_candidate':None,
    'conclusion':'Cap triangle quality improves, but the fixed 3D area invariant fails. No tetra improvement or causal proof is claimed.'})
print('STAGE 1.5 STATUS: FAIL; selected = NONE; volume and roundtrip = NOT RUN')
for r in rows:
    print(r['candidate'],{name:{'area_error':p['geometry']['relative_area_error'],'centroid_displacement_m':p['geometry']['centroid_displacement_m'],'normal_dot':p['geometry']['normal_dot_source'],'q_median':p['quality']['q_tri']['median'],'triangles':p['new_triangle_count'],'checks':p['geometry']['checks']} for name,p in r['surface_qc']['ports'].items()})
