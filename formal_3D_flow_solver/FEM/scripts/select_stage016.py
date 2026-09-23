#!/usr/bin/env python3
"""Apply the pre-run frozen policy; rejected surfaces cannot acquire a volume."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import write_json,sha256,timestamp
from fem3d.planar_port import select_candidate
R=ROOT/'reports/stage01_6';O=ROOT/'outputs/stage01_6'
policy=json.loads((R/'acceptance_policy.json').read_text()); lock=json.loads((R/'freeze_lock.json').read_text())
assert sha256(R/'acceptance_policy.json')==lock['policy_sha256']
rows=[]
for name in policy['candidates']:
 s=json.loads((O/name/'qc/surface_invariants.json').read_text())
 failures=[k for k,v in s['cap_gate']['checks'].items() if not v]
 row={'candidate':name,'surface_status':s['status'],'geometry_status':'PASS' if all(p['geometry']['status']=='PASS' for p in s['ports'].values()) else 'FAIL',
  'surface_qc':s,'rejection_reasons':failures,'volume_status':'NOT RUN: surface gate failed','volume_gate_checks':{},'proxy':None,'quality':None}
 if s['status']=='PASS':
  v=json.loads((O/name/'qc/volume_quality.json').read_text())
  row.update(volume_status=v['status'],volume_gate_checks=v['volume_gate']['checks'],proxy=v['proxy'],quality=v['quality'])
 else:
  assert not (O/name/'mesh/fluid.msh').exists() and not (O/name/'mesh/volume_mesh.npz').exists()
 rows.append(row)
selection=select_candidate(rows)
selection.update(timestamp=timestamp(),policy_sha256=lock['policy_sha256'],contract_sha256=lock['contract_sha256'],parameter_search_expanded=False,fem_solved=False,stage3_started=False)
if selection['selected_candidate'] is not None: raise RuntimeError('Selected export and roundtrip must be completed before final reporting')
selection['unperformed']=['3D tetra meshing','candidate tetra QC and P2 proxy','selected export','DOLFINx 1-rank reload','DOLFINx 2-rank reload']
write_json(R/'candidate_selection.json',selection)
b=json.loads((O/'baseline/stage1_quality.json').read_text())
write_json(R/'quality_comparison.json',{'status':'FAIL','baseline':b,'candidates':{v['candidate']:{'cap_quality':v['surface_qc']['combined_cap_quality'],'tetra_quality':v['quality'],'surface_status':v['surface_status']} for v in rows},'selected_candidate':None,'conclusion':'All new planar geometry and cap quality gates pass; density gates reject all candidates. No new tetra quality or FEM cost claim is permitted.'})
write_json(R/'cost_comparison.json',{'baseline':policy['cost']['baseline'],'baseline_cap_triangles':191,'tetra_max':200000,'P2_ratio_max':1.35,'candidates':{v['candidate']:{'cap_triangles':v['surface_qc']['total_cap_triangles'],'proxy':v['proxy'],'reason':'NOT RUN: surface density gate failed'} for v in rows},'selected_candidate':None,'FEM_space_created':False})
print('STAGE 1.6 STATUS: FAIL; selected = NONE; volume / roundtrip = NOT RUN')
