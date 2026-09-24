#!/usr/bin/env python3
"""Verify frozen sources and extract original cap rims, without changing history."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from fem3d.audit import sha256,timestamp,write_json
from fem3d.mesh_input import read_contract,FACET_NAMES
from fem3d.cap_remesh import rim_loop

policy=json.loads((ROOT/'reports/stage01_5/acceptance_policy.json').read_text())
lock=json.loads((ROOT/'reports/stage01_5/acceptance_policy_lock.json').read_text())
assert sha256(ROOT/'reports/stage01_5/acceptance_policy.json')==lock['sha256']
contract,source_lock=read_contract(ROOT)
manifest=json.loads((ROOT/'reports/stage01/medium_artifact_manifest.json').read_text())
for name,digest in manifest['files'].items(): assert sha256(ROOT/name)==digest,name
source=ROOT/'inputs/stage01/tagged_surface_si.npz'
adapter=json.loads((ROOT/'inputs/stage01/surface_adapter.json').read_text())
assert sha256(source)==adapter['adapted_surface_sha256']
qc=json.loads((ROOT/'outputs/stage01/medium/qc/geometry_qc.json').read_text())
b=policy['baseline']; quality=qc['quality']
assert qc['vertex_count']==b['vertices'] and qc['tetrahedron_count']==b['tetrahedra']
assert qc['topology']['exterior_facet_count']==b['boundary_facets']
for k in ('P1','P5','median','P95'): assert quality['gmsh_min_sicn'][k]==b[k]
assert quality['gmsh_min_sicn']['minimum']==b['minSICN_minimum']
assert quality['advisory_count']==b['total_lt_0_1']
assert quality['low_quality_nearest_boundary_counts']==b['low_quality_by_boundary']
assert quality['worst_elements'][0]['cell_index']==b['worst_cell_index']
s=np.load(source); points,tri,tags=s['points_m'],s['triangles'],s['facet_tags']
assert np.count_nonzero(tags==1)==67071
rims={}
for p in contract['inlets']+contract['outlets']:
    loop,edges=rim_loop(tri[tags==p['surface_entity_id']])
    rims[p['name']]={'entity_id':p['surface_entity_id'],'original_vertex_ids':loop.tolist(),'original_coordinates_m':points[loop].tolist(),'rim_edges':edges.tolist(),'rim_vertex_count':len(loop),'rim_edge_count':len(edges),'connected_components':1,'closed_loop_degree':2,'h_rim_m':float(np.median(np.linalg.norm(points[edges[:,0]]-points[edges[:,1]],axis=1)))}
write_json(ROOT/'outputs/stage01_5/baseline/cap_rims.json',{'status':'PASS','timestamp':timestamp(),'ports':rims,'source_sha256':sha256(source)})
record={'status':'PASS','timestamp':timestamp(),'acceptance_policy_sha256':lock['sha256'],'source_contract_sha256':source_lock['source_contract_sha256'],'source_surface_sha256':sha256(source),'stage1_medium_mesh_sha256':sha256(ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz'),'stage1_baseline_qc_sha256':sha256(ROOT/'outputs/stage01/medium/qc/geometry_qc.json'),'verified_stage1_files':manifest['files'],'fixed_baseline':b,'facet_names':FACET_NAMES,'ports':contract['inlets']+contract['outlets'],'rims':rims}
write_json(ROOT/'inputs/stage01_5/baseline_contract.json',record)
write_json(ROOT/'reports/stage01_5/baseline_verification.json',record)
print('Frozen baseline hashes and four single closed rims verified')
