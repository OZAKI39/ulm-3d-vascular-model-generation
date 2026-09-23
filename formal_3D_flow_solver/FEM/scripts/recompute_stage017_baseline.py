#!/usr/bin/env python3
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from fem3d.audit import sha256,write_json,timestamp
from fem3d.adaptive_qc import measure_volume
r=ROOT/'reports/stage01_7';o=ROOT/'outputs/stage01_7/baseline';o.mkdir(exist_ok=True)
contract_path=ROOT/'reports/stage01_6/planar_port_contract_v2.json'
ref=json.loads((r/'contract_reference.json').read_text());assert sha256(contract_path)==ref['sha256']
contract=json.loads(contract_path.read_text());gp=json.loads((ROOT/'reports/stage01_6/acceptance_policy.json').read_text())['geometry']
base=ROOT/'outputs/stage01/medium';mp=base/'mesh/volume_mesh.npz';sp=ROOT/'inputs/stage01/tagged_surface_si.npz'
history=json.loads((base/'qc/geometry_qc.json').read_text());metadata=json.loads((base/'metadata/meshing.json').read_text())
s=np.load(sp);m=np.load(mp);result=measure_volume(m,s,s,contract,{'geometry':gp})
q=result['quality'];checks={'vertices':result['proxy']['N_vertex']==history['vertex_count'],'tetra':result['proxy']['N_tetra']==history['tetrahedron_count'],
'quantiles':q['min_sicn']==history['quality']['gmsh_min_sicn'],'low_total':q['total_below_0_1']==history['quality']['advisory_count'],
'by_boundary':q['low_quality_nearest_boundary_counts']==history['quality']['low_quality_nearest_boundary_counts'],'validity':result['validity']==history['validity'],
'worst_cell':q['worst_elements'][0]['cell_index']==history['quality']['worst_elements'][0]['cell_index'],'mesh_hash':sha256(mp)==history['mesh_sha256']}
assert all(checks.values()),checks
result.update(timestamp=timestamp(),history_checks=checks,status='PASS',mesh_sha256=sha256(mp),source_surface_sha256=sha256(sp),contract_sha256=sha256(contract_path),
 h_volume_m=metadata['effective_gmsh_options']['Mesh.MeshSizeMax'],effective_gmsh_options=metadata['effective_gmsh_options'])
write_json(r/'baseline_recomputed.json',result);write_json(o/'baseline_recomputed.json',result);write_json(ROOT/'inputs/stage01_7/baseline_recomputed.json',result)
print(json.dumps({'checks':checks,'proxy':result['proxy'],'quality':q['min_sicn'],'low_total':q['total_below_0_1'],'low_cap':q['cap_adjacent_below_0_1']},indent=2))
