#!/usr/bin/env python3
"""Freeze v2 contract and acceptance policy BEFORE any sparse run."""
import json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from fem3d.audit import sha256,write_json,timestamp
from fem3d.planar_port import polygon_metrics,p2_proxy,validate_port
from fem3d.cap_remesh import check_wall,check_rim,quality_summary,project
from fem3d.mesh_input import FACET_NAMES,surface_topology
from fem3d.mesh_qc import triangle_geometry
r=ROOT/'reports/stage01_6'; inp=ROOT/'inputs/stage01_6'; out=ROOT/'outputs/stage01_6/baseline';out.mkdir(exist_ok=True)
assert (r/'PORT_CONTRACT_RATIONALE.md').exists()
assert not (r/'acceptance_policy.json').exists(),'Policy already frozen'
b=json.loads((ROOT/'inputs/stage01_5/baseline_contract.json').read_text()); planes=json.loads((ROOT/'inputs/stage01_5/cap_projection.json').read_text())['ports']
source_path=ROOT/'inputs/stage01/tagged_surface_si.npz'; mesh_path=ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz'
s=np.load(source_path); m=np.load(mesh_path)
assert sha256(source_path)==b['source_surface_sha256'] and sha256(mesh_path)==b['stage1_medium_mesh_sha256']
assert sha256(ROOT/'reports/stage00/source_contract.json')==b['source_contract_sha256']
ports={}
for old in b['ports']:
 name=old['name']; pl=planes[name];ids=pl['ccw_rim_ids'];tri=s['triangles'][s['facet_tags']==old['surface_entity_id']]
 a,c,v=triangle_geometry(s['points_m'],tri); metrics=polygon_metrics(s['points_m'][ids],old['plane_origin_m'],pl['basis'])
 xy,z=project(s['points_m'][ids],old['plane_origin_m'],pl['basis'])
 ports[name]={'name':name,'entity_id':old['surface_entity_id'],'plane_origin_m':old['plane_origin_m'],'outward_normal':old['outward_normal'],'basis':pl['basis'],
  'rim_vertex_ids':ids,'rim_coordinates_m':s['points_m'][ids].tolist(),'rim_edge_ids':b['rims'][name]['rim_edges'],
  'h_rim_m':pl['h_rim_m'],'legacy_max_nonplanarity_m':float(np.abs(z).max()),'original_cap_triangle_count':len(tri),
  'legacy_scalar_area_m2':float(a.sum()),'legacy_scalar_surface_area_m2':float(a.sum()),'legacy_scalar_area_centroid_m':(np.sum(a[:,None]*c,axis=0)/a.sum()).tolist(),
  **metrics,'relative_difference_legacy_scalar_vs_projected':float((a.sum()-metrics['formal_projected_area_m2'])/metrics['formal_projected_area_m2'])}
contract={'version':2,'timestamp':timestamp(),'source_contract_path':'reports/stage00/source_contract.json','source_contract_sha256':b['source_contract_sha256'],
 'source_surface_sha256':b['source_surface_sha256'],'facet_names':{str(k):v for k,v in FACET_NAMES.items()},'ports':ports,
 'formal_area_definition':'fixed rim projected polygon area','legacy_scalar_area_is_gate':False,
 'anatomical_wall_changed':False,'rim_changed':False,'port_plane_changed':False,'port_physical_region_changed':False,'only_cap_interior_triangulation_changed':True}
write_json(r/'planar_port_contract_v2.json',contract)
proxy=p2_proxy(m['tetra'])
options=json.loads((ROOT/'outputs/stage01/medium/metadata/meshing.json').read_text())['effective_gmsh_options']
policy={'stage':'1.6','timestamp':timestamp(),'candidates':{'sparse_A':3.5e-7,'sparse_B':4.5e-7,'sparse_C':5.5e-7},
 'sizing':{'formula':'h(d)=h_rim+(h_center-h_rim)*min(d/(0.5*sqrt(A_projected/pi)),1)','transition_R_eq':.5,
 'distance':'exact minimum Euclidean distance to frozen projected rim line segments','h_rim':'median frozen 3D rim edge length','implementation':'Gmsh setSizeCallback; two nodes per transfinite rim segment',
 'gmsh_reference':'https://gmsh.info/doc/texinfo/#gmsh_002fmodel_002fmesh_002fsetSizeCallback','algorithm_2d':6,'smoothing':0,'geometry_tolerance_m':1e-15},
 'geometry':{'wall_displacement_m_max':0.,'rim_displacement_m_max':0.,'projected_area_relative_error_max':1e-12,'vector_area_relative_error_max':1e-12,
 'normal_dot_min':1-1e-12,'projected_centroid_displacement_m_max':1e-12,'interior_plane_deviation_m_max':1e-15,'coverage_relative_roundoff_tolerance':1e-12,
 'rim_ids_and_edges_exact':True,'plane_origin_normal_basis_exact':True,'entity_mapping_exact':True,'connected_components':1,'boundary_edges':0,'nonmanifold_edges':0},
 'density':{'per_port_max':{'inlet':224,'outlet_01':176,'outlet_02':168,'outlet_03':196},'total_max':800},
 'cap_quality':{'q_below_0_1_max':0,'P5_min':.45,'median_min':.70,'positive_finite_area':True},
 'volume_meshing':{'bulk_target_m':5e-7,'effective_gmsh_options':options,'wall_time_limit_s':300,'virtual_memory_limit_gib':8,'no_local_volume_refinement':True},
 'tetra_quality':{'metric':'Gmsh minSICN','cap_adjacent_below_0_1_max':38,'total_below_0_1_max':76,'P1_min':.3359723496666575,'P5_min':.47439613749869236,'median_min':.7206460021439326,
 'zero_inverted_degenerate_nonfinite':True,'adjacency_method':'nearest exterior triangle centroid to tetrahedron centroid (same as Stage 1)','minimum_is_gate':False},
 'cost':{'tetra_max':200000,'P2_velocity_ratio_max':1.35,'baseline':proxy,'baseline_mesh_sha256':sha256(mesh_path)},
 'selection':['quality_and_cost_survivors_only','minimum_P2_velocity_proxy','minimum_tetra','minimum_cap_adjacent_low','minimum_total_low','maximum_P1'],
 'expansion_or_retuning_allowed':False,'contract_sha256':sha256(r/'planar_port_contract_v2.json')}
write_json(r/'acceptance_policy.json',policy)
lock={'timestamp':timestamp(),'policy_sha256':sha256(r/'acceptance_policy.json'),'contract_sha256':sha256(r/'planar_port_contract_v2.json'),
 'source_surface_sha256':sha256(source_path),'source_contract_sha256':b['source_contract_sha256'],'baseline_mesh_sha256':sha256(mesh_path),
 'baseline_qc_sha256':b['stage1_baseline_qc_sha256'],'frozen_before_any_sparse_candidate':True}
write_json(r/'freeze_lock.json',lock)
for name in ('planar_port_contract_v2.json','acceptance_policy.json','freeze_lock.json'): shutil.copyfile(r/name,inp/name)
# Required post-freeze baseline revalidation, independently of preparing values above.
wall=check_wall(s,s); topology=surface_topology(s['points_m'],s['triangles']); checks={}
for name,port in ports.items():
 tri=s['triangles'][s['facet_tags']==port['entity_id']]
 checks[name]=validate_port(s['points_m'],s['points_m'],tri,port,[],policy,origin=port['plane_origin_m'],normal=port['outward_normal'],basis=port['basis'])
 checks[name]['quality']=quality_summary(s['points_m'],tri)
assert p2_proxy(m['tetra'])==proxy
write_json(out/'geometry.json',{'status':'PASS','timestamp':timestamp(),'wall':wall,'topology':topology,'ports':checks,'contract_sha256':lock['contract_sha256']})
write_json(out/'cost_proxy.json',{'status':'PASS','mesh_sha256':sha256(mesh_path),'proxy':proxy,'method':'Unique vertices and undirected edges of all tetrahedra; no FEM space created'})
write_json(out/'stage1_quality.json',json.loads((ROOT/'outputs/stage01/medium/qc/geometry_qc.json').read_text()))
print(json.dumps({'freeze':lock,'baseline_proxy':proxy},indent=2))
