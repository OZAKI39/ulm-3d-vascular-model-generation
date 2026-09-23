#!/usr/bin/env python3
"""Run one predeclared cap-only candidate; reject before any volume generation."""
import argparse,json,os,platform,resource,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import gmsh
from fem3d.audit import sha256,timestamp,write_json
from fem3d.cap_remesh import check_rim,check_wall,project,lift,cross2,quality_summary,triangle_quality,coverage_check,port_geometry
from fem3d.mesh_input import surface_topology

p=argparse.ArgumentParser();p.add_argument('--candidate',required=True,choices=['candidate_A','candidate_B','candidate_C']);args=p.parse_args()
base=ROOT/'outputs/stage01_5'/args.candidate
if (base/'metadata/cap_remesh.json').exists(): raise RuntimeError('Candidate already recorded; no silent overwrite or parameter search')
if base.exists() and any(base.rglob('*')):
    import shutil
    archive=ROOT/'outputs/stage01_5/failed_attempts'/f'{args.candidate}_{time.time_ns()}'
    archive.parent.mkdir(parents=True,exist_ok=True)
    shutil.move(str(base),str(archive))
for sub in ('surface','qc','metadata'): (base/sub).mkdir(parents=True,exist_ok=True)
policy_path=ROOT/'inputs/stage01_5/acceptance_policy.json';policy=json.loads(policy_path.read_text())
b=json.loads((ROOT/'inputs/stage01_5/baseline_contract.json').read_text());planes=json.loads((ROOT/'inputs/stage01_5/cap_projection.json').read_text())['ports']
assert sha256(policy_path)==b['acceptance_policy_sha256']
source_path=ROOT/'inputs/stage01/tagged_surface_si.npz';assert sha256(source_path)==b['source_surface_sha256']
assert sha256(ROOT/'inputs/stage01/source_contract.json')==b['source_contract_sha256']
assert sha256(ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz')==b['stage1_medium_mesh_sha256']
assert sha256(ROOT/'outputs/stage01/medium/qc/geometry_qc.json')==b['stage1_baseline_qc_sha256']
s=np.load(source_path); points=s['points_m'].copy();old_tri=s['triangles'];old_tags=s['facet_tags']
triangles=[old_tri[old_tags==1]];tags=[old_tags[old_tags==1]];port_results={};quality_arrays={};start=time.perf_counter()
factor=policy['candidates'][args.candidate]
gmsh.initialize()
options={'General.NumThreads':1,'Mesh.MaxNumThreads1D':1,'Mesh.MaxNumThreads2D':1,'Mesh.Algorithm':6,'Mesh.ElementOrder':1,'Mesh.MeshSizeFromPoints':0,'Mesh.MeshSizeFromCurvature':0,'Mesh.MeshSizeExtendFromBoundary':0,'Mesh.Smoothing':0,'Geometry.Tolerance':1e-15}
try:
    for port in b['ports']:
        name=port['name']; plane=planes[name]; ids=np.asarray(plane['ccw_rim_ids']);xy=np.asarray(plane['rim_xy_m']);h=plane['h_rim_m']*factor
        gmsh.model.add(name)
        for key,value in options.items(): gmsh.option.setNumber(key,value)
        gmsh.option.setNumber('Mesh.MeshSizeMin',h);gmsh.option.setNumber('Mesh.MeshSizeMax',h)
        gp=[gmsh.model.geo.addPoint(float(v[0]),float(v[1]),0,h) for v in xy]
        lines=[gmsh.model.geo.addLine(gp[i],gp[(i+1)%len(gp)]) for i in range(len(gp))]
        loop=gmsh.model.geo.addCurveLoop(lines);face=gmsh.model.geo.addPlaneSurface([loop]);gmsh.model.geo.synchronize()
        for line in lines: gmsh.model.mesh.setTransfiniteCurve(line,2)
        gmsh.model.mesh.generate(2)
        node_to_source={};rim_node_ids=[]
        for tag,original in zip(gp,ids):
            node,coord,_=gmsh.model.mesh.getNodes(0,tag)
            assert len(node)==1
            node_to_source[int(node[0])]=int(original);rim_node_ids.append(int(node[0]))
        for i,line in enumerate(lines):
            assert len(gmsh.model.mesh.getNodes(1,line,includeBoundary=False)[0])==0,'Rim edge split'
            element_types,_,nodes=gmsh.model.mesh.getElements(1,line)
            assert list(element_types)==[1] and set(map(int,nodes[0]))=={rim_node_ids[i],rim_node_ids[(i+1)%len(ids)]}
        all_ids,coords,_=gmsh.model.mesh.getNodes();coords=coords.reshape(-1,3)
        new_ids=[]
        for node,coord in sorted(zip(all_ids,coords),key=lambda pair:pair[0]):
            node=int(node)
            if node in node_to_source: continue
            assert abs(coord[2])<=1e-20
            idx=len(points);node_to_source[node]=idx;new_ids.append(idx)
            points=np.vstack([points,lift(coord[:2],port['plane_origin_m'],np.asarray(plane['basis']))])
        types,_,nodes=gmsh.model.mesh.getElements(2,face);assert list(types)==[2]
        connectivity=np.asarray([node_to_source[int(v)] for v in nodes[0]],dtype=np.int64).reshape(-1,3)
        projected,_=project(points,port['plane_origin_m'],np.asarray(plane['basis']))
        c=projected[connectivity];flip=cross2(c[:,1]-c[:,0],c[:,2]-c[:,0])<0
        connectivity[flip]=connectivity[flip][:,[0,2,1]]
        rim=check_rim(s['points_m'],b['rims'][name]['rim_edges'],points,connectivity)
        coverage=coverage_check(projected,connectivity,ids,policy['geometry']['coverage_relative_roundoff_tolerance'])
        geometry=port_geometry(points,connectivity,port,np.asarray(new_ids,dtype=int),policy)
        quality=quality_summary(points,connectivity)
        port_results[name]={'entity_id':port['surface_entity_id'],'h_rim_m':plane['h_rim_m'],'h_cap_m':h,'old_triangle_count':int(np.count_nonzero(old_tags==port['surface_entity_id'])),'new_triangle_count':len(connectivity),'interior_vertex_count':len(new_ids),'interior_vertex_ids':new_ids,'rim_vertex_count':len(ids),'rim':rim,'coverage':coverage,'geometry':geometry,'quality':quality}
        quality_arrays[name+'_q_tri'],quality_arrays[name+'_edge_ratio']=triangle_quality(points,connectivity)
        triangles.append(connectivity);tags.append(np.full(len(connectivity),port['surface_entity_id'],dtype=np.int32))
        print(args.candidate,name,'geometry',geometry['status'],'area rel',geometry['relative_area_error'],'triangles',len(connectivity),flush=True)
        gmsh.model.remove()
finally:
    gmsh.finalize()
derived={'points_m':points,'triangles':np.concatenate(triangles),'facet_tags':np.concatenate(tags)}
wall=check_wall(s,derived);assert wall['wall_triangle_count']==67071
assert len(np.unique(derived['triangles']))==len(points),'Unused original points require explicit investigation'
topology=surface_topology(points,derived['triangles'])
status='PASS' if all(v['geometry']['status']=='PASS' for v in port_results.values()) else 'FAIL'
np.savez_compressed(base/'surface/tagged_surface_si.npz',**derived)
np.savez_compressed(base/'qc/cap_quality.npz',**quality_arrays)
record={'status':status,'candidate':args.candidate,'timestamp':timestamp(),'wall':wall,'ports':port_results,'topology':topology,'geometry_semantics':'derived CFD boundary triangulation','anatomical_wall_changed':False,'port_plane_changed':False,'port_rim_changed':False,'cap_interior_connectivity_changed':True,'volume_meshing_permitted':status=='PASS','volume_generated':False,'fem_solved':False}
write_json(base/'qc/surface_invariants.json',record)
write_json(base/'metadata/cap_remesh.json',{'run_id':os.environ.get('FEM3D_RUN_ID'),'timestamp':timestamp(),'hostname':platform.node(),'status':status,'candidate':args.candidate,'h_cap_over_h_rim':factor,'gmsh_version':gmsh.__version__,'gmsh_options':options,'rim_constraint':'Each original polygon edge is a straight CAD segment constrained to exactly 2 endpoint nodes; original SI rim coordinates copied bit for bit after 2D triangulation','source_stage0_contract_sha256':b['source_contract_sha256'],'stage1_medium_mesh_sha256':b['stage1_medium_mesh_sha256'],'stage1_baseline_qc_sha256':b['stage1_baseline_qc_sha256'],'code_source_manifest_sha256':sha256(ROOT/'remote/source_manifest.json'),'config_sha256':sha256(policy_path),'source_surface_sha256':sha256(source_path),'derived_surface_sha256':sha256(base/'surface/tagged_surface_si.npz'),'wall_time_s':time.perf_counter()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'gpu_used':False,'mpi_ranks':1,'volume_meshing_permitted':status=='PASS','port_results':port_results,'stdout':os.environ.get('FEM3D_STDOUT_LOG'),'stderr':os.environ.get('FEM3D_STDERR_LOG')})
print(args.candidate,': surface geometry',status,'; volume',('permitted' if status=='PASS' else 'BLOCKED by unchanged geometry gates'),flush=True)
