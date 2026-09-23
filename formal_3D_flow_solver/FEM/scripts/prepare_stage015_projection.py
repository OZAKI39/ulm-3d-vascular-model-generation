#!/usr/bin/env python3
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from fem3d.cap_remesh import plane_basis,project,signed_area,quality_summary,triangle_quality,coverage_check
from fem3d.mesh_qc import triangle_geometry
from fem3d.audit import write_json,timestamp,sha256
b=json.loads((ROOT/'inputs/stage01_5/baseline_contract.json').read_text())
s=np.load(ROOT/'inputs/stage01/tagged_surface_si.npz');p,t,tags=s['points_m'],s['triangles'],s['facet_tags']
records={}; arrays={}
for port in b['ports']:
    name=port['name']; ids=np.array(b['rims'][name]['original_vertex_ids']); basis=plane_basis(port['outward_normal']); origin=port['plane_origin_m']
    xy,offset=project(p,origin,basis)
    if signed_area(xy[ids])<0: ids=ids[np.r_[0,np.arange(len(ids)-1,0,-1)]]
    tri=t[tags==port['surface_entity_id']]; a,c,v=triangle_geometry(p,tri)
    coverage=coverage_check(xy,tri,ids)
    records[name]={'name':name,'entity_id':port['surface_entity_id'],'origin_m':origin,'basis':np.asarray(basis).tolist(),'ccw_rim_ids':ids.tolist(),'rim_xy_m':xy[ids].tolist(),'h_rim_m':b['rims'][name]['h_rim_m'],'q_tri':quality_summary(p,tri),'source_area_m2':float(a.sum()),'projected_polygon_area_m2':signed_area(xy[ids]),'source_area_minus_projected_relative':float((a.sum()-signed_area(xy[ids]))/port['area_m2']),'maximum_rim_plane_deviation_m':float(abs(offset[ids]).max()),'source_interior_vertex_count':len(np.setdiff1d(np.unique(tri),ids)),'coverage':coverage}
    arrays[name+'_q_tri'],arrays[name+'_edge_ratio']=triangle_quality(p,tri)
write_json(ROOT/'outputs/stage01_5/baseline/cap_projection.json',{'timestamp':timestamp(),'ports':records})
write_json(ROOT/'inputs/stage01_5/cap_projection.json',{'timestamp':timestamp(),'ports':records,'baseline_contract_sha256':sha256(ROOT/'inputs/stage01_5/baseline_contract.json')})
np.savez_compressed(ROOT/'outputs/stage01_5/baseline/cap_quality.npz',**arrays)
print({n:{'P5':d['q_tri']['q_tri']['P5'],'median':d['q_tri']['q_tri']['median'],'planarity_m':d['maximum_rim_plane_deviation_m']} for n,d in records.items()})
