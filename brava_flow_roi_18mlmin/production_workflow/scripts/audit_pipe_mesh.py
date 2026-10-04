"""Audit native cell types, geometry, layer heights and fixed analytical probe coverage."""
from pathlib import Path
import argparse,json,hashlib,itertools
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parents[1];B=R.parent
import sys
sys.path.insert(0,str(B/'vendor'))
from vascular_validation.mesh_diagnostics import tetra_sicn

def run(name):
 d=R/'straight_pipe'/name;cfg=json.loads((d/'INPUT.json').read_text());proto=json.loads((R/'straight_pipe/BENCHMARK_PROTOCOL.json').read_text());grid=pv.read(d/'volume.vtu');types,counts=np.unique(grid.celltypes,return_counts=True);rad=cfg['radius_m'];length=cfg['length_m'];points=np.asarray(grid.points,dtype=float)
 if set(types)!={10}:
  qa=dict(status='FAIL_UNSUPPORTED_NATIVE_CELL_TYPES',native_cell_types=dict(zip(map(str,types),map(int,counts))),production_particle_field_accepts='linear tetra only');(d/'MESH_QA.json').write_text(json.dumps(qa,indent=2)+'\n');return qa
 tet=grid.cells_dict[10];v=np.linalg.det(points[tet[:,1:]]-points[tet[:,:1]])/6;quality=tetra_sicn(points,tet)
 faces=np.concatenate([tet[:,ids] for ids in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]]);owner=np.tile(np.arange(len(tet)),4);keys,first,count=np.unique(np.sort(faces,axis=1),axis=0,return_index=True,return_counts=True);mask=count==1;boundary=faces[first[mask]].copy();adjacent=owner[first[mask]]
 tri=points[boundary];av=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])/2;flip=np.einsum('ni,ni->n',av,tri.mean(1)-points[tet[adjacent]].mean(1))<0;boundary[flip]=boundary[flip][:,[0,2,1]];av[flip]*=-1
 lookup={tuple(sorted(t)):i for i,t in enumerate(boundary)};tree=cKDTree(points);tags=np.zeros(len(boundary),np.int32);errors=[];mapping=0.
 for tag,label in [(1,'WALL'),(4,'INLET'),(3,'OUTLET_01')]:
  surf=pv.read(d/(label+'.vtp'));dist,ids=tree.query(surf.points);mapping=max(mapping,float(dist.max()))
  ft=ids[surf.faces.reshape(-1,4)[:,1:]]
  for face in ft:
   k=lookup.get(tuple(sorted(face)))
   if k is None:errors.append(label+' facet is not exterior');continue
   if tags[k]:errors.append('duplicate role')
   tags[k]=tag
 wallnodes=np.unique(boundary[tags==1]);allwall=np.isin(tet,wallnodes).all(1);allwallcount=int(allwall.sum());norm=av/np.linalg.norm(av,axis=1)[:,None]
 wf=boundary[tags==1];wo=adjacent[tags==1];n=norm[tags==1];cent=points[wf].mean(1);distances=abs(np.einsum('nvi,ni->nv',points[tet[wo]]-cent[:,None,:],n));height=distances.max(1);middle=(cent[:,2]>.25*length)&(cent[:,2]<.75*length)
 # Fixed physical positions inside the exact cylinder; no mesh-relative adjustment.
 theta=2*np.pi*(np.arange(proto['azimuth_count'])+proto['azimuth_offset'])/proto['azimuth_count'];levels=proto['wall_distances_m']+[f*rad for f in proto['relative_wall_distances_R']];probe=[];groups=[]
 for y in levels:
  for z in proto['analysis_z_over_L']:
   for t in theta:probe.append([(rad-y)*np.cos(t),(rad-y)*np.sin(t),z*length]);groups.append(y)
 probe=np.array(probe);cell=grid.find_containing_cell(probe);groups=np.array(groups);coverage=[dict(wall_distance_m=float(y),inside=int(np.sum((groups==y)&(cell>=0))),samples=int(np.sum(groups==y))) for y in levels]
 radius_vertices=np.linalg.norm(points[wallnodes,:2],axis=1);wall_radius=np.linalg.norm(cent[middle,:2],axis=1);volume_ref=np.pi*rad**2*length
 qa=dict(status='PASS' if not(errors or np.any(tags==0) or np.any(count>2) or np.any(v<=0) or allwallcount or mapping>1e-9 or len(tet)>proto['maximum_tetra_per_pipe'] or np.any(cell<0)) else 'FAIL_MESH_OR_ANALYTIC_PROBE_GATE',name=name,points=len(points),tetra=len(tet),native_cell_types=dict(zip(map(str,types),map(int,counts))),zero_volume=int(np.sum(v==0)),negative_Jacobian=int(np.sum(v<0)),all_WALL_fluid_tetra=allwallcount,nonmanifold_faces=int(np.sum(count>2)),unlabelled_boundary_faces=int(np.sum(tags==0)),boundary_mapping_max_m=mapping,boundary_errors=sorted(set(errors)),quality_min=float(quality.min()),quality_P1=float(np.quantile(quality,.01)),quality_median=float(np.median(quality)),volume_m3=float(v.sum()),analytic_cylinder_volume_m3=volume_ref,volume_relative_error=abs(float(v.sum())-volume_ref)/volume_ref,wall_node_radial_error_max_m=float(np.max(abs(radius_vertices-rad))),wall_facet_centroid_radial_error_max_m=float(np.max(abs(wall_radius-rad))),wall_facet_centroid_radial_error_median_m=float(np.median(abs(wall_radius-rad))),nearwall_first_tet_height_m=dict(minimum=float(height[middle].min()),median=float(np.median(height[middle])),P95=float(np.quantile(height[middle],.95))),requested_layers=cfg['policy']['layers'],layer_count_measurement='Native request only; actual radial vertex bands and cell entity IDs retained separately. Do not infer prism layers from method name.',point_arrays=list(grid.point_data),cell_arrays=list(grid.cell_data),probe_coverage=coverage,all_analytic_probes_inside=bool(np.all(cell>=0)),velocity_gradient_WSS_CFD_status='NOT_YET_RUN',mesh_sha256=hashlib.sha256((d/'volume.vtu').read_bytes()).hexdigest())
 np.savez_compressed(d/'MESH_DIAGNOSTICS.npz',points_m=points,tetra=tet,boundary_triangles=boundary,facet_tags=tags,adjacent_tetra=adjacent,first_tet_height_m=height,wall_facet_centers_m=cent,wall_normal_out=n,probes_m=probe,probe_wall_distance_m=groups,probe_cell=cell,all_wall_tetra_ids=np.flatnonzero(allwall),quality=quality)
 (d/'MESH_QA.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa,indent=2));return qa
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('case');args=ap.parse_args();run(args.case)
