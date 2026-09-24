"""Stage 1.6 tetra geometry, minSICN quality and topological cost audits."""
import numpy as np
from scipy.spatial import cKDTree
from .mesh_input import FACET_NAMES
from .mesh_qc import boundary_partition,fidelity,triangle_geometry,tetra_volumes,quantiles
from .planar_port import validate_port,p2_proxy,volume_gates


def boundary_contract(points,exterior,tags,source,original,contract,policy):
    match=fidelity(points,exterior,tags,source['points_m'],source['triangles'],source['facet_tags'],0.)
    assert match['coordinates_exactly_equal']
    used=np.unique(exterior);_,ids=cKDTree(source['points_m']).query(points[used])
    mapping=np.full(len(points),-1,int);mapping[used]=ids
    restored=np.empty_like(source['points_m']);restored[ids]=points[used]
    ports={}
    for name,port in contract['ports'].items():
        tri=mapping[exterior[tags==port['entity_id']]]
        interior=np.setdiff1d(np.unique(tri),port['rim_vertex_ids'])
        ports[name]=validate_port(original['points_m'],restored,tri,port,interior,policy,origin=port['plane_origin_m'],normal=port['outward_normal'],basis=port['basis'])
    return match,ports


def audit_volume(data,source,original,contract,policy):
    points,tetra=data['points_m'],data['tetra']
    if not np.isfinite(points).all(): raise ValueError('Nonfinite mesh coordinates')
    volume=tetra_volumes(points,tetra)
    validity={'zero_volume':int(np.count_nonzero(volume==0)),'negative_volume':int(np.count_nonzero(volume<0)),'nonfinite_volume':int(np.count_nonzero(~np.isfinite(volume)))}
    if not np.all(data['cell_tags']==100): raise ValueError('Cell tags changed')
    ext,tags,partition=boundary_partition(points,tetra,data['boundary_triangles'],data['facet_tags'])
    if partition['connected_fluid_components']!=1: raise ValueError('Disconnected volume')
    match,ports=boundary_contract(points,ext,tags,source,original,contract,policy)
    _,centers,_=triangle_geometry(points,ext)
    _,source_centers,source_vectors=triangle_geometry(source['points_m'],source['triangles'])
    origin=source['points_m'].mean(axis=0)
    enclosed=abs(float(np.einsum('ij,ij->',source_centers-origin,source_vectors)/3))
    closure=abs(float(volume.sum())-enclosed)/enclosed
    if closure>1e-10: raise ValueError('Volume closure mismatch')
    q=data['min_sicn']
    if not np.isfinite(q).all() or (q<=0).any(): raise ValueError('Invalid signed tetra quality')
    edges=np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
    lengths=np.linalg.norm(points[tetra[:,edges[:,1]]]-points[tetra[:,edges[:,0]]],axis=2)
    ratio=lengths.max(axis=1)/lengths.min(axis=1)
    cell_centers=points[tetra].mean(axis=1);tree=cKDTree(centers)
    low=np.flatnonzero(q<.1);near=tree.query(cell_centers[low])[1]
    counts={name:int(np.count_nonzero(tags[near]==tag)) for tag,name in FACET_NAMES.items()}
    order=np.argsort(q)[:20];dist,near_w=tree.query(cell_centers[order])
    worst=[{'cell_index':int(i),'gmsh_element_id':int(data['gmsh_element_ids'][i]),'min_sicn':float(q[i]),'volume_m3':float(volume[i]),'edge_ratio':float(ratio[i]),'centroid_m':cell_centers[i].tolist(),'nearest_boundary_patch':FACET_NAMES[int(tags[j])],'distance_to_nearest_boundary_triangle_center_m':float(d)} for i,j,d in zip(order,near_w,dist)]
    quality={'min_sicn':quantiles(q),'total_below_0_1':len(low),'cap_adjacent_below_0_1':len(low)-counts['WALL'],'low_quality_nearest_boundary_counts':counts,'worst_elements':worst,
      'metric':'Gmsh minSICN; same as Stage 1 medium','location_method':policy['tetra_quality']['adjacency_method'],'volume_m3':quantiles(volume),'max_to_min_edge_ratio':quantiles(ratio)}
    proxy=p2_proxy(tetra);gate=volume_gates(validity,quality,proxy,policy)
    return {'status':gate['status'],'volume_gate':gate,'validity':validity,'topology':partition,'boundary_fidelity':match,'ports':ports,
      'volume_closure':{'relative_error':closure,'tetra_volume_m3':float(volume.sum()),'source_enclosed_volume_m3':enclosed},'quality':quality,'proxy':proxy,
      'P2_velocity_ratio_to_baseline':proxy['N_P2_velocity_proxy']/policy['cost']['baseline']['N_P2_velocity_proxy'],'fem_space_created':False,'fem_solved':False}
