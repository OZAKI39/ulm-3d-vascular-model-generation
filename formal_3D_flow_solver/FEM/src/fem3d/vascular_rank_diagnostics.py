"""Topology-only pressure support audit for wall-constrained quadratic velocity."""
import numpy as np


def constrained_pressure_support(data):
    tetra=data['tetra'];points=data['points_m'];wall=data['boundary_triangles'][data['facet_tags']==1]
    nv=len(points);wall_vertices=np.unique(wall)
    pairs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    wall_edges=np.unique(np.sort(np.concatenate([wall[:,[0,1]],wall[:,[0,2]],wall[:,[1,2]]]),axis=1),axis=0)
    wall_keys=wall_edges[:,0]*np.int64(nv)+wall_edges[:,1]
    locked=np.isin(tetra,wall_vertices).all(axis=1)
    for a,b in pairs:
        edge=np.sort(tetra[:,[a,b]],axis=1)
        locked &= np.isin(edge[:,0]*np.int64(nv)+edge[:,1],wall_keys)
    total_support=np.bincount(tetra.ravel(),minlength=nv)
    free_support=np.bincount(tetra[~locked].ravel(),minlength=nv)
    vertices=np.flatnonzero((total_support>0)&(free_support==0))
    return {'wall_velocity_fully_constrained_cells':[
        {'source_cell_index':int(i),'centroid_m':points[tetra[i]].mean(axis=0).tolist(),
         'min_sicn':float(data['min_sicn'][i]),'vertex_indices':tetra[i].tolist(),
         'all_4_vertex_and_6_edge_velocity_nodes_on_wall':True} for i in np.flatnonzero(locked)],
        'uncoupled_pressure_vertices':[
        {'source_vertex_index':int(v),'coordinates_m':points[v].tolist(),
         'support_source_cells':np.flatnonzero(np.any(tetra==v,axis=1)).tolist(),
         'free_velocity_support_cells':int(free_support[v])} for v in vertices],
        'meaning':'Each listed P1 pressure basis has support only in tetrahedra whose entire P2 velocity polynomial is fixed to zero by wall no-slip. Its divergence-coupling row and column therefore vanish after wall elimination.',
        'geometry_modified':False,'pressure_constraints_added':False}
