"""Bounded geometric-band and pinned graph-cut diagnostics; never mutate bonds."""
from common import *
import csv,time
from scipy.sparse import coo_matrix,csgraph

graph=module('bridge_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')

def candidate_regions(ref):
    x=ref['X'][:int(ref['clot_count'])];u=x*1e6
    # The original notch is z=0, tip x=-50 um. The second level is midway
    # from that observed damage band to the historically loaded top surface.
    specs=[
      ('upper_notch_flap',(u[:,2]>0),(u[:,2]>=175-1e-8),np.zeros(len(x),bool)),
      ('loaded_surface_cap',(u[:,0]<-50)&(u[:,2]>200),(u[:,0]<=-125+1e-8)&(u[:,2]>=275-1e-8),(u[:,0]>=125-1e-8)|(u[:,2]<0)),
      ('left_notch_lip',(u[:,0]<-50)&(u[:,2]>0)&(u[:,2]<200),(u[:,0]<=-175+1e-8)&(u[:,2]>0)&(u[:,2]<=125+1e-8),(u[:,0]>=125-1e-8)|(u[:,2]>=275-1e-8)),
    ]
    return [dict(region=name,original_atom_ids=np.flatnonzero(region)+1,source_core_ids=np.flatnonzero(core)+1,sink_material_ids=np.flatnonzero(sink)+1,origin='Historical diagnostic volumes retained only for bridge comparisons; not loading masks or prescribed fragment boundaries') for name,region,core,sink in specs]

def pinned_cut(ref,state,region):
    n=int(ref['clot_count']);N=len(ref['X']);b=ref['bonds'];live=(state['d']<1);ids=np.flatnonzero(live);edge=b[ids]
    capacity=np.maximum(1,np.rint(1024*(1-state['d'][ids]))).astype(np.int64)
    pin=int(capacity.sum()+1);assert pin<2**31
    core=np.asarray(region['source_core_ids'])-1;sink=np.r_[np.asarray(region['sink_material_ids'],int)-1,np.arange(n,N)]
    assert len(core)>=34 and not np.intersect1d(core,sink).size
    rr=np.r_[edge[:,0],edge[:,1],np.full(len(core),N),sink]
    cc=np.r_[edge[:,1],edge[:,0],core,np.full(len(sink),N+1)]
    weights=np.r_[capacity,capacity,np.full(len(core)+len(sink),pin)].astype(np.int64)
    matrix=coo_matrix((weights,(rr,cc)),shape=(N+2,N+2)).tocsr()
    flow=csgraph.maximum_flow(matrix,N,N+1,method='dinic')
    residual=(matrix-flow.flow).tocsr();residual.data=(residual.data>0).astype(np.int64);residual.eliminate_zeros()
    reach=csgraph.breadth_first_order(residual,N,directed=True,return_predecessors=False)
    mask=np.zeros(N+2,bool);mask[reach]=True
    assert mask[core].all() and not mask[sink].any()
    crossing=mask[edge[:,0]]^mask[edge[:,1]];cut_ids=ids[crossing]
    assert int(capacity[crossing].sum())==flow.flow_value
    return dict(source_partition_atom_ids=np.flatnonzero(mask[:N])+1,original_cut_bond_ids=cut_ids+1,cut_bond_count=len(cut_ids),integer_capacity=int(flow.flow_value),integrity_capacity_proxy=float(flow.flow_value/1024),capacity_definition='max(1,round(1024*(1-d))) on live original bonds; pinned source core >=34 nodes and fixed anchors/specified far body as sink. Diagnostic proxy, not fracture energy or authorized deletions')

def live_anchor_path(ref,state,region):
    n=int(ref['clot_count']);N=len(ref['X']);b=ref['bonds'];live=state['d']<1
    # Exhibit a path entirely outside a central two-node-thick Y viewing slab.
    outside=np.abs(ref['X'][:,1])>50e-6
    keep=live&outside[b[:,0]]&outside[b[:,1]];ids=np.flatnonzero(keep);edge=b[ids]
    matrix=coo_matrix((np.ones(2*len(ids)),(np.r_[edge[:,0],edge[:,1]],np.r_[edge[:,1],edge[:,0]])),shape=(N,N)).tocsr()
    starts=np.asarray(region['source_core_ids'])-1;starts=starts[outside[starts]]
    if not len(starts):return dict(found=False)
    start=int(starts[np.argmax(ref['X'][starts,1])]);order,pred=csgraph.breadth_first_order(matrix,start,directed=False)
    anchors=order[order>=n]
    if not len(anchors):return dict(found=False,start_atom_id=start+1)
    path=[int(anchors[0])]
    while path[-1]!=start:path.append(int(pred[path[-1]]));assert len(path)<=N
    path=path[::-1];lookup={tuple(sorted(p)):int(q+1) for p,q in zip(edge,ids)};bids=[lookup[tuple(sorted(pair))] for pair in zip(path[:-1],path[1:])]
    return dict(found=True,original_atom_ids=np.array(path)+1,original_bond_ids=bids,reference_coordinates_m=ref['X'][path],actual_coordinates_m=state['x'][path],definition='An actually surviving 3D material path to a fixed anchor outside |Y|<=50 um; diagnostic restriction only, no physical deletion')

def inspect(ref,state,region,onset_by_bond=None):
    nodes=np.asarray(region['original_atom_ids'])-1;b=ref['bonds'];bt=ref['bond_types'];cut=ref['initial_cut'];boundary,_=graph.boundary(ref,nodes);boundary=boundary[~cut[boundary]];alive=boundary[state['d'][boundary]<1]
    dy=state['x'][b[alive,1]]-state['x'][b[alive,0]];stretch=np.linalg.norm(dy,axis=1)/ref['r0'][alive]-1
    mid=ref['X'][b[alive]].mean(1);onset=(S0[bt] if onset_by_bond is None else onset_by_bond)[alive];labels,_,attached=graph.topology(ref,state['d'])
    info=dict(region=region['region'],candidate_nodes=len(nodes),candidate_is_free=bool(not attached[labels[nodes]].any()),anchor_reachable_nodes=int(np.sum(attached[labels[nodes]])),original_crossing_bonds=len(boundary),surviving_crossing_bonds=len(alive),new_failed_crossing_bonds=int(np.sum(state['d'][boundary]==1)),surviving_bulk=int(np.sum(bt[alive]==1)),surviving_interface=int(np.sum(bt[alive]==2)),long_bonds_ge_2dx=int(np.sum(ref['r0'][alive]>=2*DX-1e-15)),currently_compressive=int(np.sum(stretch<0)),historically_below_onset=int(np.sum(state['kappa'][alive]<=onset)),partially_softened=int(np.sum(state['d'][alive]>0)),surviving_midpoint_mean_m=mid.mean(0) if len(mid) else None,surviving_midpoint_min_m=mid.min(0) if len(mid) else None,surviving_midpoint_max_m=mid.max(0) if len(mid) else None,right_of_notch_tip=int(np.sum(mid[:,0]>-50e-6)),original_surviving_bond_ids=alive+1,cut_proxy=pinned_cut(ref,state,region),off_midplane_anchor_path=live_anchor_path(ref,state,region))
    records=[]
    for j,bid in enumerate(alive):
        direction=(ref['X'][b[bid,1]]-ref['X'][b[bid,0]])/ref['r0'][bid]
        records.append(dict(region=region['region'],original_bond_id=int(bid+1),bond_type=int(bt[bid]),atom_i=int(b[bid,0]+1),atom_j=int(b[bid,1]+1),r0_m=ref['r0'][bid],direction_x=direction[0],direction_y=direction[1],direction_z=direction[2],midpoint_x_m=mid[j,0],midpoint_y_m=mid[j,1],midpoint_z_m=mid[j,2],current_stretch=stretch[j],kappa=state['kappa'][bid],damage=state['d'][bid],historically_below_onset=bool(state['kappa'][bid]<=onset[j]),candidate_anchor_reachable=bool(info['anchor_reachable_nodes'])))
    return info,records

