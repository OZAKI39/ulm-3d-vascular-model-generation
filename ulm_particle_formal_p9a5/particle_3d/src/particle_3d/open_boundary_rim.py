"""Exact GlobalNodeID topology gate; no geometric distance band."""
from collections import Counter, defaultdict
import numpy as np
from .particle_shapes import Capsule, EPS
from .convex_triangle import capsule_triangle_many


def perimeter_ids(triangles):
    counts=Counter(tuple(sorted((int(a),int(b)))) for tri in triangles
                   for a,b in zip(tri,np.roll(tri,-1)))
    return {edge for edge,count in counts.items() if count==1}


def build_rim_topology(wall_ids,cap_ids,source_hashes=None):
    wall_edges=perimeter_ids(wall_ids);edges=defaultdict(list);vertices=defaultdict(set)
    counts={}
    for role,ids in sorted(cap_ids.items()):
        cap_edges=perimeter_ids(ids);shared=wall_edges & cap_edges
        if shared!=cap_edges:
            raise ValueError('Open cap perimeter does not match WALL GlobalNodeID topology: '+role)
        counts[role]=len(shared)
        for edge in sorted(shared):
            edges[edge].append(role)
            for vertex in edge:vertices[vertex].add(role)
    return dict(edges=dict(edges),vertices={k:sorted(v) for k,v in vertices.items()},
        provenance=dict(method='GLOBAL_NODE_ID_PERIMETER_EDGE_INTERSECTION',
            global_node_id_base=0,source_sha256=source_hashes or {},rim_edge_counts=counts,
            distance_band_m=None))


def rim_witness(shape,wall,gap):
    topology=getattr(wall,'open_boundary_topology',None)
    if topology is None or gap.wall_feature not in ('EDGE','VERTEX'):return None
    i=int(gap.wall_triangle_id)
    # Same kernel, query center and barycentric feature tolerance as wall_gap.
    bary=capsule_triangle_many(Capsule(shape.center_m,[0,0,1],shape.radius_m,0.),wall.triangles[i:i+1])[4][0]
    active=np.flatnonzero(bary>128*EPS);gids=wall.global_node_ids[i]
    if gap.wall_feature=='EDGE' and len(active)==2:
        nodes=tuple(sorted(int(x) for x in gids[active]));roles=topology['edges'].get(nodes,[])
    elif gap.wall_feature=='VERTEX' and len(active)<=1:
        nodes=(int(gids[np.argmax(bary)]),);roles=topology['vertices'].get(nodes[0],[])
    else:return None
    return dict(feature=gap.wall_feature,global_node_ids_zero_based=list(nodes),open_boundary_roles=roles) if roles else None
