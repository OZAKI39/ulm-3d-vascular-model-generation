"""Explicit, dimensional support proxies and triangle-surface manufacturing QC.

These are design diagnostics, not a prediction of clinical flow or a certificate
of support removability. Production geometry is never edited by this module.
"""
from __future__ import annotations
import math
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from .manufacturing_roi import topology_paths


def rigid_transform(matrix):
    t = np.asarray(matrix, float)
    if t.shape != (4, 4) or not np.isfinite(t).all():
        raise ValueError('INVALID_CASTING_TRANSFORM')
    if not np.allclose(t[3], [0, 0, 0, 1],atol=1e-12,rtol=0) or not np.allclose(t[:3, :3].T @ t[:3, :3], np.eye(3), atol=1e-9,rtol=0) or not np.isclose(np.linalg.det(t[:3, :3]), 1,atol=1e-9,rtol=0):
        raise ValueError('SCALING_OR_REFLECTION_PROHIBITED')
    return t


def dimensions(points):
    xyz = np.asarray(points, float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError('INVALID_DIMENSION_POINTS')
    ext = np.ptp(xyz, axis=0)
    long, short = sorted(ext[:2], reverse=True)
    error = max(0., abs(long-60)-3) + max(0., abs(short-50)-3)
    return dict(long_side_mm=float(long), short_side_mm=float(short), Z_span_mm=float(ext[2]),
                xyz_extents_mm=ext.tolist(), size_target_met=bool(57 <= long <= 63 and 47 <= short <= 53),
                size_violation_mm=float(error), bounds_mm=[xyz.min(axis=0).tolist(), xyz.max(axis=0).tolist()])


def topology_metrics(graph):
    roots = sorted(n for n in graph if graph.in_degree(n) == 0)
    leaves = sorted(n for n in graph if graph.out_degree(n) == 0)
    paths = topology_paths(graph) if nx.is_directed_acyclic_graph(graph) else []
    lengths = [float(sum(np.linalg.norm(graph.nodes[b]['coords'][:3]-graph.nodes[a]['coords'][:3]) for a,b in zip(p,p[1:]))) for p in paths]
    twigs = [v for p,v in zip(paths,lengths) if graph.out_degree(p[-1]) == 0]
    connected = nx.number_weakly_connected_components(graph) if len(graph) else 0
    valid = connected == 1 and nx.is_arborescence(graph) and len(roots) == 1 and graph.out_degree(roots[0]) == 1 and len(leaves) in (2,3)
    return dict(root_node=roots[0] if len(roots)==1 else None, inlet_count=len(roots), outlet_nodes=leaves,
        outlet_count=len(leaves), connected_components=connected, branch_count=len(paths),
        bifurcation_count=sum(graph.out_degree(n)>1 for n in graph), centerline_length_mm=sum(lengths),
        minimum_terminal_length_mm=min(twigs,default=0.), short_terminal_count=sum(v<3 for v in twigs),
        topology_valid=bool(valid and twigs and min(twigs)>=3))


def centerline_support(graph, transform):
    t=rigid_transform(transform); paths=topology_paths(graph)
    root=next(n for n in graph if graph.in_degree(n)==0)
    depth={root:0}
    for n in nx.topological_sort(graph):
        for child in graph.successors(n): depth[child]=depth[n]+int(graph.out_degree(n)>1)
    bif=[n for n in graph if graph.out_degree(n)>1]; rows=[]
    for branch,path in enumerate(paths):
        for a,b in zip(path,path[1:]):
            pa,pb=graph.nodes[a]['coords'],graph.nodes[b]['coords']; vector=t[:3,:3]@(pb[:3]-pa[:3]); length=float(np.linalg.norm(vector))
            if length<=0: raise ValueError('ZERO_LENGTH_SEGMENT')
            # Horizontal axes maximize the downward tube surface; vertical axes
            # minimize it. This is explicitly a surrogate, not triangle area.
            severity=math.sqrt(max(0.,1-(vector[2]/length)**2))
            radius=float((pa[3]+pb[3])/2); contact=math.pi*radius*length*severity
            near=any(min(np.linalg.norm(pa[:3]-graph.nodes[n]['coords'][:3]),np.linalg.norm(pb[:3]-graph.nodes[n]['coords'][:3]))<=3 for n in bif)
            rows.append(dict(parent=a,child=b,branch=branch,branch_depth=depth[b],radius_mm=radius,
                length_mm=length,overhang_severity=severity,estimated_support_contact_mm2=contact,
                thin_branch_support_penalty_mm2=contact*max(0.,.8/radius-1),near_bifurcation=near))
    return dict(estimated_support_contact_mm2=sum(r['estimated_support_contact_mm2'] for r in rows),
        thin_branch_support_penalty=sum(r['thin_branch_support_penalty_mm2'] for r in rows),
        support_critical_branches=len({r['branch'] for r in rows if r['overhang_severity']>.5 and r['radius_mm']<.8}),
        support_near_bifurcation_count=len({r['branch'] for r in rows if r['near_bifurcation'] and r['overhang_severity']>.5}),
        segments=rows, method='TANGENT_RADIUS_LENGTH_DEPTH_SURROGATE')


def sample_companion(graph,transform,step=.2):
    t=rigid_transform(transform); xyz=[]; labels=[]; radius=[]; origins=[]; edge_labels={}
    for index,path in enumerate(topology_paths(graph)):
        label=f'B{index+1:02}'
        for a,b in zip(path,path[1:]):
            pa,pb=graph.nodes[a]['coords'],graph.nodes[b]['coords'];length=float(np.linalg.norm(pb[:3]-pa[:3]));edge_labels[a,b]=label
            for f in np.linspace(0,1,max(2,math.ceil(length/step)+1)):
                p=pa*(1-f)+pb*f;xyz.append(t[:3,:3]@p[:3]+t[:3,3]);radius.append(p[3]);labels.append(label);origins.append((a,b,f,length))
    return np.array(xyz),np.array(labels),np.array(radius),origins,edge_labels


def surface_support(mesh,graph,transform,rotation=None,threshold_deg=30):
    xyz,labels,radii,_,_=sample_companion(graph,transform)
    near=cKDTree(xyz).query(mesh.triangles_center)[1]
    rot=np.eye(3) if rotation is None else np.asarray(rotation)
    risk=(mesh.face_normals@rot.T)[:,2] < -math.cos(math.radians(threshold_deg))
    thin=radii[near]<.8
    t=np.asarray(transform);bif=[t[:3,:3]@graph.nodes[n]['coords'][:3]+t[:3,3] for n in graph if graph.out_degree(n)>1]
    bifnear=cKDTree(bif).query(mesh.triangles_center)[0]<=3 if bif else np.zeros(len(mesh.faces),bool)
    colors=np.zeros(len(mesh.faces),int);colors[risk]=1;colors[risk&(thin|bifnear)]=2
    area=mesh.area_faces
    return dict(overhang_area_mm2=float(area[risk].sum()),support_contact_candidate_area_mm2=float(area[risk].sum()),
        thin_vessel_support_area_mm2=float(area[risk&thin].sum()),bifurcation_support_area_mm2=float(area[risk&bifnear].sum()),
        support_near_bifurcation_count=len(set(labels[near[risk&bifnear]])),support_critical_branches=len(set(labels[near[risk&thin]])),
        triangle_area_used=True,threshold_from_horizontal_deg=threshold_deg,
        method='DOWNWARD_TRIANGLE_AREA_CANDIDATE_NOT_ACTUAL_SLICER_CONTACT'),colors


def vascular_ligament(mesh,graph,transform,cfg):
    from .casting_mold_qc import ligament_qc
    xyz,labels,_,origins,edge_labels=sample_companion(graph,transform)
    tree=cKDTree(xyz);face_labels=labels[tree.query(mesh.triangles_center)[1]]
    groups={label:np.flatnonzero(face_labels==label) for label in sorted(set(face_labels))}
    partition=dict(groups=groups,sample_tree=tree,sample_origins=origins)
    inputs=dict(companion=dict(graph=graph,branches=edge_labels),ports=[])
    result=ligament_qc(mesh,inputs,partition,cfg)
    result['passed']=result['minimum_ligament_mm']>=1.6
    result['status']='PASS' if result['passed'] else 'PDMS_LIGAMENT_TOO_THIN'
    return result


def reduction(old,new):
    if old is None or new is None or old<=0: return None
    return 100*(1-new/old)


def reduction_gate(old_mass,new_mass,old_contact,new_contact):
    values=[reduction(old_mass,new_mass),reduction(old_contact,new_contact)]
    return any(v is not None and v>=25 for v in values)


def interface_gate(continuity,section,loss,tolerance=.001):
    feature=next(r for r in continuity['features'] if r['angle_deg']==20)
    return bool(feature['total_length_mm']==0 and feature['circumferential_ring_count']==0
        and section['passed'] and loss<=tolerance)
