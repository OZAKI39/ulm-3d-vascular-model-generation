"""Deterministic surface-based FDM orientation heuristics; rigid transforms only."""
from __future__ import annotations

import itertools
from pathlib import Path

import networkx as nx
import numpy as np
import pyvista as pv
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from .manufacturing_roi import stable_pca, topology_paths
from .mevo_refinement import table


def cap_native_ports(surface):
    """Cap only simple boundary loops on a PRINT derivative, preserving sidewall vertices."""
    mesh = surface.extract_surface(algorithm='dataset_surface').triangulate().clean(tolerance=0.)
    faces = mesh.faces.reshape(-1,4)[:,1:]
    all_edges = np.sort(np.vstack([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
    edges,counts = np.unique(all_edges,axis=0,return_counts=True)
    if np.any(counts>2):
        raise ValueError('NATIVE_SURFACE_NONMANIFOLD')
    boundary = nx.Graph();boundary.add_edges_from(edges[counts==1].tolist())
    points = mesh.points.tolist();triangles=faces.tolist();loops=[]
    for group in sorted(nx.connected_components(boundary),key=min):
        if any(boundary.degree(n)!=2 for n in group):
            raise ValueError('NATIVE_BOUNDARY_IS_NOT_A_SIMPLE_PORT_LOOP')
        first=min(group);ordered=[first];previous=None;cursor=first
        while True:
            nxt=next(n for n in sorted(boundary.neighbors(cursor)) if n!=previous)
            if nxt==first:
                break
            ordered.append(nxt);previous,cursor=cursor,nxt
        assert len(ordered)==len(group)
        center=np.mean(mesh.points[ordered],axis=0);index=len(points);points.append(center.tolist())
        for a,b in zip(ordered,ordered[1:]+ordered[:1]):
            triangles.append([a,b,index])
        loops.append(dict(vertices=len(ordered),center_mm=center.tolist(),
                          median_radius_mm=float(np.median(np.linalg.norm(mesh.points[ordered]-center,axis=1)))))
    closed = pv.PolyData(np.asarray(points),np.c_[np.full(len(triangles),3),triangles].ravel())
    closed = closed.compute_normals(cell_normals=True,point_normals=False,consistent_normals=True,
                                    auto_orient_normals=True,split_vertices=False)
    nonmanifold=closed.extract_feature_edges(boundary_edges=False,non_manifold_edges=True,
                                             feature_edges=False,manifold_edges=False).n_cells
    if closed.n_open_edges or nonmanifold:
        raise ValueError('PRINT_PORT_CAP_QC_FAILED')
    regions=closed.connectivity(extraction_mode='all',label_regions=True)
    component_count=len(np.unique(regions.cell_data['RegionId']))
    if component_count!=1:
        raise ValueError('PRINT_SURFACE_DISCONNECTED')
    assert np.array_equal(closed.points[:mesh.n_points],mesh.points)
    return closed,dict(port_caps=loops,cap_count=len(loops),open_edges=closed.n_open_edges,
        nonmanifold_edges=nonmanifold,connected_component_count=component_count,native_sidewall_vertices_exact=True,
        solid_sacrificial_core=True,self_intersection_certified=False)


def triangle_geometry(mesh):
    triangles=mesh.points[mesh.faces.reshape(-1,4)[:,1:]]
    cross=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    norm=np.linalg.norm(cross,axis=1)
    normals=np.divide(cross,norm[:,None],out=np.zeros_like(cross),where=norm[:,None]>0)
    return norm/2,normals,triangles.mean(axis=1)


def branch_descriptors(graph):
    descriptors=[]
    for path in topology_paths(graph):
        xyzr=np.array([graph.nodes[n]['coords'] for n in path])
        vectors=np.diff(xyzr[:,:3],axis=0);lengths=np.linalg.norm(vectors,axis=1)
        descriptors.append(dict(length=float(lengths.sum()),diameter=float(np.median(2*xyzr[:,3])),
                                vectors=vectors,lengths=lengths))
    return descriptors


def transform_points(points,matrix):
    return np.asarray(points)@np.asarray(matrix)[:3,:3].T+np.asarray(matrix)[:3,3]


def bed_transform(points,rotation,volume,margin):
    rotated=points@rotation.T;low=rotated.min(axis=0);high=rotated.max(axis=0)
    center=(low+high)/2
    translation=np.array([volume[0]/2-center[0],volume[1]/2-center[1],-low[2]])
    matrix=np.eye(4);matrix[:3,:3]=rotation;matrix[:3,3]=translation
    return matrix,high-low,low+translation,high+translation


def fits_volume(extents,usable):
    return bool(np.all(np.asarray(extents)<=np.asarray(usable)+1e-8))


def orientation_rotations(points,config):
    settings=config['orientation'];pca=stable_pca(points)
    step=settings['coarse_angle_step_deg'];n=settings['tilt_steps_each_side']
    rotations=[('identity_baseline',np.eye(3),None),('pca_seed',pca,[0,0,0])]
    for flipped in [False,True]:
        base=np.diag([1,-1,-1])@pca if flipped else pca
        for x,y,z in itertools.product(range(-n,n+1),range(-n,n+1),settings['yaw_offsets_deg']):
            angles=[x*step,y*step,z]
            rotation=Rotation.from_euler('xyz',angles,degrees=True).as_matrix()@base
            rotations.append(('pca_flipped_coarse' if flipped else 'pca_coarse',rotation,angles))
    unique=[];seen=set()
    for item in rotations:
        key=tuple(np.round(item[1],10).ravel())
        if key not in seen:
            unique.append(item);seen.add(key)
    if not 1<=len(unique)<=300:
        raise ValueError('Configured orientation search exceeds 300 candidates')
    return unique


def score_orientation(points,areas,normals,centers,rotation,config,branches):
    printer=config['printer'];settings=config['orientation'];usable=np.array(config['effective']['usable_volume_mm'])
    matrix,extent,low,high=bed_transform(points,rotation,printer['build_volume_mm'],printer['bed_margin_xy_mm'])
    fits=fits_volume(extent,usable)
    rotated=transform_points(points,matrix)
    downward=normals@rotation[2] < -np.cos(np.radians(settings['overhang_angle_deg']))
    height_centers=centers@rotation[2]+matrix[2,3]
    downward &= height_centers>settings['overhang_bed_exemption_mm']
    support_area=float(areas[downward].sum());area=float(areas.sum())
    footprint=float(ConvexHull(rotated[:,:2]).volume)
    thin_length=0.;horizontal=0.
    for branch in branches:
        if branch['length']<settings['long_branch_length_mm'] or branch['diameter']>=settings['thin_branch_diameter_mm']:
            continue
        length=branch['lengths'];vertical=np.abs(branch['vectors']@rotation[2])
        angles=np.degrees(np.arcsin(np.clip(np.divide(vertical,length,out=np.zeros_like(length),where=length>0),0,1)))
        thin_length+=branch['length'];horizontal+=float(length[angles<settings['horizontal_branch_angle_deg']].sum())
    cantilever=horizontal/thin_length if thin_length else 0.
    ratio=float(extent[2]/max(np.sqrt(footprint),1e-12))
    mix=settings['cantilever_mix_in_support_score']
    components=dict(support_area=(1-mix)*support_area/area+mix*cantilever,
        z_height=float(extent[2]/usable[2]),bed_fit_margin=float(np.max(extent/usable)),
        slenderness=ratio/(1+ratio),footprint=1-min(footprint/(usable[0]*usable[1]),1.))
    score=sum(settings['score_weights'][k]*v for k,v in components.items())
    euler=Rotation.from_matrix(rotation).as_euler('xyz',degrees=True).tolist()
    return dict(fits_build_volume=fits,status='FEASIBLE' if fits else 'ROI_EXCEEDS_BUILD_VOLUME',
        score=float(score),score_components=components,transform_4x4=matrix.tolist(),euler_xyz_deg=euler,
        bbox_min_mm=low.tolist(),bbox_max_mm=high.tolist(),bbox_extents_mm=extent.tolist(),
        height_mm=float(extent[2]),footprint_mm2=footprint,height_to_footprint_ratio=ratio,
        downward_support_area_mm2=support_area,total_surface_area_mm2=area,
        support_area_fraction=support_area/area,long_thin_horizontal_fraction=cantilever,
        support_estimate_method='Downward triangle area + long thin horizontal centerlines; heuristic, not slicer support volume',
        footprint_method='Projected surface convex hull, not actual contact area',scale=1.0)


def optimize(mesh,graph,config,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    points=np.asarray(mesh.points,dtype=float);areas,normals,centers=triangle_geometry(mesh)
    # All extrema and projected convex hulls are exactly preserved by the 3-D
    # hull vertices; only area/normals calculations still use every triangle.
    hull_points=points[ConvexHull(points).vertices]
    branches=branch_descriptors(graph)
    results=[]
    with threadpool_limits(limits=1):
        for i,(method,rotation,offsets) in enumerate(orientation_rotations(points,config)):
            result=score_orientation(hull_points,areas,normals,centers,rotation,config,branches)
            result.update(candidate_id=i,method=method,pca_offset_euler_xyz_deg=offsets)
            results.append(result)
    baseline=results[0]
    baseline_mesh=mesh.copy();baseline_mesh.points=transform_points(points,baseline['transform_4x4'])
    baseline_mesh.save(output/'orientation_baseline.stl')
    feasible=sorted((r for r in results if r['fits_build_volume']),key=lambda r:(r['score'],r['candidate_id']))
    top=feasible[:config['orientation']['candidate_count_for_slicer']]
    for i,item in enumerate(top,1):
        path=output/f'orientation_candidate_{i:02d}.stl';derivative=mesh.copy()
        derivative.points=transform_points(points,item['transform_4x4']);derivative.save(path)
        item['rank']=i;item['stl']=str(path)
        reverse=np.linalg.inv(np.asarray(item['transform_4x4']))
        np.testing.assert_allclose(transform_points(derivative.points,reverse),points,atol=1e-10)
    rows=[dict(candidate_id=r['candidate_id'],method=r['method'],feasible=r['fits_build_volume'],score=r['score'],
        bbox_x_mm=r['bbox_extents_mm'][0],bbox_y_mm=r['bbox_extents_mm'][1],bbox_z_mm=r['height_mm'],
        support_area_mm2=r['downward_support_area_mm2'],footprint_mm2=r['footprint_mm2'],
        long_thin_horizontal_fraction=r['long_thin_horizontal_fraction'],
        euler_x_deg=r['euler_xyz_deg'][0],euler_y_deg=r['euler_xyz_deg'][1],euler_z_deg=r['euler_xyz_deg'][2]) for r in results]
    table(output/'orientation_scores.csv',rows)
    return dict(status='PRINT_ORIENTATION_OPTIMIZED' if top else 'ROI_EXCEEDS_BUILD_VOLUME',
        tested=len(results),feasible_count=len(feasible),baseline=baseline,top_candidates=top,all_candidates=results)


def calibration_coupon(config,path):
    settings=config['coupon'];meshes=[];rows=[]
    for i,d in enumerate(settings['diameters_mm']):
        for j,angle in enumerate(settings['angles_from_plate_deg']):
            theta=np.radians(angle)
            cylinder=pv.Cylinder(center=(0,0,0),direction=(np.cos(theta),0,np.sin(theta)),
                radius=d/2,height=settings['cylinder_length_mm'],resolution=settings['radial_segments'],capping=True).triangulate()
            cylinder.points[:,0]+=i*settings['column_spacing_mm']
            cylinder.points[:,1]+=j*settings['row_spacing_mm']
            cylinder.points[:,2]-=cylinder.bounds.z_min
            meshes.append(cylinder)
            rows.append(dict(column=i+1,row=j+1,diameter_mm=d,angle_from_plate_deg=angle,
                             length_mm=settings['cylinder_length_mm']))
    coupon=meshes[0].merge(meshes[1:],merge_points=False).clean(tolerance=0.)
    coupon.save(path)
    table(Path(path).with_suffix('.csv'),rows)
    return dict(path=str(path),cylinders=len(rows),diameters_mm=settings['diameters_mm'],
        angles_from_plate_deg=settings['angles_from_plate_deg'],
        purpose='Separate short cylinders for adhesion/stability/acetone-vapor measurement; no anatomy changes',
        bbox_extents_mm=np.ptp(coupon.points,axis=0).tolist())
