"""Exact solid/area gates and deterministic, lexicographic orientation ranking."""
import itertools
import math
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation
from shapely.geometry import box, Polygon, Point
from shapely.ops import unary_union
import trimesh
from . import sacrificial_fixture as legacy
from .sacrificial_print_frame import box_mesh,port_stem_bounds


def volume(mesh):
    return abs(float(mesh.volume)) if len(mesh.faces) else 0.


def intersection(a,b):
    return trimesh.boolean.intersection([a,b],engine='manifold')


def open_face_metrics(design,side,cfg):
    axis=design['open_axis'];other=1-axis;inner=design['inner_bounds'];outer=design['outer_bounds']
    opening=box(inner[0,other],inner[0,2],inner[1,other],inner[1,2]);blocks=[]
    # Intersect the actual boundary slab, not the projection of a remote wall.
    slab=sorted([outer[side,axis],inner[side,axis]])
    for member in design['members']:
        bounds=member['bounds']
        if min(bounds[1,axis],slab[1])-max(bounds[0,axis],slab[0])<=1e-8:continue
        projected=box(bounds[0,other],bounds[0,2],bounds[1,other],bounds[1,2]).intersection(opening)
        if projected.area>1e-12:blocks.append(projected)
    obstructed=unary_union(blocks) if blocks else Polygon();free=opening.difference(obstructed)
    xs={inner[0,other],inner[1,other]};ys={inner[0,2],inner[1,2]}
    for polygon in blocks:
        xs.update([polygon.bounds[0],polygon.bounds[2]]);ys.update([polygon.bounds[1],polygon.bounds[3]])
    best=None
    for x0,x1 in itertools.combinations(sorted(xs),2):
        for y0,y1 in itertools.combinations(sorted(ys),2):
            rectangle=box(x0,y0,x1,y1)
            if rectangle.intersection(obstructed).area>1e-8:continue
            rank=(min(x1-x0,y1-y0),rectangle.area)
            if best is None or rank>best[0]:best=(rank,(x0,y0,x1,y1))
    rectangle=best[1] if best else (0,0,0,0)
    fraction=float(free.area/opening.area)
    width,height=rectangle[2]-rectangle[0],rectangle[3]-rectangle[1]
    return dict(face=design['open_faces'][side],free_area_fraction=fraction,
        free_area_mm2=float(free.area),interior_face_area_mm2=float(opening.area),
        minimum_clear_opening_width_mm=float(width),minimum_clear_opening_height_mm=float(height),
        guaranteed_clear_rectangle=list(rectangle),
        opening_dimension_method='Largest inscribed clear rectangle, ranked by its smaller dimension then area',
        passed=bool(fraction>=cfg['frame']['minimum_open_face_fraction'] and min(width,height)>=cfg['frame']['minimum_opening_dimension_mm']))


def frame_checks(design,inputs,keep_bounds,cfg):
    frame=design['mesh'];core=inputs['core'];failures=[]
    keep=intersection(frame,box_mesh(keep_bounds));intrusion=volume(keep)
    if intrusion>cfg['pdms_keep_zone']['maximum_frame_intrusion_mm3']:failures.append('FRAME_INTRUDES_PDMS_KEEP_ZONE')
    contacts=[];masks=[]
    whole_contact=intersection(frame,core)
    for port in inputs['ports']:
        mask=box_mesh(port_stem_bounds(port,cfg['geometry']['coordinate_tolerance_mm']));masks.append(mask)
        contact=intersection(whole_contact,mask) if len(whole_contact.faces) else whole_contact.copy()
        length=float(np.ptp(contact.vertices@port['direction'])) if len(contact.vertices) else 0.
        v=volume(contact)
        passed=v>=cfg['frame']['minimum_connection_volume_mm3'] and length>=cfg['frame']['minimum_connection_length_mm']
        if not passed:failures.append('PORT_FRAME_CONNECTION_FAILED_'+port['port_id'])
        contacts.append(dict(port_id=port['port_id'],port_radius_mm=port['radius_mm'],
            frame_contact_region='EXISTING_EXTERNAL_STEM_'+port['face'],intersection_volume_mm3=v,
            intersection_length_mm=length,status='PASS' if passed else 'FAIL'))
    allowed=trimesh.boolean.union(masks,engine='manifold')
    unwanted=trimesh.boolean.difference([whole_contact,allowed],engine='manifold') if len(whole_contact.faces) else whole_contact.copy()
    unintended=volume(unwanted)
    if unintended>cfg['geometry']['volume_tolerance_mm3']:failures.append('UNINTENDED_MID_VESSEL_FRAME_CONTACT')
    native_gap=legacy.surface_distance(frame,inputs['central'])
    if native_gap<cfg['frame']['vessel_clearance_mm']:failures.append('CENTRAL_VASCULAR_CLEARANCE_TOO_SMALL')
    open_faces=[open_face_metrics(design,k,cfg) for k in (0,1)]
    if not all(r['passed'] for r in open_faces):failures.append('OPEN_FACE_TOO_BLOCKED')
    mesh=legacy.mesh_qc(frame,cfg)
    if not mesh['passed']:failures.append('FRAME_NOT_CLOSED_SINGLE_COMPONENT')
    return dict(frame_mesh=mesh,frame_volume_mm3=volume(frame),frame_intrusion_keep_zone_mm3=intrusion,
        ghost_void_volume_inside_keep_zone_mm3=intrusion,
        additional_ghost_void_volume_mm3=volume(frame)-volume(whole_contact),
        vascular_contact_volume_mm3=volume(whole_contact),port_connections=contacts,
        unintended_mid_vessel_contact_mm3=unintended,central_vessel_clearance_mm=native_gap,
        open_faces=open_faces,failures=failures,passed=not failures),keep,whole_contact


def layout_rank(design,qc,access):
    # Explicit order: access, conflicts with preferred openings, intrusion,
    # material, overhang proxy. No weighted score.
    conflicts=sum(p['face'] in design['open_faces'] for p in design['ports_for_ranking'])
    normals=design['mesh'].face_normals
    down=float(design['mesh'].area_faces[normals[:,2]<-.70710678].sum())
    return (-access['accessible_fraction'],-access['probe_accessible_fraction'],conflicts,
        qc['frame_intrusion_keep_zone_mm3'],qc['frame_volume_mm3'],down,design['layout'])


def orientation_rotations(cfg):
    yaws=cfg['orientation']['yaw_angles_deg'];tilts=cfg['orientation']['tilt_degrees']
    seeds=[('previous_print_frame',np.eye(3))]
    for a,ang in [('x',90),('x',-90),('y',90),('y',-90),('x',180)]:
        seeds.append((a+str(ang),Rotation.from_euler(a,ang,degrees=True).as_matrix()))
    choices=[]
    for name,seed in seeds:
        for yaw in yaws:
            choices.append((name+'_yaw'+str(yaw),Rotation.from_euler('z',yaw,degrees=True).as_matrix()@seed))
    for axis in ('x','y'):
        for tilt in tilts:
            for yaw in yaws:
                choices.append(('small_tilt_'+axis+str(tilt)+'_'+str(yaw),Rotation.from_euler('z',yaw,degrees=True).as_matrix()@Rotation.from_euler(axis,tilt,degrees=True).as_matrix()))
    unique=[];seen=set()
    for name,rotation in choices:
        key=tuple(np.round(rotation,10).ravel())
        if key not in seen:seen.add(key);unique.append((name,rotation))
    if not 30<=len(unique)<=100:raise ValueError('ORIENTATION_COUNT_OUTSIDE_30_TO_100')
    return unique


def bed_stability(mesh,frame,rotation,translation,cfg):
    band=cfg['orientation']['bed_contact_band_mm']
    verts=frame.vertices@rotation.T+translation
    downward=frame.face_normals@rotation[2]<-.99
    onbed=np.all(verts[frame.faces,:,][...,2]<=band+1e-8,axis=1)
    ids=np.flatnonzero(downward&onbed);area=float(frame.area_faces[ids].sum())
    centers=mesh.center_mass@rotation.T+translation
    points=verts[np.unique(frame.faces[ids]),:2] if len(ids) else np.empty((0,2))
    margin=-math.inf
    if len(points)>=3:
        try:
            poly=Polygon(points[ConvexHull(points).vertices]);point=Point(centers[:2])
            margin=float(poly.exterior.distance(point))*(1 if poly.covers(point) else -1)
        except Exception:pass
    return dict(bed_contact_area_mm2=area,support_polygon_margin_mm=margin,
        stable=bool(area>=cfg['orientation']['minimum_bed_contact_area_mm2'] and margin>=cfg['orientation']['minimum_support_polygon_margin_mm']))


def rank_orientations(mesh,frame,build_volume,access,cfg):
    settings=cfg['orientation'];build=np.array(build_volume,dtype=float);result=[]
    points=mesh.vertices[ConvexHull(mesh.vertices).vertices]
    for index,(name,rotation) in enumerate(orientation_rotations(cfg)):
        rotated=points@rotation.T;lo,hi=rotated.min(0),rotated.max(0);size=hi-lo
        translation=np.array([build[0]/2,build[1]/2,0.])-np.r_[((lo+hi)/2)[:2],lo[2]]
        transform=np.eye(4);transform[:3,:3]=rotation;transform[:3,3]=translation
        low=lo+translation;high=hi+translation
        fits=bool(np.all(size<=build-np.array([2*settings['bed_margin_mm'],2*settings['bed_margin_mm'],settings['top_margin_mm']])+1e-8))
        stability=bed_stability(mesh,frame,rotation,translation,cfg)
        h=mesh.triangles_center@rotation[2]+translation[2]
        down=(mesh.face_normals@rotation[2]<-np.cos(np.radians(settings['overhang_angle_deg'])))&(h>settings['bed_exemption_mm'])
        support=float(mesh.area_faces[down].sum())
        footprint=float(ConvexHull(rotated[:,:2]).volume)
        # Rotation cannot change physical access after removal from the bed.
        rank=(access['trapped_proxy_region_count'],-access['accessible_fraction'],support,float(size[2]),footprint,index)
        passed=fits and stability['stable'] and access['trapped_proxy_region_count']==0
        result.append(dict(candidate_id=index,method=name,transform_4x4=transform.tolist(),
            euler_xyz_deg=Rotation.from_matrix(rotation).as_euler('xyz',degrees=True).tolist(),
            bbox_min_mm=low.tolist(),bbox_max_mm=high.tolist(),bbox_extents_mm=size.tolist(),
            z_height_mm=float(size[2]),bed_footprint_mm2=footprint,downward_overhang_area_mm2=support,
            fits_build_volume=fits,**stability,lexicographic_rank=list(rank),
            passed=bool(passed),status='PASS' if passed else 'ORIENTATION_HARD_GATE_FAILED',
            connectivity_preserved=True,scale=1.,access_rotation_invariant_after_bed_removal=True))
    feasible=sorted([r for r in result if r['passed']],key=lambda r:r['lexicographic_rank'])
    top=feasible[:settings['top_count']]
    for rank,row in enumerate(top,1):row['rank']=rank
    return result,top


def interface_geometry_checks(before,after,inputs):
    import hashlib
    from .surface_continuity_qc import window,feature_edges,normal_jumps
    def geometry_hash(mesh,center,normal,radius):
        local=window(mesh.triangles_center,center,normal,radius*1.7,1.)
        rows=[]
        for triangle in mesh.triangles[local]:
            rows.append(tuple(sorted(tuple(p) for p in triangle)))
        return hashlib.sha256(np.array(sorted(rows),dtype=np.float64).tobytes()).hexdigest()
    results=[]
    for port in inputs['ports']:
        cap=port['accepted_attachment'];center=np.mean(cap['aligned_ring_print_mm'],axis=0);normal=np.array(cap['cap_normal']);r=port['radius_mm']
        oldhash=geometry_hash(before,center,normal,r);newhash=geometry_hash(after,center,normal,r)
        old,_=feature_edges(before,center,normal,r*1.7,1.,20);new,_=feature_edges(after,center,normal,r*1.7,1.,20)
        oldnormal,_=normal_jumps(before,center,normal,r*1.7,1.);newnormal,_=normal_jumps(after,center,normal,r*1.7,1.)
        passed=oldhash==newhash and old['total_length_mm']==new['total_length_mm'] and abs(oldnormal['p95_deg']-newnormal['p95_deg'])<1e-8
        results.append(dict(port_id=port['port_id'],before_geometry_hash=oldhash,after_geometry_hash=newhash,
            before_feature20_length_mm=old['total_length_mm'],after_feature20_length_mm=new['total_length_mm'],
            before_normal_p95_deg=oldnormal['p95_deg'],after_normal_p95_deg=newnormal['p95_deg'],
            status='PASS' if passed else 'FAIL'))
    return results
