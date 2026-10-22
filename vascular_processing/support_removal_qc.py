"""Top-only access to sampled support-risk patches; not removal certification."""
import itertools
import numpy as np
import trimesh
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation
from .support_access_qc import capsule,line_visible
from .print_frame_qc import bed_stability


def overhang_threshold(discovery,cfg):
    value=discovery.get('flattened',{}).get('process',{}).get('support_threshold_angle')
    try:
        angle=float(value)
        if not 0<angle<90:raise ValueError('auto/invalid threshold')
        return dict(angle_from_horizontal_deg=angle,source='ACTIVE_BAMBU_PROCESS_PROFILE',profile_value=value,
            interpretation='Downward face whose plane angle from the horizontal is below the threshold',
            reference='https://github.com/bambulab/BambuStudio/blob/master/src/libslic3r/PrintConfig.cpp')
    except (TypeError,ValueError):
        return dict(angle_from_horizontal_deg=cfg['support']['heuristic_threshold_from_horizontal_deg'],
            source='OVERHANG_THRESHOLD_HEURISTIC',profile_value=value)


def make_regions(mesh,partition,cfg):
    regions=[];length=cfg['support']['region_length_mm']
    for branch,ids in {**partition['groups'],**partition.get('inner_wall_groups',{})}.items():
        if not len(ids):continue
        # Spatial bins are grouping only. Distances/collisions use real surfaces.
        # PCA binning remains well-defined for the fixed curved port surfaces.
        points=mesh.triangles_center[ids];center=points.mean(0)
        _,_,vh=np.linalg.svd(points-center,full_matrices=False);axis=vh[0]
        if axis[np.argmax(abs(axis))]<0:axis=-axis
        coordinate=(points-center)@axis;bins=np.floor((coordinate-coordinate.min())/length).astype(int)
        for bucket in np.unique(bins):
            faces=ids[bins==bucket]
            regions.append(dict(region_id=f'{branch}_{bucket:02}',vascular_branch=branch,face_ids=faces,
                center_mm=np.average(mesh.triangles_center[faces],axis=0,weights=mesh.area_faces[faces])))
    return regions


def representative_faces(mesh,ids,limit):
    """Cover both surface-normal directions and positions deterministically."""
    if len(ids)<=limit:return ids
    directions=np.array(list(itertools.product([-1,0,1],repeat=3)),dtype=float)
    directions=directions[np.linalg.norm(directions,axis=1)>0];directions/=np.linalg.norm(directions,axis=1)[:,None]
    normals=mesh.face_normals[ids];selected=[]
    for direction in directions:
        score=normals@direction
        # Ties prefer the larger real triangle, not arbitrary tiny Boolean slivers.
        candidate=ids[np.lexsort((ids,-mesh.area_faces[ids],-np.round(score,6)))[0]]
        if candidate not in selected:selected.append(candidate)
    for face in ids[np.argsort(-mesh.area_faces[ids])]:
        if face not in selected:selected.append(face)
        if len(selected)>=limit:break
    return np.array(selected[:limit],dtype=int)


def top_path(surface,normal,diameter,box,mesh,manager,cfg,only_vertical=False):
    radius=diameter/2;gap=cfg['support_removal_probe']['surface_standoff_mm'];end=surface+normal*(radius+gap)
    top=box['inner'][1,2];height=top+radius+gap-end[2]
    if height<=0:return None
    angles=[0] if only_vertical else cfg['support']['cone_angles_deg']
    for angle in angles:
        for azimuth in ([0] if angle==0 else np.arange(0,360,cfg['support']['cone_azimuth_step_deg'])):
            shift=height*np.tan(np.radians(angle));entry=end+np.array([shift*np.cos(np.radians(azimuth)),shift*np.sin(np.radians(azimuth)),height])
            at_top=end+(entry-end)*(top-end[2])/height
            if np.any(at_top[:2]-radius<box['inner'][0,:2]) or np.any(at_top[:2]+radius>box['inner'][1,:2]):continue
            if not line_visible(box['mesh'],mesh,end,entry):continue
            shape=capsule(entry,end,radius,cfg['support_removal_probe']['capsule_sections'])
            if not manager.in_collision_single(shape):
                return dict(start=entry.tolist(),end=end.tolist(),cone_angle_deg=float(angle),diameter_mm=diameter)
    return None


def face_access(mesh,face,box,manager,cfg):
    surface=mesh.triangles_center[face];normal=mesh.face_normals[face];gap=cfg['support_removal_probe']['surface_standoff_mm']
    endpoint=surface+normal*gap;entry=endpoint.copy();entry[2]=box['inner'][1,2]+gap
    vertical=bool(line_visible(box['mesh'],mesh,endpoint,entry))
    diameter=cfg['support_removal_probe']['diameter_mm']
    path=top_path(surface,normal,diameter,box,mesh,manager,cfg)
    return dict(face_id=int(face),normal=normal.tolist(),surface_mm=surface.tolist(),vertical_visible=vertical,
        probe_accessible=path is not None,path=path)


def assess_regions(mesh,box,regions,rotation,threshold,cfg,cache=None,manager=None):
    if cache is None:cache={}
    if manager is None:
        manager=trimesh.collision.CollisionManager();manager.add_object('vascular',mesh);manager.add_object('box',box['mesh'])
    rows=[];paths=[];risk_faces=[];cosine=np.cos(np.radians(threshold['angle_from_horizontal_deg']))
    for region in regions:
        ids=np.array(region['face_ids']);selected=ids[(mesh.face_normals[ids]@rotation[2])<-cosine]
        area=float(mesh.area_faces[selected].sum())
        if area<cfg['support']['minimum_risk_region_area_mm2']:continue
        risk_faces.extend(selected.tolist());tests=[]
        for face in representative_faces(mesh,selected,cfg['support']['sample_faces_per_region']):
            if int(face) not in cache:cache[int(face)]=face_access(mesh,int(face),box,manager,cfg)
            tests.append(cache[int(face)])
            # A usable probe path is stronger than the permitted visibility
            # proxy; stop once found, keeping the full risk area in the row.
            if tests[-1]['probe_accessible']:break
        vertical=any(t['vertical_visible'] for t in tests);probe=any(t['probe_accessible'] for t in tests)
        # The task's hard definition is retained, including on synthetic traps.
        trapped=not vertical and not probe
        best=next((t for t in tests if t['probe_accessible']),next((t for t in tests if t['vertical_visible']),tests[0]))
        row=dict(region_id=region['region_id'],vascular_branch=region['vascular_branch'],overhang_area_mm2=area,
            vertical_visible_to_top=vertical,probe_accessible_to_top=probe,
            minimum_access_width_mm=cfg['support_removal_probe']['diameter_mm'] if probe else None,
            access_width_interpretation='Verified lower bound along an accepted capsule path; null means no width established',
            tested_surface_samples=len(tests),representative_surface_mm=best['surface_mm'],
            status='TRAPPED_SUPPORT_RISK' if trapped else 'GEOMETRY_PROXY_SUPPORT_REMOVAL_PASS')
        rows.append(row)
        if best['path']:paths.append(dict(region_id=region['region_id'],**best['path']))
    result=dict(candidate_regions=len(rows),trapped_risk_count=sum(r['status']=='TRAPPED_SUPPORT_RISK' for r in rows),
        top_accessible_count=sum(r['vertical_visible_to_top'] or r['probe_accessible_to_top'] for r in rows),
        vertical_visible_count=sum(r['vertical_visible_to_top'] for r in rows),probe_accessible_count=sum(r['probe_accessible_to_top'] for r in rows),
        internal_support_risk_area_mm2=sum(r['overhang_area_mm2'] for r in rows),analysis_level='GEOMETRIC_PROXY',
        removal_exit='TOP_OPENING_ONLY_IN_CASTING_COORDINATES',probe_diameter_mm=cfg['support_removal_probe']['diameter_mm'],
        threshold=threshold,rows=rows,paths=paths,risk_face_ids=risk_faces,
        disclaimer='Probe accessibility is a geometry proxy, not a guarantee that real support can be removed.')
    result['status']='TRAPPED_SUPPORT_RISK' if result['trapped_risk_count'] else 'GEOMETRY_PROXY_SUPPORT_REMOVAL_PASS'
    return result


def rotations(cfg):
    choices=[]
    for axis,angle in [('z',0),('x',90),('x',-90),('y',90),('y',-90)]:
        seed=Rotation.from_euler(axis,angle,degrees=True).as_matrix()
        for yaw in cfg['orientation']['yaw_angles_deg']:
            choices.append((('TOP_OPENING_UP' if angle==0 else f'SIDE_{axis}{angle}')+f'_yaw{yaw}',Rotation.from_euler('z',yaw,degrees=True).as_matrix()@seed))
    for axis in ('x','y'):
        for angle in cfg['orientation']['tilt_degrees']:
            for yaw in cfg['orientation']['yaw_angles_deg']:
                choices.append((f'TILT_{axis}{angle}_yaw{yaw}',Rotation.from_euler('z',yaw,degrees=True).as_matrix()@Rotation.from_euler(axis,angle,degrees=True).as_matrix()))
    if not 30<=len(choices)<=80:raise ValueError('ORIENTATION_COUNT_OUTSIDE_REQUESTED_RANGE')
    return choices


def orientation_row(index,name,rotation,mesh,box,build,access,cfg):
    settings=cfg['orientation'];vertices=mesh.vertices@rotation.T;lo=vertices.min(0);hi=vertices.max(0);size=hi-lo;build=np.array(build)
    translation=np.r_[build[:2]/2-(lo+hi)[:2]/2,-lo[2]];transform=np.eye(4);transform[:3,:3]=rotation;transform[:3,3]=translation
    fits=bool(np.all(size<=build-[2*settings['bed_margin_mm'],2*settings['bed_margin_mm'],settings['top_margin_mm']]))
    stability=bed_stability(mesh,box['mesh'],rotation,translation,cfg)
    center=mesh.triangles_center@rotation.T+translation
    downward=mesh.face_normals@rotation[2]<-np.cos(np.radians(access['threshold']['angle_from_horizontal_deg']))
    above_bed=center[:,2]>settings['bed_contact_band_mm'];total=float(mesh.area_faces[downward&above_bed].sum())
    internal=access['internal_support_risk_area_mm2'];external=max(0,total-internal)
    footprint=float(ConvexHull(vertices[:,:2]).volume)
    rank=[access['trapped_risk_count'],round(internal,6),-access['probe_accessible_count'],round(external,6),round(total,6),round(float(size[2]),6),round(footprint,6),index]
    return dict(candidate_id=index,method=name,transform_4x4=transform.tolist(),casting_restore_transform=np.linalg.inv(transform).tolist(),
        euler_xyz_deg=Rotation.from_matrix(rotation).as_euler('xyz',degrees=True).tolist(),bbox_min_mm=(lo+translation).tolist(),bbox_max_mm=(hi+translation).tolist(),
        bbox_extents_mm=size.tolist(),z_height_mm=float(size[2]),bed_footprint_mm2=footprint,**stability,fits_build_volume=fits,
        trapped_risk_count=access['trapped_risk_count'],internal_support_risk_area_mm2=internal,external_support_area_mm2=external,
        downward_overhang_area_mm2=total,top_accessible_count=access['top_accessible_count'],probe_accessible_count=access['probe_accessible_count'],
        lexicographic_rank=rank,geometry_placement_pass=bool(fits and stability['stable']),
        passed=bool(fits and stability['stable'] and access['trapped_risk_count']==0),support=access,scale=1.)
