"""V2 exact Boolean diagnostics, continuous assembly sweeps, and trim audits."""
from __future__ import annotations
import itertools
import numpy as np
import trimesh
import manifold3d as md
from scipy.spatial.transform import Rotation
from shapely.geometry import Polygon
from .casting_fixture_geometry import shape_mesh,bounds_box
from .casting_fixture_v2_geometry import world,port_transform,fuse
from .casting_fixture_audit import sweep_translation
from .support_removal import swept_path
from .casting_mold_qc import fcl_distance
from .support_inputs import save,table


def manifold(mesh):
    m=md.Manifold(md.Mesh64(np.ascontiguousarray(mesh.vertices,dtype=np.float64),np.ascontiguousarray(mesh.faces,dtype=np.uint64)))
    if m.status()!=md.Error.NoError:raise ValueError('BOOLEAN_INPUT_NOT_MANIFOLD '+str(m.status()))
    return m


def mesh(m):
    if m.status()!=md.Error.NoError:raise ValueError('BOOLEAN_FAILED '+str(m.status()))
    x=m.to_mesh64();return trimesh.Trimesh(x.vert_properties[:,:3],x.tri_verts,process=False)


def cad_manifold(shape):
    # CAD tessellation can produce 1e-14 mm differences on nominally shared
    # planes. Those generate sliver faces that collapse when STL is float32.
    # Snap only CAD operands to 1e-8 mm before Boolean operations. The frozen
    # ABS core vertices are NOT snapped, refitted, or otherwise changed.
    m=shape_mesh(shape);m.vertices=np.round(m.vertices,8);return manifold(m)


def choose_strategy(local_volume,strip_volume,port_count,strip_valid,limit=.1):
    if port_count>1 and strip_valid and strip_volume<=local_volume*(1+limit):
        return 'SIDE_TRIM_STRIP','One common cut; additional removed volume is within the explicit 10% allowance.'
    return 'LOCAL_PORT_NOTCH_CUT','Local notches remove less material; the common strip is ineligible or exceeds the explicit 10% allowance. Each notch needs three vertical cuts.'


def contact_zone(point,p,tol=.06,distal_stop=False):
    q=np.asarray(point)-np.asarray(p['straight_start']);n=np.asarray(p['axis']);s=float(q@n)
    radial=np.linalg.norm(q-s*n);limit=p['straight_length_mm']
    # Distal flat face is the single unilateral master axial datum; it is not
    # the upstream cap-to-cerebral-tree attachment.
    start=float((np.asarray(p['approved_contact_start'])-np.asarray(p['straight_start']))@n)
    end=float((np.asarray(p['approved_contact_end'])-np.asarray(p['straight_start']))@n)
    return bool(radial<=p['radius_mm']+tol and
                (start-tol<=s<=end+tol or (distal_stop and abs(s-limit)<=tol)))


def trim(l,g,cfg,out):
    out.mkdir(parents=True,exist_ok=True);H=l['high'][2];lo=l['low'];hi=l['high'];cav=g['cavity'];cvm=cad_manifold(cav);corem=manifold(l['core']);centralm=manifold(l['central'])
    localcuts={};strategies=[];byface={}
    for p in l['ports']:
        half=cfg['trim_notch_half_width_mm'];z=-p['trim_origin'][2]
        localcuts[p['port_id']]=world(bounds_box([-half,0,z-1],[half,p['wall_distance_from_trim_mm']+3,H+z+1]),p).intersect(cav)
        byface.setdefault(p['face'],[]).append(p)
    selected=[];choice={};cutgroups=[]
    for face,ps in byface.items():
        local=fuse([localcuts[p['port_id']] for p in ps]);p=ps[0];n=np.array(p['axis']);ax=int(np.argmax(abs(n)));a=lo.copy();b=hi.copy()
        if n[ax]>0:a[ax]=p['trim_origin'][ax]
        else:b[ax]=p['trim_origin'][ax]
        strip=bounds_box(a,b);sm=manifold(shape_mesh(strip));same_plane=all(abs(np.dot(np.array(v['trim_origin'])-p['trim_origin'],n))<1e-6 for v in ps)
        central_collision=float((sm^centralm).volume());valid=bool(same_plane and central_collision<1e-8)
        name,reason=choose_strategy(local.Volume(),strip.Volume(),len(ps),valid,cfg['side_strip_extra_volume_limit'])
        selected_shape=strip if name=='SIDE_TRIM_STRIP' else local;selected.append(selected_shape)
        cutgroups.append(dict(face=face,ports=[v['port_id'] for v in ps],strategy=name,cut_count=1 if name=='SIDE_TRIM_STRIP' else 3*len(ps)))
        for v in ps:choice[v['port_id']]=name
        strategies.append(dict(wall=face,ports=[v['port_id'] for v in ps],local_removed_envelope_mm3=local.Volume(),strip_removed_envelope_mm3=strip.Volume(),strip_valid=valid,strip_central_intersection_mm3=central_collision,selected=name,reason=reason,volume_tolerance_fraction=cfg['side_strip_extra_volume_limit']))
    cutter=fuse(selected);cutm=cad_manifold(cutter);final_region=cvm-cutm
    frame_m=cad_manifold(g['frame']);pla=frame_m
    for r in g['retainers'].values():pla=pla+cad_manifold(r)
    before=cvm-corem-pla;after=before-cutm;without_imprint=final_region-corem
    residual=float((pla^final_region).volume());missing=float((without_imprint-after).volume());extra=float((after-without_imprint).volume())
    # The final block is a vertical prism apart from vascular holes. An
    # extended prism detects a wall/guide overhang which would mechanically
    # hook it during upward release. Material adhesion remains unvalidated.
    # Extend this constant-XY diagnostic prism in Z directly. Going via the
    # 2D CrossSection/Clipper grid would re-quantize X/Y and invent micrometre
    # wall overlaps. No transform is applied to the vascular core here.
    release_prism=final_region.scale((1,1,(H+20)/H))
    interlock=float((frame_m^release_prism).volume())
    central_removed=float((cutm^centralm).volume());outside=float((corem-cvm).volume());contactrows=[];sections=[];planemesh=[];plane_records=[];wrong_port_volumes=[]
    for p in l['ports']:
        k=p['port_id'];n=np.array(p['axis']);pt=np.array(p['trim_origin']);pm=l['source_ports'][k];section=pm.section(plane_normal=n,plane_origin=pt)
        if section is None:raise ValueError('TRIM_PLANE_MISSES_EXTENSION_'+k)
        mat=port_transform(p);coords=(section.discrete[0]-pt)@mat[:3,:3];area=abs(Polygon(coords[:,[0,2]]).area);diam=2*np.sqrt(area/np.pi)
        # Finitely bounded cut surface: an infinite plane would cross an
        # unrelated parallel port outside a local notch. Audit the actual cut.
        half=cfg['trim_notch_half_width_mm'];z0=-pt[2]
        patch=world(bounds_box([-half,-1e-5,z0],[half,1e-5,H+z0]),p).intersect(cav);patchmesh=shape_mesh(patch)
        planes_for_port=[patchmesh]
        if choice[k]=='LOCAL_PORT_NOTCH_CUT':
            planemesh.append(patchmesh);plane_records.append(dict(port=k,kind='PORT_EXPOSURE',origin=pt.tolist(),normal=n.tolist()))
            for sign in [-1,1]:
                side=world(bounds_box([sign*half-1e-5,0,z0],[sign*half+1e-5,p['wall_distance_from_trim_mm']+.01,H+z0]),p).intersect(cav)
                # A corner notch already ends at the neighbouring wall: do
                # not count a fictional cut outside the actual PDMS cavity.
                if side.Volume()>1e-7:
                    smesh=shape_mesh(side);planes_for_port.append(smesh);planemesh.append(smesh)
                    plane_records.append(dict(port=k,kind='NOTCH_SIDE',origin=(pt+mat[:3,0]*sign*half).tolist(),normal=(mat[:3,0]*sign).tolist()))
        elif not any(x.get('wall')==p['face'] for x in plane_records):
            # Side-strip is one common full-width exposure plane.
            ax=int(np.argmax(abs(n)));aa=lo.copy();bb=hi.copy();aa[ax]=pt[ax]-1e-5;bb[ax]=pt[ax]+1e-5
            shared=shape_mesh(bounds_box(aa,bb));planemesh.append(shared);planes_for_port=[shared]
            plane_records.append(dict(port=k,wall=p['face'],kind='SHARED_EXPOSURE',origin=pt.tolist(),normal=n.tolist()))
        distances=[fcl_distance(l['central'],patchmesh)[0]]
        intended={v['port_id'] for v in byface[p['face']]} if choice[k]=='SIDE_TRIM_STRIP' else {k}
        for other,othermesh in l['source_ports'].items():
            if other not in intended:
                distances.extend(fcl_distance(othermesh,patchpart)[0] for patchpart in planes_for_port)
                wrong_port_volumes.append(float((manifold(shape_mesh(localcuts[k]))^manifold(othermesh)).volume()))
        nonport_gap=float(min(distances));cut_vol=float((manifold(shape_mesh(localcuts[k]))^manifold(pm)).volume())
        # Verify real final PDMS cut face surrounds the *actual* faceted hole.
        u=mat[:3,0];v=mat[:3,2];angles=np.linspace(0,2*np.pi,24,endpoint=False);rad=p['radius_mm']+.08
        ring=pt+rad*(np.cos(angles)[:,None]*u+np.sin(angles)[:,None]*v)
        finalmesh=mesh(after);inside=finalmesh.contains(ring-n*.025);outside_ring=finalmesh.contains(ring+n*.025)
        exposed=bool(inside.all() and not outside_ring.any() and cut_vol>0 and len(section.discrete)==1)
        row=dict(port=k,wall=p['face'],trim_strategy=choice[k],trim_plane_origin=pt.tolist(),trim_plane_normal=n.tolist(),extension_axis=n.tolist(),axis_plane_angle=90.,distance_to_bifurcation=p['distance_to_bifurcation_mm'],distance_to_curvature_mm=p['distance_to_curvature_mm'],distance_to_nonport_vascular_surface_mm=nonport_gap,usable_straight_length_after_trim=p['usable_straight_length_after_trim_mm'],fixture_imprint_final_side=residual,exposed_diameter_area_equivalent_mm=float(diam),nominal_diameter_mm=p['diameter_mm'],section_area_mm2=float(area),actual_surface_exposed=exposed,section_closed_contours=len(section.discrete),status='TRIM_GUIDE_GEOMETRY_VALIDATED' if exposed and nonport_gap>0 and residual<1e-7 else 'NEEDS_ADJUSTMENT')
        contactrows.append(row);sections.append(row.copy())
    proof=dict(fixture_imprint_remaining_after_trim_mm3=residual,ghost_channel_missing_volume_mm3=missing,ghost_channel_extra_volume_mm3=extra,central_vascular_removed_mm3=central_removed,core_outside_cavity_mm3=outside,all_ports_exposed=all(r['actual_surface_exposed'] for r in contactrows),minimum_nonport_vascular_trim_clearance_mm=min(r['distance_to_nonport_vascular_surface_mm'] for r in contactrows),before_pdms_volume_mm3=before.volume(),after_pdms_volume_mm3=after.volume(),removed_pdms_volume_mm3=before.volume()-after.volume(),final_pdms_components=len(after.decompose()),trim_plane_count=sum(x['cut_count'] for x in cutgroups),port_exposure_plane_count=len(contactrows),alignment_rib_count=len(g['ribs']),strategy_per_port=choice,strategies=strategies,cut_groups=cutgroups,cut_with_frame_attached=True,peripheral_to_final_boolean_intersection_mm3=float(((before-after)^after).volume()),physical_trim_validated=False,finite_patch_definition='Local notches: one transverse port-exposure plane plus two wall-to-port side cuts; no infinite-plane claim.')
    proof.update(trim_plane_count=len(plane_records),actual_cut_planes=plane_records,maximum_unrelated_port_cut_volume_mm3=max(wrong_port_volumes,default=0.),vertical_release_prism_PLA_intrusion_mm3=interlock,no_obvious_vertical_interlock=interlock<1e-7,cad_boolean_snap_mm=1e-8,frozen_core_snapped=False,fixture_imprint_remaining_after_trim=residual,minimum_nonport_vascular_trim_clearance=proof['minimum_nonport_vascular_trim_clearance_mm'])
    for group in proof['cut_groups']:group['cut_count']=sum(x['port'] in group['ports'] for x in plane_records)
    diagnostic=out.parent/'diagnostic';diagnostic.mkdir(exist_ok=True)
    diagnostic_qc=[]
    for name,obj in [('PDMS_before_trim',before),('PDMS_after_trim',after),('diagnostic_final_trimmed_PDMS',after),('sacrificial_trim_regions',cutm),('FINAL_PDMS_SIDE',final_region)]:
        path=diagnostic/(name+'.stl');mesh(obj).export(path);readback=trimesh.load_mesh(path)
        diagnostic_qc.append(dict(name=name,watertight=bool(readback.is_watertight),volume_error_mm3=abs(float(readback.volume)-obj.volume()),passed=bool(readback.is_watertight and abs(float(readback.volume)-obj.volume())<.05)))
    proof['diagnostic_stl_roundtrip']=diagnostic_qc;proof['diagnostic_mesh_qc_passed']=all(x['passed'] for x in diagnostic_qc)
    import pyvista as pv
    from .casting_fixture_review import poly
    poly(trimesh.util.concatenate(planemesh)).save(out/'trim_planes.vtp')
    trimesh.util.concatenate([shape_mesh(s) for s in g['ribs'].values()]).export(out/'trim_alignment_ribs.stl')
    table(out/'trimmed_port_sections.csv',sections);table(out.parent/'tables/port_trim_design.csv',contactrows);save(out/'trim_strategy.json',proof)
    return proof,contactrows,dict(before=mesh(before),after=mesh(after),cuts=mesh(cutm),final_region=mesh(final_region),planes=planemesh)


def audit(l,g,meshes,cfg,out):
    core=l['core'];central=l['central'];frame=meshes['frame'];ports=l['ports'];H=l['high'][2];contacts=[]
    print('V2 audit: continuous core insertion',flush=True)
    top=sweep_translation(core,[frame],[0,0,H+5],cfg['sweep_tolerance_mm']);top['status']='COLLISION_FREE_TOP_INSERTION_PASS' if top['passed'] else top['status']
    for name,part,p in [('FRAME',frame,None)]+[(kind+'_'+p['port_id'],shape_mesh(g[kind][p['port_id']]),p) for kind in ['supports','retainers','ribs'] for p in ports]:
        d,ca,cb=fcl_distance(core,part);centralgap,_,_=fcl_distance(central,part)
        pp=p or min(ports,key=lambda x:np.linalg.norm(ca-np.array(x['trim_origin'])))
        # Ribs and walls intentionally have no contact, even under nominal
        # clearance; near-contact cradle/retainer/stop points must be approved.
        intended=name.startswith(('supports','retainers')) or name=='FRAME'
        approved=contact_zone(ca,pp,distal_stop=pp['role']=='MASTER') if intended else d>.05
        contacts.append(dict(part=name,port_id=pp['port_id'],design_contact=intended,nominal_hard_contact=d<1e-6,minimum_core_clearance_mm=d,central_clearance_mm=centralgap,central_contact=centralgap<1e-6,closest_ABS_point=ca.tolist(),closest_PLA_point=cb.tolist(),approved_extension_zone=bool(approved)))
    present=[frame,core];insertion=[];locks=[]
    manager=trimesh.collision.CollisionManager();manager.add_object('FRAME',frame)
    for p in ports:
        k=p['port_id'];part=meshes['retainers'][k];slide=g['local'][k]['retainer_slide_mm'];a=np.eye(4);a[:3,3]=-np.array(p['axis'])*slide+[0,0,H+5];b=np.eye(4);b[:3,3]=-np.array(p['axis'])*slide
        r=swept_path(part,present,[a,b,np.eye(4)],tolerance=cfg['sweep_tolerance_mm']);insertion.append(dict(port_id=k,sequence=['TOP_DOWN','SHORT_WALLWARD_SLIDE'],wallward_slide_mm=slide,**r));present.append(part)
        moved=np.eye(4);moved[2,3]=.8;blocked=bool(manager.in_collision_single(part,transform=moved));locks.append(dict(port_id=k,up_0p8mm_blocked_by_inner_guides=blocked))
    cm=trimesh.collision.CollisionManager()
    for i,m in enumerate([frame,*meshes['retainers'].values()]):cm.add_object(str(i),m)
    pivot=np.mean([p['trim_origin'] for p in ports],axis=0);perturb=[]
    for name,ax,ang,dz in [('Z+',None,0,1.),('pitch+','y',2,0),('pitch-','y',-2,0),('roll+','x',2,0),('roll-','x',-2,0)]:
        m=np.eye(4)
        if ax:m[:3,:3]=Rotation.from_euler(ax,ang,degrees=True).as_matrix()
        m[:3,3]=pivot-m[:3,:3]@pivot+[0,0,dz];perturb.append(dict(perturbation=name,degrees=ang,translation_z_mm=dz,blocked=bool(cm.in_collision_single(core,transform=m))))
    # Analytic CAD difference certifies all four continuous full-height wall
    # slabs, not just a few test rays. Fins are allowed at the bottom perimeter;
    # there must be no floor in the retained central cavity.
    wall_missing=float(g['walls'].cut(g['frame']).Volume());lo=l['low'];hi=l['high'];small=bounds_box(central.bounds[0]*[1,1,0]+[0,0,.0001],np.r_[central.bounds[1,:2],.02]);bottom_intrusion=float(g['frame'].intersect(small).Volume())
    origins=np.array([[0,0,z] for z in np.linspace(.1,H-.1,41)]);wall_rays=[]
    for direction in [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]]:wall_rays.extend(frame.ray.intersects_any(origins,np.tile(direction,(len(origins),1))).tolist())
    result=dict(top_insertion=top,retainer_insertions=insertion,retainer_insertion_pass=all(x['passed'] for x in insertion),contacts=contacts,central_contact_count=sum(r['central_contact'] for r in contacts),all_design_contacts_on_extensions=all(r['approved_extension_zone'] for r in contacts if r['design_contact']),minimum_central_PLA_clearance_mm=min(r['central_clearance_mm'] for r in contacts),nominal_hard_contact_count=sum(r['nominal_hard_contact'] for r in contacts if r['design_contact']),support_contact_location_count=len(ports),four_walls_closed=wall_missing<1e-8 and all(wall_rays),wall_missing_volume_mm3=wall_missing,lateral_wall_ray_count=len(wall_rays),wall_ray_hits=sum(wall_rays),bottom_open=bottom_intrusion<1e-8,bottom_test_intrusion_mm3=bottom_intrusion,top_open=top['passed'],anti_float=dict(retainer_count=len(ports),geometric_Z_constraint=all(x['blocked'] for x in perturb) and all(x['up_0p8mm_blocked_by_inner_guides'] for x in locks),finite_pose_checks=perturb,retainer_lock_checks=locks,physical_force_validated=False,status='ANTI_FLOAT_FORCE_NOT_VALIDATED'))
    table(out/'fixture_contact_audit.csv',contacts);table(out/'core_insertion_clearance.csv',[top]);table(out/'retainer_insertion_clearance.csv',insertion);save(out/'assembly_audit.json',result)
    return result
