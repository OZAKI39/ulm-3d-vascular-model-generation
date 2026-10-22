"""Bottomless Y/straight molds with external stopped dovetails and cross keys.

A separate flat casting substrate defines z=0; no substrate is exported as a
printable component. CAD channel dimensions are reused, not rescaled. Each
external sliding post captures its core collar and is locked by a removable
cross key. Geometric capture and assembly access are tested, not inferred from
watertightness. Physical fit/load capacity still require a printed trial.
"""
from __future__ import annotations
import itertools
import numpy as np
import trimesh
import pyvista as pv
from shapely.geometry import Polygon, box as polygon_box
from scipy.spatial.transform import Rotation
from . import y_channel_casting as y
from . import straight_channel_casting as s
from . import y_uniform_channel as yu

MODE = 'bottomless_external_dovetail'
PORT_NAMES = {'Y': ['INLET', 'OUTLET_A', 'OUTLET_B'], 'Straight': ['INLET', 'OUTLET1']}


def part_names(kind):
    return [kind+'_lumen_core', kind+'_bottomless_box_wall2mm',
            *[kind+'_dovetail_retainer_'+n for n in PORT_NAMES[kind]],
            *[kind+'_dovetail_key_'+n for n in PORT_NAMES[kind]]]


def validate_config(cfg):
    if cfg.get('fixture_mode') != MODE or cfg['base_mm'] != 0:
        raise ValueError('BOTTOMLESS_REQUIRES_ZERO_INTERNAL_BASE')
    if cfg.get('geometry_mode', 'uniform_diameter') != 'uniform_diameter':
        raise ValueError('BOTTOMLESS_REQUIRES_UNIFORM_CHANNEL')
    if cfg.get('chamber_mode', 'closed_wall_internal_end') not in ['closed_wall_internal_end','closed_wall_two_vertical_ports']:
        raise ValueError('BOTTOMLESS_REQUIRES_NO_SIDE_HOLES')
    # Unlike the retired wall clips, these arms have no feet below the rim.
    # Keep a real liquid-to-arm gap while reducing only unused headroom.
    if not np.isfinite(cfg['freeboard_mm']) or cfg['freeboard_mm'] < .8-1e-9:
        raise ValueError('BOTTOMLESS_FREEBOARD_AT_LEAST_0P8_MM')
    gap=cfg.get('retainer_above_rim_mm',0.)
    if not np.isfinite(gap) or gap<0 or ('retainer_above_rim_mm' in cfg and gap==0):
        raise ValueError('INVALID_RETAINER_ABOVE_RIM')
    if limits := cfg.get('height_limits'):
        for key in ['fixture_max_mm', 'cad_reserve_mm', 'core_max_mm',
                    *(['retainer_max_mm'] if 'retainer_max_mm' in limits else [])]:
            if not np.isfinite(limits[key]) or limits[key] <= 0:
                raise ValueError('INVALID_HEIGHT_LIMIT: '+key)
        if not limits['cad_reserve_mm'] < limits['fixture_max_mm'] <= limits['core_max_mm']:
            raise ValueError('INVALID_HEIGHT_BUDGET')
        if 'retainer_max_mm' in limits and not limits['cad_reserve_mm']<limits['retainer_max_mm']<=limits['core_max_mm']:
            raise ValueError('INVALID_RETAINER_HEIGHT_BUDGET')
    p = cfg['dovetail']
    for key in ['external_foot_mm', 'fit_clearance_mm', 'collar_capture_clearance_mm', 'upper_fork_mm']:
        if not np.isfinite(p[key]) or p[key] <= 0: raise ValueError('INVALID_DOVETAIL_DIMENSION: '+key)
    if not .15 <= p['fit_clearance_mm'] <= .3: raise ValueError('DOVETAIL_CLEARANCE_OUT_OF_DESIGN_RANGE')
    if p['external_foot_mm'] != 2.: raise ValueError('THIS_FIXTURE_REQUIRES_2MM_EXTERNAL_FOOT')
    if p['collar_capture_clearance_mm']+p['upper_fork_mm']+cfg['collar_thickness_mm'] >= cfg['port_above_retainer_mm']:
        raise ValueError('UPPER_FORK_BLOCKS_PORT_TIP')
    if cfg.get('magnet_slots'):
        from .channel_magnet_slots import validate_config as validate_magnets
        validate_magnets(cfg)


def prism(points, z0, z1):
    poly = points if isinstance(points, Polygon) else Polygon(points)
    m = trimesh.creation.extrude_polygon(poly, z1-z0)
    m.apply_translation([0,0,z0])
    return m


def moved(mesh, transform):
    m=mesh.copy();m.apply_transform(transform);return m


def shifted(mesh, delta):
    m=mesh.copy();m.apply_translation(delta);return m


def part_height_limit(cfg, name):
    limits=cfg.get('height_limits')
    if limits is None:return None
    if 'core' in name:return limits['core_max_mm']
    maximum=limits.get('retainer_max_mm',limits['fixture_max_mm']) if 'retainer' in name else limits['fixture_max_mm']
    return maximum-limits['cad_reserve_mm']


def casting_height_qc(cfg, parts):
    """Measure assembled tops from the external substrate, not part extents.

    Reserve is a design allowance without taking credit for ABS shrinkage.
    An external casting plate's thickness is NOT included in this Z0 datum.
    """
    limits=cfg.get('height_limits')
    rows={name:dict(bottom_z_mm=float(m.bounds[0,2]),top_z_mm=float(m.bounds[1,2]),
                   separate_part_height_mm=float(m.extents[2]),core_exception='core' in name,
                   raised_retainer_exception='retainer' in name and 'retainer_max_mm' in (limits or {}),
                   cad_top_limit_mm=part_height_limit(cfg,name))
          for name,m in parts.items()}
    fixture_top=max(r['top_z_mm'] for r in rows.values() if not r['core_exception'])
    core_top=max(r['top_z_mm'] for r in rows.values() if r['core_exception'])
    budget=None if limits is None else limits['fixture_max_mm']-limits['cad_reserve_mm']
    passed=limits is None or all(r['top_z_mm']<=r['cad_top_limit_mm']+1e-5 and r['bottom_z_mm']>=-1e-5 for r in rows.values())
    return dict(enforced=limits is not None,reference='External casting substrate top = Z0',
        includes_external_substrate_thickness=False,limits=limits,fixture_cad_limit_mm=budget,
        fixture_limit_applies_to='Frame and keys; retainers use retainer_max_mm when configured',
        frame_and_key_top_z_mm=max(r['top_z_mm'] for n,r in rows.items() if 'core' not in n and 'retainer' not in n),
        retainer_top_z_mm=max(r['top_z_mm'] for n,r in rows.items() if 'retainer' in n),
        fixture_top_z_mm=fixture_top,core_top_z_mm=core_top,parts=rows,passed=passed,
        physical_height_validated=False)


def outward(name):
    return np.array([-1.,0.]) if name=='INLET' else np.array([0.,-1.]) if name=='OUTLET_B' else np.array([0.,1.])


def joint(port, inner, rim, cfg):
    """Local x points out of the chamber; z is the dovetail insertion axis."""
    center=np.asarray(port['top_center_mm'][:2]);n=outward(port['name']);axis=int(np.argmax(abs(n)))
    inner_wall=inner[0 if n[axis]<0 else 1,axis]
    d=abs(inner_wall-center[axis]);A=d+cfg['wall_mm']
    p=cfg['dovetail'];foot=p['external_foot_mm'];fit=p['fit_clearance_mm']
    # A 4.4 mm deep dovetail has 8.0 mm back width and 5.0 mm throat width.
    # Its 1.4 mm solid backing is OUTSIDE the uninterrupted 2 mm mold wall.
    male=Polygon([(A+1.4,-4),(A+5.8,-2.5),(A+5.8,2.5),(A+1.4,4)])
    cavity=male.buffer(fit,join_style=2).union(polygon_box(A+5.79,-2.5-fit,A+8,2.5+fit))
    socket=y.union([y.box([[A-.15,-6,0],[A+5.8,6,rim]]),
                    y.box([[A-.15,-7,0],[A+8,7,foot]])])
    # One blind lower stop, open upper entry, no hole into the PDMS chamber.
    key_z=5.;slot=y.box([[A+2.4,-8,key_z-1.2],[A+4.8,8,key_z+1.2]])
    socket=trimesh.boolean.difference([socket, prism(cavity,foot,rim+1),slot],engine='manifold')
    arm_bottom=y.retainer_arm_bottom(cfg,rim)
    collar_z=arm_bottom+cfg['retainer_thickness_mm'];capture_top=collar_z+cfg['collar_thickness_mm']+p['collar_capture_clearance_mm']
    upper=capture_top+p['upper_fork_mm']
    arm=y.union([y.box([[-3.5,-3.5,arm_bottom],[A+7.2,3.5,upper]]),
                 y.box([[A+1.2,-4.3,arm_bottom],[A+7.2,4.3,upper]])])
    stem_radius=port['radius_mm']+cfg['retainer_fit_clearance_mm']
    collar_radius=port['radius_mm']+cfg['collar_radial_width_mm']+p['collar_capture_clearance_mm']
    stem_cut=y.union([y.cylinder([0,0],stem_radius,rim-.1,upper+.1),
                     y.box([[-4,-stem_radius,rim-.1],[0,stem_radius,upper+.1]])])
    collar_cut=y.union([y.cylinder([0,0],collar_radius,collar_z,capture_top),
                       y.box([[-4,-collar_radius,collar_z],[0,collar_radius,capture_top]])])
    arm=trimesh.boolean.difference([arm,stem_cut,collar_cut],engine='manifold')
    # Extend the external column, keeping rail engagement, lower stop and key
    # at their existing heights. The arm floats above the full box rim.
    post=y.union([prism(male,foot,arm_bottom+.1),y.box([[A+5.5,-2,foot],[A+7.2,2,arm_bottom+.1]])])
    retainer=trimesh.boolean.difference([y.union([post,arm]),slot],engine='manifold')
    # A headed tapered cross key. Its taper contacts the entry roof when seated;
    # it is not modeled as a spring or as an unverified interference snap fit.
    key_profile=Polygon([(-6.05,key_z-1.1),(7,key_z-1.1),(7,key_z+.7),(-6.05,key_z+1.2)])
    key=trimesh.creation.extrude_polygon(key_profile,2.)
    key.vertices=key.vertices[:,[2,0,1]]+[A+2.6,0,0]
    key=y.union([key,y.box([[A+1.8,-7.55,key_z-1.1],[A+5.4,-6.05,key_z+2.1]])])
    T=np.eye(4);T[:2,0]=n;T[:2,1]=[-n[1],n[0]];T[:2,3]=center
    record=dict(port=port['name'],local_to_casting=T,outward_axis=np.r_[n,0],
        key_insertion_axis=T[:3,1],inner_wall_distance_mm=d,wall_outer_local_x_mm=A,
        dovetail_back_width_mm=8.,dovetail_throat_width_mm=5.,dovetail_depth_mm=4.4,
        dovetail_engagement_mm=rim-foot,normal_fit_clearance_mm=fit,stop_z_mm=foot,
        collar_support_z_mm=collar_z,collar_capture_top_z_mm=capture_top,
        box_rim_z_mm=rim,arm_bottom_z_mm=arm_bottom,arm_to_rim_clearance_mm=arm_bottom-rim,
        upper_fork_top_z_mm=upper,key_center_z_mm=key_z,
        external_foot_thickness_mm=foot,external_foot_extension_mm=8.,
        key_slot_size_xz_mm=[2.4,2.4],key_nominal_width_mm=2.,
        local_key_head_center_mm=[A+3.6,-6.8,key_z+.5])
    return moved(socket,T),moved(retainer,T),moved(key,T),record


def frame_and_joints(cfg,inner,rim,ports,kind):
    outer=inner+[[-cfg['wall_mm']]*2,[cfg['wall_mm']]*2]
    frame=trimesh.boolean.difference([y.box([np.r_[outer[0],0],np.r_[outer[1],rim]]),
        y.box([np.r_[inner[0],-1],np.r_[inner[1],rim+1]])],engine='manifold')
    mounts=[];retainers={};keys={};records=[]
    for p in ports:
        socket,retainer,key,record=joint(p,inner,rim,cfg)
        mounts.append(socket);retainers[kind+'_dovetail_retainer_'+p['name']]=retainer
        keys[kind+'_dovetail_key_'+p['name']]=key;records.append(record)
    return y.union([frame,*mounts]),retainers,keys,records


def retainer_rim_clearance_qc(cfg,parts,inner,rim,records,kind):
    """Measure the underside over each wall and the entire air-gap volume.

    External dovetail columns are intentionally excluded from the wall ring;
    their contact with the guide/stop is the required mechanical attachment.
    """
    gap=cfg.get('retainer_above_rim_mm',0.)
    outer=inner+[[-cfg['wall_mm']]*2,[cfg['wall_mm']]*2]
    air=None
    if gap>2e-4:
        air=trimesh.boolean.difference([
            y.box([np.r_[outer[0],rim+1e-4],np.r_[outer[1],rim+gap-1e-4]]),
            y.box([np.r_[inner[0],rim],np.r_[inner[1],rim+gap]])],engine='manifold')
    rows=[]
    for r in records:
        post=parts[kind+'_dovetail_retainer_'+r['port']];T=np.asarray(r['local_to_casting'])
        wall_middle=r['inner_wall_distance_mm']+cfg['wall_mm']/2
        origins=np.array([(T@np.array([wall_middle,t,rim-.01,1]))[:3] for t in [-3.,0.,3.]])
        hits,ids,_=post.ray.intersects_location(origins,np.tile([0,0,1],(len(origins),1)))
        measured=[]
        for i in range(len(origins)):
            z=hits[ids==i,2];z=z[z>=rim-1e-5]
            measured.append(None if len(z)==0 else float(z.min()-rim))
        occupied=0. if air is None else s.volume_intersection(post,air)
        passed=all(v is not None and abs(v-gap)<1e-4 for v in measured) and occupied<1e-4
        rows.append(dict(port=r['port'],measured_wall_to_arm_gaps_mm=measured,
            retainer_in_required_air_gap_mm3=occupied,passed=passed))
    return dict(requested_gap_mm=gap,box_rim_z_mm=rim,ray_checks_per_port=3,
        scope='Arm above complete four-wall rim; external dovetail connection remains engaged',
        parts=rows,passed=all(r['passed'] for r in rows))


def fixture_qc(cfg,parts,inner,fill_z,rim,records,kind):
    frame=parts[kind+'_bottomless_box_wall2mm'];core=parts[kind+'_lumen_core']
    # Test the entire central prism, not just the presence of a hole at centre.
    interior=y.box([np.r_[inner[0]+1e-4,-.5],np.r_[inner[1]-1e-4,rim+.5]])
    internal_abs=s.volume_intersection(frame,interior)
    outer=inner+[[-cfg['wall_mm']]*2,[cfg['wall_mm']]*2]
    missing={}
    for axis,side in itertools.product([0,1],[0,1]):
        b=np.array([np.r_[outer[0],0],np.r_[outer[1],rim]])
        if side==0:b[1,axis]=inner[0,axis]
        else:b[0,axis]=inner[1,axis]
        slab=y.box(b);missing[str(axis)+'_'+str(side)]=max(0.,float(slab.volume)-s.volume_intersection(frame,slab))
    probes=np.array([[x,v,-.5] for x in np.linspace(inner[0,0]+.1,inner[1,0]-.1,9)
                    for v in np.linspace(inner[0,1]+.1,inner[1,1]-.1,9)])
    bottom_hits=frame.ray.intersects_location(probes,np.tile([0,0,1],(len(probes),1)))[0]
    fill=y.box([np.r_[inner[0],0],np.r_[inner[1],fill_z]])
    fixture_pdms={name:s.volume_intersection(m,fill) for name,m in parts.items() if 'core' not in name}
    collisions={a+' / '+b:s.volume_intersection(parts[a],parts[b]) for a,b in itertools.combinations(parts,2)}
    joints=[]
    for r in records:
        name=r['port'];post=parts[kind+'_dovetail_retainer_'+name];key=parts[kind+'_dovetail_key_'+name]
        n=np.array(r['outward_axis']);t=np.array(r['key_insertion_axis'])
        # Retainers first slide onto their collars OUTSIDE the mold. Then the
        # whole core/retainer group lowers together through the rail entries.
        collar_insertion=[s.volume_intersection(shifted(post,n*travel),core) for travel in np.linspace(0,12,13)]
        vertical_insertion=[s.volume_intersection(shifted(post,[0,0,lift]),frame) for lift in np.linspace(0,rim+2,15)]
        key_insertion=[s.volume_intersection(shifted(key,-t*travel),y.union([frame,post])) for travel in [0,1,2,4,8,16]]
        lateral_capture=[s.volume_intersection(shifted(post,delta),frame) for delta in [n*.8,-n*.8,t*.8,-t*.8]]
        lifted=shifted(post,[0,0,.7])
        joints.append(dict(port=name,collar_lateral_insertion_intersection_mm3=collar_insertion,
            vertical_rail_insertion_intersection_mm3=vertical_insertion,key_insertion_intersection_mm3=key_insertion,
            lateral_capture_intersection_mm3=lateral_capture,
            locked_lift_key_intersection_mm3=s.volume_intersection(lifted,key),
            unlocked_lift_frame_intersection_mm3=s.volume_intersection(lifted,frame),
            downward_stop_intersection_mm3=s.volume_intersection(shifted(post,[0,0,-.5]),frame),
            collar_upward_capture_intersection_mm3=s.volume_intersection(shifted(core,[0,0,.5]),post)))
    core_path=[s.volume_intersection(shifted(core,[0,0,lift]),frame) for lift in np.linspace(0,rim+2,15)]
    rim_clearance=retainer_rim_clearance_qc(cfg,parts,inner,rim,records,kind)
    gates=dict(no_internal_abs_bottom=internal_abs<1e-4 and len(bottom_hits)==0,
        raised_arms_clear_of_box_rim=rim_clearance['passed'],
        four_unbroken_walls=max(missing.values())<1e-3,
        fixtures_outside_pdms=max(fixture_pdms.values())<1e-4,
        no_assembly_penetration=max(collisions.values())<1e-4,
        core_lowering_clear=max(core_path)<1e-4,
        collar_side_assembly_clear=all(max(r['collar_lateral_insertion_intersection_mm3'])<1e-4 for r in joints),
        dovetail_vertical_assembly_clear=all(max(r['vertical_rail_insertion_intersection_mm3'])<1e-4 for r in joints),
        cross_key_assembly_clear=all(max(r['key_insertion_intersection_mm3'])<1e-4 for r in joints),
        dovetail_lateral_capture=all(min(r['lateral_capture_intersection_mm3'])>.01 for r in joints),
        cross_key_blocks_lift=all(r['locked_lift_key_intersection_mm3']>.01 for r in joints),
        rail_releases_after_key_removal=all(r['unlocked_lift_frame_intersection_mm3']<1e-4 for r in joints),
        positive_downward_stops=all(r['downward_stop_intersection_mm3']>.01 for r in joints),
        collar_upward_capture=all(r['collar_upward_capture_intersection_mm3']>.01 for r in joints))
    return dict(gates=gates,internal_frame_volume_mm3=internal_abs,open_bottom_ray_count=len(probes),
        retainer_rim_clearance_qc=rim_clearance,
        open_bottom_ray_hits=len(bottom_hits),wall_missing_volume_mm3=missing,
        fixtures_inside_pdms_mm3=fixture_pdms,assembly_intersections_mm3=collisions,
        joints=joints,core_lowering_intersection_mm3=core_path,
        physical_fit_validated=False,physical_retention_load_validated=False,
        external_substrate_required=True,casting_seal_validated=False)


def geometry(cfg,out,kind):
    validate_config(cfg)
    if kind=='Y':
        source,_,planes=y.import_source(cfg)
        core_cad,channel_cad,ports,heights,samples=yu.build(cfg,planes)
        core,channel=s.cad_mesh(core_cad,cfg),s.cad_mesh(channel_cad,cfg)
        inner=y.enforce_pdms_width(channel.bounds[:,:2]+[[-cfg['side_clearance_mm']]*2,[cfg['side_clearance_mm']]*2],cfg)
        rim,fill_z=heights['rim_z_mm'],heights['fill_z_mm'];design=None
        mounting_ports=ports
    else:
        design=s.design(cfg);core_cad,channel_cad=s.build_core(cfg,design)
        core,channel=s.cad_mesh(core_cad,cfg),s.cad_mesh(channel_cad,cfg)
        inner=design['inner'];rim,fill_z=design['rim'],design['fill']
        heights=dict(fill_z_mm=fill_z,rim_z_mm=rim,port_top_z_mm=design['top'],pdms_floor_z_mm=0.)
        mounting_ports=[dict(name=name,radius_mm=design['radius'],top_center_mm=[x,0,design['top']])
                        for name,x in [('INLET',design['inlet_x']),('OUTLET1',design['junction_x'])]]
        ports=[dict(name=p['name'],label_id=i+2,top_center_mm=p['top_center_mm'],center_mm=p['top_center_mm'],outward_axis=[0,0,1])
               for i,p in enumerate(mounting_ports)]
        if s.two_top_ports(cfg):
            from .straight_two_top_ports import section_samples
            samples=section_samples(design)
        else:
            ports.append(dict(name='OUTLET2',label_id=4,center_mm=[design['horizontal_end'],0,design['center_z']],outward_axis=[1,0,0]))
            samples=[dict(region='straight',center_mm=[x,0,design['center_z']],normal=[1,0,0])
                     for x in [-13,-7,0,6,9,15,design['horizontal_end']-.5]]
            for angle in np.linspace(.12,np.pi/2-.12,5):
                samples.append(dict(region='elbow',center_mm=[design['bend_x']-design['bend_radius']*np.sin(angle),0,
                    design['center_z']+design['bend_radius']*(1-np.cos(angle))],normal=[-np.cos(angle),0,np.sin(angle)]))
            for p in mounting_ports:
                for z in [design['center_z']+design['bend_radius']+.5,fill_z-.5,design['top']-.4]:
                    samples.append(dict(region='riser',center_mm=[*p['top_center_mm'][:2],z],normal=[0,0,1]))
    frame,retainers,keys,records=frame_and_joints(cfg,inner,rim,mounting_ports,kind)
    magnet_slots=None
    if cfg.get('magnet_slots'):
        from .channel_magnet_slots import add_slots
        frame,magnet_slots=add_slots(cfg,frame,{kind+'_lumen_core':core,**retainers,**keys},inner,out)
    parts={kind+'_lumen_core':core,kind+'_bottomless_box_wall2mm':frame,**retainers,**keys}
    if list(parts)!=part_names(kind):raise ValueError('PART_ORDER_MISMATCH')
    assembly=out/'assembly';assembly.mkdir(parents=True,exist_ok=True)
    for name,m in parts.items():m.export(assembly/(name+'.stl'))
    import cadquery as cq
    cq.exporters.export(core_cad,str(assembly/(kind+'_lumen_core.step')))
    ref=kind+('_channel_without_positioning_collars_reference.stl' if kind=='Y' else '_channel_without_collars_reference.stl')
    channel.export(assembly/ref)
    if kind=='Y':
        source.apply_translation([0,0,heights['source_z_translation_mm']]);source.export(assembly/'Y_original_body_mm_reference.stl')
    for r in samples:r['equivalent_diameter_mm']=yu.local_section_diameter(channel,r['center_mm'],r['normal'])
    labels=np.ones(len(channel.faces),dtype=np.int32)
    for p in ports:
        center=np.asarray(p.get('center_mm',p.get('top_center_mm')));normal=np.asarray(p.get('outward_axis',[0,0,1]))
        cap=(abs((channel.triangles-center)@normal).max(axis=1)<1e-5)&(np.linalg.norm(channel.triangles_center-center,axis=1)<.81)
        if not cap.any():raise ValueError('MISSING_PORT_CAP')
        p['measured_cap_diameter_mm']=float(2*np.sqrt(channel.area_faces[cap].sum()/np.pi))
        p['minimum_cap_normal_dot']=float((channel.face_normals[cap]@normal).min());labels[cap]=p['label_id']
    poly=pv.wrap(channel);poly.cell_data['ModelFaceID']=labels;poly.field_data['coordinate_scale_to_m']=[.001]
    poly.field_data['ModelFaceNames']=['WALL']+[p['name'] for p in ports]
    poly.save(out/('labelled_y_D1p6_lumen_mm.vtp' if kind=='Y' else 'labelled_straight_lumen_mm.vtp'))
    qc=fixture_qc(cfg,parts,inner,fill_z,rim,records,kind)
    width_qc=y.pdms_width_qc(frame,inner,cfg,fill_z)
    measured_h=float(frame.bounds[1,2]-cfg['freeboard_mm']-frame.bounds[0,2])
    height_budget=y.dimension_budget(cfg,'height',measured_h)
    height_qc=dict(axis='Z',external_substrate_top_z_mm=0.,measured_internal_abs_base_mm=0.,
        measured_inner_height_mm=measured_h,allowance=height_budget,
        measured_from='Coplanar lower edges of walls on external flat substrate; no internal printed floor',
        passed=height_budget['passed'] and abs(measured_h-fill_z)<1e-5 and qc['gates']['no_internal_abs_bottom'])
    meshes={n:y.mesh_info(m) for n,m in parts.items()}
    manufacturing_height=casting_height_qc(cfg,parts)
    gates={**qc['gates'], 'single_closed_parts':all(m.is_volume and len(m.split())==1 for m in parts.values()),
        'casting_height_limit':manufacturing_height['passed'],
        'pdms_width_minimum':width_qc['passed'],'pdms_height_minimum':height_qc['passed'],
        'diameter_1p6_all_sampled_sections':all(abs(r['equivalent_diameter_mm']-1.6)<.003 for r in samples),
        'diameter_1p6_all_caps':all(abs(p['measured_cap_diameter_mm']-1.6)<.003 for p in ports),
        'correct_port_normals':all(p['minimum_cap_normal_dot']>.99999 for p in ports),
        'requested_floor':abs(core.bounds[0,2]-cfg['pdms_floor_mm'])<.003}
    if magnet_slots:gates['external_magnet_slots']=magnet_slots['qc']['passed']
    chamber_qc=None
    if kind=='Straight' and s.two_top_ports(cfg):
        from .straight_two_top_ports import channel_qc
        chamber_qc=channel_qc(design,channel,ports)
        gates['two_top_ports_no_horizontal_stub']=chamber_qc['passed']
    elif kind=='Straight':
        extension=trimesh.intersections.slice_mesh_plane(channel,[1,0,0],
            [design['junction_x']+design['radius']+cfg['tee_fillet_mm'],0,0],cap=True)
        fill=y.box([np.r_[inner[0],0],np.r_[inner[1],fill_z]])
        outside=max(0.,float(extension.volume)-s.volume_intersection(extension,fill))
        clearance=float(inner[1,0]-extension.bounds[1,0])
        gates['horizontal_extension_inside_pdms']=outside<1e-4 and clearance>=cfg['horizontal_end_clearance_mm']-1e-5
        chamber_qc=dict(horizontal_extension_outside_pdms_mm3=outside,horizontal_end_to_inner_wall_mm=clearance,
            manual_cut_reference=dict(plane_x_mm=design['horizontal_end']-cfg['manual_trim_overlap_mm'],
                pdms_to_remove_from_end_face_mm=clearance+cfg['manual_trim_overlap_mm']))
    else:gates['source_file_unchanged']=y.digest(cfg['source'])==cfg['source_sha256']
    fill=y.box([np.r_[inner[0],0],np.r_[inner[1],fill_z]])
    summary=dict(status='GEOMETRY_PASSED' if all(gates.values()) else 'GEOMETRY_FAILED',units='mm',
        variant=kind,channel_topology=cfg.get('channel_topology','three_ports'),
        fixture_mode=MODE,diameter_mm=1.6,wall_mm=cfg['wall_mm'],base_mm=0.,
        external_foot_mm=cfg['dovetail']['external_foot_mm'],heights=heights,design=design,
        cavity_xy_bounds_mm=inner,cavity_xy_size_mm=np.ptp(inner,axis=0),mold_outer_size_mm=frame.extents,
        pdms_fill_height_mm=fill_z,pdms_minimum_floor_mm=float(core.bounds[0,2]),
        pdms_net_fill_ml=(fill.volume-s.volume_intersection(core,fill))/1000,
        pdms_width_qc=width_qc,pdms_height_qc=height_qc,ports=ports,section_checks=samples,
        fixture_qc=qc,dovetail_joints=records,closed_chamber_qc=chamber_qc,
        manufacturing_height_qc=manufacturing_height,
        magnet_slots=magnet_slots,
        part_materials=cfg.get('part_materials'),
        source_sha256=cfg.get('source_sha256'),source_wall_preserved=False,global_geometric_scaling=False,
        production_part_names=list(parts),meshes=meshes,gates=gates,physical_fit_validated=False,
        physical_casting_validated=False,external_flat_substrate_required=True,side_gate_external_seal_required=False,
        assembly_sha256={p.name:y.digest(p) for p in sorted(assembly.iterdir())})
    if kind=='Y':
        summary['planes_sha256']=y.digest(cfg['planes'])
        summary['source_to_mm_scale']=1000.
    y.save(out/'effective_config.json',cfg);y.save(out/'geometry_qc.json',summary)
    if not all(gates.values()):raise ValueError('BOTTOMLESS_GEOMETRY_GATES_FAILED: '+str({k:v for k,v in gates.items() if not v}))
    return parts,summary


def arrange(parts,out):
    """Separated native objects; temporary supports on core and capture forks."""
    records=[];meshes=[];ni=nk=0
    for name,part in parts.items():
        T=np.eye(4)
        if 'core' in name:xy=(65.,65.)
        elif 'box' in name:xy=(150.,65.)
        elif 'retainer' in name:
            xy=(55.+55*ni,130.);ni+=1;T[:3,:3]=Rotation.from_euler('x',180,degrees=True).as_matrix()
        else:xy=(60.+45*nk,185.);nk+=1
        m=moved(part,T)
        delta=np.r_[xy,0.]-np.r_[m.bounds.mean(0)[:2],m.bounds[0,2]];m.apply_translation(delta);T[:3,3]=delta
        records.append(dict(name=name,support_enabled='core' in name or 'retainer' in name,
            casting_to_plate_transform=T,plate_to_casting_transform=np.linalg.inv(T),plate_bounds_mm=m.bounds))
        meshes.append(m)
    envelopes=[m.bounds[:,:2]+[[-7,-7],[7,7]] for m in meshes]
    for i,a in enumerate(envelopes):
        if not ((a[0]>=[20,25]).all() and (a[1]<=[245,245]).all()):raise ValueError('BOTTOMLESS_PLATE_BOUNDS')
        for b in envelopes[i+1:]:
            if not ((a[1]<b[0]).any() or (b[1]<a[0]).any()):raise ValueError('BOTTOMLESS_PLATE_OVERLAP')
    (out/'print_parts').mkdir(exist_ok=True)
    for m,r in zip(meshes,records):m.export(out/'print_parts'/(r['name']+'.stl'))
    y.save(out/'print_layout.json',dict(passed=True,envelope_margin_mm=7,records=records))
    return meshes,records


def render(parts,summary,out,only_views=None):
    (out/'QC').mkdir(exist_ok=True);kind=summary['variant'];inner=np.asarray(summary['cavity_xy_bounds_mm'])
    views=['assembly','top','bottom','side','core','exploded','dovetail_detail','dovetail_section']
    if summary.get('magnet_slots'):views+=['magnet_detail','magnet_section']
    if kind=='Straight' and summary.get('channel_topology')!='two_vertical_ports':views.append('internal_end')
    if only_views is not None:views=[v for v in views if v in only_views]
    for name in views:
        plot=pv.Plotter(off_screen=True,window_size=(1600,1100));plot.set_background('white')
        for key,part in parts.items():
            if name=='core' and 'core' not in key:continue
            m=part.copy()
            if name=='magnet_section':
                if 'box' not in key:continue
                center=summary['magnet_slots']['placements'][1]['pocket_floor_center_mm']
                m=trimesh.intersections.slice_mesh_plane(m,[1,0,0],center,cap=True)
            if name=='dovetail_section':
                if 'core' in key or 'key' in key:continue
                m=trimesh.intersections.slice_mesh_plane(m,[0,0,-1],[0,0,8],cap=True)
                if not len(m.faces):continue
            if name=='exploded':
                if 'key' in key:
                    r=next(r for r in summary['dovetail_joints'] if key.endswith('_'+r['port']))
                    m.apply_translation(-np.asarray(r['key_insertion_axis'])*15)
                elif 'retainer' in key or 'core' in key:m.apply_translation([0,0,16])
            color='#c45142' if 'core' in key else '#84a9bb' if 'box' in key else '#3f8065' if 'key' in key else '#d4a64c'
            plot.add_mesh(pv.wrap(m),color=color,opacity=.35 if 'box' in key and name not in ['bottom','dovetail_section','magnet_detail','magnet_section'] else 1.,smooth_shading='core' in key)
        if name in ['assembly','side','top']:
            plane=pv.Plane(center=[*inner.mean(0),summary['heights']['fill_z_mm']],direction=[0,0,1],
                i_size=np.ptp(inner[:,0]),j_size=np.ptp(inner[:,1]));plot.add_mesh(plane,color='#64c6d1',opacity=.12)
        variant_title=kind+(' (two TOP ports only)' if summary.get('channel_topology')=='two_vertical_ports' else '')
        materials='ABS core + PLA fixture' if summary.get('part_materials') else 'ABS'
        plot.add_text(f"{variant_title} | BOTTOMLESS + external dovetails + cross keys | P1S / 0.4 / {materials}\n"
            f"PDMS CAD height 9.60 mm | channel D1.60 mm | floor clearance 1.50 mm\n"
            f"Fixture top {summary['manufacturing_height_qc']['fixture_top_z_mm']:.2f} mm | "
            f"Core top {summary['manufacturing_height_qc']['core_top_z_mm']:.2f} mm | wall 2.00 mm | NO internal base",
            font_size=13,color='#263238')
        if name in ['side','dovetail_detail']:
            gap=summary['fixture_qc']['retainer_rim_clearance_qc']['requested_gap_mm']
            plot.add_text(f"Clear air gap: arm underside to box rim = {gap:.2f} mm",position='lower_right',font_size=13,color='#263238')
        plot.add_axes();plot.show_bounds(xtitle='X (mm)',ytitle='Y (mm)',ztitle='Z (mm)',fmt='%.2f',font_size=12,color='#33434b',location='outer')
        if name=='top':plot.view_xy()
        elif name=='bottom':plot.view_xy(negative=True)
        elif name=='side':plot.view_xz()
        else:plot.view_isometric()
        plot.camera.zoom(.85)
        if name=='dovetail_detail':
            r=summary['dovetail_joints'][0];T=np.asarray(r['local_to_casting']);A=r['wall_outer_local_x_mm']
            focus=(T@np.array([A+3,0,7,1]))[:3];eye=focus+T[:3,:3]@np.array([22,-25,17])
            plot.camera_position=(eye,focus,[0,0,1]);plot.camera.parallel_projection=True;plot.camera.parallel_scale=12
        if name=='dovetail_section':
            r=summary['dovetail_joints'][0];T=np.asarray(r['local_to_casting']);A=r['wall_outer_local_x_mm']
            focus=(T@np.array([A+3,0,8,1]))[:3]
            plot.camera_position=(focus+[0,0,30],focus,T[:3,0]);plot.camera.parallel_projection=True;plot.camera.parallel_scale=10
        if name=='internal_end':
            d=summary['design'];focus=[d['inner'][1][0],0,d['center_z']+2]
            plot.camera_position=([focus[0]+15,-18,focus[2]+13],focus,[0,0,1])
        if name in ['magnet_detail','magnet_section']:
            center=np.array(summary['magnet_slots']['placements'][1]['pocket_floor_center_mm'])
            eye=[-18,23,22] if name=='magnet_section' else [18,23,22]
            plot.camera_position=(center+eye,center,[0,0,1]);plot.camera.parallel_projection=True;plot.camera.parallel_scale=9
            plot.add_text('Source-sized blind pocket: D6.40 x 2.00 mm | floor 1.00 mm\nTwo opposite exterior slots; continuous bottomless chamber walls',
                position='lower_left',font_size=13,color='#263238')
        plot.screenshot(out/'QC'/(name+'.png'));plot.close()


def toolpath_qc(objects, model_segments, core_raft_lift):
    """Check each exported fixture against real positive model extrusion."""
    from .bambu_casting_mold import sample_segments
    from scipy.spatial import cKDTree
    mid=model_segments.mean(1);rows=[]
    for obj in objects:
        name=obj['metadata']['name'];m=obj['mesh'];lo,hi=m.bounds[:,:2]
        segments=model_segments[((mid[:,:2]>=lo-.5)&(mid[:,:2]<=hi+.5)).all(1)]
        lift=core_raft_lift if 'core' in name else 0.
        max_distance=None
        if len(segments):
            distances=cKDTree(sample_segments(segments,.15)).query(m.triangles_center+[0,0,lift])[0]
            max_distance=float(distances.max())
        rows.append(dict(part=name,positive_model_segments=len(segments),raft_lift_mm=lift,
            surface_max_distance_to_model_toolpath_mm=max_distance,
            passed=max_distance is not None and max_distance<.65))
    return dict(parts=rows,passed=all(r['passed'] for r in rows),
        scope='Digital model/toolpath proximity, not a physical strength or sealing test')
