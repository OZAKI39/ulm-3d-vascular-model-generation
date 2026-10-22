"""New endpoints and side-wall routes, reusing the accepted real-cap pipeline."""
from __future__ import annotations
import copy
import json
import math
import hashlib
import shutil
from pathlib import Path
import numpy as np
import trimesh
import yaml
from . import sacrificial_fixture as f
from . import abs_casting_mold as mold
from . import casting_mold_qc as qc
from .port_attachment_alignment import align_known_cap,clean_boolean_attachment_slivers
from .all_port_attachment_run import staged_union
from .surface_continuity_qc import port_continuity
from .sacrificial_fixture_review import save_json
from .support_aware_roi import ROOT,read_csv,table
from .swc_export import read_source
from .topbrain_qc import sha256


def candidate_configuration(directory,cfg):
    """Apply traceable routing feedback tied to this exact vascular surface."""
    path=Path(directory)/'port_route_review.json'
    if not path.exists():return cfg
    review=json.loads(path.read_text())
    if review['vascular_sha256']!=sha256(Path(directory)/'vascular_core.stl'):
        raise ValueError('PORT_ROUTE_REVIEW_SOURCE_MISMATCH')
    for faces in review['allowed_faces_by_endpoint'].values():
        if not faces or not set(faces)<={'-X','+X','-Y','+Y'}:raise ValueError('INVALID_REVIEW_SIDEWALL')
    cfg=copy.deepcopy(cfg);cfg['route_review']=review
    return cfg


def ensure_current_design(directory,transform,cfg):
    """Retire geometry AND dependent slices when routing/box inputs change."""
    directory=Path(directory);cfg=candidate_configuration(directory,cfg);path=directory/'manufacturing_inputs.json'
    inputs=dict(algorithm_revision=2,vascular_sha256=sha256(directory/'vascular_core.stl'),
        transform=np.asarray(transform).tolist(),configuration=cfg)
    signature=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    previous=json.loads(path.read_text()) if path.exists() else {}
    if previous.get('signature')==signature:return
    names=['manufacturing','orientations','evaluation.json','evaluation_error.log','manufacturing_inputs.json']
    present=[directory/n for n in names if (directory/n).exists()]
    if present:
        history=directory/'history';history.mkdir(exist_ok=True)
        archive=history/f'revision_{len(list(history.iterdir()))+1:02}';archive.mkdir()
        relocated={str(p):str(archive/p.name) for p in present}
        for p in present:shutil.move(str(p),str(archive/p.name))
        save_json(archive/'relocation.json',dict(reason='ROUTING_OR_BOX_INPUTS_CHANGED',paths=relocated,
            note='Historical JSON keeps original absolute references; use this relocation map.'))
    save_json(path,dict(signature=signature,inputs=inputs))


def configuration():
    stage=yaml.safe_load((ROOT/'config/support_aware_roi_60x50.yaml').read_text())
    cfg=yaml.safe_load((ROOT/stage['baseline_config']).read_text())
    route=yaml.safe_load((ROOT/'config/sacrificial_box_BG001.yaml').read_text())
    for key in ('ports','routing','attachment_alignment','surface_continuity'):cfg[key]=copy.deepcopy(route[key])
    cfg['box'].update(route['box']);cfg['box']['policy']='BOX_DIMENSIONS_REUSED'
    cfg['orientation'].pop('candidate_ids',None)
    cfg['evaluation_printer']=stage['evaluation_printer'];cfg['printer'].update(stage['printer'])
    cfg['target']=stage['target'];cfg['box'].update(stage['box']);cfg['baseline_config']=stage['baseline_config']
    cfg['ports'].update(stage.get('ports',{}))
    if cfg['target']['scale_factor']!=1. or min(stage['clearance_levels_mm'])<1.6:raise ValueError('MANUFACTURING_HARD_GATE_CONFIG_CHANGED')
    cfg['production_name']='BG001_RMCA_MeVO_ROI_60x50_support_optimized'
    # A failed preferred 3 mm clearance may be relaxed locally, never below1.6.
    cfg['clearance_levels_mm']=stage['clearance_levels_mm']
    return cfg


def load_companion(directory,transform,cfg,align=True):
    import pyvista as pv
    from .print_orientation import cap_native_ports
    directory=Path(directory);result=json.loads((directory/'model_result.json').read_text())
    graph=read_source(Path(result['VascularMD']['outputs']['swc'])).graph
    endpoints,branches=f.decode_endpoints(graph,read_csv(directory/'compensated_mapping.csv'),read_csv(directory/'compensated_edges.csv'),transform,cfg)
    precise=directory/'vascular_core.vtp'
    if not precise.exists():
        closed,_=cap_native_ports(pv.read(result['VascularMD']['outputs']['surface_vtk']))
        closed.transform(transform,inplace=True);closed.save(precise)
    surface=pv.read(precise)
    mesh=trimesh.Trimesh(surface.points,surface.faces.reshape(-1,4)[:,1:],process=True)
    delivered=trimesh.load_mesh(directory/'vascular_core.stl')
    if len(delivered.vertices)!=len(mesh.vertices):raise ValueError('VTP_STL_TOPOLOGY_MISMATCH')
    from scipy.spatial import cKDTree
    if cKDTree(delivered.vertices).query(mesh.vertices)[0].max()>1e-4:raise ValueError('VTP_STL_GEOMETRY_MISMATCH')
    if align:
        for e in endpoints:
            coord=graph.nodes[e.swc_id]['coords']
            caps=[c for c in result['cap_qc']['port_caps'] if np.linalg.norm(coord[:3]-c['center_mm'])<1e-4 and abs(coord[3]-c['median_radius_mm'])<1e-4]
            if len(caps)!=1 or caps[0]['vertices']!=24:raise ValueError('CAP_PROVENANCE_UNRESOLVED_'+e.endpoint_id)
            align_cap_local(e,mesh,caps[0],cfg)
    return dict(graph=graph,branches=branches,endpoints=endpoints,transform=transform,mesh=mesh)


def align_cap_local(endpoint,source,cap,cfg):
    """The accepted algorithm in a translated frame, before STL quantization.

    No radius, topology, frame rotation or tolerance changes. Float32 Boolean
    kernels lose several micrometres at world coordinates around 150 mm;
    carrying the native double-precision cap and translating the calculation
    avoids mistaking that loss for a protruding buried collar.
    """
    local=copy.deepcopy(endpoint);origin=endpoint.position.copy();local.position-=origin
    mesh=source.copy();mesh.apply_translation(-origin)
    record=align_known_cap(local,mesh,cap,cfg)
    endpoint.position=local.position+origin;endpoint.tangent=local.tangent;endpoint.cap_ring=local.cap_ring+origin
    endpoint.attachment_collar=local.attachment_collar.translate(tuple(origin))
    for key in ('original_ring_print_mm','aligned_ring_print_mm'):record[key]=(np.asarray(record[key])+origin).tolist()
    record.update(numerical_frame='CAP_CENTER_TRANSLATION_ONLY',calculation_origin_mm=origin.tolist(),
        source_precision='NATIVE_VASCULARMD_DOUBLE_PRECISION_VTP',tolerances_changed=False)
    endpoint.attachment_qc=record
    return record


def new_box(bounds,cfg):
    fitted=f.compute_box(bounds,cfg)
    # Only the stable five-wall constructor is reused. These bounds were just
    # derived from this NEW vascular mesh, not read from a legacy mold.
    record=dict(inner_bounds_mm=fitted.inner,outer_bounds_mm=fitted.outer,
        wall_thickness_mm=cfg['box']['wall_thickness_mm'],bottom_thickness_mm=cfg['box']['bottom_thickness_mm'])
    box=mold.build_intact_box(dict(old_box=record),cfg)
    box['dimensions_policy']='NEW_VASCULAR_SURFACE_BOUNDS_PLUS_MARGINS'
    return box,fitted


def adaptive_box(core,ports,seed,cfg):
    """Fit unconnected faces to the complete new vascular/port envelope.

    Assigned wall planes are fixed by the new routes. Other planes approach
    actual core bounds, with explicit edge material retained around every port.
    No vertex, vessel radius, or port trajectory is changed by this operation.
    """
    if not cfg['box'].get('adaptive_unconnected_faces'):return seed
    inner=seed['inner'].copy();assigned={p['face'] for p in ports};changes=[]
    gap=cfg['box']['unconnected_face_gap_mm'];edge=cfg['ports']['min_wall_edge_clearance_mm']
    bounds=core.bounds
    for axis,side,name in [(0,0,'-X'),(0,1,'+X'),(1,0,'-Y'),(1,1,'+Y'),(2,0,'BOTTOM'),(2,1,'TOP_OPENING')]:
        previous=float(inner[side,axis]);clearance=cfg['box'].get('top_headroom_mm',gap) if name=='TOP_OPENING' else gap
        if name in assigned:
            reason='ROUTED_PORT_WALL_RETAINED';proposed=previous
        else:
            sign=-1 if side==0 else 1;proposed=float(bounds[side,axis]+sign*clearance)
            # All ports on OTHER axes need solid edge material. Along the
            # assigned axis, this is an intentional wall crossing instead.
            relevant=[p['target'][axis] for p in ports if f.face_axis(p['face'])[0]!=axis]
            unconstrained=proposed
            if relevant:proposed=min(proposed,min(relevant)-edge) if side==0 else max(proposed,max(relevant)+edge)
            reason='ACTUAL_CORE_ENVELOPE_PLUS_GAP' if abs(proposed-unconstrained)<1e-9 else 'PORT_EDGE_MATERIAL_LIMIT'
            limited=max(previous,proposed) if side==0 else min(previous,proposed)
            if abs(limited-proposed)>1e-9:reason='SEED_ENVELOPE_LIMIT'
            proposed=limited
            inner[side,axis]=proposed
        actual=float(bounds[0,axis]-proposed) if side==0 else float(proposed-bounds[1,axis])
        changes.append(dict(face=name,port_connected=name in assigned,original_plane_mm=previous,adaptive_plane_mm=proposed,
            movement_toward_core_mm=abs(proposed-previous),requested_gap_mm=clearance if name not in assigned else None,
            core_envelope_gap_mm=actual if name not in assigned else None,reason=reason))
    outer=inner.copy();wall=cfg['box']['wall_thickness_mm'];bottom=cfg['box']['bottom_thickness_mm']
    outer[0]-=[wall,wall,bottom];outer[1]+=[wall,wall,0.]
    result=box_from_bounds(inner,outer,cfg);result['dimensions_policy']='ADAPTIVE_UNCONNECTED_FACES_FROM_NEW_CORE_AND_ROUTES'
    from .sacrificial_print_frame import box_mesh
    for row in changes:
        if not row['port_connected'] and row['face']!='TOP_OPENING':
            row['fcl_surface_distance_mm']=float(f.surface_distance(core,box_mesh(result['members'][row['face']])))
            row['absolute_clearance_pass']=row['fcl_surface_distance_mm']>=1.6
    result['adaptive_faces']=changes;result['seed_inner_bounds_mm']=seed['inner'];result['seed_outer_bounds_mm']=seed['outer']
    return result


def box_from_bounds(inner,outer,cfg):
    return mold.build_intact_box(dict(old_box=dict(inner_bounds_mm=inner,outer_bounds_mm=outer,
        wall_thickness_mm=cfg['box']['wall_thickness_mm'],bottom_thickness_mm=cfg['box']['bottom_thickness_mm'])),cfg)


def interface_qc(core,source,endpoints,routes,cfg):
    lost=trimesh.boolean.difference([source,core],engine='manifold');loss=abs(float(lost.volume)) if len(lost.faces) else 0.
    rows=[]
    for e,route in zip(endpoints,routes):
        section=f.attachment_section_qc(route,cfg);continuity=port_continuity(core,e,cfg)
        area_error=abs(section['measured_area_mm2']/section['expected_area_mm2']-1)
        feature=continuity['features'][1]
        rows.append(dict(port_id=e.endpoint_id,actual_cap=e.attachment_qc,section=section,continuity=continuity,
            feature_edge_20_length_mm=feature['total_length_mm'],new_closed_sharp_ring=feature['circumferential_ring_count'],
            section_area_change_fraction=area_error,vascular_material_loss_mm3=loss,
            passed=bool(feature['total_length_mm']==0 and feature['circumferential_ring_count']==0 and area_error<.01 and section['passed'] and loss<=cfg['geometry']['volume_tolerance_mm3'])))
    return rows


def delivered_interface_qc(union,companion,interfaces,cfg):
    """Recheck delivered final STL rather than relying on the intermediate union."""
    from .surface_continuity_qc import cross_section
    rows=[]
    for endpoint in companion['endpoints']:
        reference=next(r for r in interfaces if r['port_id']==endpoint.endpoint_id)
        detail=port_continuity(union,endpoint,cfg);feature=detail['features'][1]
        section=cross_section(union,endpoint.position+.02*endpoint.tangent,endpoint.tangent,endpoint.radius)
        error=abs(section['area_mm2']/reference['section']['expected_area_mm2']-1)
        rows.append(dict(port_id=endpoint.endpoint_id,feature_edge_20_length_mm=feature['total_length_mm'],
            new_closed_sharp_ring=feature['circumferential_ring_count'],section_area_change_fraction=error,
            passed=bool(feature['total_length_mm']==0 and feature['circumferential_ring_count']==0 and error<.01)))
    return rows


def select_new_layout(options,cfg,minimum):
    """Shortest feasible new layout across all four legal side walls.

    The old four-port fixture restricted outlets to at most two adjacent walls.
    That historical preference is not an anatomical or manufacturing constraint
    on a freshly selected subtree. Real FCL collisions and spacing remain hard.
    """
    tolerance=2*cfg['geometry']['mesh_tolerance_mm']
    options=[[r for r in rows if r.vessel_clearance>=minimum+tolerance] for rows in options]
    review=cfg.get('route_review',{}).get('allowed_faces_by_endpoint',{})
    if review:
        options=[[r for r in rows if r.face in review.get(r.endpoint.endpoint_id,('-X','+X','-Y','+Y'))] for rows in options]
    if not options or any(not rows for rows in options):return [],dict(status='NO_ENDPOINT_ROUTE')
    ordered=[sorted(rows,key=lambda r:(r.length,r.turn,r.option_id)) for rows in options]
    minima=[min(r.length for r in rows) for rows in ordered]
    preferred=cfg['ports'].get('preferred_maximum_extension_mm',35.)
    cache={};best=[];best_rank=(math.inf,math.inf,math.inf);tol=2*cfg['geometry']['mesh_tolerance_mm']
    def compatible(a,b):
        if a.face==b.face and np.linalg.norm(a.target-b.target)<f.required_spacing(a.endpoint,b.endpoint,cfg):return False
        key=tuple(sorted((a.option_id,b.option_id)))
        if key not in cache:cache[key]=f.surface_distance(a.mesh,b.mesh)
        return cache[key]>=minimum+tol
    def visit(routes,length):
        nonlocal best,best_rank
        j=len(routes)
        lengths=[r.length for r in routes]+minima[j:]
        bound=(sum(v>preferred for v in lengths),max(lengths),sum(lengths))
        if bound>=best_rank:return
        if j==len(ordered):
            best=routes[:];best_rank=bound
            return
        for route in ordered[j]:
            if all(compatible(route,r) for r in routes):visit(routes+[route],length+route.length)
    visit([],0.)
    for route in best:
        others=[r for r in best if r is not route]
        route.nearest_other_port=min(f.surface_distance(route.mesh,r.mesh) for r in others)
        route.wall_spacing=min((np.linalg.norm(route.target-r.target) for r in others if r.face==route.face),default=math.inf)
    return best,dict(status='PASS' if best else 'NO_COLLISION_FREE_LAYOUT',method='FOUR_SIDE_WALLS_SHORTEST_FEASIBLE',
        minimum_port_gap_mm=minimum,total_length_mm=sum(r.length for r in best) if best else None,
        maximum_extension_mm=max((r.length for r in best),default=None),preferred_maximum_extension_mm=preferred,
        rank_order=['count over preferred length','maximum extension','total extension'],evaluated_port_pairs=len(cache))


def design(directory,transform,cfg,logger):
    directory=Path(directory);cfg=candidate_configuration(directory,cfg);ensure_current_design(directory,transform,cfg)
    out=directory/'manufacturing';out.mkdir(exist_ok=True)
    if (out/'geometry.json').exists():return load_design(directory,transform,cfg)
    companion=load_companion(directory,transform,cfg);source=companion['mesh'];box,bounds=new_box(source.bounds,cfg)
    collision=f.VesselCollision(source,companion['graph'],companion['branches'],transform,companion['endpoints'],cfg)
    options=[];records=[];cap_records={};failures=[]
    for e in companion['endpoints']:
        cap_records[e.endpoint_id]=e.attachment_qc;feasible=[]
        for level in cfg['clearance_levels_mm']:
            feasible,rows=f.generate_options(e,bounds,cfg,collision,clearance_override=level,stage=f'GAP_{level}')
            records.extend(rows)
            if feasible:break
        options.append(feasible)
        logger.info('%s %s feasible routes=%s',directory.name,e.endpoint_id,len(feasible))
        f.write_csv(out/'port_route_candidates.csv',records)
    routes=[];layout={}
    for level in cfg['clearance_levels_mm']:
        routes,layout=select_new_layout(options,cfg,level)
        if routes:break
    if not routes:
        # A few individually feasible targets can all conflict as a set. Expand
        # only the NEW wall targets; retain the real caps, radii and bend gates.
        expanded=copy.deepcopy(cfg)
        expanded['routing'].update(target_lateral_offsets_mm=[0,-10,10,-20,20],
            target_vertical_offsets_mm=[0,-4,4,-8,8,-12,12],options_per_face=24)
        options=[]
        for endpoint in companion['endpoints']:
            feasible,rows=f.generate_options(endpoint,bounds,expanded,collision,clearance_override=1.6,stage='EXPANDED_LEGAL_SIDEWALLS')
            records.extend(rows);options.append(feasible)
            logger.info('%s expanded %s routes=%s',directory.name,endpoint.endpoint_id,len(feasible))
            f.write_csv(out/'port_route_candidates.csv',records)
        for level in cfg['clearance_levels_mm']:
            routes,layout=select_new_layout(options,expanded,level)
            if routes:break
    save_json(out/'routing_audit.json',dict(layout=layout,actual_caps=cap_records,start_clearances=collision.start_checks,
        expanded_search_performed=any(r['search_stage']=='EXPANDED_LEGAL_SIDEWALLS' for r in records)))
    if not routes:raise ValueError('NEW_ENDPOINT_SIDEWALL_LAYOUT_FAILED: '+repr(layout))
    # Store route meshes and endpoint metadata before union for recoverable QC.
    for route in routes:f.export_checked_stl(route.mesh,out/(route.endpoint.endpoint_id+'.stl'),cfg)
    core,stages=staged_union(source,routes,cfg)
    core,cleanup=clean_boolean_attachment_slivers(core,source,companion['endpoints'],cfg)
    f.export_checked_stl(core,out/'core_with_ports.stl',cfg);core=trimesh.load_mesh(out/'core_with_ports.stl')
    interfaces=interface_qc(core,source,companion['endpoints'],routes,cfg)
    if not all(r['passed'] for r in interfaces):failures.append('NEW_PORT_INTERFACE_QC_FAILED')
    ports=[]
    for route in routes:
        e=route.endpoint
        ports.append(dict(port_id=e.endpoint_id,radius_mm=e.radius,face=route.face,target=route.target,end=route.end,
            direction=route.direction,parent_branch=e.branch_id,swc_id=e.swc_id,reference_mesh=route.mesh,
            accepted_attachment=e.attachment_qc,length_mm=route.length,vessel_clearance_mm=route.vessel_clearance,
            outside_stub_mm=cfg['ports']['outside_stub_mm'],preferred_maximum_extension_mm=cfg['ports']['preferred_maximum_extension_mm'],
            moderate_length_target_met=route.length<=cfg['ports']['preferred_maximum_extension_mm'],
            minimum_port_clearance_mm=route.nearest_other_port,local_clearance_relaxed=route.vessel_clearance<3 or route.nearest_other_port<3))
    inputs=dict(core=core,central=source,ports=ports,companion=companion)
    box=adaptive_box(core,ports,box,cfg)
    if any(not r.get('absolute_clearance_pass',True) for r in box.get('adaptive_faces',[])):failures.append('ADAPTIVE_FREE_WALL_TOO_CLOSE')
    boxqc=qc.intact_box_qc(box,cfg);overlaps,contacts=qc.port_overlaps(inputs,box,cfg)
    if not all(boxqc[k] for k in ('four_side_walls_closed','bottom_closed','top_fully_open','no_assembly_clearance_holes')):failures.append('FIVE_WALL_BOX_STRUCTURE_FAILED')
    if any(min(p['vessel_clearance_mm'],p['minimum_port_clearance_mm'])<1.6 for p in ports):failures.append('PORT_CLEARANCE_BELOW_ABSOLUTE_MINIMUM')
    if any(r['status']!='PASS' for r in overlaps):failures.append('PORT_WALL_OVERLAP_FAILED')
    union,boolean=mold.combine(core,box['mesh'],cfg)
    f.export_checked_stl(box['mesh'],out/'intact_open_top_box.stl',cfg)
    f.export_checked_stl(union,out/'final_abs_casting_mold.stl',cfg);union=trimesh.load_mesh(out/'final_abs_casting_mold.stl')
    audit=qc.union_audit(source,box['mesh'],union,cfg)
    # Full core (including all new ports) must also be retained by the box union.
    core_audit=qc.union_audit(core,box['mesh'],union,cfg)
    failures.extend(audit['failures']+core_audit['failures'])
    final_interfaces=delivered_interface_qc(union,companion,interfaces,cfg)
    if not all(r['passed'] for r in final_interfaces):failures.append('FINAL_UNION_PORT_INTERFACE_QC_FAILED')
    partition=qc.partition_vascular(union,inputs,box,cfg)
    pdms=qc.ligament_qc(union,inputs,partition,cfg);walls=qc.wall_clearance(union,partition,inputs,box,cfg)
    if pdms['minimum_ligament_mm']<1.6:failures.append('PDMS_LIGAMENT_TOO_THIN')
    if min(r['distance_mm'] for r in walls if not r['vascular_branch'].startswith(('I','O')))<1.6:failures.append('CORE_WALL_TOO_CLOSE')
    record=dict(failures=sorted(set(failures)),route_review=cfg.get('route_review'),interfaces=interfaces,final_interfaces=final_interfaces,ports=[{k:v for k,v in p.items() if k!='reference_mesh'} for p in ports],
        box={**boxqc,'inner_bounds_mm':box['inner'],'outer_bounds_mm':box['outer'],
            'adaptive_faces':box.get('adaptive_faces',[]),'seed_inner_bounds_mm':box.get('seed_inner_bounds_mm'),
            'seed_outer_bounds_mm':box.get('seed_outer_bounds_mm')},union=core_audit,
        vascular_only_preservation=audit,staged_union=stages,cleanup=cleanup,pdms=pdms,wall_clearance=walls,port_overlaps=overlaps)
    save_json(out/'geometry.json',record)
    export_geometry_tables(out,record)
    return dict(inputs=inputs,box=box,union=union,partition=partition,record=record)


def export_geometry_tables(out,record):
    # Connected faces intentionally have no free-wall distance measurement.
    # Use the union of columns rather than assuming the first row has all keys.
    for name,rows in [('port_interfaces',record['interfaces']),('final_port_interfaces',record['final_interfaces']),
        ('port_design',record['ports']),('port_wall_union',record['port_overlaps']),
        ('pdms_ligament_top10',record['pdms']['top10']),('wall_clearance',record['wall_clearance']),
        ('adaptive_box_faces',record['box'].get('adaptive_faces',[]))]:
        if rows:table(Path(out)/(name+'.csv'),rows)


def load_design(directory,transform,cfg):
    directory=Path(directory);out=directory/'manufacturing';record=json.loads((out/'geometry.json').read_text())
    companion=load_companion(directory,transform,cfg,align=False);source=companion['mesh']
    box=box_from_bounds(np.array(record['box']['inner_bounds_mm']),np.array(record['box']['outer_bounds_mm']),cfg)
    box['dimensions_policy']=record['box']['dimensions_policy']
    ports=[]
    for saved in record['ports']:
        p=copy.deepcopy(saved)
        for k in ('target','end','direction'):p[k]=np.array(p[k])
        p['reference_mesh']=trimesh.load_mesh(out/(p['port_id']+'.stl'));ports.append(p)
    inputs=dict(core=trimesh.load_mesh(out/'core_with_ports.stl'),central=source,ports=ports,companion=companion)
    union=trimesh.load_mesh(out/'final_abs_casting_mold.stl');partition=qc.partition_vascular(union,inputs,box,cfg)
    if 'final_interfaces' not in record:
        aligned=load_companion(directory,transform,cfg,align=True)
        record['final_interfaces']=delivered_interface_qc(union,aligned,record['interfaces'],cfg)
        if not all(r['passed'] for r in record['final_interfaces']):record['failures'].append('FINAL_UNION_PORT_INTERFACE_QC_FAILED')
        save_json(out/'geometry.json',record);f.write_csv(out/'final_port_interfaces.csv',record['final_interfaces'])
    export_geometry_tables(out,record)
    return dict(inputs=inputs,box=box,union=union,partition=partition,record=record)
