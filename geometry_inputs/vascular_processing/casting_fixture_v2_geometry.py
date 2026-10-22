"""Closed-wall V2 fixture. Only rigid placement of the immutable ABS mesh."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import cadquery as cq
import trimesh
from scipy.spatial.transform import Rotation
from .support_inputs import transform, swc_branches
from .casting_fixture_geometry import rigid, bounds_box, shape_mesh, mesh_qc


def port_transform(p):
    n=np.array(p['axis']); u=np.cross(n,[0.,0.,1.]); m=np.eye(4)
    m[:3,:3]=np.column_stack((u,n,[0,0,1.]));m[:3,3]=p['trim_origin'];return rigid(m)


def world(shape,p):
    m=port_transform(p); angle=np.degrees(np.arctan2(m[1,0],m[0,0]))
    return shape.rotate((0,0,0),(0,0,1),float(angle)).translate(tuple(m[:3,3]))


def straight_port(raw,cfg):
    n=np.asarray(raw['direction'],float);n/=np.linalg.norm(n);a=np.array(raw['target']);b=np.array(raw['end']);L=float((b-a)@n);k=raw['port_id']
    if abs(n[2])>1e-8 or np.linalg.norm(b-a-L*n)>1e-5:
        raise ValueError('PORT_STRAIGHT_FIXTURE_ZONE_CURVED_'+k)
    # A midpoint cut must leave room for two thin cradle webs, axial float and
    # the conservative transition/end margin; no anatomy repair is attempted.
    if L/2 < 2.35+.33+cfg['secondary_axial_float_mm']+cfg['straight_zone_margin_mm']:
        raise ValueError('PORT_STRAIGHT_FIXTURE_ZONE_TOO_SHORT_'+k)
    return dict(port_id=k,face=raw['face'],axis=n,straight_start=a,end=b,straight_length_mm=L,
                radius_mm=float(raw['radius_mm']),diameter_mm=2*raw['radius_mm'],attachment=raw['accepted_attachment']['calculation_origin_mm'],
                parent_branch=raw['parent_branch'],swc_id=raw['swc_id'],total_route_length_mm=raw['length_mm'])


def layout(inputs,cfg,line):
    if cfg['scale_factor']!=1:raise ValueError('SCALING_FORBIDDEN')
    original=trimesh.load_mesh(inputs['core']);central0=trimesh.load_mesh(inputs['vascular'])
    raw=json.loads(Path(inputs['geometry']).read_text())['ports'];ports0=[straight_port(p,cfg) for p in raw]
    placements=[]
    for yaw in cfg['yaw_candidates_degrees']:
        r=Rotation.from_euler('z',yaw,degrees=True).as_matrix();points=original.vertices@r.T;central=central0.vertices@r.T
        for setback in cfg['setback_candidates_mm']:
            low=points.min(0);high=points.max(0)
            for ax in (0,1):
                for side,sgn in ((0,-1),(1,1)):
                    group=[p for p in ports0 if (r@p['axis'])[ax]*sgn>.99]
                    margin=setback if group else cfg['unported_wall_margin_mm']
                    (low if side==0 else high)[ax]+=sgn*margin
                    if side==0:low[ax]=min(low[ax],central[:,ax].min()-cfg['central_wall_margin_mm'])
                    else:high[ax]=max(high[ax],central[:,ax].max()+cfg['central_wall_margin_mm'])
            low[2]=points[:,2].min()-cfg['bottom_clearance_mm'];high[2]=points[:,2].max()+cfg['top_clearance_mm']
            spans=[];tips=[]
            for p in ports0:
                n=r@p['axis'];ax=int(np.argmax(abs(n)));bound=high[ax] if n[ax]>0 else low[ax];end=r@p['end'];tip=abs(bound-end[ax]);tips.append(tip);spans.append(tip/abs(n[ax])+p['straight_length_mm']/4)
            placements.append(dict(yaw_degrees=yaw,setback_candidate_mm=setback,low=low,high=high,rotation=r,max_span_mm=max(spans),sum_span_mm=sum(spans),minimum_tip_setback_mm=min(tips)))
    # Explicit lexicographic objective: shortest worst span, then total span,
    # then least rotation. No opaque multi-objective score.
    chosen=min(placements,key=lambda x:(round(x['max_span_mm'],6),round(x['sum_span_mm'],6),abs(x['yaw_degrees'])))
    r=chosen['rotation'];lo=chosen['low'];hi=chosen['high'];m=np.eye(4);m[:3,:3]=r;m[:3,3]=[-(lo[0]+hi[0])/2,-(lo[1]+hi[1])/2,-lo[2]];rigid(m)
    core=original.copy();core.apply_transform(m);central=central0.copy();central.apply_transform(m);lo=lo+m[:3,3];hi=hi+m[:3,3]
    old=rigid(json.loads(Path(inputs['casting_transform']).read_text())['transform_4x4']);_,top=swc_branches(inputs['smoothed_swc'],m@old,inputs['branch_mapping']);bif=top['bifurcations']
    master=max(ports0,key=lambda p:(round(p['straight_length_mm'],4),p['radius_mm']))['port_id'];ports=[];source_ports={}
    for p0 in ports0:
        p=dict(p0);p['axis']=r@p0['axis']
        for k in ['straight_start','end','attachment']:p[k]=transform([p[k]],m)[0]
        n=p['axis'];ax=int(np.argmax(abs(n)));bound=hi[ax] if n[ax]>0 else lo[ax]
        p['trim_origin']=(p['straight_start']+p['end'])/2;p['trim_normal']=n.copy();p['axis_plane_angle_degrees']=90.
        p['wall_distance_from_trim_mm']=float((bound-p['trim_origin'][ax])/n[ax]);p['extension_wall_setback_mm']=abs(bound-p['end'][ax]);p['role']='MASTER' if p['port_id']==master else 'SECONDARY'
        p['cradle_clearance_mm']=line*cfg['selected_master_clearance_lines' if p['role']=='MASTER' else 'selected_secondary_clearance_lines'];p['anti_float_gap_mm']=line
        p['axial_float_mm']=cfg['secondary_axial_float_mm'] if p['role']=='SECONDARY' else 0.;p['usable_straight_length_after_trim_mm']=p['straight_length_mm']/2
        p['distance_to_bifurcation_mm']=float(np.linalg.norm(bif-p['trim_origin'],axis=1).min()) if len(bif) else None
        p['distance_to_cap_attachment_mm']=float(np.linalg.norm(p['trim_origin']-p['attachment']));p['distance_to_curvature_mm']=p['straight_length_mm']/2
        p['approved_contact_start']=p['straight_start']+n*cfg['straight_zone_margin_mm'];p['approved_contact_end']=p['end']-n*cfg['straight_zone_margin_mm']
        mp=trimesh.load_mesh(Path(inputs['candidate'])/'manufacturing'/f"{p['port_id']}.stl");mp.apply_transform(m);source_ports[p['port_id']]=mp
        sections=[]
        for f in [.25,.5,.75]:
            pt=p['straight_start']+f*p['straight_length_mm']*n;s=mp.section(plane_origin=pt,plane_normal=n)
            if s is None:raise ValueError('PORT_STRAIGHT_SECTION_MISSING_'+p['port_id'])
            rr=np.linalg.norm(s.vertices-pt,axis=1);ok=len(s.discrete)==1 and rr.max()<=p['radius_mm']+5e-5 and rr.min()>=p['radius_mm']*np.cos(np.pi/24)-5e-5
            if not ok:raise ValueError('PORT_STRAIGHT_SECTION_NOT_CONFIRMED_'+p['port_id'])
            sections.append(dict(fraction=f,radius_min_mm=float(rr.min()),radius_max_mm=float(rr.max()),passed=bool(ok)))
        p['straight_sections']=sections;ports.append(p)
    return dict(core=core,central=central,source_ports=source_ports,ports=ports,master_port=master,low=lo,high=hi,transform=m,placements=placements,selected_placement=chosen,bifurcations=bif,source_qc=mesh_qc(original))


def fuse(parts):
    result=parts[0]
    for part in parts[1:]:result=result.fuse(part)
    return result.clean()


def build(l,cfg,line,lines=2):
    lo=l['low'];hi=l['high'];H=hi[2];wall=cfg['wall_thickness_mm'];fw=cfg['flange_width_mm'];t=line*lines
    outerlo=lo-[wall,wall,0];outerhi=hi+[wall,wall,0]
    cavity=bounds_box(lo-[0,0,1],hi+[0,0,1]);walls=bounds_box(outerlo,outerhi).cut(cavity)
    fl=bounds_box(outerlo-[fw,fw,0],np.array([outerhi[0]+fw,outerhi[1]+fw,cfg['flange_thickness_mm']])).cut(cavity)
    parts=[walls,fl];supports={};retainers={};ribs={};local={};markers=[];designs=[]
    for p in l['ports']:
        k=p['port_id'];r=p['radius_mm'];R=r+p['cradle_clearance_mm'];w=p['wall_distance_from_trim_mm'];z0=-p['trim_origin'][2];top=H-p['trim_origin'][2];inner=R+.55
        # Twin vertical radial fins connect the transverse cradle webs to the
        # inner wall. Every fin starts on the dish, not on a suspended roof.
        cradle_ys=[2.2] if p['role']=='MASTER' else [1.65,2.35]
        fin_front=min(cradle_ys)-t/2
        retz=r+p['anti_float_gap_mm'];rett=3*line;retback=w-.24;retfront=1.05;slide=1.6;guide_start=w-1.1;rear_rise=.7
        solids=[]
        for sgn in [-1,1]:
            a,b=sorted([sgn*inner,sgn*(inner+t)])
            solids.append(bounds_box([a,fin_front,z0],[b,w+.2,retz+rett+rear_rise+1.1]))
            # 45-degree inner ledge: supported as it grows, no long roof.
            x=sgn*inner;z=retz+rett+rear_rise+.12
            triangle=[(x,z),(x,z+.55),(x-sgn*.55,z+.55)]
            ledge=cq.Workplane('XZ').polyline(triangle).close().extrude(w+.2-guide_start).val()
            # XZ workplane extrudes toward -Y.
            ledge=ledge.translate((0,w+.2,0));solids.append(ledge)
        cutter=cq.Solid.makeCylinder(R,w+2,cq.Vector(0,-.5,0),cq.Vector(0,1,0))
        for yy in cradle_ys:
            cradle=bounds_box([-inner-t,yy-t/2,z0],[inner+t,yy+t/2,0]).cut(cutter);solids.append(cradle)
        if p['role']=='MASTER':
            # Unilateral distal reference, with a positive end-face gap; the
            # attachment to the cerebral tree is far upstream, never touched.
            stop_y=p['straight_length_mm']/2+line
            solids.append(bounds_box([-t/2,stop_y,z0],[t/2,w+.2,r*.5]))
        support=fuse(solids)
        # Separate flat plate. A short wallward slide engages both ledges;
        # reverse slide is needed before lifting it out. No wall penetration.
        # Raise the rear tongue over the distal cap. Only the front ceiling
        # approaches the approved straight contact zone; print on its planar
        # YZ face, so this modest step requires no print support.
        outline=[(retfront,retz),(3.,retz),(3.5,retz+rear_rise),(retback,retz+rear_rise),
                 (retback,retz+rear_rise+rett),(3.5,retz+rear_rise+rett),(3.,retz+rett),(retfront,retz+rett)]
        ret=cq.Workplane('YZ').polyline(outline).close().extrude(2*(inner-.12)).val().translate((-inner+.12,0,0))
        # Tall L-shaped alignment rib; its inner FACE, not its centre, is on
        # the trim plane, so no half-thickness imprint leaks to the final side.
        rib_outer=cfg['trim_notch_half_width_mm']-.3
        rib=fuse([bounds_box([inner+t+.15,0,z0],[rib_outer,t,top-.5]),
                  bounds_box([rib_outer-t,0,z0],[rib_outer,w+.2,top-.5])])
        marker=bounds_box([rib_outer-t,w-.05,top],[rib_outer,w+wall,top+t])
        support_world=world(support,p);rib_world=world(rib,p);ret_world=world(ret,p);mark_world=world(marker,p)
        parts.extend([support_world,rib_world,mark_world]);supports[k]=support_world;retainers[k]=ret_world;ribs[k]=rib_world;markers.append(mark_world)
        local[k]=dict(support=support,retainer=ret,rib=rib,retainer_slide_mm=slide,retainer_thickness_mm=rett,cradle_ys=cradle_ys,rib_outer_mm=rib_outer)
        designs.append(dict(port_id=k,support_blade_thickness_mm=t,support_blade_line_count=lines,wall_to_extension_span_mm=w-min(cradle_ys),cradle_type='OPEN_LOWER_SEMICIRCULAR_ARC',cradle_clearance_mm=p['cradle_clearance_mm'],contact_arc_degrees=180,contact_length_mm=t*len(cradle_ys),cradle_web_centres_mm=cradle_ys,trim_side='SACRIFICIAL_TRIM_SIDE',axial_float_plus_minus_mm=p['axial_float_mm'],retainer_top_gap_mm=p['anti_float_gap_mm'],retainer_slide_mm=slide,rib_thickness_mm=t,rib_inner_face_on_trim_plane=True))
    for d in designs:
        d.update(support_blade_thickness=d['support_blade_thickness_mm'],wall_to_extension_span=d['wall_to_extension_span_mm'],cradle_clearance=d['cradle_clearance_mm'])
    return dict(frame=fuse(parts),walls=walls,flange=fl,supports=supports,retainers=retainers,ribs=ribs,markers=markers,local=local,designs=designs,cavity=bounds_box(lo,hi),flange_bounds=np.array([outerlo-[fw,fw,0],outerhi+[fw,fw,0]]),lines=lines,line_width=line)


def clearance_study(l,cfg,line,out):
    """Export actual alternative cradles; clearances are engineering candidates,
    not a measured PLA shrinkage distribution or a force calculation.
    """
    import copy
    from .support_inputs import table
    from .casting_mold_qc import fcl_distance
    rows=[]
    for f in cfg['radial_clearance_line_multipliers']:
        alt=copy.deepcopy(l)
        for p in alt['ports']:p['cradle_clearance_mm']=line*f
        g=build(alt,cfg,line,2);folder=Path(out)/f'clearance_{line*f:.2f}mm';folder.mkdir(parents=True,exist_ok=True)
        for p in alt['ports']:
            k=p['port_id'];part=shape_mesh(g['supports'][k]);part.export(folder/(k+'_lower_support.stl'));cq.exporters.export(g['supports'][k],str(folder/(k+'_lower_support.step')))
            distance,_,_=fcl_distance(l['core'],part)
            rows.append(dict(port_id=k,port_diameter_mm=p['diameter_mm'],clearance_mm=line*f,clearance_lines=f,minimum_actual_support_core_distance_mm=distance,collision_free=distance>0,selected=abs(next(x for x in l['ports'] if x['port_id']==k)['cradle_clearance_mm']-line*f)<1e-8,physical_tolerance_calibrated=False))
    table(Path(out)/'comparison.csv',rows);return rows
