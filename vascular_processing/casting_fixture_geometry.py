"""BCF: read-only frozen core, parameterized bottomless CAD and separate retainers."""
from __future__ import annotations
import copy,json,math
from pathlib import Path
import numpy as np
import trimesh
import cadquery as cq
from .support_inputs import ROOT,discover,digest,save,transform,swc_branches


def rigid(matrix):
    m=np.asarray(matrix,float)
    if m.shape!=(4,4) or not np.allclose(m[3],[0,0,0,1]) or not np.allclose(m[:3,:3].T@m[:3,:3],np.eye(3),atol=1e-10) or not np.isclose(np.linalg.det(m[:3,:3]),1):raise ValueError('SCALING_OR_NONRIGID_PLACEMENT')
    return m


def mesh_qc(mesh):
    counts=np.bincount(mesh.edges_unique_inverse)
    r=dict(watertight=bool(mesh.is_watertight),manifold=bool(np.all(counts==2)),positive_volume=bool(mesh.volume>0),volume_mm3=float(mesh.volume),connected_components=len(mesh.split(only_watertight=False)),degenerate_faces=int(np.count_nonzero(mesh.area_faces<1e-12)),winding_consistent=bool(mesh.is_winding_consistent))
    r['passed']=bool(r['watertight'] and r['manifold'] and r['positive_volume'] and r['connected_components']==1 and not r['degenerate_faces'] and r['winding_consistent']);return r


def bounds_box(low,high):
    low=np.array(low,float);high=np.array(high,float)
    if np.any(high<=low):raise ValueError('NONPOSITIVE_BOX')
    return cq.Workplane('XY').box(*(high-low)).translate(tuple((low+high)/2)).val()


def local_shape(shape,port):
    n=np.asarray(port['axis']);u=np.cross(n,[0,0,1.]);angle=math.degrees(math.atan2(u[1],u[0]))
    return shape.rotate((0,0,0),(0,0,1),angle).translate(tuple(port['wall_intersection']))


def cylinder_y(radius,start,end):
    return cq.Solid.makeCylinder(radius,end-start,cq.Vector(0,start,0),cq.Vector(0,1,0))


def shape_mesh(shape):
    vertices,faces=shape.tessellate(.01,.12)
    return trimesh.Trimesh([v.toTuple() for v in vertices],faces,process=True)


def export_shape(shape,path,cfg):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    cq.exporters.export(shape,str(path.with_suffix('.step')))
    cq.exporters.export(shape,str(path.with_suffix('.stl')),tolerance=cfg['cad_tolerance_mm'],angularTolerance=cfg['cad_angular_tolerance'])
    mesh=trimesh.load_mesh(path.with_suffix('.stl'));qc=mesh_qc(mesh)
    if not qc['passed']:raise ValueError('CAD_STL_QC_FAILED '+str(path)+' '+str(qc))
    return mesh,qc


def source_layout(inputs,cfg,line_width):
    if cfg['scale_factor']!=1:raise ValueError('SCALING_FORBIDDEN')
    core=trimesh.load_mesh(inputs['core']);central=trimesh.load_mesh(inputs['vascular'])
    meta=json.loads(inputs['geometry'].read_text());ports=[];wall=float(cfg['wall_thickness_mm']);fl=float(cfg['selected_axial_float_mm']);safe=cfg['safe_straight_margin_mm']
    c=cfg['selected_radial_clearance_line_multiplier']*line_width
    # Conservative, directly saved straight segment target -> end. Curved route
    # reconstruction and all vascular fitting are intentionally unnecessary.
    for p in meta['ports']:
        n=np.asarray(p['direction'],float);n/=np.linalg.norm(n);a=np.asarray(p['target']);b=np.asarray(p['end']);length=float((b-a)@n)
        if abs(n[2])>1e-8 or np.count_nonzero(abs(n)>.999999)!=1 or np.linalg.norm((b-a)-length*n)>1e-5:raise ValueError('PORT_STRAIGHT_METADATA_INVALID_'+p['port_id'])
        if length < wall+2*(safe+fl):raise ValueError('PORT_EXTENSION_TOO_SHORT_FOR_FIXTURE_'+p['port_id'])
        ports.append(dict(port_id=p['port_id'],face=p['face'],axis=n.tolist(),straight_start=a.tolist(),end=b.tolist(),diameter_mm=2*p['radius_mm'],radius_mm=p['radius_mm'],straight_length_mm=length,attachment=p['accepted_attachment']['calculation_origin_mm'],total_route_length_mm=p['length_mm'],parent_branch=p['parent_branch'],swc_id=p['swc_id'],radial_clearance_mm=c,axial_float_mm=fl))
    low=central.bounds[0].copy();high=central.bounds[1].copy();low[:2]-=cfg['side_margin_mm'];high[:2]+=cfg['side_margin_mm']
    for axis in (0,1):
        for side,sign in ((0,-1),(1,1)):
            group=[p for p in ports if np.asarray(p['axis'])[axis]*sign>.99]
            if group:
                # Feasible wall positions come from real frozen straight lengths.
                lower=max(sign*np.asarray(p['straight_start'])[axis]+safe+fl for p in group)
                upper=min(sign*np.asarray(p['end'])[axis]-wall-safe-fl for p in group)
                if lower>upper:raise ValueError('PORT_EXTENSION_TOO_SHORT_FOR_SHARED_WALL')
                proposed=sign*(low[axis] if side==0 else high[axis]);value=sign*np.clip(proposed,lower,upper)
            else:value=core.bounds[side,axis]+sign*cfg['unconnected_side_margin_mm']
            (low if side==0 else high)[axis]=value
    low[2]=core.bounds[0,2]-cfg['bottom_clearance_mm'];high[2]=core.bounds[1,2]+cfg['top_clearance_mm']
    matrix=np.eye(4);matrix[:3,3]=[-(low[0]+high[0])/2,-(low[1]+high[1])/2,-low[2]];rigid(matrix)
    core.apply_transform(matrix);central.apply_transform(matrix);low+=matrix[:3,3];high+=matrix[:3,3]
    casting=rigid(json.loads(inputs['casting_transform'].read_text())['transform_4x4']);_,top=swc_branches(inputs['smoothed_swc'],matrix@casting,inputs['branch_mapping'])
    for p in ports:
        for k in ['straight_start','end','attachment']:p[k]=transform([p[k]],matrix)[0].tolist()
        n=np.asarray(p['axis']);axis=int(np.argmax(abs(n)));side=int(n[axis]>0);center=np.asarray(p['straight_start']).copy();center[axis]=(high if side else low)[axis]
        p['wall_intersection']=center.tolist();p['capture_center']=(center+n*wall/2).tolist()
        p['inner_distance_from_transition_mm']=float((center-n*fl-n*safe-np.asarray(p['straight_start']))@n)
        p['end_reserve_mm']=float((np.asarray(p['end'])-center-n*(wall+fl+safe))@n)
        p['distance_to_attachment_mm']=float(np.linalg.norm(center-np.asarray(p['attachment'])))
        p['distance_to_nearest_bifurcation_mm']=float(np.min(np.linalg.norm(top['bifurcations']-center,axis=1))) if len(top['bifurcations']) else None
        source_port=trimesh.load_mesh(inputs['candidate']/'manufacturing'/f"{p['port_id']}.stl");source_port.apply_transform(matrix);sections=[]
        for offset in [-fl,wall/2,wall+fl]:
            plane=center+n*offset;section=source_port.section(plane_normal=n,plane_origin=plane)
            if section is None:raise ValueError('PORT_STRAIGHT_SECTION_MISSING_'+p['port_id'])
            radius=np.linalg.norm(section.vertices-plane,axis=1);contours=len(section.discrete)
            good=bool(contours==1 and radius.max()<=p['radius_mm']+5e-5 and radius.min()>=p['radius_mm']*np.cos(np.pi/24)-5e-5)
            sections.append(dict(wall_normal_offset_mm=offset,min_radius_mm=float(radius.min()),max_radius_mm=float(radius.max()),closed_contours=contours,passed=good))
            if not good:raise ValueError('FROZEN_STRAIGHT_GEOMETRY_MISMATCH_'+p['port_id'])
        p['straight_section_checks']=sections
        if p['distance_to_attachment_mm']<safe or (p['distance_to_nearest_bifurcation_mm'] is not None and p['distance_to_nearest_bifurcation_mm']<safe):raise ValueError('UNSAFE_FIXTURE_NEAR_VASCULAR_TRANSITION_'+p['port_id'])

    # Longest known straight segment, then diameter, then bifurcation clearance.
    master=max(ports,key=lambda p:(round(p['straight_length_mm'],5),p['diameter_mm'],p['distance_to_nearest_bifurcation_mm'] or 0))['port_id']
    for p in ports:p['role']='MASTER_DATUM' if p['port_id']==master else 'FLOATING_SECONDARY';p['axial_float_mm']=0. if p['role']=='MASTER_DATUM' else fl
    return dict(core=core,central=central,ports=ports,low=low,high=high,transform=matrix,master_port=master,source_qc=mesh_qc(trimesh.load_mesh(inputs['core'])))


def build_frame(layout,cfg,line_width):
    lo=np.asarray(layout['low']);hi=np.asarray(layout['high']);H=hi[2];wall=cfg['wall_thickness_mm'];flange=min(cfg['flange_width_candidates_mm']);fheight=cfg['flange_thickness_mm'];slide=line_width*cfg['slide_clearance_line_multiplier'];key=cfg['key_width_mm'];pin_size=cfg['locking_pin_size_mm'];pw=pin_size+2*slide
    outerlow=lo-np.array([wall,wall,0]);outerhigh=hi+np.array([wall,wall,0]);flo=outerlow-np.array([flange,flange,0]);fhi=outerhigh+np.array([flange,flange,0]);fhi[2]=fheight
    cavity=bounds_box(lo-[0,0,1],hi+[0,0,1]);frame=bounds_box(outerlow,outerhigh).fuse(bounds_box(flo,fhi)).cut(cavity)
    # Continuous rectangular annular channel strictly outside cavity.
    midlo=outerlow.copy();midhi=outerhigh.copy();midlo[:2]-=flange/2;midhi[:2]+=flange/2;midlo[2]=0;midhi[2]=0
    sw=cfg['seal_channel_width_mm'];depth=cfg['seal_channel_depth_mm']
    seal=bounds_box(midlo-[sw/2,sw/2,.05],midhi+[sw/2,sw/2,depth]).cut(bounds_box(midlo+[sw/2,sw/2,-.1],midhi+[-sw/2,-sw/2,depth+.1]))
    frame=frame.cut(seal);closures={};pins={};port_seals={};local_parts={}
    for p in layout['ports']:
        R=p['radius_mm']+p['radial_clearance_mm'];h=H-p['wall_intersection'][2];span=R+key;pinz=h-3.0;pinv=wall/2
        if h<6 or wall<pw+.3:raise ValueError('FIXTURE_WALL_OR_HEADROOM_TOO_SMALL')
        # Narrow inner mouth + wider outer keyway, both open all the way to top.
        slot_profile=[(-R,-.05),(R,-.05),(span+slide,wall+.05),(-span-slide,wall+.05)]
        mouth=cq.Workplane('XY').polyline(slot_profile).close().extrude(h+1).val().fuse(cylinder_y(R,-.05,wall+.05))
        frame=frame.cut(local_shape(mouth,p))
        from shapely.geometry import Polygon,box as polygon_box
        inset=Polygon(slot_profile).buffer(-slide,join_style=2).intersection(polygon_box(-100,slide,100,wall-slide))
        key_profile=list(inset.exterior.coords)[:-1]
        closure=cq.Workplane('XY').polyline(key_profile).close().extrude(h-slide).translate((0,0,slide)).val().cut(cylinder_y(R,-1,wall+1))
        # Removable outward-facing locking pin engages the outer T key AND
        # the solid inner rail. A blind inner socket keeps PLA outside cavity.
        pinu=R+key/2
        skin=2*line_width
        bore=bounds_box([pinu-pw/2,skin,pinz-pw/2],[pinu+pw/2,wall+2,pinz+pw/2])
        frame=frame.cut(local_shape(bore,p));closure=closure.cut(bore)
        pin=bounds_box([pinu-pin_size/2,skin+slide,pinz-pin_size/2],[pinu+pin_size/2,wall+.6,pinz+pin_size/2])
        head=bounds_box([pinu-pin_size,wall+slide,pinz-pin_size],[pinu+pin_size,wall+slide+1,pinz+pin_size]);pin=pin.fuse(head)
        if p['role']=='MASTER_DATUM':
            end=float((np.asarray(p['end'])-p['wall_intersection'])@np.asarray(p['axis']));stop=end+p['radial_clearance_mm']
            arm=bounds_box([-span,wall+slide,R+1],[span,stop+1,R+3])
            connector=bounds_box([-R+slide,wall/2+slide,R+1],[R-slide,wall+slide+.1,R+3]);arm=arm.fuse(connector)
            stopbody=bounds_box([-span,stop,-R-.5],[span,stop+1,R+3]);closure=closure.fuse(arm).fuse(stopbody)
            p['axial_stop_gap_mm']=p['radial_clearance_mm'];p['axial_stop_normal']=p['axis'];p['end_stop_position']=(np.asarray(p['end'])+np.asarray(p['axis'])*p['radial_clearance_mm']).tolist()
        # External rectangular moat plus split annulus around the port. The
        # same cutter acts on both halves; assembly gaps need physical sealant.
        su=span+.8;zb=-R-1.3;zt=h-.7;w=sw/2
        moat=bounds_box([-su-w,wall-depth,zb-w],[su+w,wall+.05,zt+w]).cut(bounds_box([-su+w,wall-depth-.1,zb+w],[su-w,wall+.1,zt-w]))
        annulus=cylinder_y(R+.65,wall-depth,wall+.05).cut(cylinder_y(R+.25,wall-depth-.1,wall+.1))
        cutters=moat.fuse(annulus);frame=frame.cut(local_shape(cutters,p));closure=closure.cut(cutters)
        closures[p['port_id']]=local_shape(closure,p);pins[p['port_id']]=local_shape(pin,p);port_seals[p['port_id']]=local_shape(cutters,p)
        local_parts[p['port_id']]=dict(closure=closure,pin=pin)
        p.update(blind_pin_socket_inner_skin_mm=skin,slot_width_mm=2*R,slot_key_width_mm=2*(span+slide),slide_clearance_mm=slide,closure_top_mm=H,locking_pin_center=transform([[pinu,wall/2,pinz]],port_matrix(p))[0].tolist(),locking_pin_direction=p['axis'])
    return dict(frame=frame,closures=closures,pins=pins,local_parts=local_parts,seal_shape=seal,seal_loop=np.array([[midlo[0],midlo[1],0],[midhi[0],midlo[1],0],[midhi[0],midhi[1],0],[midlo[0],midhi[1],0],[midlo[0],midlo[1],0]]),port_seals=port_seals,flange_width=flange,flange_bounds=np.array([flo,fhi]),minimum_dish_diameter=float(np.linalg.norm(fhi[:2]-flo[:2])),height=H)


def port_matrix(port):
    n=np.asarray(port['axis']);m=np.eye(4);m[:3,0]=np.cross(n,[0,0,1.]);m[:3,1]=n;m[:3,3]=port['wall_intersection'];return rigid(m)
