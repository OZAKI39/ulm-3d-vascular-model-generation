"""s1-10 geometry: closed ABS mold, separate retainers, real solvent lands.

The s1-9 files are immutable history. Shared read-only placement/CAD kernels
are reused, but this module builds a new closed-bottom production solid.
"""
from __future__ import annotations
import numpy as np
import cadquery as cq
from shapely.geometry import Polygon,box
from shapely.ops import unary_union
from .casting_fixture_v2_geometry import layout as fit_frozen_core,port_transform,world,fuse
from .casting_fixture_geometry import bounds_box,shape_mesh,rigid


def layout(inputs,cfg,line):
    try:l=fit_frozen_core(inputs,cfg,line)
    except ValueError as e:
        if 'PORT_STRAIGHT_FIXTURE_ZONE_TOO_SHORT_' in str(e):raise ValueError(str(e).replace('PORT_STRAIGHT_FIXTURE_ZONE_TOO_SHORT_','PORT_EXTENSION_FIXTURE_ZONE_TOO_SHORT_')) from e
        raise
    shift=np.array([0.,0.,cfg['bottom_thickness_mm']]);t=np.eye(4);t[:3,3]=shift
    l['transform']=rigid(t@l['transform'])
    for m in [l['core'],l['central'],*l['source_ports'].values()]:m.apply_translation(shift)
    for p in l['ports']:
        for k in ['straight_start','end','attachment','trim_origin','approved_contact_start','approved_contact_end']:p[k]=np.asarray(p[k])+shift
    l['low']+=shift;l['high']+=shift;l['bifurcations']+=shift
    return l


def conformal_land(core,p,width,length,center):
    """Loft vertical prisms to the actual downward-facing frozen STL planes.

    This is a separate auxiliary ABS part touching a finite patch; no offset,
    refit, cutter penetration or modification of the vascular mesh is used.
    """
    m=port_transform(p);local=core.copy();local.apply_transform(np.linalg.inv(m));tri=local.triangles;normals=local.face_normals
    domain=box(-width/2,center-length/2,width/2,center+length/2);lo=np.array(domain.bounds[:2]);hi=np.array(domain.bounds[2:]);r=p['radius_mm']
    select=(normals[:,2]<-.5)&(tri[:,:,:2].min(1)<=hi).all(1)&(tri[:,:,:2].max(1)>=lo).all(1)&(tri[:,:,2].max(1)>-r-.02)&(tri[:,:,2].min(1)<-r+.25)
    pieces=[];patches=[];witnesses=[];areas=[];face_ids=[];zbase=-r-p['cradle_clearance_mm']-.12
    for face_id in np.flatnonzero(select):
        xyz=tri[face_id];normal=normals[face_id];polygon=Polygon(xyz[:,:2]).intersection(domain)
        if polygon.is_empty or polygon.area<1e-9:continue
        xy=np.asarray(polygon.exterior.coords)[:-1];z=xyz[0,2]-((xy-xyz[0,:2])@normal[:2])/normal[2];top=np.c_[xy,z];bottom=np.c_[xy,np.full(len(xy),zbase)]
        if np.min(z)<=zbase:raise ValueError('ACCESS_LAND_BASE_ABOVE_CORE')
        bw=cq.Wire.makePolygon([cq.Vector(*q) for q in bottom],close=True);tw=cq.Wire.makePolygon([cq.Vector(*q) for q in top],close=True)
        pieces.append(cq.Solid.makeLoft([bw,tw],ruled=True));patches.append(polygon);areas.append(polygon.area/abs(normal[2]));face_ids.append(int(face_id))
        xy0=np.array(polygon.representative_point().coords[0]);z0=xyz[0,2]-((xy0-xyz[0,:2])@normal[:2])/normal[2];witnesses.append(np.r_[xy0,z0]@m[:3,:3].T+m[:3,3])
    covered=unary_union(patches).area if patches else 0.
    if abs(covered-domain.area)>1e-6:raise ValueError('SOLVENT_ACCESS_LAND_NOT_ON_FROZEN_SURFACE_'+p['port_id'])
    return fuse(pieces),dict(port_id=p['port_id'],contact_area_mm2=float(sum(areas)),projected_contact_area_mm2=float(covered),width_mm=width,length_mm=length,source_face_ids=face_ids,contact_witnesses=witnesses,center_from_trim_mm=center,nominal_gap_mm=0.,interference_required=False,contact_role='FINITE_AREA_SOLVENT_ACCESS_LAND_ON_APPROVED_STRAIGHT_ZONE',source_geometry_changed=False)


def build(l,cfg,line):
    lo=l['low'];hi=l['high'];wall=cfg['wall_thickness_mm'];base=cfg['bottom_thickness_mm'];H=hi[2];t=line*cfg['blade_line_count'];depth=cfg['trim_guide_depth_mm'];half=cfg['trim_notch_half_width_mm']
    outerlo=lo-[wall,wall,base];outerhi=hi+[wall,wall,0]
    shell=bounds_box(outerlo,outerhi).cut(bounds_box(lo,hi+[0,0,1]))
    wall_slab=bounds_box(np.r_[outerlo[:2],base],outerhi).cut(bounds_box(lo-[0,0,1],hi+[0,0,1]))
    bottom=bounds_box(outerlo,np.r_[outerhi[:2],base]);parts=[shell];baseparts=[shell];fixturebases={};fixtures={};retainers={};lands={};guides={};localparts={};access=[];design=[]
    for p in l['ports']:
        k=p['port_id'];r=p['radius_mm'];R=r+p['cradle_clearance_mm'];w=p['wall_distance_from_trim_mm'];zfloor=base-p['trim_origin'][2];ztop=H-p['trim_origin'][2];inner=R+.55;ys=[2.2] if p['role']=='MASTER' else [1.65,2.35]
        retz=r+p['anti_float_gap_mm'];rett=3*line;rear=.7;front=1.05;back=w-.24;slide=1.6;guide_start=w-1.1;solids=[]
        for sgn in [-1,1]:
            a,b=sorted([sgn*inner,sgn*(inner+t)]);solids.append(bounds_box([a,min(ys)-t/2,zfloor-.05],[b,w+.2,retz+rett+rear+1.1]))
            x=sgn*inner;z=retz+rett+rear+.12
            ledge=cq.Workplane('XZ').polyline([(x,z),(x,z+.55),(x-sgn*.55,z+.55)]).close().extrude(w+.2-guide_start).val().translate((0,w+.2,0));solids.append(ledge)
        bore=cq.Solid.makeCylinder(R,w+2,cq.Vector(0,-.5,0),cq.Vector(0,1,0))
        for y in ys:solids.append(bounds_box([-inner-t,y-t/2,zfloor-.05],[inner+t,y+t/2,0]).cut(bore))
        if p['role']=='MASTER':solids.append(bounds_box([-t/2,p['straight_length_mm']/2+line,zfloor-.05],[t/2,w+.2,r*.5]))
        # Narrow vertical access web grows from bottom and inner wall. The
        # conformal land is the only zero-gap contact; the cradle stays loose.
        y=cfg['access_land_center_from_trim_mm'];length=cfg['access_land_length_mm'];solids.append(bounds_box([-t/2,y-length/2,zfloor-.05],[t/2,w+.2,-r-p['cradle_clearance_mm']-.08]))
        land,record=conformal_land(l['core'],p,t,length,y);access.append(record);land_world=world(land,p);lands[k]=land_world
        fixturebases[k]=world(fuse(solids),p);baseparts.append(fixturebases[k])
        fixture=world(fuse(solids+[land]),p);fixtures[k]=fixture;parts.append(fixture)
        outline=[(front,retz),(3.,retz),(3.5,retz+rear),(back,retz+rear),(back,retz+rear+rett),(3.5,retz+rear+rett),(3.,retz+rett),(front,retz+rett)]
        ret=cq.Workplane('YZ').polyline(outline).close().extrude(2*(inner-.12)).val().translate((-inner+.12,0,0))
        # A top-exposed sacrificial stem makes each otherwise isolated
        # retainer directly bath-accessible, independent of its guide gaps.
        stem_back=guide_start-.25
        chimney=bounds_box([-inner+.12,stem_back-3*line,retz+rear],[inner-.12,stem_back,ztop+cfg['retainer_exposure_above_fill_mm']]);ret=fuse([ret,chimney]);retainers[k]=world(ret,p)
        # Shallow relief on the INNER BOTTOM transfers a U-shaped guide into
        # PDMS. All relief lies on the sacrificial side. No persistent tall rib.
        marks=[bounds_box([-half,0,zfloor-.05],[half,t,zfloor+depth])]
        for sign in [-1,1]:
            xx=[-half,-half+t] if sign<0 else [half-t,half]
            marks.append(bounds_box([xx[0],0,zfloor-.05],[xx[1],w+.1,zfloor+depth]))
            marks.append(bounds_box([xx[0],w-depth,zfloor-.05],[xx[1],w+.1,ztop]))
        mark=world(fuse(marks),p);guides[k]=mark;parts.append(mark);baseparts.append(mark)
        localparts[k]=dict(retainer=ret,cradle_ys=ys,retainer_slide_mm=slide,inner_guide_x=inner,guide_top_local_z=retz+rett+rear+1.1,stem_center_y_mm=stem_back-1.5*line)
        design.append(dict(port_id=k,role=p['role'],cradle_clearance_mm=p['cradle_clearance_mm'],cradle_diameter_mm=2*R,retainer_guide_clearance_mm=.12,retainer_top_gap_mm=p['anti_float_gap_mm'],blade_thickness_mm=t,blade_line_count=cfg['blade_line_count'],access_channel_minimum_width_mm=t,access_land_contact_area_mm2=record['contact_area_mm2'],secondary_axial_float_mm=p['axial_float_mm'],wall_to_cradle_span_mm=w-min(ys),retainer_exposure_above_fill_mm=cfg['retainer_exposure_above_fill_mm'],shallow_trim_guide_depth_mm=depth))
    return dict(mold=fuse(parts),mold_without_lands=fuse(baseparts),fixture_bases=fixturebases,shell=shell,walls=wall_slab,bottom=bottom,cavity=bounds_box(lo,hi),fixtures=fixtures,retainers=retainers,lands=lands,guides=guides,access_records=access,design=design,local_parts=localparts,outer_bounds=np.array([outerlo,outerhi]),line_width_mm=line,blade_thickness_mm=t)
