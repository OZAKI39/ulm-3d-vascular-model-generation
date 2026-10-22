"""New independent support geometry: continuous ridges over segmented curtains."""
from __future__ import annotations
import copy,pickle
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import save,table
from .support_rfus_geometry import lower_points,resample,union,boolean,strip_transition
from .support_v3_geometry import sweep,window_cutter
from .support_v2_mesh import support_stl
from .core_support_inputs import place_branches


def overlap(a,b):
    if not len(a.faces) or not len(b.faces):return 0.
    if np.any(a.bounds[1]<b.bounds[0]) or np.any(b.bounds[1]<a.bounds[0]):return 0.
    m=boolean([a,b],'intersection');return abs(float(m.volume)) if len(m.faces) else 0.


def bed_roots(support):
    """Conservative audit of positive material shells, not internal voids.

    Closed, negatively oriented cavity boundaries are not floating pieces.
    Positive disconnected shells remain failures even when very small.
    """
    return [dict(component_id=i,minimum_z_mm=float(m.bounds[0,2]),
                 signed_volume_mm3=float(m.volume),
                 bed_rooted=abs(m.bounds[0,2])<2e-5)
            for i,m in enumerate(support.split()) if m.volume>0.]


def envelope(mesh,branches,profile,gap):
    result=[]
    for branch in branches:
        b=resample(branch,min(profile['line_width'],profile['support_line_width'])/2)
        low=lower_points(mesh,b['xyz'])
        # At a downward-facing terminal cap, the centreline endpoint itself
        # is already on the lower surface. Starting below it skips that cap.
        # Retry from a float-precision distance ABOVE the endpoint and keep
        # only real ray intersections at/below the original endpoint.
        for i in np.flatnonzero(~np.isfinite(low).all(1)):
            q=b['xyz'][i];eps=4*float(np.spacing(np.float32(max(1.,np.max(np.abs(q))))))
            hits,_,_=mesh.ray.intersects_location([q+[0,0,eps]],[ [0,0,-1.] ],multiple_hits=True)
            hits=hits[hits[:,2]<=q[2]+eps] if len(hits) else hits
            if len(hits):low[i]=hits[np.argmax(hits[:,2])]
        b['lower']=low;b['valid']=np.isfinite(low).all(1)
        tangent=b['tangent'][:,:2];normal=np.c_[-tangent[:,1],tangent[:,0]];norm=np.linalg.norm(normal,axis=1);normal[norm<1e-8]=[1,0];normal/=np.linalg.norm(normal,axis=1)[:,None];b['lateral']=normal
        # Bambu samples mesh at layer mid-height but extrudes at layer top.
        # Half-layer subtraction gives an actual extrusion gap in [gap,gap+h].
        b['roof_top']=low[:,2]-gap-profile['layer_height']/2-1e-5
        result.append(b)
    return result


def intervals(b):
    edges=[0];i=0;s=b['arc'];n=len(s)
    while i<n-1:
        length=5. if b['curvature'][i]>.08 or b['roof_top'][i]>18 else 8. if b['curvature'][i]>.025 else 10.
        # Each anatomical junction belongs to a supported endpoint; shorter
        # terminal segments keep its removal piece compact.
        if not b['is_port'] and (s[i]<3 or s[-1]-s[i]<7):length=5.
        j=min(n-1,max(i+1,int(np.searchsorted(s,s[i]+length))));edges.append(j);i=j
    return list(zip(edges[:-1],edges[1:]))


def ruled(xy,top,w,rw,h,offset=0.,twin=False):
    """Bed-rooted walls, at most 0.5 mm/mm local convergence under the ridge."""
    xy=np.asarray(xy);top=np.asarray(top);t=np.gradient(xy,axis=0);norm=np.linalg.norm(t,axis=1)
    if np.any(norm<1e-8):raise ValueError('VERTICAL_OR_FOLDED_PROJECTION')
    side=np.c_[-t[:,1],t[:,0]]/norm[:,None];roofbottom=top-3*h
    rise=max(2*h,abs(offset)/.45+2*h,abs(w-rw)/.9+2*h)
    bodytop=roofbottom-rise
    if bodytop.min()<.4:raise ValueError('INSUFFICIENT_LOCAL_TRANSITION_HEIGHT')
    roof=sweep(xy,roofbottom,top,rw);walls=[]
    for off in ([abs(offset),-abs(offset)] if twin else [offset]):
        path=xy+side*off
        walls.extend([sweep(path,0,bodytop+.001,w),strip_transition(path,xy,bodytop,roofbottom+.001,w,rw)])
    curtain=union(walls);mesh=union([curtain,roof])
    return dict(mesh=mesh,roofs=roof,curtains=curtain,xy=xy,top=top,lateral=side,offset_mm=offset,twin=twin,width_mm=w,roof_width_mm=rw)


def route(xy,top,vascular,profile,lines,roof_lines):
    w=profile['support_line_width'];h=profile['layer_height'];attempts=[];direct=None
    trials=[('DIRECT',0.,False)]+[('LATERAL_OFFSET',s*w,False) for s in [1,-1,2,-2]]+ [('TWIN_CURTAINS',s*w,True) for s in range(3,10)]
    for kind,off,twin in trials:
        try:
            part=ruled(xy,top,lines*w,roof_lines*w,h,off,twin)
            if direct is None:direct=part
            v=overlap(part['mesh'],vascular);attempts.append(dict(route=kind,offset_mm=off,intersection_mm3=v))
            if v<1e-7:part.update(route=kind,attempts=attempts);return part
        except ValueError as e:attempts.append(dict(route=kind,offset_mm=off,error=str(e)))
    if direct is not None:
        cross=boolean([direct['mesh'],vascular],'intersection')
        if len(cross.faces):
            t=np.asarray(xy[-1])-xy[0];t/=max(np.linalg.norm(t),1e-12);proj=cross.vertices[:,:2]@t;width=float(np.ptp(proj)+2*(profile['roof_gap_mm']+w));center=cross.bounds.mean(0);center[:2]+=t*((proj.min()+proj.max())/2-center[:2]@t)
            apex=cross.bounds[1,2]+profile['roof_gap_mm']+width/.9
            if apex<min(top)-3*h and width<6.:
                cut,_=window_cutter(center,t,width,-1,apex,max(5,4*lines*w),slope=.45)
                p=copy.deepcopy(direct)
                for key in ['mesh','roofs','curtains']:p[key]=boolean([p[key],cut],'difference')
                v=overlap(p['mesh'],vascular);attempts.append(dict(route='TWO_PILLARS_LOCAL_ROOF',intersection_mm3=v,span_mm=width))
                if v<1e-7:p.update(route='TWO_PILLARS_LOCAL_ROOF',attempts=attempts);return p
    # Keep the defect in diagnostics. Never delete an unsupported target.
    if direct is None:return dict(mesh=trimesh.Trimesh(),roofs=trimesh.Trimesh(),curtains=trimesh.Trimesh(),route='UNRESOLVED',attempts=attempts)
    direct.update(route='UNRESOLVED',attempts=attempts);return direct


def layer_reinforcement(vascular,profile,gap,lines):
    """Narrow extra ridges beneath required layer regions, all bed rooted.

    The primary SWC ridge remains continuous; each extra region is widened to
    an integer-line-width strip before slicing, rather than left as a sliver.
    """
    from .support_detection import section_polygon
    h=profile['layer_height'];first=profile['first_layer_height'];w=profile['support_line_width'];tops=np.arange(first,vascular.bounds[1,2]+h,h);cuts=vascular.section_multiplane([0,0,0],[0,0,1],tops-np.where(tops<=first+.001,first,h)/2)
    from shapely.geometry import GeometryCollection
    import shapely
    previous=GeometryCollection();solids=[];rows=[]
    for z,path in zip(tops,cuts):
        poly=section_polygon(path);required=poly.difference(previous.buffer(min(h,w/2)));previous=poly
        if required.is_empty or z<=first+.001:continue
        top=z-gap-h/2-1e-5
        if top<.2:continue
        for region in shapely.get_parts(required):
            tol=region.length*4*np.spacing(np.float32(256.))
            if region.area<=tol:continue
            footprint=region.buffer(lines*w/2,quad_segs=4).simplify(.002,preserve_topology=True)
            for part in shapely.get_parts(footprint):
                if part.geom_type!='Polygon':continue
                mesh=trimesh.creation.extrude_polygon(part,top,engine='triangle');solids.append(mesh);rows.append(dict(layer_top_mm=float(z),ridge_width_lines=lines,required_area_mm2=float(region.area),support_top_mm=float(top)))
    return union(solids),rows


def base_foot(support,profile,expansion=.3):
    from .support_detection import section_polygon
    import shapely
    z=profile['first_layer_height']/2;path=support.section_multiplane([0,0,0],[0,0,1],[z])[0];foot=section_polygon(path).buffer(expansion,quad_segs=6);parts=[]
    for p in shapely.get_parts(foot):
        if p.geom_type=='Polygon':parts.append(trimesh.creation.extrude_polygon(p,profile['first_layer_height'],engine='triangle'))
    return union([support,*parts])


def generate(source,branches,orientation,profile,out,gap=.2,lines=2,roof_lines=2,reinforce=False,wide_foot=False):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);M=np.asarray(orientation['transform']);vascular=source.copy();vascular.apply_transform(M)
    bs=place_branches(branches,M);bs=envelope(vascular,bs,profile,gap);parts=[];notches=[];segments=[];errors=[];p={**profile,'roof_gap_mm':gap}
    for bi,b in enumerate(bs):
        if not b['valid'].all():errors.append(dict(branch=b['branch_id'],reason='LOWER_RAY_MISSING',indices=np.flatnonzero(~b['valid']).tolist()));continue
        # Construct the whole ridge without segmentation gaps. Local pieces
        # share endpoints; cuts below are placed in the body only.
        ranges=intervals(b);bp=[]
        for si,(a,z) in enumerate(ranges):
            # Adjacent ruled solids need finite shared material, not merely
            # coincident cap edges that become nonmanifold in binary STL.
            aa=max(0,a-1);zz=min(len(b['arc'])-1,z+1)
            xy=b['xyz'][aa:zz+1,:2].copy();top=b['roof_top'][aa:zz+1].copy()
            # Extend branch endpoints by half an extrusion width: the slicer's
            # end perimeter otherwise retreats from a centreline cap target.
            if a==0 and np.linalg.norm(xy[1]-xy[0])>1e-8:
                t=(xy[1]-xy[0])/np.linalg.norm(xy[1]-xy[0]);xy=np.vstack([xy[0]-t*profile['support_line_width']/2,xy]);top=np.r_[top[0],top]
            if z==len(b['arc'])-1 and np.linalg.norm(xy[-1]-xy[-2])>1e-8:
                t=(xy[-1]-xy[-2])/np.linalg.norm(xy[-1]-xy[-2]);xy=np.vstack([xy,xy[-1]+t*profile['support_line_width']/2]);top=np.r_[top,top[-1]]
            s=route(xy,top,vascular,p,lines,roof_lines);s.update(segment_id=f'B{bi:02}_S{si:02}',branch_id=b['branch_id'],is_port=b['is_port'],start_index=a,end_index=z,arc_start_mm=float(b['arc'][a]),arc_end_mm=float(b['arc'][z]),length_mm=float(b['arc'][z]-b['arc'][a]))
            if s['route']=='UNRESOLVED':errors.append(dict(branch=b['branch_id'],segment=s['segment_id'],reason='PROJECTION_COLLISION_UNRESOLVED'))
            bp.append(s);segments.append(s)
        valid=[s for s in bp if len(s['mesh'].faces)]
        if not valid:continue
        part=union([s['mesh'] for s in valid]);roof=union([s['roofs'] for s in valid]);curtain=union([s['curtains'] for s in valid])
        cutters=[]
        for ni,(_,i) in enumerate(ranges[:-1]):
            t=b['tangent'][i,:2];norm=np.linalg.norm(t)
            if norm<1e-8:continue
            t=t/norm;z=b['roof_top'][i]-3*profile['layer_height']-.15
            cut,_=window_cutter(np.r_[b['xyz'][i,:2],0],t,profile['support_line_width']*2,-1,z,6.,slope=.45);cutters.append(cut)
            notches.append(dict(notch_id=f'N{bi:02}_{ni:02}',branch_id=b['branch_id'],arc_mm=float(b['arc'][i]),point=np.r_[b['xyz'][i,:2],max(.8,z-2)],plane_origin=np.r_[b['xyz'][i,:2],z],tangent=np.r_[t,0],lateral=np.r_[-t[1],t[0],0],roof_preserved=True,cut_width_mm=profile['support_line_width']*2))
        if cutters:part=boolean([part,union(cutters)],'difference');curtain=boolean([curtain,union(cutters)],'difference')
        parts.append(dict(branch_id=b['branch_id'],mesh=part,roofs=roof,curtains=curtain))
    raw_support=union([p['mesh'] for p in parts]);extra=[]
    if reinforce:
        reinforcement,extra=layer_reinforcement(vascular,p,gap,roof_lines)
        raw_support=union([raw_support,reinforcement])
        # Clip NEW support only, keeping the actual STL lower surface and
        # its requested gap. No Boolean result replaces the vascular input.
        cutters=[]
        for dz in np.arange(-gap-profile['layer_height']/2,.051,.05):
            shifted=vascular.copy();shifted.apply_translation([0,0,dz]);cutters.append(shifted)
        raw_support=boolean([raw_support,union(cutters)],'difference')
    if wide_foot:raw_support=base_foot(raw_support,p)
    with (out/'pre_export_support.pkl').open('wb') as f:pickle.dump(raw_support,f)
    support=support_stl(raw_support);roofs=union([p['roofs'] for p in parts]);curtains=union([p['curtains'] for p in parts])
    if reinforce:
        shifted=support.copy();shifted.apply_translation([0,0,-3*profile['layer_height']])
        roofs=boolean([support,shifted],'difference');curtains=boolean([support,roofs],'difference')
    vascular.export(out/'vascular_core.stl');support.export(out/'custom_support.stl');roofs.export(out/'support_roofs.stl');curtains.export(out/'support_curtains.stl');trimesh.util.concatenate([vascular,support]).export(out/'vascular_core_plus_support_review.stl')
    save(out/'vascular_print_transform.json',dict(transform_4x4=M,scale_factor=1.));save(out/'inverse_transform.json',dict(transform_4x4=np.linalg.inv(M),scale_factor=1.))
    result=dict(candidate_id=out.name,folder=str(out),orientation=orientation,gap_mm=gap,roof_width_lines=roof_lines,curtain_lines=lines,roof_width_mm=roof_lines*profile['support_line_width'],curtain_width_mm=lines*profile['support_line_width'],support_volume_mm3=float(support.volume),intersection_mm3=overlap(support,vascular),segment_count=len(segments),notch_count=len(notches),errors=errors,layer_reinforcement_count=len(extra),layer_reinforcement=extra,wide_base_foot=wide_foot,gap_status='PHYSICAL_REMOVAL_GAP_UNCALIBRATED',geometry_ready=not errors)
    result['bed_roots']=bed_roots(support)
    save(out/'geometry.json',result);save(out/'notches.json',notches)
    with (out/'geometry.pkl').open('wb') as f:pickle.dump(dict(vascular=vascular,support=support,roofs=roofs,curtains=curtains,branches=bs,parts=parts,segments=segments,notches=notches,profile=p),f)
    print('GEOMETRY',out.name,'segments',len(segments),'volume',round(support.volume,2),'intersection',result['intersection_mm3'],'errors',len(errors),flush=True)
    return result
