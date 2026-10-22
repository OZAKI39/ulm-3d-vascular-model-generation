"""Full-branch curtains. No sparse detector participates in support selection.

Only new support is constructed. Surface rays use the frozen vascular STL;
the SWC defines arc coordinates and local frames, never a replacement vessel.
"""
from __future__ import annotations
import copy
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from .support_v3_geometry import sweep, cut_section
from .support_v2_mesh import support_stl
from .support_custom_qc import intersection_volume


def boolean(meshes,operation='union'):
    """Use the installed Manifold double-precision API for thin support joins."""
    import manifold3d as md
    meshes=list(meshes)
    if not meshes or (operation=='difference' and not len(meshes[0].faces)) or (operation=='intersection' and any(not len(m.faces) for m in meshes)):
        return trimesh.Trimesh()
    solids=[md.Manifold(md.Mesh64(np.array(m.vertices,dtype=np.float64,order='C',copy=True),np.array(m.faces,dtype=np.uint64,order='C',copy=True))) for m in meshes if len(m.faces)]
    if not solids:return trimesh.Trimesh()
    if any(s.status()!=md.Error.NoError for s in solids):raise ValueError('RFUS_BOOLEAN_INPUT_NOT_MANIFOLD')
    op={'union':md.OpType.Add,'difference':md.OpType.Subtract,'intersection':md.OpType.Intersect}[operation]
    solid=md.Manifold.batch_boolean(solids,op);m=solid.to_mesh64()
    return trimesh.Trimesh(m.vert_properties[:,:3],m.tri_verts,process=False)


def union(meshes):return boolean(meshes)


def resample(branch, step):
    b=copy.deepcopy(branch);s=np.asarray(b['arc']);arc=np.linspace(0,s[-1],max(3,int(np.ceil(s[-1]/step))+1))
    for key in ('xyz','radius'):
        a=np.asarray(b[key]);b[key]=np.column_stack([np.interp(arc,s,a[:,j]) for j in range(a.shape[1])]) if a.ndim==2 else np.interp(arc,s,a)
    b['arc']=arc;t=np.gradient(b['xyz'],arc,axis=0);t/=np.maximum(np.linalg.norm(t,axis=1)[:,None],1e-12)
    b['tangent']=t;b['curvature']=np.linalg.norm(np.gradient(t,arc,axis=0),axis=1)
    return b


def lower_points(mesh, xyz):
    """First frozen surface below each centreline point (not lowest vessel)."""
    xyz=np.asarray(xyz);direction=np.tile([0.,0.,-1.],(len(xyz),1))
    loc,ids,_=mesh.ray.intersects_location(xyz+direction*1e-6,direction,multiple_hits=True)
    result=np.full_like(xyz,np.nan,dtype=float)
    if len(ids):
        order=np.lexsort((-loc[:,2],ids));ids=ids[order];loc=loc[order]
        unique,index=np.unique(ids,return_index=True);result[unique]=loc[index]
    # Rays exactly in an end-cap plane can be ignored by the triangle kernel.
    # Retry at sub-float32 XY offsets, still against the frozen surface. No
    # radius estimate or unrelated lower vessel is substituted for a hit.
    for i in np.flatnonzero(~np.isfinite(result).all(1)):
        eps=float(np.spacing(np.float32(max(1.,np.max(np.abs(xyz[i]))))))
        origins=xyz[i]+np.array([[eps,0,-1e-6],[-eps,0,-1e-6],[0,eps,-1e-6],[0,-eps,-1e-6]])
        hits,_,_=mesh.ray.intersects_location(origins,np.tile([0.,0.,-1.],(4,1)),multiple_hits=True)
        if len(hits):result[i]=hits[np.argmax(hits[:,2])]
    return result


def corridor(branch, mesh, inner, profile, cfg):
    b=resample(branch,profile['line_width']*cfg['sample_spacing_lines'])
    low=lower_points(mesh,b['xyz']);valid=np.isfinite(low).all(1)
    # No radius-based invented surface is allowed when the frozen ray misses.
    b.update(lower=low,surface_valid=valid,gravity=np.tile([0.,0.,-1.],(len(low),1)))
    t=b['tangent'];normal=np.c_[-t[:,1],t[:,0]];n=np.linalg.norm(normal,axis=1)
    normal[n<1e-8]=[1.,0.];normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-12)
    b['lateral']=normal
    lo,hi=np.asarray(inner);xy=b['xyz'][:,:2];inside=np.all((xy>lo[:2])&(xy<hi[:2]),axis=1)
    b['root_domain']=np.where(inside,'BOX_INNER_BOTTOM','BUILD_PLATE')
    b['base']=np.where(inside,lo[2],0.)
    h=profile['layer_height'];first=profile['first_layer_height']
    b['roof_top']=first+np.floor((low[:,2]-profile['roof_gap_mm']-first)/h)*h
    return b


def segment_ranges(b, topology, cfg):
    """Adaptive arc segmentation; no branch or interval is sparsified."""
    arc=b['arc'];curv=b['curvature'];bif=np.asarray(topology.get('bifurcations',[])).reshape(-1,3)
    distance=cKDTree(bif).query(b['xyz'])[0] if len(bif) else np.full(len(arc),np.inf)
    changes=np.flatnonzero(b['root_domain'][1:]!=b['root_domain'][:-1])+1
    edges=[0];start=0
    while start<len(arc)-1:
        if distance[start]<max(2.,3*b['radius'][start]):length=cfg['bifurcation_segment_length_mm']
        elif curv[start]>.12 or b['xyz'][start,2]-b['base'][start]>18:length=cfg['segment_length_mm'][0]
        elif curv[start]>.04:length=cfg['segment_length_mm'][1]
        else:length=cfg['segment_length_mm'][2]
        end=min(len(arc)-1,max(start+1,int(np.searchsorted(arc,arc[start]+length))))
        crossings=changes[(changes>start)&(changes<end)]
        if len(crossings):end=int(crossings[0])
        edges.append(end);start=end
    return list(zip(edges[:-1],edges[1:]))


def strip_transition(xy0, xy1, bottom, top, width0, width1):
    """Ruled strip with independently displaced lower and upper stations."""
    xy0=np.asarray(xy0);xy1=np.asarray(xy1);t=np.gradient(xy1,axis=0)
    normal=np.c_[-t[:,1],t[:,0]]/np.maximum(np.linalg.norm(t,axis=1)[:,None],1e-12)
    bottom=np.broadcast_to(bottom,(len(t),));top=np.broadcast_to(top,(len(t),));vertices=[]
    for i in range(len(t)):
        vertices.extend([np.r_[xy0[i]-normal[i]*width0/2,bottom[i]],np.r_[xy0[i]+normal[i]*width0/2,bottom[i]],np.r_[xy1[i]-normal[i]*width1/2,top[i]],np.r_[xy1[i]+normal[i]*width1/2,top[i]]])
    n=len(t);faces=[[0,1,3],[0,3,2],[4*n-4,4*n-2,4*n-1],[4*n-4,4*n-1,4*n-3]]
    for i in range(n-1):
        a,b=4*i,4*(i+1)
        for j,k in [(0,2),(2,3),(3,1),(1,0)]:faces.extend([[a+j,b+j,b+k],[a+j,b+k,a+k]])
    base=trimesh.Trimesh(vertices,faces,process=True)
    trimesh.repair.fix_normals(base)
    if not base.is_volume:raise ValueError('TRANSITION_NOT_CLOSED')
    return base


def stacked_curtain(xy,lateral,base,top,profile,width,roofwidth,neckwidth,offset):
    """One closed ruled solid from foot to ridge; avoid coplanar Boolean seams."""
    from shapely.geometry import Polygon
    h=profile['layer_height'];necktop=base+6*h;flared=necktop+np.ceil((max(width-neckwidth,0)/2/.5+h)/h)*h
    roofbottom=top-3*h;bodytop=roofbottom-max(2*h,np.ceil((abs(offset)/.5+2*h)/h)*h)
    z=np.array([np.maximum(base-.03,0),base+2*h,base+2*h,necktop,flared,bodytop,roofbottom,top])
    widths=np.array([width,width,neckwidth,neckwidth,width,width,roofwidth,roofwidth])
    offsets=np.array([offset]*6+[0.,0.])
    keep=[0]+[i for i in range(1,len(z)) if not (np.allclose(z[i],z[i-1],atol=1e-12) and widths[i]==widths[i-1] and offsets[i]==offsets[i-1])]
    z=z[keep];widths=widths[keep];offsets=offsets[keep];levels=len(z);ring=2*levels;vertices=[]
    for station in range(len(xy)):
        for side,sequence in [(1,range(levels)),(-1,range(levels-1,-1,-1))]:
            for j in sequence:vertices.append(np.r_[xy[station]+lateral[station]*(offsets[j]+side*widths[j]/2),z[j,station]])
    vertices=np.asarray(vertices);faces=[]
    for station in range(len(xy)-1):
        a=station*ring;b=a+ring
        for j in range(ring):k=(j+1)%ring;faces.extend([[a+j,b+j,b+k],[a+j,b+k,a+k]])
    for station in [0,len(xy)-1]:
        ids=np.arange(station*ring,(station+1)*ring);xyz=vertices[ids];coords=np.c_[(xyz[:,:2]-xy[station])@lateral[station],xyz[:,2]]
        coords=np.c_[np.r_[offsets+widths/2,(offsets-widths/2)[::-1]],np.r_[z[:,station],z[::-1,station]]]
        vv,ff=trimesh.creation.triangulate_polygon(Polygon(coords),engine='triangle',triangle_args='pY');indices=cKDTree(coords).query(vv)[1]
        faces.extend(ids[indices[ff]].tolist())
    mesh=trimesh.Trimesh(vertices,faces,process=True);mesh.update_faces(mesh.area_faces>1e-14);mesh.remove_unreferenced_vertices();trimesh.repair.fix_normals(mesh)
    if not mesh.is_volume:raise ValueError('RFUS_STACKED_CURTAIN_NOT_CLOSED')
    return mesh


def make_segment(xy, top, base, lateral, profile, lines=2, roof_lines=2, neck_lines=1, offset=0., twin=False):
    w,h=profile['line_width'],profile['layer_height'];xy=np.asarray(xy);top=np.asarray(top);base=np.broadcast_to(base,top.shape)
    tangent=np.gradient(xy,axis=0);length=np.linalg.norm(tangent,axis=1)
    if np.any(length<1e-8):raise ValueError('FOLDED_OR_VERTICAL_RFUS_PROJECTION')
    lateral=np.c_[-tangent[:,1],tangent[:,0]]/length[:,None]
    width=lines*w;roofwidth=roof_lines*w;neckwidth=neck_lines*w
    if lines<2 or neck_lines not in (1,2):raise ValueError('RFUS_REQUIRES_TWO_LINE_CURTAIN')
    if np.any(top-base<1.5):raise ValueError('INSUFFICIENT_CURTAIN_HEIGHT')
    roofbottom=top-3*h;transition_height=max(2*h,np.ceil((abs(offset)/.5+2*h)/h)*h)
    bodytop=roofbottom-transition_height;necktop=base+6*h;flared=necktop+np.ceil((max(width-neckwidth,0)/2/.5+h)/h)*h
    if np.any(bodytop<=flared+h):raise ValueError('OFFSET_TRANSITION_TOO_TALL')
    paths=[xy+lateral*offset] if not twin else [xy+lateral*abs(offset),xy-lateral*abs(offset)]
    feet=[];necks=[];bodies=[];transitions=[]
    for path in paths:
        feet.append(sweep(path,np.maximum(base-.03,0),base+2*h,width))
        necks.append(sweep(path,base+2*h-.01,necktop,neckwidth))
        flare=sweep(path,necktop-.01,flared,width,neckwidth)
        body=sweep(path,flared-.01,bodytop+.01,width)
        bodies.extend([flare,body]);transitions.append(strip_transition(path,xy,bodytop,roofbottom+.01,width,roofwidth))
    roof=sweep(xy,roofbottom,top,roofwidth)
    parts=dict(feet=union(feet),necks=union(necks),curtains=union(bodies+transitions),roofs=roof,handles=trimesh.Trimesh())
    solids=[stacked_curtain(xy,lateral,base,top,profile,width,roofwidth,neckwidth,sign*abs(offset) if twin else offset) for sign in ([1,-1] if twin else [1])]
    parts['mesh']=union(solids);parts.update(path_xy=xy,base=base,top=top,lateral=lateral,offset_mm=float(offset),twin=bool(twin),width_mm=width,roof_width_mm=roofwidth,neck_width_mm=neckwidth,neck_cut_z=float(np.mean(base)+4*h))
    return parts


def collision_clear(segment, vascular, box):
    hit=intersection_volume(segment['mesh'],vascular)
    overlap=boolean([segment['mesh'],box],'intersection')
    bad=boolean([overlap,segment['feet']],'difference') if len(overlap.faces) else overlap
    boxhit=abs(float(bad.volume)) if len(bad.faces) else 0.
    return dict(passed=hit<1e-5 and boxhit<1e-4,vascular_port_intersection_mm3=hit,illegal_box_intersection_mm3=boxhit)


def route_segment(xy,top,base,lateral,profile,cfg,vascular,box,lines,roof_lines,neck_lines):
    w=profile['line_width'];attempts=[]
    trials=[('DIRECT',0.,False)]+[('LATERAL_OFFSET',j*w,False) for j in [1,-1,2,-2]]
    trials += [('TWIN_CURTAINS',j*w,True) for j in [3,4,5,6] if j*w<=cfg['maximum_local_transition_mm']]
    trials += [('LOCAL_PILLAR_TRANSITION',j*w,False) for j in [3,-3,4,-4,5,-5,6,-6] if abs(j*w)<=cfg['maximum_local_transition_mm']]
    direct=None
    for kind,offset,twin in trials:
        try:
            r=make_segment(xy,top,base,lateral,profile,lines,roof_lines,neck_lines,offset,twin)
            if kind=='DIRECT':direct=r
            check=collision_clear(r,vascular,box);attempts.append(dict(route=kind,offset_mm=offset,**check))
            if check['passed']:r.update(route=kind,route_attempts=attempts);return r,attempts
        except (ValueError,RuntimeError) as exc:attempts.append(dict(route=kind,offset_mm=offset,passed=False,error=str(exc)))
    # A perpendicular crossing cannot be avoided by moving an entire curtain
    # parallel to the lower vessel. Keep the full roof and split its lower web
    # into two legs joined by a pointed, self-growing transition above it.
    if direct is not None:
        from .support_v3_geometry import window_cutter
        crossing=boolean([direct['mesh'],vascular],'intersection')
        if len(crossing.faces):
            tangent=np.asarray(xy[-1])-np.asarray(xy[0]);tangent/=max(np.linalg.norm(tangent),1e-12)
            projections=crossing.vertices[:,:2]@tangent;lo,hi=projections.min(),projections.max();width=float(hi-lo+2*(w+profile['roof_gap_mm']))
            apex=float(crossing.bounds[1,2]+profile['roof_gap_mm']+width)
            if width/2<=cfg['maximum_local_transition_mm'] and apex<np.min(top)-3*profile['layer_height']:
                center=crossing.bounds.mean(0);center[:2]+=tangent*((lo+hi)/2-center[:2]@tangent)
                cutter,_=window_cutter(center,tangent,width,float(np.min(base)-1),apex,max(3.,roof_lines*w*6),slope=.5)
                candidate=copy.deepcopy(direct)
                for key in ['mesh','feet','necks','curtains','roofs']:candidate[key]=boolean([candidate[key],cutter],'difference')
                check=collision_clear(candidate,vascular,box);attempts.append(dict(route='TWIN_CURTAINS_AROUND_CROSSING',transition_half_span_mm=width/2,apex_mm=apex,**check))
                if check['passed']:
                    candidate.update(route='TWIN_CURTAINS_AROUND_CROSSING',route_attempts=attempts,twin=True,crossing_transition_width_mm=width)
                    return candidate,attempts
    return None,attempts


def add_handle(segment,vascular,box,others,inner,profile,cfg):
    """Root handles in their own curtain above its cut, never in the box."""
    w,h=profile['line_width'],profile['layer_height'];lo,hi=np.asarray(inner)
    access_z=hi[2]-(hi[2]-lo[2])*cfg['handle_access_fraction'];width=cfg['handle_lines']*w
    mesh=segment['mesh'];xy=segment['path_xy'];normals=segment['lateral'];attempts=[]
    targetz=max(access_z+.5,float(np.max(segment['top'])+.8))
    for fraction in [.5,.25,.75,.1,.9]:
        i=int(round(fraction*(len(xy)-1)));normal=normals[i]
        for sign in [1.,-1.]:
            for lines in [5,7,9,12,15,18]:
                distance=lines*w;targetxy=xy[i]+sign*normal*distance
                if not np.all((targetxy>lo[:2]+width)&(targetxy<hi[:2]-width)):continue
                startz=max(float(segment['base'][i])+1.,float(segment['top'][i])-distance/.45-1.)
                endz=startz+distance/.45
                if endz>=targetz-h:continue
                # Ascending, supported transition, with a vertical upper grip.
                from .support_v2_primitives import loft
                stem=loft(np.r_[xy[i]+normal*segment['offset_mm'],startz],np.r_[targetxy,endz],[width,width])
                post=loft(np.r_[targetxy,endz-.02],np.r_[targetxy,targetz],[width,width])
                handle=union([stem,post]);collision=intersection_volume(handle,vascular)>1e-5 or intersection_volume(handle,box)>1e-5
                collision=collision or any(intersection_volume(handle,m)>1e-5 for m in others)
                attempts.append(dict(fraction=fraction,side=sign,offset_mm=distance,collision=collision))
                if collision:continue
                combined=union([mesh,handle])
                if len(combined.split())!=len(mesh.split()):continue
                segment.update(handles=handle,mesh=combined,handle_target=np.r_[targetxy,targetz-.25].tolist(),handle_width_mm=width,handle_access_zone_min_z_mm=float(access_z),handle_attempts=attempts)
                return True
    segment.update(handle_target=None,handle_attempts=attempts);return False


def neck_audit(segment):
    from shapely.geometry import box
    lo,hi=segment['mesh'].bounds;z=segment['neck_cut_z'];origin=np.r_[(lo[:2]+hi[:2])/2,z]
    aperture=box(*(np.r_[lo[:2]-origin[:2]-.01,hi[:2]-origin[:2]+.01]))
    return cut_section(segment['mesh'],segment['necks'],origin,[0,0,1],aperture)


def disconnect_base(segment):
    """Only material below the designed cut plane is retained at the root."""
    m=segment['mesh'];lo,hi=m.bounds;z=segment['neck_cut_z']+.03
    cutter=trimesh.creation.box(np.r_[hi[:2]-lo[:2]+2,z-lo[2]+1],transform=trimesh.transformations.translation_matrix(np.r_[(hi[:2]+lo[:2])/2,(z+lo[2]-1)/2]))
    return boolean([m,cutter],'difference')


def exported_support(mesh):
    """Resolve support-only binary-STL rounding seams, with an explicit bound."""
    from .support_v2_mesh import exact_stl
    try:return exact_stl(mesh)
    except ValueError:pass
    import manifold3d as md
    solid=md.Manifold(md.Mesh64(np.array(mesh.vertices,dtype=np.float64,order='C',copy=True),np.array(mesh.faces,dtype=np.uint64,order='C',copy=True)))
    for tol in [2e-5,5e-5,.0001,.0005,.001]:
        m=solid.simplify(tol).to_mesh64();candidate=trimesh.Trimesh(m.vert_properties[:,:3],m.tri_verts,process=False)
        try:
            result=exact_stl(candidate);result.metadata['support_only_simplify_bound_mm']=tol;return result
        except ValueError:pass
    for grid in [1e-5,5e-5,1e-4,5e-4,.001]:
        m=trimesh.Trimesh(np.round(mesh.vertices/grid)*grid,mesh.faces,process=True)
        m.update_faces(m.unique_faces());m.update_faces(m.area_faces>0);m.remove_unreferenced_vertices()
        trimesh.repair.fix_normals(m,multibody=True)
        if not m.is_volume:continue
        try:
            result=exact_stl(m);result.metadata['support_export_grid_mm']=grid;result.metadata['maximum_vertex_snap_mm']=float(grid*np.sqrt(3)/2)
            return result
        except ValueError:continue
    raise ValueError('RFUS_SUPPORT_EXPORT_NOT_CLOSED_AFTER_BOUNDED_NUMERICAL_REPAIR')
