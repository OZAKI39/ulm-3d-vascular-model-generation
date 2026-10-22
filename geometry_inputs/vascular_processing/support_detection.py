"""Surface, cross-section, local-minimum and SWC support requirements."""
from __future__ import annotations
import time
import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import Polygon,GeometryCollection,Point
from shapely.ops import unary_union
from .support_inputs import save,table
from . import support_freefloating as ff

def overhang(mesh,angle=45.):
    normals=mesh.face_normals;theta=np.degrees(np.arccos(np.clip(-normals[:,2],-1,1)))
    return dict(centroids=mesh.triangles_center,area=mesh.area_faces,angle_from_down=theta,required=theta<angle)

def section_polygon(path):
    if path is None:return GeometryCollection()
    polygons=[]
    for loop in path.discrete:
        if len(loop)>=4 and np.linalg.norm(loop[0]-loop[-1])<1e-5:
            p=Polygon(loop[:,:2]);p=p if p.is_valid else p.buffer(0)
            if not p.is_empty:polygons.append(p)
    result=GeometryCollection()
    for p in sorted(polygons,key=lambda p:-p.area):result=result.symmetric_difference(p)
    return result

def unsupported_layer(current,previous,height,width,angle=45.):
    # At most half a line width may project over empty space; also bound the
    # horizontal step by the configured angle from vertical.
    reach=min(width*.5,height*np.tan(np.radians(angle)))
    return current.difference(previous.buffer(reach)),reach

def layers(mesh,profile,min_area=.005,supporting_mesh=None):
    h=profile['layer_height'];first=profile['first_layer_height']
    levels=np.arange(first/2,mesh.bounds[1,2]+h,h)
    # path vertices use the same XY axes for all horizontal cuts.
    cuts=mesh.section_multiplane([0,0,0],[0,0,1],levels)
    othercuts=supporting_mesh.section_multiplane([0,0,0],[0,0,1],levels) if supporting_mesh is not None else [None]*len(levels)
    previous=GeometryCollection();rows=[];polys=[];unsupported=[];risk=[]
    for n,(z,path) in enumerate(zip(levels,cuts)):
        current=section_polygon(path);bad,reach=unsupported_layer(current,previous,h,profile['line_width'])
        if n==0:bad=GeometryCollection() # bed supports the actual first layer
        pieces=list(bad.geoms) if hasattr(bad,'geoms') else [bad]
        for p in pieces:
            if p.area>=min_area:
                center=p.representative_point();risk.append([center.x,center.y,z]);rows.append(dict(layer=n,z_mm=float(z),unsupported_area_mm2=float(p.area),x=center.x,y=center.y))
        polys.append(current);unsupported.append(bad);previous=current.union(section_polygon(othercuts[n]))
    return dict(levels=levels,polygons=polys,unsupported_polygons=unsupported,points=np.asarray(risk).reshape(-1,3),rows=rows,reach_mm=reach)

def layer_field_distance(points,radii,layer):
    """Distance to entire unsupported footprints, NOT their representative point.

    A long straight tube can create one long skinny unsupported polygon. Using
    only its centroid would silently miss both ends, including external ports.
    """
    distances=[];levels=layer['levels'];polygons=layer['unsupported_polygons']
    for p,r in zip(points,radii):
        bound=max(.6,1.5*r);lo=np.searchsorted(levels,p[2]-bound);hi=np.searchsorted(levels,p[2]+bound,side='right');point=Point(*p[:2]);best=np.inf
        for index in range(lo,hi):
            polygon=polygons[index]
            if polygon.is_empty:continue
            dz=abs(levels[index]-p[2])
            if dz>=best:continue
            best=min(best,float(np.hypot(point.distance(polygon),dz)))
        distances.append(best)
    return np.array(distances)

def bottom_points(mesh,xyz):
    hits,rays,_=mesh.ray.intersects_location(xyz,np.tile([0,0,-1],(len(xyz),1)),multiple_hits=True)
    bottom=np.full_like(xyz,np.nan);best=np.full(len(xyz),np.inf)
    for p,r in zip(hits,rays):
        d=xyz[r,2]-p[2]
        if d>=-1e-5 and d<best[r]:bottom[r]=p;best[r]=d
    return bottom

def contiguous_runs(mask):
    a=np.flatnonzero(np.diff(np.r_[False,mask,False]));return [(int(i),int(j-1)) for i,j in zip(a[::2],a[1::2])]

def long_spans(branch,bottom,layer_points,profile,cfg,already_supported=None,layer_distances=None):
    # Ray to the next lower material: self interior was already exited at bottom.
    angle=np.degrees(np.arcsin(np.clip(np.abs(branch['tangent'][:,2]),0,1)))
    if layer_distances is not None:distance=np.asarray(layer_distances)
    elif len(layer_points):distance=cKDTree(layer_points).query(bottom)[0]
    else:distance=np.full(len(bottom),np.inf)
    need=(angle<=cfg['horizontal_tangent_deg']) & (distance<=np.maximum(1.5*branch['radius'],.6))
    if already_supported is not None:need &= ~already_supported
    runs=[]
    for start,end in contiguous_runs(need):
        length=float(branch['arc'][end]-branch['arc'][start])
        runs.append(dict(branch=branch['branch_id'],arc_start=float(branch['arc'][start]),arc_end=float(branch['arc'][end]),length_mm=length,
            mean_angle=float(angle[start:end+1].mean()),radius=float(branch['radius'][start:end+1].mean()),
            risk='CRITICAL' if length>cfg['critical_max_span_mm'] else 'ATTENTION',coverage_status='UNASSESSED',start_index=start,end_index=end))
    return need,angle,distance,runs

def detect(mesh,branches,topology,profile,cfg,out,box=None):
    start=time.monotonic();surf=overhang(mesh,cfg['surface_overhang_from_down_deg'])
    layer=layers(mesh,profile,cfg['minimum_layer_unsupported_area_mm2'],box)
    floating,ffinfo=ff.detect(mesh,profile['layer_height'],cfg['freefloating_pixel_mm'],cfg['freefloating_effective_radius_mm'])
    minima=ff.local_minima(mesh);raw=[];regions=[];critical_keys=set()
    down_tree=cKDTree(surf['centroids'][surf['required']]);float_tree=cKDTree(floating) if len(floating) else None
    bif=topology['bifurcations'];bif_tree=cKDTree(bif) if len(bif) else None
    for b in branches:
        bottom=bottom_points(mesh,b['xyz']);valid=np.isfinite(bottom).all(axis=1)
        if not valid.all():bottom[~valid]=b['xyz'][~valid]-np.column_stack([np.zeros((~valid).sum()),np.zeros((~valid).sum()),b['radius'][~valid]])
        already=box.contains(bottom-[0,0,.02]) if box is not None else np.zeros(len(bottom),bool)
        ld=layer_field_distance(bottom,b['radius'],layer)
        need,angle,distance,runs=long_spans(b,bottom,layer['points'],profile,cfg,already,ld);regions.extend(runs)
        sd=down_tree.query(bottom)[0];fd=float_tree.query(bottom)[0] if float_tree else np.full(len(bottom),np.inf)
        for r in runs:
            if r['risk']=='CRITICAL':critical_keys.update((b['branch_id'],i) for i in range(r['start_index'],r['end_index']+1))
        for i,(p,s,r,k) in enumerate(zip(bottom,b['arc'],b['radius'],b['curvature'])):
            if already[i]:continue
            sources=[]
            if sd[i]<.4:sources.append('SURFACE_OVERHANG')
            if distance[i]<max(.6,1.5*r):sources.append('LAYER_UNSUPPORTED')
            if fd[i]<cfg['sample_step_mm']*1.5:sources.append('FREEFLOATING')
            if need[i]:sources.append('SWC_LONG_SPAN')
            if not sources:continue
            critical=(b['branch_id'],i) in critical_keys
            # Start of a new floating island must have a real support anchor.
            local_min=0<i<len(bottom)-1 and p[2]<=min(bottom[i-1,2],bottom[i+1,2]) and p[2]<max(bottom[i-1,2],bottom[i+1,2])
            if local_min and 'LAYER_UNSUPPORTED' in sources:critical=True;sources.append('LOCAL_MINIMUM');critical_keys.add((b['branch_id'],i))
            row=dict(point_id=f"{b['branch_id']}:{i}",branch_id=b['branch_id'],sample_index=i,xyz=p.tolist(),arc_position=float(s),local_radius=float(r),local_curvature=float(k),
                detectors_triggered=sources,unsupported_span=max([a['length_mm'] for a in runs if a['start_index']<=i<=a['end_index']],default=0),
                risk_level='CRITICAL' if critical else 'ATTENTION',distance_to_bifurcation=float(bif_tree.query(p)[0]) if bif_tree else 1e6,
                layer_distance_mm=float(distance[i]) if np.isfinite(distance[i]) else None,surface_hit_valid=bool(valid[i]),is_port=b.get('is_port',False))
            raw.append(row)
        b.update(bottom=bottom,needs_support=need)
    # Surface-only risk on external artificial ports remains separately visible.
    # These points are never silently labelled as anatomical SWC branches.
    save(out/'support_detection_config.json',dict(**cfg,profile=profile,lower_layer_reach_mm=layer['reach_mm'],formula='min(line_width / 2, layer_height * tan(45 deg))',
        cross_section_method='Triangle plane intersections; polygon XOR preserves holes; previous section buffered by extrusion reach',engineering_span_not_physical_law=True))
    table(out/'tables/layer_unsupported.csv',layer['rows']);table(out/'tables/all_candidate_support_points.csv',raw)
    table(out/'tables/unsupported_spans.csv',regions)
    save(out/'freefloating_audit.json',ffinfo)
    np.savez_compressed(out/'support_requirement_surface.npz',**surf)
    np.savez_compressed(out/'detector_points.npz',freefloating=floating,local_minima=minima,layer=layer['points'])
    result=dict(raw=raw,regions=regions,freefloating=ffinfo,runtime_seconds=time.monotonic()-start,critical_region_count=sum(r['risk']=='CRITICAL' for r in regions),
        longest_dangerous_span_mm=max([r['length_mm'] for r in regions if r['risk']=='CRITICAL'],default=0))
    save(out/'detection.json',result)
    return result,surf,layer,floating
