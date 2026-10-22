"""Reuse V1 section polygons; classify and audit actual layer requirements."""
from __future__ import annotations
import numpy as np
from shapely.geometry import Point,LineString,GeometryCollection
from shapely.ops import unary_union
from .support_detection import layers
from .support_inputs import table,save

def dimensions(poly):
    rect=poly.minimum_rotated_rectangle
    if rect.geom_type!='Polygon':return 0.,0.,np.array([1.,0.])
    p=np.asarray(rect.exterior.coords)[:4];v=np.roll(p,-1,axis=0)-p;lengths=np.linalg.norm(v,axis=1);i=int(np.argmax(lengths))
    return float(lengths[i]),float(lengths[(i+1)%4]),v[i]/max(lengths[i],1e-12)

def bridge_evidence(poly,previous,limit,width):
    """A short dimension alone is insufficient: require two opposing anchors."""
    if previous.is_empty:return dict(bridgeable=False,span_mm=None,direction=[1.,0.],both_ends_anchored=False)
    length,narrow,axis=dimensions(poly);origin=np.array(poly.representative_point().coords[0]);best=None
    for direction in [axis,np.array([-axis[1],axis[0]])]:
        across=np.array([-direction[1],direction[0]]);verts=np.array(poly.exterior.coords);projections=(verts-origin)@across;lo,hi=projections.min(),projections.max()
        offsets=np.linspace(lo+1e-5,hi-1e-5,max(3,int(np.ceil((hi-lo)/(width/2)))+1));anchors=previous.buffer(width*.1);spans=[];valid=True
        for shift in offsets:
            center=origin+across*shift;line=LineString([center-direction*(length+2*width),center+direction*(length+2*width)]);cut=poly.intersection(line)
            for part in (list(cut.geoms) if hasattr(cut,'geoms') else [cut]):
                if part.geom_type!='LineString' or part.length<.00001:continue
                ends=np.array(part.coords)[[0,-1]];span=float(np.linalg.norm(ends[1]-ends[0]));unit=(ends[1]-ends[0])/max(span,1e-12)
                anchored=all(anchors.intersects(Point(p+sign*unit*width*.6)) for p,sign in zip(ends,[-1,1]))
                if not anchored or span>limit+1e-6:valid=False;break
                spans.append(span)
            if not valid:break
        if valid and spans and (best is None or max(spans)<best['span_mm']):best=dict(bridgeable=True,span_mm=max(spans),direction=direction.tolist(),both_ends_anchored=True,cross_section_checks=len(spans),bridge_state='VALIDATED_BRIDGE_CANDIDATE',validation_scope='Software two-anchor geometry only; material calibration remains separate')
    return best or dict(bridgeable=False,span_mm=None,direction=axis.tolist(),both_ends_anchored=False)

def analyse(mesh,profile,kind,foundation=None,out=None,bridge_limit=None):
    limit=profile['bridge_limit_mm'] if bridge_limit is None else bridge_limit
    field=layers(mesh,profile,.005,foundation);rows=[];polygons=[]
    for i,(z,bad) in enumerate(zip(field['levels'],field['unsupported_polygons'])):
        previous=field['polygons'][i-1] if i else GeometryCollection()
        for poly in (list(bad.geoms) if hasattr(bad,'geoms') else [bad]):
            if poly.area<.005:continue
            length,width,direction=dimensions(poly);bridge=bridge_evidence(poly,previous,limit,profile['line_width'])
            status='BRIDGEABLE' if bridge['bridgeable'] else 'SUPPORT_REQUIRED'
            row=dict(region_id=f'{kind}_{len(rows):05}',geometry_class=kind,layer_index=i,z_mm=float(z),area_mm2=float(poly.area),
                width_mm=width,span_length_mm=length,span_direction=direction.tolist(),classification=status,bridge=bridge,
                bridge_calibration=profile['bridge_calibration_status'],bounds=list(poly.bounds),representative_xy=list(poly.representative_point().coords[0]))
            rows.append(row);polygons.append(poly)
    if out:
        table(out/f'{kind}_layer_regions.csv',rows);save(out/f'{kind}_layer_regions.json',rows)
    return dict(rows=rows,polygons=polygons,levels=field['levels'],layer_polygons=field['polygons'],self_supported_layers=len(field['levels'])-len({r['layer_index'] for r in rows}),
        source_model=mesh,foundation=foundation)

def residual_support_area(poly,anchors,profile,previous=None):
    remaining=poly.difference(anchors);bad=[];bridges=[]
    material=anchors.union(previous) if previous is not None else anchors
    for part in (list(remaining.geoms) if hasattr(remaining,'geoms') else [remaining]):
        if part.is_empty or part.area<.000001:continue
        evidence=bridge_evidence(part,material,profile['bridge_limit_mm'],profile['line_width'])
        if evidence['bridgeable']:bridges.append((part,evidence))
        else:bad.append(part)
    return unary_union(bad),bridges

def coverage(field,roofs,profile,actual_ids=None,bridge_segments=None):
    rows=[];area=0.;miss=0.
    for row,poly in zip(field['rows'],field['polygons']):
        if row['classification']!='SUPPORT_REQUIRED':continue
        eligible=[]
        for roof in roofs:
            if actual_ids is not None and roof['roof_id'] not in actual_ids:continue
            dz=row['z_mm']-roof['top_z_mm']
            if -profile['layer_height']*.5<=dz<=roof['gap_mm']+profile['layer_height']*2.0:
                eligible.append(roof['_footprint'] if actual_ids is not None else roof['_footprint'].buffer(profile['line_width']/2))
        cover=unary_union(eligible) if eligible else GeometryCollection()
        index=row.get('layer_index',0);previous=field.get('layer_polygons',[])[index-1] if index>0 and 'layer_polygons' in field else GeometryCollection()
        remaining,bridges=residual_support_area(poly,cover,profile,previous)
        if actual_ids is not None:
            residuals=[remaining]
            if bridge_segments is not None and len(bridge_segments):
                seg=bridge_segments[np.abs(bridge_segments.mean(axis=1)[:,2]-row['z_mm'])<=profile['layer_height']*1.5]
            else:seg=np.empty((0,2,3))
            delta=seg[:,1,:2]-seg[:,0,:2];length=np.linalg.norm(delta,axis=1)
            for part,evidence in bridges:
                aligned=np.abs(delta@np.array(evidence['direction']))/np.maximum(length,1e-9)>.85
                actual=unary_union([LineString(s[:,:2]).buffer(profile['line_width']/2,cap_style=2) for s in seg[aligned]])
                residuals.append(part.difference(actual))
            remaining=unary_union(residuals)
        uncovered=float(remaining.area);area+=poly.area;miss+=uncovered
        rows.append(dict(**row,covered_area_mm2=float(poly.area-uncovered),uncovered_area_mm2=uncovered,covered=uncovered<.005,bridge_evidence=[b for _,b in bridges]))
    return dict(percent=100*(1-miss/area) if area else 100.,required_area_sum_mm2=area,uncovered_area_sum_mm2=miss,
        region_count=len(rows),uncovered_region_count=sum(not r['covered'] for r in rows),all_required_covered=all(r['covered'] for r in rows),rows=rows,
        method='Complete layer polygons: roof coverage plus multi-chord two-anchor short bridges; actual mode requires direction-aligned actual span extrusion coverage (bridge infill or perimeter); no 5 mm arc substitution')
