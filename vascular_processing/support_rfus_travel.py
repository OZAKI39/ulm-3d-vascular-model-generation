"""Chronological real-G-code near-miss proxy; never a nozzle collision proof."""
from __future__ import annotations
import re,zipfile,json,pickle
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import LineString
from scipy.spatial import cKDTree
from .bambu_casting_mold import extruder_coordinate_frame
from .bambu_slice_adapter import arc_polyline
from .support_inputs import table,save


def events(archive):
    """Retractions, E modes, G92, coordinate offsets and true XY arc travels."""
    extr=[];travel=[];features=[];exids=[];trids=[];widths=[];heights=[];filament=[];clock=0;arcs=0
    with zipfile.ZipFile(archive) as z:
        offset,frame=extruder_coordinate_frame(z)
        from .support_inputs import numeric
        settings=json.loads(z.read('Metadata/project_settings.config'))
        default_width=numeric(settings.get('line_width',.22),numeric(settings.get('nozzle_diameter',.2)))
        default_height=numeric(settings.get('layer_height',.1))
        for name in z.namelist():
            if not name.endswith('.gcode'):continue
            position=np.zeros(3);eprev=0.;relative_e=True;absolute=True;feature='';plane='G17';started=False;width=default_width;height=default_height
            for raw in z.read(name).decode(errors='replace').splitlines():
                if raw.startswith('; FEATURE:'):feature=raw.split(':',1)[1].strip()
                if raw.startswith('; LINE_WIDTH:'):width=float(raw.split(':',1)[1])
                if raw.startswith('; LAYER_HEIGHT:'):height=float(raw.split(':',1)[1])
                if raw.startswith('; CHANGE_LAYER'):started=True
                line=raw.split(';',1)[0].strip();clock+=1
                if line=='M83':relative_e=True
                elif line=='M82':relative_e=False
                elif line=='G90':absolute=True
                elif line=='G91':absolute=False
                elif line in ('G17','G18','G19'):plane=line
                match=re.match(r'^(G0|G1|G2|G3|G92)\s+(.*)',line)
                if not match:continue
                vals={k:float(v) for k,v in re.findall(r'([XYZEIJR])([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)',match[2])};code=match[1]
                if code=='G92':
                    eprev=vals.get('E',eprev)
                    for k,a in enumerate('XYZ'):
                        if a in vals:position[k]=vals[a]
                    continue
                new=position.copy()
                for k,a in enumerate('XYZ'):
                    if a in vals:new[k]=vals[a]+(0 if absolute else position[k])
                e=vals.get('E',0.) if relative_e else vals.get('E',eprev)-eprev
                if 'E' in vals:eprev=eprev+vals['E'] if relative_e else vals['E']
                if code in ('G2','G3'):
                    if plane!='G17':raise ValueError('UNSUPPORTED_TRAVEL_ARC_PLANE')
                    points=arc_polyline(position,new,vals,code=='G2');arcs+=1
                else:points=np.array([position,new])
                pieces=np.stack([points[:-1]+offset,points[1:]+offset],axis=1)
                pieces=pieces[np.linalg.norm(pieces[:,1]-pieces[:,0],axis=1)>1e-8]
                if e>1e-8 and feature:
                    extr.extend(pieces);features.extend([feature]*len(pieces));exids.extend([clock]*len(pieces));widths.extend([width]*len(pieces));heights.extend([height]*len(pieces));started=True
                    lengths=np.linalg.norm(pieces[:,1]-pieces[:,0],axis=1);filament.extend(e*lengths/max(lengths.sum(),1e-12))
                elif e<=1e-8 and started and len(pieces):travel.extend(pieces);trids.extend([clock]*len(pieces))
                position=new
    return dict(extrusion=np.asarray(extr).reshape(-1,2,3),travel=np.asarray(travel).reshape(-1,2,3),extrusion_ids=np.asarray(exids),travel_ids=np.asarray(trids),extrusion_widths=np.asarray(widths),extrusion_heights=np.asarray(heights),extrusion_filament_mm=np.asarray(filament),features=features,coordinate_frame=frame,arc_command_count=arcs)


def actual_sections(parsed,mesh,profile,native_support=False,with_filament=False):
    """Use archived per-move LINE_WIDTH, including Arachne variable-width walls."""
    from .support_detection import section_polygon
    tags=parsed['features']
    select=np.array([('support' in t.lower()) if native_support else not any(s in t.lower() for s in ('support','brim','skirt','wipe','prime','custom')) for t in tags])
    segments=parsed['extrusion'][select];widths=parsed['extrusion_widths'][select]
    z=np.round(segments[:,:,2].mean(1),5);levels=np.unique(z);h=profile['layer_height'];first=profile['first_layer_height']
    cuts=mesh.section_multiplane([0,0,0],[0,0,1],levels-np.where(levels<=first+.001,first,h)/2);result={};filament=0.
    for level,path in zip(levels,cuts):
        p=section_polygon(path)
        if p.is_empty:continue
        indices=np.flatnonzero(z==level);lines=shapely.linestrings(segments[indices,:,:2]);clipped=shapely.intersection(lines,p.buffer(.025))
        valid=~shapely.is_empty(clipped)
        if valid.any():result[float(level)]=shapely.union_all(shapely.buffer(clipped[valid],widths[indices][valid]/2,cap_style='flat'))
        if with_filament:
            fraction=shapely.length(clipped)/np.maximum(np.linalg.norm(segments[indices,1]-segments[indices,0],axis=1),1e-12)
            filament+=float(np.sum(parsed['extrusion_filament_mm'][select][indices]*np.minimum(fraction,1.)))
    return (result,filament) if with_filament else result


def audit(archive,c,record):
    folder=Path(record['folder']);cache=folder/'travel_near_miss.json'
    if cache.exists():return json.loads(cache.read_text())
    e=events(archive);ex=e['extrusion'];trav=e['travel'];ids=e['extrusion_ids'];cfg=c['cfg'];w=c['profile']['line_width'];h=c['profile']['layer_height']
    with (folder/'actual_footprints.pkl').open('rb') as f:fp=pickle.load(f)
    import trimesh
    fp['roof']=actual_sections(e,trimesh.load_mesh(folder/'geometry/support_roofs.stl'),c['profile'])
    fp['handle']=actual_sections(e,trimesh.load_mesh(folder/'geometry/support_handles.stl'),c['profile'])
    # Keep only real extrusions attributed to anatomy or explicit support.
    z=np.round(ex[:,:,2].mean(1),5);middle=ex.mean(1);labels=np.full(len(ex),'OTHER',dtype=object)
    with (folder/'segment_geometry.pkl').open('rb') as f:parts=pickle.load(f)
    handles=[s['handle_target'] for s in parts if s.get('handle_target')]
    for zz in np.unique(z):
        ix=np.flatnonzero(z==zz);xy=middle[ix,:2]
        for key,label in [('support','CURTAIN_TOP_OR_BODY'),('roof','ROOF_RIDGE'),('handle','REMOVAL_HANDLE'),('vascular','VASCULAR_OR_PORT')]:
            p=fp[key].get(float(zz))
            if p is not None:labels[ix[shapely.intersects_xy(p,xy[:,0],xy[:,1])]]=label
    keep=labels!='OTHER';widths=e['extrusion_widths'][keep];ex=ex[keep];ids=ids[keep];labels=labels[keep];z=z[keep]
    # Segment queries use exact XY line distances. Chronology filters out all
    # future extrusion, including material not yet printed in the same layer.
    lines=shapely.linestrings(ex[:,:,:2]);tree=shapely.STRtree(lines);rows=[];r=cfg['travel_body_radius_mm'];tip=cfg['travel_tip_clearance_mm'];candidate_min={'nozzle_tip':None,'body':None};by_feature={}
    for i,(segment,event_id) in enumerate(zip(trav,e['travel_ids'])):
        if segment[:,2].min()<0 or np.linalg.norm(segment[1]-segment[0])<w:continue
        line=LineString(segment[:,:2]);inds=tree.query(line.buffer(r+w/2),predicate='intersects')
        if not len(inds):continue
        # Conservative Z interval for sloping travel, exact temporal order.
        inds=inds[(ids[inds]<event_id)&(z[inds]>=segment[:,2].min()-tip-h)&(z[inds]<=segment[:,2].max()+cfg['travel_body_height_mm']+h)]
        if not len(inds):continue
        length=line.length
        if length>w*2:
            a=np.asarray(line.interpolate(w).coords[0]);b=np.asarray(line.interpolate(length-w).coords[0]);line=LineString([a,b])
        else:continue # departure/arrival contact is not a traversing near-miss
        distance=shapely.distance(lines[inds],line)-widths[inds]/2
        # Printed bead is [layer top-h, layer top]; use swept travel Z range.
        below=np.maximum(0.,segment[:,2].min()-z[inds]);above=np.maximum(0.,z[inds]-h-segment[:,2].max())
        tipdistance=np.hypot(np.maximum(distance,0),np.maximum(below,above));j=int(np.argmin(tipdistance));minimum=float(tipdistance[j]);idx=int(inds[j])
        bodybottom=segment[:,2].min()+cfg['travel_body_start_height_mm'];bodytop=segment[:,2].max()+cfg['travel_body_height_mm']
        vertical=(z[inds]>=bodybottom)&(z[inds]-h<=bodytop);bodydistance=float(np.min(distance[vertical]-r)) if vertical.any() else None
        risk=minimum<tip or (bodydistance is not None and bodydistance<0)
        for name in np.unique(labels[inds]):
            value=float(np.min(tipdistance[labels[inds]==name]));by_feature[str(name)]=min(value,by_feature.get(str(name),value))
        for name,value in [('nozzle_tip',minimum),('body',bodydistance)]:
            if value is not None:candidate_min[name]=value if candidate_min[name] is None else min(value,candidate_min[name])
        if risk:rows.append(dict(travel_id=i,gcode_event=int(event_id),start=segment[0].tolist(),end=segment[1].tolist(),nearest_feature=str(labels[idx]),nearest_extrusion_event=int(ids[idx]),nozzle_tip_min_clearance_mm=minimum,body_clearance_mm=bodydistance,status='TRAVEL_NEAR_MISS',chronology_verified=bool(ids[idx]<event_id)))
    result=dict(status='GEOMETRIC_TRAVEL_PROXY',travel_segments=len(trav),audited_extrusion_segments=len(ex),near_miss_segments=len(rows),minimum_clearance_mm=candidate_min,minimum_tip_clearance_by_feature_mm=by_feature,nozzle_tip_clearance_mm=tip,body_radius_mm=r,body_height_mm=cfg['travel_body_height_mm'],body_start_height_mm=cfg['travel_body_start_height_mm'],endpoint_exclusion_mm=w,method='Chronological actual extrusion, archived variable width, XY segment distance and conservative swept Z interval; no future material, no full nozzle-body collision proof',physical_collision_validated=False,coordinate_frame=e['coordinate_frame'],rows=rows)
    table(folder/'travel_near_miss.csv',rows,fields=None if rows else ['travel_id','status','start','end','nozzle_tip_min_clearance_mm']);save(cache,result)
    return result
