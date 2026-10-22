"""Chronological nozzle-tip and cylinder proxy on actual non-extrusion moves."""
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import LineString
from .support_rfus_travel import events
from .support_inputs import save


def audit(archive,out,line_width=.22):
    e=events(archive);mask=np.array([not any(s in t.lower() for s in ['prime','wipe','custom','brim','skirt']) for t in e['features']]);ex=e['extrusion'][mask];ids=e['extrusion_ids'][mask];width=e['extrusion_widths'][mask];height=e['extrusion_heights'][mask];z=ex[:,:,2].mean(1);lines=shapely.linestrings(ex[:,:,:2]);tree=shapely.STRtree(lines);rows=[];tip=.2;radius=1.5
    for i,(seg,event) in enumerate(zip(e['travel'],e['travel_ids'])):
        line=LineString(seg[:,:2]);length=line.length
        if length<=2*line_width or seg[:,2].min()<0:continue
        # Departure/arrival are intentional print contacts, not crossings.
        line=LineString([line.interpolate(line_width),line.interpolate(length-line_width)])
        inds=tree.query(line.buffer(radius+width.max()/2),predicate='intersects')
        inds=inds[(ids[inds]<event)&(z[inds]>=seg[:,2].min()-tip)&(z[inds]-height[inds]<=seg[:,2].max()+3)]
        if not len(inds):continue
        dxy=shapely.distance(lines[inds],line)-width[inds]/2;below=np.maximum(0,seg[:,2].min()-z[inds]);above=np.maximum(0,z[inds]-height[inds]-seg[:,2].max());distance=np.hypot(np.maximum(0,dxy),np.maximum(below,above));best=int(np.argmin(distance));minimum=float(distance[best]);body=(z[inds]>=seg[:,2].min()+.2)&(z[inds]-height[inds]<=seg[:,2].max()+3);bodygap=float(np.min(dxy[body]-radius)) if body.any() else None
        if minimum<tip or (bodygap is not None and bodygap<0):rows.append(dict(travel_id=i,event_id=int(event),start=seg[0],end=seg[1],nearest_extrusion_event=int(ids[inds[best]]),tip_clearance_mm=minimum,body_clearance_mm=bodygap))
    result=dict(status='GEOMETRIC_TRAVEL_PROXY',travel_segments=len(e['travel']),near_miss_segments=len(rows),tip_clearance_mm=tip,body_radius_mm=radius,body_bottom_above_tip_mm=.2,body_height_mm=3.,endpoint_exclusion_mm=line_width,chronological=True,physical_collision_validated=False,rows=rows,suggested_profile_change='Review avoidance of crossing walls and travel/retraction settings on a profile copy; this run did not change them')
    save(Path(out)/'travel_near_miss.json',result);return result
