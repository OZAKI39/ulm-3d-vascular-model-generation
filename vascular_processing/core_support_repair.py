"""Explicit localized geometry corrections driven by saved actual deficits."""
from pathlib import Path
import copy,json,pickle
import numpy as np
import trimesh,shapely
from shapely.geometry import Point
from .support_inputs import save
from .core_support_geometry import union,overlap,base_foot,boolean,bed_roots
from .support_rfus_geometry import lower_points
from .support_v2_mesh import support_stl
from .support_v3_geometry import window_cutter


def repair(source_folder,out):
    source_folder=Path(source_folder);out=Path(out);out.mkdir(parents=True,exist_ok=True);g=pickle.load((source_folder/'geometry.pkl').open('rb'));p=g['profile'];vascular=g['vascular'];w=p['support_line_width'];h=p['layer_height'];first=p['first_layer_height'];gap=p['roof_gap_mm'];tiles=[];records=[];blocked=[]
    cov=json.loads((source_folder/'coverage.json').read_text());layer=json.loads((source_folder/'layer_coverage.json').read_text());selfcheck=json.loads((source_folder/'support_self.json').read_text())
    def add(xy,top,kind,identity,width=3):
        if top<=0:return
        size=width*w;m=trimesh.creation.box([size,size,top+.049]);m.apply_translation([*xy,(top+.049)/2]);v=overlap(m,vascular)
        row=dict(kind=kind,target_id=identity,xy=np.asarray(xy),top_layer_mm=top,roof_width_lines=width,intersection_mm3=v)
        if v>=1e-7:blocked.append(row);return
        tiles.append(m);records.append(row)
    for r in cov['rows']:
        if r['supported'] or r['lower_surface'] is None:continue
        q=np.asarray(r['lower_surface']);top=first+np.floor((q[2]-gap-first)/h)*h;add(q[:2],float(top),'LOWER_ENVELOPE',f"{r['branch_id']}:{r['sample_index']}")
    for ri,r in enumerate(layer['rows']):
        if r['passed']:continue
        poly=shapely.from_wkt(r['uncovered_wkt']);lo=np.array(poly.bounds[:2]);hi=np.array(poly.bounds[2:]);xs=np.arange(lo[0],hi[0]+w/2,w/2);ys=np.arange(lo[1],hi[1]+w/2,w/2)
        grid=np.array([[x,y] for x in xs for y in ys if poly.buffer(w/4).covers(Point(x,y))])
        if not len(grid):grid=np.asarray(poly.representative_point().coords)
        xyz=np.c_[grid,np.full(len(grid),r['z_mm']-h/2+1e-5)];low=lower_points(vascular,xyz)
        for i,q in enumerate(low):
            if np.isfinite(q).all():add(q[:2],float(first+np.floor((q[2]-gap-first)/h)*h),'LAYER_DOMAIN',f'{ri}:{i}',2)
    for ri,r in enumerate(selfcheck['rows']):
        poly=shapely.from_wkt(r['wkt']);pt=np.asarray(poly.representative_point().coords[0]);add(pt,float(r['z_mm']-h),'SUPPORT_PREVIOUS_LAYER',str(ri),3)
    support=union([g['support'],*tiles]);support=base_foot(support,p,expansion=4*w)
    cuts=[]
    for n in g['notches']:
        cut,_=window_cutter(n['plane_origin'],np.asarray(n['tangent'])[:2],n['cut_width_mm'],-1,n['plane_origin'][2],6.,slope=.45);cuts.append(cut)
    if cuts:support=boolean([support,union(cuts)],'difference')
    support=support_stl(support);shifted=support.copy();shifted.apply_translation([0,0,-3*h]);roofs=boolean([support,shifted],'difference');curtains=boolean([support,roofs],'difference');g.update(support=support,roofs=roofs,curtains=curtains)
    for name,mesh in [('vascular_core',vascular),('custom_support',support),('support_roofs',roofs),('support_curtains',curtains)]:mesh.export(out/(name+'.stl'))
    trimesh.util.concatenate([vascular,support]).export(out/'vascular_core_plus_support_review.stl')
    for name in ['vascular_print_transform.json','inverse_transform.json','notches.json']:(out/name).write_bytes((source_folder/name).read_bytes())
    with (out/'geometry.pkl').open('wb') as f:pickle.dump(g,f)
    r=json.loads((source_folder/'geometry.json').read_text());r.update(candidate_id=out.name,folder=str(out),repair_parent=source_folder.name,support_volume_mm3=float(support.volume),intersection_mm3=overlap(support,vascular),correction_count=len(records),blocked_corrections=len(blocked),cut_notches_reapplied=True)
    r['bed_roots']=bed_roots(support)
    save(out/'geometry.json',r);save(out/'actual_deficit_repair.json',dict(accepted=records,blocked=blocked,geometry_only_changed='CUSTOM_SUPPORT',coverage_tolerance_unchanged=True));print('REPAIR',out.name,len(records),'blocked',len(blocked),'intersection',r['intersection_mm3'],flush=True);return r
