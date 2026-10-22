"""V3 ribs retained; add separately gapped inverted-box ceiling curtains."""
from __future__ import annotations
import pickle, json, copy
from pathlib import Path
import numpy as np
import trimesh
import shapely
from shapely.geometry import Point
from .support_inputs import save, table
from .support_detection import bottom_points, contiguous_runs
from .support_v3_geometry import curtain, union, cut_section
from .support_v3_planner import plan
from .support_v2_layers import analyse
from .support_v2_mesh import support_stl, exact_stl
from .support_custom_qc import intersection_volume


def box_ribs(boxfield, box, vascular, profile):
    """Vertical ruled curtains sampled only inside frozen REQUIRED polygons.

    No union to the target ceiling. Rows are separated by one roof width;
    clear rays are only a preliminary filter, and every whole rib is checked
    against blood and box. Incomplete ceiling coverage remains a hard fail.
    """
    w,h=profile['line_width'],profile['layer_height']; pitch=9*w
    points={}
    for row,poly in zip(boxfield['rows'],boxfield['polygons']):
        if row['classification']!='SUPPORT_REQUIRED':continue
        x0,y0,x1,y1=poly.bounds
        for iy in range(int(np.floor(y0/pitch)),int(np.ceil(y1/pitch))+1):
            y=iy*pitch
            xs=np.arange(np.floor(x0/(2*w)),np.ceil(x1/(2*w))+1)*2*w
            inside=shapely.contains_xy(poly,xs,np.full(len(xs),y))
            for x in xs[inside]:
                key=(iy,round(float(x),6));z=float(row['z_mm'])
                points[key]=min(points.get(key,z),z)
    ribs=[]; rejected=[]
    for iy in sorted({p[0] for p in points}):
        line=sorted((x,z) for (j,x),z in points.items() if j==iy)
        xyz=np.array([[x,iy*pitch,z] for x,z in line]);
        if len(xyz)<2:continue
        tops=np.floor((xyz[:,2]-profile['roof_gap_mm'])/h)*h
        roots=bottom_points(box,np.c_[xyz[:,:2],tops-.1])
        bases=np.where(np.isfinite(roots[:,2]),roots[:,2],0.)
        # First priority is the plate; a box hit is a separately audited root.
        active=tops-bases>3.
        for j in np.flatnonzero((np.diff(xyz[:,0])>2*w+.001)|(abs(np.diff(bases))>h+.001)):
            active[j]=False;active[j+1]=False
        for number,(a,b) in enumerate(contiguous_runs(active)):
            if b<=a:continue
            rid=f'CEILING_Y{iy}_{number}'
            try:
                rib=curtain(xyz[a:b+1,:2],tops[a:b+1],bases[a:b+1],profile,3,pitch,5.,3.)
                hit=intersection_volume(rib['mesh'],vascular)
                allowed=rib['feet']; bh=trimesh.boolean.intersection([rib['mesh'],box],engine='manifold')
                illegal=trimesh.boolean.difference([bh,allowed],engine='manifold') if len(bh.faces) else bh
                if hit>1e-5 or (len(illegal.faces) and abs(illegal.volume)>1e-4):
                    # Split a blocked curtain into short local runs and test
                    # them independently; do not silently trim frozen anatomy.
                    subruns=[(i,min(i+5,b)) for i in range(a,b,5)]
                    for k,(c,d) in enumerate(subruns):
                        if d<=c:continue
                        sub=curtain(xyz[c:d+1,:2],tops[c:d+1],bases[c:d+1],profile,3,pitch,5.,3.)
                        overlap=intersection_volume(sub['mesh'],vascular)
                        bh=trimesh.boolean.intersection([sub['mesh'],box],engine='manifold')
                        extra=trimesh.boolean.difference([bh,sub['feet']],engine='manifold') if len(bh.faces) else bh
                        if overlap>1e-5 or (len(extra.faces) and abs(extra.volume)>1e-4):
                            rejected.append(dict(rib_id=f'{rid}_{k}',reason='CEILING_PATH_BLOCKED_BY_FROZEN_MODEL'));continue
                        sub.update(rib_id=f'{rid}_{k}',interval_id=rid,branch_id='BOX',is_port=False,component='BOX_REQUIRED')
                        ribs.append(sub)
                    continue
                rib.update(rib_id=rid,interval_id=rid,branch_id='BOX',is_port=False,component='BOX_REQUIRED');ribs.append(rib)
            except (ValueError,RuntimeError) as exc:rejected.append(dict(rib_id=rid,reason=str(exc)))
    return ribs,rejected


def write_geometry(record,ribs,profile,baseline_ribs=None):
    folder=Path(record['folder'])
    record.pop('post_export_reverted_ribs',None)
    if baseline_ribs is not None:
        if {r['rib_id'] for r in baseline_ribs}!={r['rib_id'] for r in ribs}:raise ValueError('BASELINE_RIB_ID_SET_CHANGED')
        record['baseline_ribs_preserved']=True
    try:support=support_stl(union([r['mesh'] for r in ribs]))
    except ValueError:
        if baseline_ribs is not None:
            # Phase B may reject a NEW edit, but must never discard an
            # unchanged baseline rib to make a problematic union exportable.
            baseline_by_id={r['rib_id']:r for r in baseline_ribs}
            current=[copy.deepcopy(r) for r in baseline_ribs]
            support=support_stl(union([r['mesh'] for r in current]));rejections=[]
            for candidate in ribs:
                i=next(j for j,r in enumerate(current) if r['rib_id']==candidate['rib_id'])
                proposed=current.copy();proposed[i]=candidate
                try:new=support_stl(union([r['mesh'] for r in proposed]))
                except ValueError as exc:
                    rejections.append(dict(rib_id=candidate['rib_id'],reason='NEW_EDIT_REVERTED_FOR_EXPORT_PRECISION',detail=str(exc)));continue
                current=proposed;support=new
            ribs=current;record['post_export_reverted_ribs']=rejections;record['baseline_ribs_preserved']=True
            table(folder/'post_export_reverted_edits.csv',rejections)
        else:
        # A failing local support seam must not terminate the orientation
        # experiment. Retain exportable curtains, explicitly reject the rest;
        # their absent coverage remains visible in the actual hard gates.
            accepted=[];rejections=[];support=trimesh.Trimesh()
            for rib in ribs:
                try:new=support_stl(union([support,rib['mesh']]))
                except ValueError as exc:
                    rejections.append(dict(rib_id=rib['rib_id'],reason=str(exc)));continue
                accepted.append(rib);support=new
            if not accepted:raise ValueError('NO_EXPORTABLE_V3_1_SUPPORT')
            ribs=accepted;table(folder/'support_precision_rejections.csv',rejections)
            record['support_precision_rejections']=rejections
    model=trimesh.load_mesh(folder/'mold_without_support.stl');v=trimesh.load_mesh(folder/'vascular_with_ports.stl');box=trimesh.load_mesh(folder/'box.stl')
    raw=union([model,support]);error=None
    try:assembly=exact_stl(raw)
    except ValueError as exc:assembly=raw;error=str(exc);raw.export(folder/'assembly_indexed.ply')
    support.export(folder/'custom_support_v3.stl');support.export(folder/'custom_support.stl');assembly.export(folder/'mold_with_custom_support_v3.stl')
    for name,key in [('support_ribs.stl','ribs'),('support_roofs.stl','roofs'),('support_cut_necks.stl','necks')]:union([r[key] for r in ribs]).export(folder/name)
    cuts=[];intended=union([c['intended'] for r in ribs for c in r['cuts']])
    for rib in ribs:
        # Plate-rooted legs travel out via plate release; their old V3 notch
        # geometry stays printable but is not claimed as a required cut.
        rib['root_domain']='BUILD_PLATE' if np.max(rib['base'])<.001 else 'SACRIFICIAL_BOX_SURFACE'
        for j,c in enumerate(rib['cuts']):
            c['cut_id']=f"{rib['rib_id']}_C{j:03}"
            check=cut_section(support,intended,c['origin'],c['normal'],c['aperture'])
            cuts.append(dict(cut_id=c['cut_id'],rib_id=rib['rib_id'],kind=c['kind'],origin=np.asarray(c['origin']).tolist(),**check))
    with (folder/'rib_geometry_cache.pkl').open('wb') as f:pickle.dump(ribs,f)
    table(folder/'cut_points.csv',cuts)
    sf=analyse(support,profile,'support',box,folder,bridge_limit=0.)
    sr=[r for r in sf['rows'] if r['classification']=='SUPPORT_REQUIRED']
    record.update(support_volume_mm3=float(support.volume),component_count=len(support.split()),
                  vascular_intersection_mm3=intersection_volume(support,v),assembly_export_error=error,
                  self_buildability=dict(passed=not sr,unsupported_region_count=len(sr),rows=sr,bridge_credit=False,area_tolerance_mm2=.005),
                  cuts_valid=all(r['passed'] for r in cuts),cut_count=len(cuts),
                  rib_segments=sum(r['segment_count'] for r in ribs),roof_ridges=len(ribs),
                  box_ceiling_rib_count=sum(r.get('component')=='BOX_REQUIRED' for r in ribs))
    from .support_v2_orientation import fits
    record['build_volume_fit']=fits(assembly.vertices,np.eye(4),profile['build_bounds_mm'])
    save(folder/'candidate.json',record)
    return record


def plan_inverted(folder,row,inputs,profile,cfg,mold,v,box,branches,topology):
    record=plan(folder,row,inputs,profile,cfg,mold,v,box,branches,topology,3)
    with (folder/'frozen_detection_cache.pkl').open('rb') as f:_,_,bf,_,_=pickle.load(f)
    with (folder/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
    extra,rejected=box_ribs(bf,trimesh.load_mesh(folder/'box.stl'),trimesh.load_mesh(folder/'vascular_with_ports.stl'),profile)
    table(folder/'box_support_rejections.csv',rejected)
    record['box_ceiling_planning']=dict(accepted=len(extra),rejected=len(rejected),gap_mm=profile['roof_gap_mm'],target_union_allowed=False)
    return write_geometry(record,ribs+extra,profile)
