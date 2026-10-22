"""Open-space cut access, angular wrap and swept fragment release proxies."""
from __future__ import annotations
import pickle
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from .support_inputs import save,table
from .core_support_geometry import union,boolean,overlap
from .support_removal import swept_path


def directions():
    d=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(1,1,0),(1,-1,0),(-1,1,0),(-1,-1,0),(1,0,1),(-1,0,1),(0,1,1),(0,-1,1)]
    return [np.array(v)/np.linalg.norm(v) for v in d]


def notch_access(notch,vascular,support,allow_lower=True):
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular',vascular);manager.add_object('support',support);target=np.asarray(notch['point']);side=np.asarray(notch['lateral']);attempts=[]
    for direction in [side,-side,*directions()]:
        # Tool axis approaches a notch side. The nose stays outside the
        # material; its 3 mm diameter must have a clear full approach.
        # The centre of the open slit is air. Probe the two actual cut lips,
        # not a ray that incorrectly requires material in the empty notch.
        for sign in [-1,1]:
            lip=target+np.asarray(notch['tangent'])*sign*(notch.get('cut_width_mm',.44)/2+.05)
            origins=np.array([lip+direction*200]);hits,_,_=support.ray.intersects_location(origins,[-direction],multiple_hits=True)
            if not len(hits):continue
            distance=(hits-lip)@direction;near=hits[np.argmin(np.abs(distance))]
            if np.linalg.norm(near-lip)>.8:continue
            end=near+direction*.05;start=end+direction*200
            matrix=trimesh.geometry.align_vectors([0,0,1],direction);matrix[:3,3]=(start+end)/2
            probe=trimesh.creation.cylinder(radius=1.5,height=200,sections=32,transform=matrix)
            hit=bool(manager.in_collision_single(probe));attempts.append(dict(direction=direction.tolist(),lip_sign=sign,end=end.tolist(),collision=hit))
            if not hit:return dict(notch_id=notch['notch_id'],passed=True,point=target,cut_lip=near,direction=direction,probe_diameter_mm=3.,attempts=attempts,claim='GEOMETRIC_TOOL_ACCESS_PROXY')
    if allow_lower and target[2]>1.6:
        # The same open body notch continues down to the bed. After lifting
        # the print off the plate, its lower lips are also valid cut access.
        lower=dict(notch);lower['point']=np.r_[target[:2],1.5]
        retry=notch_access(lower,vascular,support,allow_lower=False)
        if retry['passed']:retry['access_strategy']='LOWER_LIP_AFTER_BED_RELEASE';return retry
        attempts.extend(retry['attempts'])
    return dict(notch_id=notch['notch_id'],passed=False,point=target,probe_diameter_mm=3.,attempts=attempts,claim='GEOMETRIC_TOOL_ACCESS_PROXY')


def wrap_audit(vascular,support,branches,profile):
    """Angular extent of support near each local vessel cross-section.

    Only a finite radial band adjacent to the actual lower surface is counted;
    the long bed-rooted wall is not a fictitious encircling contact patch.
    """
    sv=support.vertices;tree=cKDTree(sv);rows=[]
    for b in branches:
        for i in range(0,len(b['xyz']),2):
            q=b['xyz'][i];t=b['tangent'][i];low=b['lower'][i]
            if not np.isfinite(low).all():continue
            v=low-q;down=v-t*(v@t);norm=np.linalg.norm(down)
            if norm<1e-8:continue
            down/=norm;side=np.cross(t,down);radius=np.linalg.norm(v);ids=tree.query_ball_point(q,radius+profile['roof_gap_mm']+profile['line_width']+.15)
            points=sv[ids]-q if ids else np.empty((0,3));along=points@t;planar=points-along[:,None]*t;rad=np.linalg.norm(planar,axis=1)
            select=(np.abs(along)<=profile['line_width']/2)&(rad<=radius+profile['roof_gap_mm']+profile['line_width']+.05)&(rad>1e-8)
            pts=planar[select]
            if len(pts):
                a=np.degrees(np.arctan2(pts@side,pts@down));angle=float(a.max()-a.min())
            else:angle=0.
            rows.append(dict(branch_id=b['branch_id'],sample_index=i,angular_coverage_deg=angle,pass_120=angle<=120,preferred_below_90=angle<90))
    angles=[r['angular_coverage_deg'] for r in rows]
    return dict(passed=bool(rows) and all(r['pass_120'] for r in rows),wrap_traps=sum(not r['pass_120'] for r in rows),maximum_wrap_angle_deg=max(angles,default=0),wrap_angle_p95=float(np.percentile(angles,95)) if angles else 0,rows=rows,method='Near-surface support angular coverage projected onto SWC tangent-normal cross-sections; 120 degree hard bound, not physical release force')


def audit(record):
    folder=Path(record['folder']);g=pickle.load((folder/'geometry.pkl').open('rb'));vascular=g['vascular'];support=g['support'];notches=g['notches'];cutters=[];cutrows=[]
    for n in notches:
        t=np.asarray(n['tangent']);l=np.asarray(n['lateral']);z=n['plane_origin'][2]+.46;center=np.r_[n['point'][:2],(z-1)/2]
        m=np.eye(4);m[:3,0]=t;m[:3,1]=l;m[:3,2]=[0,0,1];m[:3,3]=center
        cutter=trimesh.creation.box([.10,6.,z+1],transform=m);hit=overlap(cutter,vascular)
        cutrows.append(dict(notch_id=n['notch_id'],vascular_cut_intersection_mm3=hit,passed=hit<1e-7));cutters.append(cutter)
    remaining=boolean([support,union(cutters)],'difference') if cutters else support
    fragments=[m for m in remaining.split(only_watertight=False) if m.volume>1e-7];access=[notch_access(n,vascular,support) for n in notches]
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular',vascular)
    for i,f in enumerate(fragments):manager.add_object(f'F{i:03}',f)
    center=vascular.bounds.mean(0);bif=np.array([b['xyz'][0] for b in g['branches'] if not b['is_port']][1:]).reshape(-1,3)
    near_bif=[bool(len(bif) and cKDTree(f.vertices).query(bif)[0].min()<2.) for f in fragments]
    order=sorted(range(len(fragments)),key=lambda i:(near_bif[i],-np.linalg.norm(fragments[i].centroid[:2]-center[:2])));rows=[];sequence=[]
    for i in order:
        part=fragments[i];name=f'F{i:03}';manager.remove_object(name);pivot=part.centroid;local=part.copy();local.apply_translation(-pivot);initial=np.eye(4);initial[:3,3]=pivot;attempts=[];accepted=None
        for direction in [np.array([0.,0.,-1.]),*directions()]:
            final=initial.copy();final[:3,3]+=direction*float(np.linalg.norm(vascular.extents)+np.linalg.norm(part.extents)+5)
            check=swept_path(local,[],[initial,final],static_manager=manager);attempts.append(dict(direction=direction,**check))
            if check['passed']:accepted=dict(direction=direction,transforms=[initial,final]);break
        if accepted is None:
            from scipy.spatial.transform import Rotation
            for axis in ['x','y']:
                for angle in [-10,10]:
                    lowered=initial.copy();lowered[2,3]-=.5
                    tilted=lowered.copy();tilted[:3,:3]=Rotation.from_euler(axis,angle,degrees=True).as_matrix();final=tilted.copy();final[2,3]-=float(np.linalg.norm(vascular.extents)+np.linalg.norm(part.extents)+5)
                    check=swept_path(local,[],[initial,lowered,tilted,final],static_manager=manager);attempts.append(dict(rotation_axis=axis,rotation_deg=angle,**check))
                    if check['passed']:accepted=dict(direction=[0,0,-1],rotation_axis=axis,rotation_deg=angle,transforms=[initial,lowered,tilted,final]);break
                if accepted is not None:break
        if accepted:sequence.append(name)
        else:manager.add_object(name,part)
        rows.append(dict(fragment_id=name,passed=accepted is not None,near_bifurcation=near_bif[i],volume_mm3=float(part.volume),centroid=part.centroid,path=accepted,attempts=attempts,suggested_order=len(sequence) if accepted else None))
    (folder/'removal_fragments').mkdir(exist_ok=True)
    for i,f in enumerate(fragments):f.export(folder/'removal_fragments'/f'F{i:03}.stl')
    wrap=wrap_audit(vascular,support,g['branches'],g['profile']);save(folder/'wrap_audit.json',wrap)
    lengths=[s['length_mm'] for s in g['segments']];result=dict(passed=all(a['passed'] for a in access) and all(r['passed'] for r in rows) and all(r['passed'] for r in cutrows) and wrap['passed'],segment_count=len(lengths),mean_segment_length_mm=float(np.mean(lengths)),median_segment_length_mm=float(np.median(lengths)),max_segment_length_mm=max(lengths),cut_notches=len(notches),accessible_notches=sum(a['passed'] for a in access),accessible_notch_fraction=sum(a['passed'] for a in access)/len(access) if access else 1.,fragment_count=len(fragments),removable_segments=sum(r['passed'] for r in rows),removable_segment_fraction=sum(r['passed'] for r in rows)/len(rows) if rows else 0.,wrap_traps=wrap['wrap_traps'],wrap_angle_p95=wrap['wrap_angle_p95'],maximum_wrap_angle_deg=wrap['maximum_wrap_angle_deg'],sequence=sequence,cut_checks=cutrows,access=access,fragments=rows,physical_validation=False,method='Finite virtual cuts; 3 mm cylinder approach; full-fragment FCL bounded swept translation and optional small rotation against core and all remaining fragments')
    save(folder/'removal.json',result);return result
