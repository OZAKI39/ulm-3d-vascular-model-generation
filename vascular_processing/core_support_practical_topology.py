"""Strategic material-graph cuts and actual C3 necks, within three rounds."""
from pathlib import Path
import json,pickle
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from shapely.geometry import Point
from .core_support_geometry import boolean,union,overlap
from .core_support_practical_geometry import material_components,write_candidate
from .core_support_closure import local_box
from .core_support_safe_removal import actual_surface_normal,lever_arm,tool_access,release_path,closed_cage,preserved_feature
from .core_support_removal import wrap_audit
from .support_inputs import save,table


def contact_measure(parts,g,exact=False):
    samples=np.vstack([b['lower'] for b in g['branches']]);arc=np.concatenate([b['arc'] for b in g['branches']]);names=np.concatenate([[b['branch_id']]*len(b['arc']) for b in g['branches']]);tree=cKDTree(samples);rows=[]
    reach=g['profile']['roof_gap_mm']+max(3*g['profile']['layer_height'],g['profile']['support_line_width'])
    for i,part in enumerate(parts):
        if exact:
            candidate=np.flatnonzero(np.all(samples>=part.bounds[0]-reach,axis=1)&np.all(samples<=part.bounds[1]+reach,axis=1));contacts=[]
            for batch in np.array_split(candidate,max(1,int(np.ceil(len(candidate)/64)))):
                if not len(batch):continue
                distance=part.nearest.on_surface(samples[batch])[1];contacts.extend(batch[distance<=reach].tolist())
            ids=np.array(contacts,dtype=int)
        else:
            ds,ix=tree.query(part.vertices);ids=np.unique(ix[ds<=reach])
        spans={str(n):float(np.ptp(arc[ids][names[ids]==n])) for n in np.unique(names[ids])}
        rows.append(dict(component_id=f'M{i:02}',volume_mm3=float(part.volume),bbox=part.bounds.tolist(),contact_arc_mm=sum(spans.values()),branch_contacts=spans,sample_indices=ids.tolist(),maximum_height_mm=float(part.extents[2])))
    return rows


def cut_effect(before,after,g):
    a=contact_measure(before,g);b=contact_measure(after,g)
    old=max((r['contact_arc_mm'] for r in a),default=0.);new=max((r['contact_arc_mm'] for r in b),default=0.)
    # Reject extra crumbs or a decorative notch that leaves the same large
    # vascular-contact beam. Significant reduction is one actual bead length.
    old_energy=sum(r['contact_arc_mm']**2 for r in a);new_energy=sum(r['contact_arc_mm']**2 for r in b)
    # Also permit a useful split of the second-largest group. A crumb or
    # decorative hole has no meaningful reduction of contact-span energy.
    effective=(len(after)>len(before) and (new<old-g['profile']['support_line_width'] or old_energy-new_energy>2*g['profile']['support_line_width']*old))
    return dict(effective=effective,before_components=len(before),after_components=len(after),old_max_contact_arc_mm=old,new_max_contact_arc_mm=new,old_contact_arc_squared=old_energy,new_contact_arc_squared=new_energy)


def strategic_candidates(support,g,round_number=1):
    """Read actual narrow cross-sections at existing seams and branch middles."""
    p=g['profile'];w=p['support_line_width'];h=p['layer_height'];v=g['vascular'];candidates=[];rejected=[]
    forks=np.unique(np.asarray([b['xyz'][0] for b in g['branches'] if not b['is_port']][1:]),axis=0)
    sources=[]
    for b in g['branches']:
        ids=[int(np.argmin(abs(b['arc']-n['arc_mm']))) for n in g['notches'] if n['branch_id']==b['branch_id']]
        ids+= [len(b['arc'])//2]
        if round_number>=2 and b['arc'][-1]>sum(x['arc'][-1] for x in g['branches'])/len(g['branches']):
            # Only bisect the long existing branch spans left by the first
            # graph cut. This is not an orientation or parameter sweep.
            ids += [len(b['arc'])//4,3*len(b['arc'])//4]
        if round_number==3 and b['arc'][-1]>sum(x['arc'][-1] for x in g['branches'])/len(g['branches']):
            ids += [len(b['arc'])//8,3*len(b['arc'])//8,5*len(b['arc'])//8,7*len(b['arc'])//8]
        for i in sorted(set(ids)):
            if i==0 or i==len(b['arc'])-1:continue
            tangent=b['tangent'][i,:2];norm=np.linalg.norm(tangent)
            if norm<1e-8:continue
            t=tangent/norm;side=np.array([-t[1],t[0]]);q=np.r_[b['xyz'][i,:2],b['roof_top'][i]-h]
            sources.append((b,i,q,t,side))
    for k,(b,i,q,t,side) in enumerate(sources):
        key=f'SC{k:02}';path=support.section(plane_origin=q,plane_normal=np.r_[t,0.])
        if path is None:continue
        flat,to3d=path.to_2D();qq=np.linalg.inv(to3d)@np.r_[q,1.];polys=list(flat.polygons_full)
        if not polys:continue
        poly=min(polys,key=lambda x:x.distance(Point(qq[:2])))
        if poly.distance(Point(qq[:2]))>2*w:continue
        coords=np.asarray(poly.exterior.coords);pts=np.c_[coords,np.zeros(len(coords)),np.ones(len(coords))]@to3d.T;pts=pts[:,:3]
        lateral=pts[:,:2]@side;lo,hi=lateral.min(),lateral.max();zlo,zhi=pts[:,2].min(),pts[:,2].max()
        center_xy=q[:2]+side*((lo+hi)/2-q[:2]@side)
        # A finite blade envelope. If it intersects the frozen vessel, the
        # proposed cut is rejected; clipping the blade to support would hide it.
        cutter=local_box(center_xy,zhi+.01*h,hi-lo+2*w,w,t,zlo-h)
        iv=overlap(cutter,v)
        cm=trimesh.collision.CollisionManager();cm.add_object('v',v);distance=float(cm.min_distance_single(cutter))
        row=dict(cut_id=key,branch=b['branch_id'],arc_mm=float(b['arc'][i]),point=q.tolist(),tangent=np.r_[t,0.].tolist(),section_area_mm2=float(poly.area),cut_kerf_mm=w,lateral_width_mm=float(hi-lo+2*w),cut_top_mm=float(zhi+.01*h),cut_bottom_mm=float(zlo-h),cut_center_xy=center_xy.tolist(),distance_to_bifurcation_mm=float(np.min(np.linalg.norm(forks-q,axis=1))) if len(forks) else None,nearest_vascular_surface_mm=distance,vascular_intersection_mm3=iv)
        if iv>=1e-7 or distance<=0:rejected.append({**row,'reason':'BLADE_ENVELOPE_INTERSECTS_VASCULAR'});continue
        kerf=boolean([support,cutter],'intersection')
        if not len(kerf.faces):rejected.append({**row,'reason':'NO_MATERIAL_CONNECTION'});continue
        # Actual material point, preferentially away from vessel, for the
        # bulk-cut approach. The complete finite blade is checked above.
        distances=cKDTree(v.vertices).query(kerf.vertices)[0];point=kerf.vertices[int(np.argmax(distances))]
        clearance=float(v.nearest.on_surface([point])[1][0])
        row.update(bulk_cut_grip_point=point.tolist(),bulk_point_vascular_clearance_mm=clearance)
        candidates.append(dict(row=row,cutter=cutter,kerf=kerf))
    return candidates,rejected


def select_strategic_cuts(g,round_number=1):
    if round_number not in [1,2,3]:raise ValueError('THREE_TOPOLOGY_ROUNDS_MAXIMUM')
    support=g['support'];parts=material_components(support);before=contact_measure(parts,g);initial=max(r['contact_arc_mm'] for r in before)
    candidates,rejected=strategic_candidates(support,g,round_number);selected=[];remaining=support.copy()
    print('GRAPH_CUT_SECTIONS','eligible',len(candidates),'rejected',len(rejected),flush=True)
    anatomical=sum(not b['is_port'] for b in g['branches']);target=initial/max(2,anatomical+round_number-1)
    # A graph cut is scored by its effect, not by a desired number of chunks.
    while candidates and max(r['contact_arc_mm'] for r in contact_measure(parts,g))>target:
        options=[]
        for candidate in candidates:
            after=boolean([remaining,candidate['cutter']],'difference');newparts=material_components(after);effect=cut_effect(parts,newparts,g)
            if round_number>=2 and effect['new_max_contact_arc_mm']>=effect['old_max_contact_arc_mm']-g['profile']['support_line_width'] and effect['old_contact_arc_squared']-effect['new_contact_arc_squared']<.02*effect['old_contact_arc_squared']:
                # Avoid adding cuts to already-small pieces while leaving the
                # dominant beam intact. Component count is not the objective.
                continue
            if effect['effective']:
                options.append((effect['new_max_contact_arc_mm'],effect['new_contact_arc_squared'],candidate['row']['section_area_mm2'],candidate,after,newparts,effect))
        if not options:break
        _,_,_,item,remaining,parts,effect=min(options,key=lambda x:x[:3]);item['row'].update(effect);selected.append(item);candidates.remove(item)
        print('GRAPH_CUT_ACCEPT',item['row']['cut_id'],effect,flush=True)
    for c in candidates:rejected.append({**c['row'],'reason':'INEFFECTIVE_OR_NOT_NEEDED_FOR_CONTACT_SPAN_REDUCTION'})
    return selected,rejected,parts,dict(initial_components=len(before),initial_max_contact_arc_mm=initial,target_max_contact_arc_mm=target,after_components=len(parts),after_max_contact_arc_mm=max(r['contact_arc_mm'] for r in contact_measure(parts,g)))


def topology_candidate(parent,folder,round_number=1):
    parent=Path(parent);folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);g=pickle.load((parent/'geometry.pkl').open('rb'));p=g['profile'];w=p['support_line_width'];h=p['layer_height'];v=g['vascular'];support=g['support'];selected,rejected,virtual,graph=select_strategic_cuts(g,round_number)
    additions=[];notches=[];links=[]
    for item in selected:
        row=item['row'];t=np.asarray(row['tangent'])[:2];xy=np.asarray(row['cut_center_xy']);roof_bottom=row['cut_top_mm']-.01*h-3*h;web_bottom=roof_bottom-3*h
        # R3 clears the complete intended neck band before inserting C3.
        # Earlier diagnostics left parent curtain material beside the web;
        # slicer preservation alone must not mislabel that as an isolated C3.
        notch=local_box(xy,roof_bottom if round_number==3 else web_bottom,row['lateral_width_mm'],2*w,t,max(0.,row['cut_bottom_mm']))
        web=local_box(xy,roof_bottom,3*w,6*w,t,web_bottom)
        if overlap(notch,v)>=1e-7 or overlap(web,v)>=1e-7:row['thin_web_status']='GEOMETRIC_C3_PLACEMENT_REJECTED';continue
        notches.append(notch);additions.append(web)
        levels=np.arange(p['first_layer_height'],roof_bottom+h,h);levels=levels[(levels>=web_bottom+h/2)&(levels<=roof_bottom+h/2)]
        links.append(dict(link_id=row['cut_id'],width_mm=3*w,thickness_mm=3*h,source_fixture='actual_thin_web_C3',center=[*xy,(roof_bottom+web_bottom)/2],witness_points=[[*xy,float(z)] for z in levels],neck_band_bottom_mm=web_bottom,neck_band_top_mm=roof_bottom,parent_curtain_cleared=round_number==3))
        row['thin_web_status']='C3_PENDING_ACTUAL_SLICE'
    changed=boolean([support,union(notches)],'difference') if notches else support.copy()
    if additions:changed=union([changed,*additions])
    metadata=dict(stage='STRATEGIC_TOPOLOGY',round=round_number,strategic_cuts=[x['row'] for x in selected],rejected_cuts=rejected,graph=graph,weak_links=links,roof_entirely_regenerated=False)
    r=write_candidate(parent,folder,changed,metadata)
    save(folder/'strategic_cut_plan.json',metadata);table(folder/'strategic_cut_seams.csv',metadata['strategic_cuts'])
    with (folder/'strategic_cut_solids.pkl').open('wb') as f:pickle.dump(selected,f)
    print('STRATEGIC_TOPOLOGY',round_number,'cuts',len(selected),'graph',graph,flush=True)
    return r
