"""Bounded, actual-G-code driven local closure of the existing s1-11 support.

No global orientation search and no source/architecture regeneration. Removal
geometry is locked until all three coverage domains and self-support pass.
"""
from __future__ import annotations
import copy,json,pickle,shutil
import io
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import Point
import trimesh
from .core_support_inputs import default_output,discover_source,freeze,current_profiles
from .core_support_geometry import union,boolean,overlap,bed_roots
from .core_support_slicing import evaluate
from .support_v2_mesh import support_stl
from .support_rfus_geometry import lower_points
from .support_inputs import save,table,digest,verify,ROOT


def canonical_support(mesh):
    """Keep double-precision joins via ASCII STL if binary STL collapses them.

    This changes encoding, not dimensions or coverage tolerances. Both routes
    are round-tripped and required to remain closed positive volumes.
    """
    try:return support_stl(mesh),'stl'
    except ValueError as error:
        raw=mesh.export(file_type='stl_ascii').encode()
        restored=trimesh.load_mesh(io.BytesIO(raw),file_type='stl',process=False)
        restored.merge_vertices(digits_vertex=14)
        valid=restored.area_faces>0
        restored.update_faces(valid);restored.remove_unreferenced_vertices()
        if not restored.is_volume:raise ValueError('SUPPORT_ASCII_PRECISION_ROUNDTRIP_FAILED') from error
        return restored,'stl_ascii'


def coverage_pass(e):
    c=e.get('coverage',{});l=e.get('layer_coverage',{})
    return (c.get('centerline_underbody_percent')==100. and
            c.get('port_underbody_percent')==100. and c.get('max_uncovered_arc_mm')==0.
            and c.get('passed',False) and l.get('passed',False)
            and l.get('percent')==100.)


def closure_pass(e):
    return coverage_pass(e) and e.get('support_self',{}).get('self_unsupported_regions')==0


def loss(e):
    if 'coverage' not in e:return (10**9,float('inf'),10**9)
    return (e['coverage']['unsupported_samples'],e['layer_coverage']['uncovered_area_sum_mm2'],e['support_self']['self_unsupported_regions'])


def diagnostic_rank(e):
    """Keep all three required coverage domains visible in diagnostic ranking."""
    deficits=[100.-e['coverage']['centerline_underbody_percent'],100.-e['coverage']['port_underbody_percent'],100.-e['layer_coverage']['percent']]
    return (not closure_pass(e),max(deficits),sum(deficits),e['support_self']['self_unsupported_regions'])


def stagnated(history):
    """Two consecutive actual comparisons with neither coverage metric improved."""
    if len(history)<3:return False
    return all(b[0]>=a[0] and b[1]>=a[1] for a,b in zip(history[-3:-1],history[-2:]))


def choose_production(candidates):
    valid=[x for x in candidates if closure_pass(x[1]) and not x[2]]
    if not valid:return None
    # Removal/material are considered only after exact, real coverage gates.
    return min(valid,key=lambda x:(x[0].get('safe_removal_risk',0),x[0]['support_volume_mm3']))


def quantized_top(surface_z,gap,profile):
    h=profile['layer_height'];first=profile['first_layer_height']
    return float(first+np.floor((surface_z-gap-first)/h)*h)


def local_box(xy,top,width,length,tangent=(1.,0.),bottom=0.):
    t=np.asarray(tangent,dtype=float);t/=max(np.linalg.norm(t),1e-12)
    if np.linalg.norm(t)<.5:t=np.array([1.,0.])
    M=np.eye(4);M[:2,0]=t;M[:2,1]=[-t[1],t[0]];M[:3,3]=[*xy,(top+bottom)/2]
    return trimesh.creation.box([length,width,top-bottom],transform=M)


def repair_local(parent,out,iteration):
    """Repair only measured deficits; every modification has target provenance."""
    parent=Path(parent);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    g=pickle.load((parent/'geometry.pkl').open('rb'));p=g['profile'];w=p['support_line_width'];h=p['layer_height'];gap=p['roof_gap_mm'];v=g['vascular'];s=g['support']
    cov=json.loads((parent/'coverage.json').read_text());layers=json.loads((parent/'layer_coverage.json').read_text());sc=json.loads((parent/'support_self.json').read_text())
    stage=min(iteration,6);width_lines=3 if stage>=4 else 2;extend_lines=1 if stage==1 else 2
    pads=[];cuts=[];operations=[];blocked=[];seen=set();byid={b['branch_id']:b for b in g['branches']}
    all_lower=np.vstack([b['lower'] for b in g['branches']]);all_lower=all_lower[np.isfinite(all_lower).all(1)]
    # Mesh top just below the NEXT slicing plane. This is not a relaxed gap:
    # the intended extrusion top is the exact profile layer top recorded below.
    eps=4*float(np.spacing(np.float32(max(1.,abs(v.bounds).max()))))
    def add(xy,top,kind,identity,tangent=(1.,0.),too_close=False):
        key=(*np.round(xy,5),round(top,5),kind)
        if key in seen or top<=0:return
        seen.add(key);mesh_top=top+h/2-eps
        length=(2+extend_lines)*w;wid=width_lines*w
        offsets=[0.] if stage<3 else [0.,w,-w]
        t=np.asarray(tangent);t=t/max(np.linalg.norm(t),1e-12);side=np.array([-t[1],t[0]])
        row=dict(target_kind=kind,target_id=identity,point=[*np.asarray(xy).tolist(),top],length_line_multiples=2+extend_lines,width_line_multiples=width_lines,layer_top_mm=top,slicer_half_layer_mm=h/2,float_epsilon_mm=eps,strategy=['ROOF_EXTEND_1_LINE','ROOF_EXTEND_2_LINES','ROOF_LATERAL_SHIFT','ROOF_THICKEN_3_LINES','CURTAIN_THICKEN_3_LINES','SHORT_VERTICAL_CURTAIN'][stage-1])
        for shift in offsets:
            center=np.asarray(xy)+side*shift
            delta=all_lower[:,:2]-center
            near=(np.abs(delta@t)<=length/2+w/2)&(np.abs(delta@side)<=wid/2+w/2)
            if np.any(all_lower[near,2]-top<gap-eps):
                row['preflight_rejection']='PATCH_WOULD_INTRODUCE_TOO_CLOSE_ROOF_AT_OTHER_SAMPLE'
                continue
            if too_close:
                clipped=(np.abs(delta@t)<=w/2+w/2)&(np.abs(delta@side)<=w+w/2)
                if np.any(all_lower[clipped,2]-top>gap+h+eps):
                    row['preflight_rejection']='LOCAL_CAP_WOULD_REMOVE_OTHER_REQUIRED_ROOF'
                    continue
            # Coverage repair starts near the roof, joining existing curtain.
            # Self-support repair can extend a short vertical web to the bed.
            bottom=0. if kind=='SUPPORT_SELF' or stage>=5 else max(0.,mesh_top-4*h)
            patch=local_box(center,mesh_top,wid,length,t,bottom)
            iv=overlap(patch,v)
            if iv>=1e-7:continue
            pads.append(patch);row.update(lateral_shift_line_multiples=shift/w,intersection_mm3=iv,bottom_mm=bottom)
            if too_close:
                # Actual upper extrusion invalidates the nominal gap. Lower
                # only this local roof, by complete actual layer increments.
                high=min(float(s.bounds[1,2]+h),mesh_top+3*h)
                # Do not flatten an entire sloping neighbourhood: cap just
                # one line-width along the arc, with two lines across it.
                cuts.append(local_box(center,high,2*w,w,t,mesh_top))
                row['remove_too_close_upper_layers']=True
            operations.append(row);return
        blocked.append({**row,'reason':row.get('preflight_rejection','ALL_LOCAL_PATCH_POSITIONS_INTERSECT_FROZEN_VASCULAR')})
    for r in cov['rows']:
        if r['supported'] or r['lower_surface'] is None:continue
        q=np.asarray(r['lower_surface']);b=byid[r['branch_id']];i=r['sample_index']
        add(q[:2],quantized_top(q[2],gap,p),'UNDERBODY',f"{r['branch_id']}:{i}",b['tangent'][i,:2],bool(r.get('too_close_support_layer_tops')))
    for ri,r in enumerate(layers['rows']):
        if r['passed']:continue
        poly=shapely.from_wkt(r['uncovered_wkt'])
        for pi,part in enumerate(shapely.get_parts(poly)):
            if part.is_empty or part.area<=0:continue
            lo=np.asarray(part.bounds[:2]);hi=np.asarray(part.bounds[2:]);xy=[np.asarray(part.representative_point().coords[0])]
            for x in np.arange(lo[0],hi[0]+w/2,w/2):
                for y in np.arange(lo[1],hi[1]+w/2,w/2):
                    if part.covers(Point(x,y)):xy.append([x,y])
            xy=np.asarray(xy);origin=np.c_[xy,np.full(len(xy),r['z_mm']-h/2+eps)]
            lows=lower_points(v,origin)
            for j,q in enumerate(lows):
                if np.isfinite(q).all():add(q[:2],quantized_top(q[2],gap,p),'LAYER_DOMAIN',f'{ri}:{pi}:{j}')
    for ri,r in enumerate(sc['rows']):
        part=shapely.from_wkt(r['wkt']);lo=np.asarray(part.bounds[:2]);hi=np.asarray(part.bounds[2:]);xy=[np.asarray(part.representative_point().coords[0])]
        for x in np.arange(lo[0],hi[0]+w,w):
            for y in np.arange(lo[1],hi[1]+w,w):
                if part.covers(Point(x,y)):xy.append([x,y])
        for j,q in enumerate(xy):add(q,float(r['z_mm']-h),'SUPPORT_SELF',f'{ri}:{j}')
    changed=union([s,*pads]) if pads else s.copy()
    if cuts:changed=boolean([changed,union(cuts)],'difference')
    result,encoding=canonical_support(changed)
    # Preserve source triangle array and transformed coordinates byte-for-byte.
    for name in ['vascular_core.stl','vascular_print_transform.json','inverse_transform.json','notches.json']:
        if (parent/name).exists():shutil.copy2(parent/name,out/name)
    result.export(out/'custom_support.stl',file_type=encoding);review=trimesh.util.concatenate([v,result]);review.export(out/'vascular_core_plus_support_review.stl',file_type=encoding);review.export(out/'vascular_plus_support_review.stl',file_type=encoding)
    shifted=result.copy();shifted.apply_translation([0,0,-3*h]);g['roofs']=boolean([result,shifted],'difference');g['curtains']=boolean([result,g['roofs']],'difference');g['support']=result
    g['roofs'].export(out/'support_roofs.stl');g['curtains'].export(out/'support_curtains.stl')
    with (out/'geometry.pkl').open('wb') as f:pickle.dump(g,f)
    record=json.loads((parent/'geometry.json').read_text());record.update(candidate_id=out.name,folder=str(out),closure_parent=str(parent),closure_iteration=iteration,support_volume_mm3=float(result.volume),intersection_mm3=overlap(result,v),bed_roots=bed_roots(result),support_stl_encoding=encoding)
    save(out/'geometry.json',record);save(out/'local_repairs.json',dict(operations=operations,blocked=blocked,stage=stage,source_sha256=digest(out/'vascular_core.stl'),parent_source_sha256=digest(parent/'vascular_core.stl'),quantization=dict(line_width=w,layer_height=h),removed_target_count=0,architecture_regenerated=False))
    print('LOCAL_REPAIR',out.name,'operations',len(operations),'blocked',len(blocked),'intersection',record['intersection_mm3'],flush=True)
    return record


def baseline_jobs(base):
    summary=json.loads((base/'vascular_core_support_summary.json').read_text())
    comparisons=json.loads((base/'tables/candidate_comparison.json').read_text())
    jobs=[base/'candidates'/summary['selected_candidate']];orientations=['O002']
    for row in comparisons:
        if row['orientation'] not in orientations:orientations.append(row['orientation'])
        if len(orientations)==3:break
    # Use the already generated candidates. No architecture or 48-way search.
    for ori in orientations:
        for gap in [.2,.25]:
            choices=[r for r in comparisons if r['orientation']==ori and r['gap_mm']==gap]
            if not choices:continue
            best=min(choices,key=lambda r:r['layer_coverage']['uncovered_area_sum_mm2'])
            path=base/'candidates'/best['candidate']
            if path not in jobs:jobs.append(path)
    return jobs


def run(base=None,max_iterations=8,review_only=False):
    base=Path(base or default_output()).resolve();out=base/'closure_v2'
    for n in ['reports','validation','candidates','selected','QC_SAFE','tables']:(out/n).mkdir(parents=True,exist_ok=True)
    inputs=discover_source();hashes=freeze(inputs,base);save(out/'validation/protected_before.json',hashes)
    if not (out/'reports/baseline').exists():
        snap=out/'reports/baseline';snap.mkdir()
        summary=json.loads((base/'vascular_core_support_summary.json').read_text());folder=base/'candidates'/summary['selected_candidate']
        for p in [base/'vascular_core_support_report.md',base/'vascular_core_support_summary.json',base/'tables/candidate_comparison.json',*[folder/(n+'.json') for n in ['coverage','layer_coverage','removal','wrap_audit','travel_near_miss']]]:shutil.copy2(p,snap/p.name)
    profile,actual,engine=current_profiles(inputs,out/'bambu/profiles');rows=[];results=[];stop=[]
    for job in baseline_jobs(base):
        parent=job;e=json.loads((parent/'evaluation.json').read_text());r=json.loads((parent/'geometry.json').read_text());hist=[loss(e)];results.append((r,e));print('BASELINE',job.name,loss(e),flush=True)
        for iteration in range(1,max_iterations+1):
            if closure_pass(e):break
            folder=out/'candidates'/f'{job.name}_I{iteration:02}'
            if (folder/'geometry.json').exists():
                prior_record=json.loads((folder/'geometry.json').read_text())
                if prior_record.get('closure_parent')!=str(parent):
                    backup=folder.with_name(folder.name+'_rejected_parent')
                    if not backup.exists():folder.rename(backup)
                    else:raise RuntimeError('AMBIGUOUS_CLOSURE_RESUME_PARENT')
            if (folder/'evaluation.json').exists():
                rr=json.loads((folder/'geometry.json').read_text());ee=json.loads((folder/'evaluation.json').read_text())
            elif review_only:
                error_file=folder/'runtime_error.json'
                if error_file.exists():
                    error=json.loads(error_file.read_text())
                    stop.append(dict(baseline=job.name,reason='LOCAL_REPAIR_RUNTIME_FAILURE_NEXT_LOCAL_VARIANT',iteration=iteration,error=error['error']))
                    continue
                break
            else:
                try:
                    rr=json.loads((folder/'geometry.json').read_text()) if (folder/'geometry.json').exists() else repair_local(parent,folder,iteration)
                    ee=evaluate(rr,actual,engine)
                except Exception as error:
                    import traceback
                    save(folder/'runtime_error.json',dict(error=str(error),traceback=traceback.format_exc()))
                    print('REPAIR_ERROR',folder.name,str(error),flush=True)
                    # A geometry/runtime failure is not an actual coverage
                    # measurement; never count it as two-round stagnation.
                    stop.append(dict(baseline=job.name,reason='LOCAL_REPAIR_RUNTIME_FAILURE_NEXT_LOCAL_VARIANT',iteration=iteration,error=str(error)));continue
            if not ee.get('coverage'):
                stop.append(dict(baseline=job.name,reason='NATIVE_SLICE_FAILED_NO_ACTUAL_MEASUREMENT',iteration=iteration,candidate=folder.name))
                continue
            current=loss(ee);prior=loss(e);previous_trial=hist[-1];results.append((rr,ee))
            accepted=current[0]<=prior[0] and current[1]<=prior[1] and (current[0]<prior[0] or current[1]<prior[1] or current[2]<prior[2])
            hist.append(current)
            rows.append(dict(baseline=job.name,iteration=iteration,candidate=folder.name,uncovered_samples=current[0],uncovered_area_mm2=current[1],self_unsupported_regions=current[2],sample_improved=current[0]<prior[0],area_improved=current[1]<prior[1],accepted=accepted,parent=str(parent),centerline_percent=ee.get('coverage',{}).get('centerline_underbody_percent'),port_percent=ee.get('coverage',{}).get('port_underbody_percent'),layer_percent=ee.get('layer_coverage',{}).get('percent'),coverage_pass=coverage_pass(ee),closure_pass=closure_pass(ee),native_slice_success=ee.get('slice',{}).get('slicer_success',False)))
            rows[-1].update(previous_actual_trial_uncovered_samples=previous_trial[0],previous_actual_trial_uncovered_area_mm2=previous_trial[1],sample_improved_over_previous_trial=current[0]<previous_trial[0],area_improved_over_previous_trial=current[1]<previous_trial[1])
            table(out/'support_closure_iterations.csv',rows);save(out/'reports/closure_progress.json',dict(rows=rows,stop=stop,production_selected=None))
            if accepted:parent=folder;r,e=rr,ee
            if stagnated(hist):
                stop.append(dict(baseline=job.name,reason='TWO_CONSECUTIVE_REAL_ROUNDS_WITHOUT_SAMPLE_OR_AREA_IMPROVEMENT',iteration=iteration,history=hist));break
        else:stop.append(dict(baseline=job.name,reason='ITERATION_LIMIT',iteration=max_iterations,history=hist))
    valid=[x for x in results if 'coverage' in x[1]]
    # Diagnostic priority is coverage first. Never rank it as production.
    final_r,final_e=min(valid,key=lambda x:diagnostic_rank(x[1]))
    save(out/'reports/closure_progress.json',dict(rows=rows,stop=stop,production_selected=None,diagnostic_selected=final_r['candidate_id']))
    table(out/'support_closure_iterations.csv',rows)
    from .core_support_safe_removal import finish
    return finish(base,out,final_r,final_e,rows,stop,hashes,actual,engine,alternatives=[x for x in valid if closure_pass(x[1])])
