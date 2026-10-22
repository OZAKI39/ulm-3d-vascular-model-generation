"""RFUS production generator, O000 only; all historical experiments read-only."""
from __future__ import annotations
import json,pickle,shutil,time,copy
from pathlib import Path
import numpy as np
import trimesh
import yaml
from .support_inputs import ROOT,save,table,digest,verify,discover,transform
from .support_v2_inputs import load,oriented
from .support_rfus_geometry import union,boolean
from .support_v2_mesh import support_stl,exact_stl
from .support_rfus_geometry import corridor,segment_ranges,route_segment,add_handle,neck_audit


def preparation(source=None,output=None,config=None):
    inputs=discover(source=source);base=inputs['root'].parent/'vascular_support_planner';out=base/'robust_full_underbody'
    if output and Path(output).resolve()!=out.resolve():raise ValueError('RFUS_WRITES_ONLY_ROBUST_FULL_UNDERBODY')
    out.mkdir(parents=True,exist_ok=True);validation=out/'validation';validation.mkdir(exist_ok=True)
    hp=validation/'protected_before.json'
    if not hp.exists():
        hashes=json.loads((base/'inverted_study/validation/protected_before.json').read_text())
        for p in base.rglob('*'):
            if p.is_file() and 'robust_full_underbody' not in p.relative_to(base).parts:hashes[str(p)]=digest(p)
        # Freeze historical implementation as well as all historical data.
        for p in (ROOT/'vascular_processing').glob('support_v*.py'):hashes[str(p)]=digest(p)
        save(hp,hashes)
        shutil.copy2(ROOT/'tests/test_s1_8_vascular_support_planner.py',validation/'historical_tests_before.py.txt')
    protection=verify(json.loads(hp.read_text()))
    reads={}
    for name,path in [('v3_report',base/'v3/vascular_support_v3_report.md'),('v3_summary',base/'v3/vascular_support_v3_summary.json'),('v31_report',base/'v3_1/orientation_then_v3_1_report.md'),('v31_summary',base/'v3_1/support_orientation_and_closure_summary.json'),('v31_coverage',base/'v3_1/iteration_02/actual_coverage.csv')]:
        text=path.read_text();reads[name]=dict(path=str(path),sha256=digest(path),content=json.loads(text) if path.suffix=='.json' else text)
    save(validation/'baseline_read.json',reads)
    inputs,profile,discovery,mold,v,box,branches,top=load(validation/'inputs',source)
    cfg=yaml.safe_load(Path(config or ROOT/'config/robust_underbody.yaml').read_text())
    if cfg['orientation']!='O000' or cfg['scale_factor']!=1.:raise ValueError('RFUS_O000_SCALE_ONE_ONLY')
    if not 0<cfg['sample_spacing_lines']<=1:raise ValueError('RFUS_SAMPLING_MUST_NOT_EXCEED_LINE_WIDTH')
    old=reads['v3_summary']['content']['diagnostic_selected']['orientation'];matrix=np.asarray(old['print_transform'])
    if not np.allclose(matrix[:3,:3],np.eye(3),atol=1e-10):raise ValueError('BASELINE_IS_NOT_O000')
    profile['roof_gap_mm']=cfg['roof_gap_mm'];profile['gap_status']='PHYSICAL_GAP_UNCALIBRATED'
    for mesh in (mold,v,box):mesh.apply_transform(matrix)
    branches,top=oriented(branches,top,matrix);inner=transform(np.asarray(inputs['inner_bounds']),matrix)
    save(out/'robust_underbody_config.json',{**cfg,'profile':profile,'orientation_record':old,'source_unchanged':protection['all_unchanged']})
    save(out/'print_transform.json',dict(transform_4x4=matrix,scale_factor=1.))
    save(out/'casting_restore_transform.json',dict(transform_4x4=np.linalg.inv(matrix),scale_factor=1.))
    cache=out/'corridors.pkl'
    if cache.exists():
        with cache.open('rb') as f:corridors=pickle.load(f)
        from .support_rfus_geometry import lower_points
        changed=False
        for b in corridors:
            missing=np.flatnonzero(~b['surface_valid'])
            if not len(missing):continue
            hits=lower_points(v,b['xyz'][missing]);ok=np.isfinite(hits).all(1);ids=missing[ok]
            if len(ids):
                b['lower'][ids]=hits[ok];b['surface_valid'][ids]=True
                first=profile['first_layer_height'];h=profile['layer_height']
                b['roof_top'][ids]=first+np.floor((hits[ok,2]-profile['roof_gap_mm']-first)/h)*h
                changed=True
        if changed:
            (out/'history').mkdir(exist_ok=True)
            if not (out/'history/corridors_before_cap_ray_fix.pkl').exists():shutil.copy2(cache,out/'history/corridors_before_cap_ray_fix.pkl')
            with cache.open('wb') as f:pickle.dump(corridors,f)
    else:
        corridors=[corridor(b,v,inner,profile,cfg) for b in branches]
        with cache.open('wb') as f:pickle.dump(corridors,f)
    rows=[]
    for b in corridors:
        for i,s in enumerate(b['arc']):rows.append(dict(branch_id=b['branch_id'],is_port=bool(b['is_port']),sample_index=i,arc_mm=float(s),xyz=b['xyz'][i].tolist(),tangent=b['tangent'][i].tolist(),gravity=[0,0,-1],lower_surface=b['lower'][i].tolist() if b['surface_valid'][i] else None,lower_surface_valid=bool(b['surface_valid'][i]),radius_mm=float(b['radius'][i]),curvature=float(b['curvature'][i]),root_domain=b['root_domain'][i],roof_top_mm=float(b['roof_top'][i]) if b['surface_valid'][i] else None))
    table(out/'underbody_corridors.csv',rows)
    print('RFUS INPUTS',len(corridors),'branches/ports',len(rows),'full-arc samples','profile',profile['line_width'],flush=True)
    return dict(base=base,out=out,inputs=inputs,profile=profile,discovery=discovery,cfg=cfg,mold=mold,vascular=v,box=box,corridors=corridors,topology=top,inner=inner,orientation=old,corridor_digest=digest(cache))


def trim_overlap(segment,previous,profile):
    """Partition overlapping support volumes; previously assigned support stays.

    A narrow explicit seam belongs to the support, never to the frozen model.
    All original underbody samples remain in the final actual coverage audit.
    """
    if not previous:return segment
    from .support_custom_qc import intersection_volume
    margin=profile['line_width']/2
    neighbours=[m for m in previous if np.all(m.bounds[1]+margin>=segment['mesh'].bounds[0]) and np.all(m.bounds[0]-margin<=segment['mesh'].bounds[1])]
    if not neighbours:return segment
    blocker=union(neighbours)
    manager=trimesh.collision.CollisionManager();manager.add_object('previous',blocker)
    if manager.min_distance_single(segment['mesh'])>margin:return segment
    # An exact coincident interface becomes a four-face nonmanifold edge when
    # STL welds it. Open an actual narrow XY seam in the NEW support only.
    copies=[blocker]
    dz=profile['layer_height']/10
    for delta in [(margin,0,0),(-margin,0,0),(0,margin,0),(0,-margin,0),(0,0,dz),(0,0,-dz)]:
        shifted=blocker.copy();shifted.apply_translation(delta);copies.append(shifted)
    blocker=union(copies)
    for key in ['mesh','curtains','roofs','necks','feet','handles']:
        if not len(segment[key].faces):continue
        segment[key]=boolean([segment[key],blocker],'difference')
    segment['shared_projection_ownership']=True;segment['ownership_seam_margin_mm']=margin
    return segment


def plan(context,lines,neck_lines,repair_of=None,feature_upgrades=None):
    c=context;out=c['out'];profile=c['profile'];cfg=c['cfg'];w=profile['line_width'];h=profile['layer_height']
    identity=f'RFUS_C{lines}_R{lines}_N{neck_lines}'+('_RMS' if repair_of and feature_upgrades and feature_upgrades.get('removal_segments') else '_P1' if repair_of else '')
    folder=out/'candidates'/identity;folder.mkdir(parents=True,exist_ok=True);geometry=folder/'geometry';geometry.mkdir(exist_ok=True)
    cached=folder/'candidate.json'
    if cached.exists() and json.loads(cached.read_text()).get('planning_complete'):return json.loads(cached.read_text())
    segments=[];unresolved=[];exceptions=[];seams=[];attempts=[];sample_rows=[]
    for bi,b in enumerate(c['corridors']):
        ranges=segment_ranges(b,c['topology'],cfg);jobs=[]
        requests=(feature_upgrades or {}).get('removal_segments',[])
        for ri,(a,z) in enumerate(ranges):
            sid=f'B{bi:02}_S{ri:03}'
            if any(r['branch_id']==b['branch_id'] and r['arc_end_mm']>b['arc'][a]+.01 and r['arc_start_mm']<b['arc'][z]-.01 for r in requests) and z-a>=4:
                mid=(a+z)//2;jobs.extend([(sid+'_A',a,mid,sid),(sid+'_B',mid,z,sid)])
            else:jobs.append((sid,a,z,sid))
        for sid,a,z,parent_sid in jobs:
            s=b['arc'][a:z+1].copy();xyz=b['xyz'][a:z+1].copy();tops=b['roof_top'][a:z+1].copy();norm=b['lateral'][a:z+1].copy();base=b['base'][a:z+1].copy()
            # Keep a complete target set. Samples embedded in box material are
            # omitted geometrically only provisionally, pending actual paths.
            valid=b['surface_valid'][a:z+1].copy();low=b['lower'][a:z+1]
            contact=np.zeros(len(low),bool)
            if valid.any():contact[valid]=c['box'].contains(low[valid]-[0,0,.02])
            build=valid&~contact&(tops-base>1.5)
            from .support_detection import contiguous_runs
            for i in np.flatnonzero(~build):exceptions.append(dict(segment_id=sid,branch_id=b['branch_id'],arc_mm=float(s[i]),reason='MODEL_CONTACT_REQUIRES_ACTUAL_PATH_PROOF' if contact[i] else 'NO_FROZEN_LOWER_RAY' if not valid[i] else 'INSUFFICIENT_HEIGHT_REQUIRES_ACTUAL_MODEL_PROOF'))
            for run,(u,t) in enumerate(contiguous_runs(build)):
                if t<=u:
                    unresolved.append(dict(segment_id=sid,branch_id=b['branch_id'],arc_start=float(s[u]),arc_end=float(s[t]),reason='ISOLATED_TARGET_SAMPLE'));continue
                coords=xyz[u:t+1,:2].copy();tt=tops[u:t+1].copy();bb=base[u:t+1].copy();nn=norm[u:t+1].copy();ss=s[u:t+1].copy()
                # Full-height slit, half an actual line wide. Its sub-line
                # length is reported, never silently credited as supported.
                slit=cfg['seam_width_lines']*w
                if a>0 and u==0:
                    d=coords[1]-coords[0];length=np.linalg.norm(d)
                    if length>slit/2:coords[0]+=d/length*slit/2
                if z<len(b['arc'])-1 and t==len(s)-1:
                    d=coords[-1]-coords[-2];length=np.linalg.norm(d)
                    if length>slit/2:coords[-1]-=d/length*slit/2
                if a==0 and u==0:
                    d=coords[1]-coords[0];coords[0]-=d/max(np.linalg.norm(d),1e-9)*w
                if z==len(b['arc'])-1 and t==len(s)-1:
                    d=coords[-1]-coords[-2];coords[-1]+=d/max(np.linalg.norm(d),1e-9)*w
                key=f'{sid}_{run}';parent_key=f'{parent_sid}_{run}'
                def upgrade(kind):return feature_upgrades and any(k in feature_upgrades.get(kind,[]) for k in (key,parent_key))
                nl=2 if upgrade('necks') else neck_lines
                cl=3 if upgrade('curtains') else lines
                rl=3 if upgrade('roofs') else lines
                r,tries=route_segment(coords,tt,bb,nn,profile,cfg,c['vascular'],c['box'],cl,rl,nl)
                attempts.extend([dict(segment_id=key,**x) for x in tries])
                if r is None:
                    # Adaptive segmentation at a blocked projection preserves
                    # the targets; each shorter piece retries the same routing.
                    pieces=[]
                    for k in range(0,len(coords)-1,max(2,int(round(1.5/(w*cfg['sample_spacing_lines']))))):
                        end=min(len(coords),k+max(2,int(round(1.5/(w*cfg['sample_spacing_lines']))))+1)
                        sub,subtries=route_segment(coords[k:end],tt[k:end],bb[k:end],nn[k:end],profile,cfg,c['vascular'],c['box'],cl,rl,nl)
                        subid=key+f'_L{k:03}';attempts.extend([dict(segment_id=subid,**x) for x in subtries])
                        if sub is None:unresolved.append(dict(segment_id=subid,branch_id=b['branch_id'],arc_start=float(ss[k]),arc_end=float(ss[end-1]),reason='ALL_OFFSET_TWIN_PILLAR_ROUTES_FAILED'));continue
                        pieces.append((subid,sub,ss[k],ss[end-1]))
                else:pieces=[(key,r,ss[0],ss[-1])]
                for subid,r,start,end in pieces:
                    r=trim_overlap(r,[p['mesh'] for p in segments],profile)
                    if not len(r['mesh'].faces):
                        exceptions.append(dict(segment_id=subid,branch_id=b['branch_id'],arc_mm=float(start),reason='SHARED_PROJECTION_ASSIGNED_TO_PREVIOUS_CURTAIN; ACTUAL_COVERAGE_STILL_REQUIRED'));continue
                    r.update(segment_id=subid,branch_id=b['branch_id'],branch_index=bi,is_port=bool(b['is_port']),arc_start_mm=float(start),arc_end_mm=float(end),arc_length_mm=float(end-start),internal=bool(np.all(bb>0)),curtain_lines=cl,roof_lines=rl,neck_lines=nl,root_domain='BOX_INNER_BOTTOM' if np.all(bb>0) else 'BUILD_PLATE')
                    r['neck_audit']=neck_audit(r);segments.append(r)
                seams.append(dict(branch_id=b['branch_id'],arc_mm=float(s[-1]),position=xyz[-1].tolist(),slit_width_mm=slit,style='FULL_HEIGHT_VERTICAL_SLIT_WITH_SUB_LINE_ROOF_BREAK',actual_coverage_required=True))
    print('RFUS CURTAINS',identity,len(segments),'unresolved',len(unresolved),flush=True)
    for i,r in enumerate(segments):
        if r['internal']:add_handle(r,c['vascular'],c['box'],[p['mesh'] for j,p in enumerate(segments) if i!=j],c['inner'],profile,cfg)
        else:r['handle_target']=None
        r['neck_audit']=neck_audit(r)
        if i%10==0:print('RFUS HANDLES',identity,i,'/',len(segments),flush=True)
    return write_plan(c,folder,segments,identity,lines,neck_lines,unresolved,exceptions,seams,attempts,repair_of)


def write_plan(c,folder,segments,identity,lines,neck_lines,unresolved,exceptions,seams,attempts,repair_of=None):
    with (folder/'pre_export_segments.pkl').open('wb') as f:pickle.dump(segments,f)
    g=folder/'geometry';parts={key:union([s[key] for s in segments if len(s[key].faces)]) for key in ['mesh','roofs','curtains','handles','necks']}
    parts['mesh'].export(g/'support_indexed_before_export.ply')
    from .support_rfus_geometry import exported_support
    support_export_error=None
    try:support=exported_support(parts['mesh'])
    except ValueError as exc:
        # Keep the exact diagnostic shape and continue real slicing. This is
        # an explicit hard failure, never a silently repaired production mesh.
        support=parts['mesh'];support_export_error=str(exc)
    support.export(g/'custom_full_underbody_support.stl')
    for key in ['roofs','curtains','handles','necks']:parts[key].export(g/f'support_{key}.stl')
    for key,name in [('mold','mold_without_support'),('vascular','vascular_with_ports'),('box','box')]:c[key].export(g/f'{name}.stl')
    combined=union([c['mold'],support]);assembly_error=None
    try:combined=exact_stl(combined)
    except ValueError as exc:assembly_error=str(exc);combined.export(g/'assembly_indexed.ply')
    combined.export(g/'mold_with_robust_support.stl')
    sd=g/'segments';sd.mkdir(exist_ok=True)
    part_export_notes=[]
    for s in segments:
        try:support_stl(s['mesh']).export(sd/(s['segment_id']+'.stl'))
        except ValueError as exc:
            # The unioned production support has already passed binary STL
            # verification. Preserve exact indexed individual topology for
            # decomposition/QC if a zero-width ownership edge cannot round-trip.
            s['mesh'].export(sd/(s['segment_id']+'.ply'))
            part_export_notes.append(dict(segment_id=s['segment_id'],format='INDEXED_PLY',reason=str(exc),source_mesh_closed=bool(s['mesh'].is_volume)))
    with (folder/'segment_geometry.pkl').open('wb') as f:pickle.dump(segments,f)
    def row(s):return {k:v for k,v in s.items() if k not in ['mesh','roofs','curtains','handles','necks','feet','path_xy','top','base','lateral','route_attempts','handle_attempts']}
    table(folder/'curtain_segments.csv',[row(s) for s in segments]);table(folder/'roof_ridges.csv',[dict(segment_id=s['segment_id'],width_mm=s['roof_width_mm'],lines=s['roof_lines'],gap_mm=c['profile']['roof_gap_mm'],arc_start_mm=s['arc_start_mm'],arc_end_mm=s['arc_end_mm']) for s in segments])
    table(folder/'removal_handles.csv',[dict(segment_id=s['segment_id'],target=s.get('handle_target'),width_mm=s.get('handle_width_mm'),access_zone_min_z_mm=s.get('handle_access_zone_min_z_mm'),planned=bool(s.get('handle_target'))) for s in segments if s['internal']])
    table(folder/'base_necks.csv',[dict(segment_id=s['segment_id'],lines=s['neck_lines'],width_mm=s['neck_width_mm'],cut_z_mm=s['neck_cut_z'],**s['neck_audit']) for s in segments])
    table(folder/'segment_seams.csv',seams);table(folder/'route_attempts.csv',attempts);table(folder/'unresolved_corridors.csv',unresolved);table(folder/'provisional_model_contacts.csv',exceptions)
    from .support_custom_qc import intersection_volume
    hits=intersection_volume(support,c['vascular']);boxhit=trimesh.boolean.intersection([support,c['box']],engine='manifold');feet=union([s['feet'] for s in segments]);illegal=trimesh.boolean.difference([boxhit,feet],engine='manifold') if len(boxhit.faces) else boxhit
    volume=abs(float(illegal.volume)) if len(illegal.faces) else 0.
    record=dict(mode='ROBUST_FULL_UNDERBODY_SUPPORT',candidate_id=identity,folder=str(folder),planning_complete=True,source_unchanged=True,orientation=c['orientation'],scale_factor=1.,curtain_lines=lines,neck_lines=neck_lines,segment_count=len(segments),internal_segment_count=sum(s['internal'] for s in segments),handle_count=sum(bool(s.get('handle_target')) for s in segments),mean_segment_length_mm=float(np.mean([s['arc_length_mm'] for s in segments])),support_volume_mm3=float(support.volume),vascular_port_intersection_mm3=hits,illegal_box_intersection_mm3=volume,all_neck_planes_valid=all(s['neck_audit']['passed'] for s in segments),unresolved_corridors=unresolved,provisional_model_contacts=len(exceptions),seams=seams,assembly_export_error=assembly_error,repair_of=repair_of,roof_wraps_vessel=False,wrap_rule='Only narrow lower ridges; no saddle/circumferential primitives',slice={},removal={},physical_status=['PHYSICAL_GAP_UNCALIBRATED','CUTTABLE_BASE_NECK_NOT_PHYSICALLY_CALIBRATED','CURTAIN_STRENGTH_UNCALIBRATED','HANDLE_STRENGTH_UNCALIBRATED','MANUAL_REMOVAL_FORCE_UNCALIBRATED'])
    record['support_export_error']=support_export_error
    record['individual_segment_export_notes']=part_export_notes
    save(folder/'candidate.json',record);return record


def run(source=None,output=None,stage='all',config=None):
    began=time.monotonic();c=preparation(source,output,config);records=[]
    if stage=='report':
        records=[json.loads(p.read_text()) for p in sorted((c['out']/'candidates').glob('*/candidate.json')) if json.loads(p.read_text()).get('evaluation_complete')]
        if not records:raise ValueError('NO_COMPLETED_RFUS_RUN')
    else:
        for lines in c['cfg']['curtain_lines']:
            for neck in c['cfg']['neck_lines']:
                r=plan(c,lines,neck);records.append(r)
                if stage in ('plan','geometry'):continue
                from .support_rfus_validation import evaluate
                r=evaluate(c,r);records[-1]=r
                adjustments=copy.deepcopy(r.get('feature_upgrades',{}))
                removal_segments=[{k:s[k] for k in ('segment_id','branch_id','arc_start_mm','arc_end_mm')} for s in r['removal']['segments'] if s['status']!='REMOVAL_PATH_FOUND']
                if removal_segments:adjustments['removal_segments']=removal_segments
                if adjustments and c['cfg']['maximum_preservation_repairs']:
                    repaired=plan(c,lines,neck,repair_of=r['candidate_id'],feature_upgrades=adjustments);repaired=evaluate(c,repaired);records.append(repaired)
    verify(json.loads((c['out']/'validation/protected_before.json').read_text()))
    if stage in ('plan','geometry'):return dict(status='RFUS_GEOMETRY_COMPLETE',output=str(c['out']))
    from .support_rfus_review import generate
    return generate(c,records,time.monotonic()-began)
