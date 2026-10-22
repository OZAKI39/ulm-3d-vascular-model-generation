"""RFUS evidence from actual saved X1C toolpaths, not support STL coverage."""
from __future__ import annotations
import json,pickle,shutil,zipfile
from pathlib import Path
import numpy as np
import trimesh
import shapely
from shapely.geometry import Point,GeometryCollection
from shapely.ops import unary_union
from scipy.spatial import cKDTree
from .support_inputs import save,table,digest,numeric
from .support_v3_validation import extrusion_sections,actual_self_support,coverage
from .support_v2_layers import analyse
from .support_qc import parse_actual_paths
from .support_rfus_geometry import union


def underbody_coverage(corridors,footprints,model_footprints,profile):
    rows=[];w,h=profile['line_width'],profile['layer_height'];gap=profile['roof_gap_mm'];zs=np.array(sorted(footprints));mz=np.array(sorted(model_footprints));maxarc=0.;totalarc=missingarc=0.
    support_layers={};model_layers={}
    for b in corridors:
        flags=[]
        for i,(arc,low) in enumerate(zip(b['arc'],b['lower'])):
            valid=bool(np.isfinite(low).all());supported=False;direct=False;eligible=[];nearest=None
            if valid:
                ulp=float(np.spacing(np.float32(max(1.,np.max(np.abs(low))))));tol=4*ulp
                eligible=zs[(zs>=low[2]-gap-h-1e-6)&(zs<=low[2]-gap+tol)]
                p=Point(low[:2]);key=tuple(eligible)
                if key not in support_layers:support_layers[key]=unary_union([footprints[float(z)] for z in eligible])
                actual=support_layers[key]
                nearest=float(p.distance(actual)) if not actual.is_empty else None
                supported=bool(not actual.is_empty and nearest<=tol)
                # Already printed material must be truly below the first
                # target layer. The target's own first extrusion gets no credit.
                lower=mz[(mz>=low[2]-1.5*h-tol)&(mz<=low[2]-.5*h+tol)]
                key=tuple(lower)
                if key not in model_layers:model_layers[key]=unary_union([model_footprints[float(z)] for z in lower])
                existing=model_layers[key]
                direct=bool(not existing.is_empty and existing.distance(p)<=tol)
            else:tol=0.
            ok=supported or direct;flags.append(ok)
            rows.append(dict(branch_id=b['branch_id'],is_port=bool(b['is_port']),sample_index=i,arc_mm=float(arc),xyz=b['xyz'][i].tolist(),lower_surface=low.tolist() if valid else None,supported=ok,actual_custom_support=supported,actual_existing_model=direct,nearest_actual_xy_mm=nearest,actual_support_layer_top_mm=eligible.tolist() if valid else [],xy_numeric_tolerance_mm=tol,allowed_gap_mm=gap,layer_quantization_mm=h,surface_valid=valid))
        bad=~np.asarray(flags);arc=b['arc'];weights=np.r_[np.diff(arc),0.]/2+np.r_[0.,np.diff(arc)]/2
        totalarc+=arc[-1];missingarc+=weights[bad].sum()
        from .support_detection import contiguous_runs
        for a,z in contiguous_runs(bad):
            # Include half cells at run boundaries; a single missing sample
            # has nonzero extent rather than falsely reporting zero arc.
            lo=arc[a] if a==0 else (arc[a-1]+arc[a])/2;hi=arc[z] if z==len(arc)-1 else (arc[z]+arc[z+1])/2
            maxarc=max(maxarc,float(hi-lo))
    return dict(percent=100*sum(r['supported'] for r in rows)/len(rows) if rows else 0.,arc_weighted_percent=100*(1-missingarc/totalarc) if totalarc else 0.,sample_count=len(rows),unsupported_samples=sum(not r['supported'] for r in rows),all_samples_supported=bool(rows) and all(r['supported'] for r in rows),max_uncovered_underbody_arc_mm=maxarc,arc_total_mm=float(totalarc),uncovered_arc_mm=float(missingarc),sampling_max_spacing_mm=max(np.diff(b['arc']).max() for b in corridors),support_source='ACTUAL_GCODE_EXTRUSION_FOOTPRINT',rows=rows)


def features(segments,footprints,profile):
    from .support_detection import section_polygon
    zs=np.array(sorted(footprints));h=profile['layer_height'];rows=[];upgrades={}
    for s in segments:
        for key in ['curtains','roofs','necks','handles']:
            mesh=s[key]
            if not len(mesh.faces):continue
            lo,hi=mesh.bounds[:,2];levels=np.linspace(lo+(hi-lo)*.2,hi-(hi-lo)*.2,3)
            for i,z in enumerate(levels):
                near=zs[np.abs(zs-h/2-z)<=h/2+1e-5]
                cuts=mesh.section_multiplane([0,0,0],[0,0,1],near-h/2);nominal=unary_union([section_polygon(p) for p in cuts]);actual=unary_union([footprints[float(t)] for t in near]);ratio=nominal.intersection(actual).area/nominal.area if nominal.area else 0.
                ok=ratio>.5;rows.append(dict(segment_id=s['segment_id'],feature=key,station=i,z_mm=float(z),nominal_area_mm2=float(nominal.area),actual_area_fraction=float(ratio),preserved=ok,criterion='Feature existence >50% nominal cross-section; full coverage is a separate stricter gate'))
                current=s.get({'curtains':'curtain_lines','roofs':'roof_lines','necks':'neck_lines','handles':'handle_lines'}[key],3)
                if not ok and current<(2 if key=='necks' else 3):upgrades.setdefault(key,set()).add(s['segment_id'])
    return dict(all_preserved=bool(rows) and all(r['preserved'] for r in rows),rows=rows),{k:sorted(v) for k,v in upgrades.items()}


def vascular_features(corridors,footprints,profile):
    """Sample every branch and port against actual isolated model extrusion."""
    points=[]
    for z,p in footprints.items():
        for line in shapely.get_parts(shapely.segmentize(p.boundary,profile['line_width']/2)):
            if line.geom_type not in ('LineString','LinearRing'):continue
            xy=np.asarray(line.coords);points.append(np.c_[xy,np.full(len(xy),z)])
    tree=cKDTree(np.vstack(points)) if points else None;rows=[]
    for b in corridors:
        d=tree.query(b['xyz'])[0] if tree is not None else np.full(len(b['xyz']),1e9);limit=b['radius']+profile['line_width']+profile['layer_height']
        rows.append(dict(branch_id=b['branch_id'],is_port=bool(b['is_port']),samples=len(d),missing_samples=int((d>limit).sum()),maximum_nearest_actual_path_mm=float(d.max()),preserved=bool(np.all(d<=limit)),method='Every arc sample versus actual vascular-model extrusion; radius + one line + one layer search bound'))
    return dict(all_preserved=all(r['preserved'] for r in rows),rows=rows)


def gates(record):
    s=record.get('slice',{});u=s.get('underbody',{});r=record.get('removal',{})
    checks={
        'SOURCE_CHANGED':record.get('source_unchanged',False),
        'SCALE_CHANGED':record['scale_factor']==1.,
        'ORIENTATION_NOT_O000':record['orientation']['orientation_id']=='O000' and np.allclose(np.asarray(record['orientation']['print_transform'])[:3,:3],np.eye(3)),
        'AUTO_SUPPORT_NOT_OFF':s.get('auto_support_off',False) and s.get('native_support_toolpaths',-1)==0,
        'CENTERLINE_UNDERBODY_INCOMPLETE':u.get('all_samples_supported',False),
        'LONG_UNCOVERED_UNDERBODY_ARC':u.get('max_uncovered_underbody_arc_mm',1e9)<s.get('line_width_mm',0.),
        'ACTUAL_LAYER_COVERAGE_INCOMPLETE':s.get('required_layer_coverage',{}).get('all_required_covered',False),
        'ACTUAL_SUPPORT_SELF_UNSUPPORTED':s.get('self_printability',{}).get('unsupported_region_count',-1)==0,
        'VASCULAR_OR_PORT_INTERSECTION':record.get('vascular_port_intersection_mm3',1.)<1e-5,
        'ILLEGAL_BOX_INTERSECTION':record.get('illegal_box_intersection_mm3',1.)<1e-4,
        'ROOF_WRAPS_VESSEL':record.get('roof_wraps_vessel',True) is False,
        'BASE_NECK_BYPASSED':record.get('all_neck_planes_valid',False),
        'INACCESSIBLE_INTERNAL_SEGMENT':r.get('inaccessible_segments',-1)==0,
        'TRAPPED_SEGMENT':r.get('all_segments_removable',False),
        'SEGMENT_NOT_INDEPENDENT':record.get('segment_topology',{}).get('all_independent',False),
        'REAL_X1C_SLICE_FAILED':s.get('slicer_success',False) and s.get('support_only_slicer_success',False),
        'VASCULAR_TOOLPATH_DROPPED':s.get('vascular_features',{}).get('all_preserved',False),
        'CUSTOM_FEATURE_DROPPED':s.get('support_features',{}).get('all_preserved',False),
        'UNRESOLVED_CRITICAL_WARNINGS':s.get('unresolved_critical_warning_count',-1)==0,
        'SUPPORT_EXPORT_FAILED':record.get('support_export_error') is None,
        'ASSEMBLY_EXPORT_FAILED':record.get('assembly_export_error') is None}
    return [k for k,v in checks.items() if not v]


def refresh_underbody(c,record):
    """Re-audit saved toolpaths if numerical end-cap ray recovery changed targets."""
    fingerprint=c.get('corridor_digest')
    if not fingerprint or record['slice'].get('target_fingerprint')==fingerprint:return record
    folder=Path(record['folder'])
    with (folder/'actual_footprints.pkl').open('rb') as f:fp=pickle.load(f)
    model={z:fp['box'].get(z,GeometryCollection()).union(fp['vascular'].get(z,GeometryCollection())) for z in fp['box'].keys()|fp['vascular'].keys()}
    u=underbody_coverage(c['corridors'],fp['support'],model,c['profile'])
    prior=folder/'underbody_before_cap_ray_recovery.json'
    if not prior.exists():shutil.copy2(folder/'underbody_coverage.json',prior)
    table(folder/'underbody_coverage.csv',u['rows']);save(folder/'underbody_coverage.json',u)
    record['slice']['underbody']={k:v for k,v in u.items() if k!='rows'}
    record['slice']['target_fingerprint']=fingerprint
    save(folder/'bambu_validation.json',record['slice'])
    return record


def evaluate(c,record):
    folder=Path(record['folder']);g=folder/'geometry'
    if record.get('evaluation_complete'):return refresh_underbody(c,record)
    from .support_slicer_adapter import discover_engine,engine_profiles
    from .vascular_support_pipeline import sliced
    engine=discover_engine('Bambu',c['discovery'],folder);profiles=engine_profiles(c['discovery'],'Bambu',{'enable_support':'0'})
    source=g/'mold_with_robust_support.stl';result=sliced(engine,profiles,source,source,folder/'bambu',c['cfg']['slice_timeout_seconds'])
    save(folder/'bambu_raw_result.json',result)
    if not result.get('slicer_success'):raise RuntimeError('REAL_BAMBU_SLICE_FAILED: '+str(result))
    archive=folder/'X1C_diagnostic_full_support.3mf';shutil.copy2(result['archive'],archive)
    profile=c['profile'].copy()
    with zipfile.ZipFile(archive) as z:
        settings=json.loads(z.read('Metadata/project_settings.config'));crc=z.testzip()
    profile['initial_layer_line_width']=numeric(settings.get('initial_layer_line_width',profile['line_width']),profile['nozzle_diameter'])
    paths=parse_actual_paths(archive);extrusions=paths['segments'][paths['kind']==0]
    support=trimesh.load_mesh(g/'custom_full_underbody_support.stl')
    from .support_rfus_travel import events,actual_sections
    parsed=events(archive)
    fp,filament=actual_sections(parsed,support,profile,with_filament=True);boxfp=actual_sections(parsed,c['box'],profile);vfp=actual_sections(parsed,c['vascular'],profile)
    model_fp={z:boxfp.get(z,GeometryCollection()).union(vfp.get(z,GeometryCollection())) for z in boxfp.keys()|vfp.keys()}
    with (folder/'actual_footprints.pkl').open('wb') as f:pickle.dump(dict(support=fp,box=boxfp,vascular=vfp),f)
    u=underbody_coverage(c['corridors'],fp,model_fp,profile);table(folder/'underbody_coverage.csv',u['rows']);save(folder/'underbody_coverage.json',u)
    sf=actual_self_support(fp,boxfp,profile);table(folder/'support_self_printability.csv',sf['rows']);save(folder/'support_self_printability.json',sf)
    # Complete layer geometry is an independent validation target; it never
    # selects or removes production support corridors.
    field=analyse(c['vascular'],profile,'vascular',c['box'],bridge_limit=0.);lc,hit,miss=coverage(field,fp,profile)
    table(folder/'actual_required_layer_coverage.csv',lc['rows'])
    with (folder/'actual_layer_geometry.pkl').open('wb') as f:pickle.dump((hit,miss),f)
    with (folder/'segment_geometry.pkl').open('rb') as f:segments=pickle.load(f)
    features_result,upgrades=features(segments,fp,profile);vf=vascular_features(c['corridors'],vfp,profile)
    table(folder/'feature_preservation.csv',features_result['rows']);table(folder/'vascular_feature_preservation.csv',vf['rows'])
    foundation=folder/'support_with_box_foundation.stl';union([support,c['box']]).export(foundation)
    solo=sliced(engine,profiles,foundation,foundation,folder/'bambu_support_only',c['cfg']['slice_timeout_seconds'])
    native=result.get('warning_audit',{}).get('critical',[]);warnings=[]
    for row in lc['rows']:
        if not row['covered']:warnings.append(dict(source='ACTUAL_LAYER_RESIDUAL',z_mm=row['z_mm'],polygon_wkt=row['uncovered_wkt']))
    for row in sf['rows']:warnings.append(dict(source='ACTUAL_SUPPORT_SELF',**row))
    # Native floating warnings remain unresolved if any actual evidence fails.
    # Do not dismiss an unlocated slicer warning as an intentional gap.
    warnings.extend(dict(source='NATIVE_BAMBU_WARNING',detail=x) for x in native)
    table(folder/'warning_localization.csv',warnings)
    s=dict(slicer_success=True,support_only_slicer_success=bool(solo.get('slicer_success')),auto_support_off=str(settings.get('enable_support')) in ['0','false','False'],native_support_toolpaths=int((paths['kind']==1).sum()),archive_crc_passed=crc is None,line_width_mm=profile['line_width'],underbody={k:v for k,v in u.items() if k!='rows'},self_printability=sf,required_layer_coverage={k:v for k,v in lc.items() if k!='rows'},support_features=features_result,vascular_features=vf,native_warnings=native,unresolved_critical_warning_count=len(warnings),estimated_print_time_seconds=result.get('estimated_print_time_seconds'),actual_profile={k:settings.get(k) for k in ['layer_height','line_width','nozzle_diameter','enable_support','filament_type']})
    filament_diameter=numeric(settings['filament_diameter'])
    s['actual_support_extrusion_volume_mm3']=float(filament*np.pi*(filament_diameter/2)**2)
    s['extrusion_width_source']='PER_MOVE_ARCHIVED_LINE_WIDTH; no uniform-width approximation'
    record['slice']=s;record['feature_upgrades']=upgrades;save(folder/'bambu_validation.json',s);save(folder/'candidate.json',record)
    print('RFUS ACTUAL',record['candidate_id'],'underbody',u['percent'],'layer',lc['percent'],'self',sf['unsupported_region_count'],flush=True)
    from .support_rfus_removal import run as removal
    record['removal']=removal(c,record)
    from .support_rfus_removal import topology_audit
    record=topology_audit(record)
    from .support_rfus_travel import audit
    record['travel']=audit(archive,c,record)
    record['failed_gates']=gates(record);record['status']='NEEDS_ADJUSTMENT' if record['failed_gates'] else 'READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED';record['evaluation_complete']=True
    save(folder/'candidate.json',record);return record
