"""Real Bambu model extrusion, clipped by support sections, never nominal roofs."""
from __future__ import annotations
import json, pickle, shutil, zipfile,hashlib
from pathlib import Path
import numpy as np
import trimesh
import shapely
from shapely.geometry import GeometryCollection, Point, LineString
from shapely.ops import unary_union
from .support_inputs import save, table,digest
from .support_qc import parse_actual_paths
from .support_detection import section_polygon
from .support_slicer_adapter import discover_engine, engine_profiles
from .vascular_support_pipeline import sliced
from .support_v3_geometry import union, curtain
from .support_v3_planner import POLICY


def extrusion_sections(segments, mesh, profile):
    """Clip actual line centres to explicit-support sections before buffering.

    This excludes ordinary anatomical/box paths without sampling millions of
    points or counting the whole assembly as support. Extrusion widths come
    from the resolved profile, not the STL footprint. G-code Z is layer TOP.
    """
    h,w=profile['layer_height'],profile['line_width']
    z=np.round(segments[:,:,2].mean(axis=1),5)
    levels=np.unique(z);heights=np.where(levels<=profile['first_layer_height']+.001,profile['first_layer_height'],h)
    cuts=mesh.section_multiplane([0,0,0],[0,0,1],levels-heights/2)
    order=np.argsort(z);sorted_z=z[order];result={};kept={}
    for level,path in zip(levels,cuts):
        poly=section_polygon(path)
        if poly.is_empty:continue
        ids=order[np.searchsorted(sorted_z,level):np.searchsorted(sorted_z,level,side='right')]
        lines=shapely.linestrings(segments[ids,:,:2]);mask=shapely.intersects(lines,poly.buffer(.025))
        if not mask.any():continue
        clipped=shapely.intersection(lines[mask],poly.buffer(.025))
        width=profile.get('initial_layer_line_width',w) if level<=profile['first_layer_height']+.001 else w
        result[float(level)]=shapely.union_all(shapely.buffer(clipped,width/2,cap_style='flat'))
        kept[float(level)]=clipped
    return result,kept


def coverage(field,footprints,profile):
    rows=[];supported=[];missing=[];area=miss=0.;zs=np.array(sorted(footprints));h=profile['layer_height'];gap=profile['roof_gap_mm']
    tol=POLICY['coverage_area_tolerance_mm2']
    for r,p in zip(field['rows'],field['polygons']):
        if r['classification']!='SUPPORT_REQUIRED':continue
        # Frozen detector samples at first/2+n*h; include its +/-half-layer
        # phase and quantised gap. No free-bridge credit in V3 production.
        eligible=zs[(zs>=r['z_mm']-gap-2*h-1e-6)&(zs<=r['z_mm']-gap+h/2+1e-6)]
        actual=unary_union([footprints[float(z)] for z in eligible])
        hit=p.intersection(actual);bad=p.difference(actual);u=float(bad.area);area+=p.area;miss+=u
        supported.append((r['z_mm'],hit));missing.append((r['z_mm'],bad))
        rows.append(dict(**r,covered_area_mm2=float(hit.area),uncovered_area_mm2=u,covered=u<=tol,
                         actual_layer_top_z_mm=eligible.tolist(),uncovered_wkt=bad.wkt,supported_wkt=hit.wkt))
    return dict(percent=100*(1-miss/area) if area else 100.,required_area_sum_mm2=area,uncovered_area_sum_mm2=miss,
        uncovered_region_count=sum(not r['covered'] for r in rows),region_count=len(rows),all_required_covered=all(r['covered'] for r in rows),
        area_tolerance_mm2=tol,bridge_credit=False,method='ACTUAL_GCODE_LINE_CENTRES_CLIPPED_TO_EXPLICIT_SUPPORT_THEN_EXTRUSION_WIDTH_BUFFERED',rows=rows),supported,missing


def preserved_features(ribs,footprints,profile):
    rows=[];h=profile['layer_height'];zs=np.array(sorted(footprints))
    for rib in ribs:
        for name,part in [('RIB',rib['ribs']),('ROOF_RIDGE',rib['roofs']),('CUT_NECK',rib['necks'])]:
            levels=np.linspace(part.bounds[0,2]+h,part.bounds[1,2]-h,3) if name=='RIB' else np.array([part.bounds.mean(0)[2]])
            for j,level in enumerate(levels):
                nearby=zs[np.abs(zs-h/2-level)<h]
                cuts=part.section_multiplane([0,0,0],[0,0,1],nearby-h/2)
                nominal=unary_union([section_polygon(p) for p in cuts]);actual=unary_union([footprints[float(z)] for z in nearby])
                ratio=float(nominal.intersection(actual).area/nominal.area) if nominal.area else 0.
                rows.append(dict(feature_id=f"{rib['rib_id']}_{name}_{j}",kind=name,nominal_area_mm2=float(nominal.area),actual_overlap_ratio=ratio,
                    preserved=ratio>.5,check='Cross-section preservation; coverage gate separately requires complete required polygons'))
    return dict(passed=all(r['preserved'] for r in rows),rows=rows)


def actual_self_support(footprints,box_footprints,profile):
    """Conservative layer-growth audit of actual extruded support footprints."""
    rows=[];h=profile['layer_height'];reach=min(profile['line_width']/2,h)
    for z,current in sorted(footprints.items()):
        if z<=profile['first_layer_height']+.001:continue
        previous=unary_union([p for zz,p in footprints.items() if abs(zz-(z-h))<.001]+[p for zz,p in box_footprints.items() if abs(zz-(z-h))<.001])
        bad=current.difference(previous.buffer(reach))
        for polygon in (bad.geoms if hasattr(bad,'geoms') else [bad]):
            if polygon.area>.005:
                rows.append(dict(z_mm=z,area_mm2=float(polygon.area),representative_xy=list(polygon.representative_point().coords[0]),polygon_wkt=polygon.wkt))
    return dict(passed=not rows,unsupported_region_count=len(rows),rows=rows,area_tolerance_mm2=.005,
                method='Actual support extrusion footprints minus previous actual support and approved box footprints buffered by min(width/2, layer height)')


def contact_burden(kept,vascular,profile):
    """Same downward-face metric as V2, with the same V3 path isolation."""
    from scipy.spatial import cKDTree
    from .bambu_casting_mold import sample_segments
    from .support_qc import contact_face_metrics
    pieces=[];lo,hi=vascular.bounds
    for z,lines in kept.items():
        if z<lo[2]-.5 or z>hi[2]:continue
        for line in shapely.get_parts(lines):
            if line.geom_type!='LineString' or line.is_empty:continue
            xy=np.asarray(line.coords)
            if np.any(xy.max(0)<lo[:2]-.5) or np.any(xy.min(0)>hi[:2]+.5):continue
            xyz=np.c_[xy,np.full(len(xy),z)]
            pieces.extend(np.stack([xyz[:-1],xyz[1:]],axis=1))
    if not pieces:return dict(contact_area_mm2=0.,method='Actual isolated support paths / downward vascular triangle area')
    points=sample_segments(np.asarray(pieces),.12)
    metrics,_=contact_face_metrics(vascular,cKDTree(points),profile)
    return {k:v for k,v in metrics.items() if k not in ['contact_face_ids','patch_areas_mm2']}


def audit(record,result,profile):
    folder=Path(record['folder']);out=folder/'bambu'
    if not result.get('slicer_success'):return dict(slicer_success=False,status=result.get('status'),error=result.get('error'))
    with zipfile.ZipFile(result['archive']) as z:settings=json.loads(z.read('Metadata/project_settings.config'))
    from .support_inputs import numeric
    profile={**profile,'initial_layer_line_width':numeric(settings.get('initial_layer_line_width',profile['line_width']),profile['nozzle_diameter'])}
    paths=parse_actual_paths(result['archive']);segments=paths['segments'][paths['kind']==0]
    support=trimesh.load_mesh(folder/'custom_support_v3.stl')
    footprints,kept=extrusion_sections(segments,support,profile)
    boxpaths,_=extrusion_sections(segments,trimesh.load_mesh(folder/'box.stl'),profile)
    actual_self=actual_self_support(footprints,boxpaths,profile);table(out/'actual_support_self_regions.csv',actual_self['rows'])
    with (folder/'frozen_detection_cache.pkl').open('rb') as f:det,field,boxfield,intervals,bs=pickle.load(f)
    with (folder/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
    vc,supported,missing=coverage(field,footprints,profile);bc,_,_=coverage(boxfield,footprints,profile)
    features=preserved_features(ribs,footprints,profile)
    table(out/'actual_coverage.csv',vc['rows']);table(out/'actual_box_coverage.csv',bc['rows']);table(out/'feature_preservation.csv',features['rows'])
    with (out/'actual_support_footprints.pkl').open('wb') as f:pickle.dump(footprints,f)
    with (out/'actual_coverage_geometry.pkl').open('wb') as f:pickle.dump((supported,missing),f)
    with zipfile.ZipFile(result['archive']) as z:
        settings=json.loads(z.read('Metadata/project_settings.config'));crc=z.testzip()
    from scipy.spatial import cKDTree
    raw=det['raw'];tree=cKDTree([p['xyz'] for p in raw]);warnings=[]
    for row in vc['rows']:
        nearest=raw[int(tree.query([*row['representative_xy'],row['z_mm']])[1])]
        warnings.append(dict(region_id=row['region_id'],xyz=[*row['representative_xy'],row['z_mm']],
            classification='INTENTIONAL_AIR_GAP' if row['covered'] else 'PORT_UNSUPPORTED' if nearest['is_port'] else 'VASCULAR_UNSUPPORTED',
            uncovered_area_mm2=row['uncovered_area_mm2'],source='Independent frozen layer detector / actual G-code; native warning has no coordinates'))
    for row in record['self_buildability']['rows']:
        warnings.append(dict(region_id=row['region_id'],xyz=[*row['representative_xy'],row['z_mm']],classification='SUPPORT_UNSUPPORTED',uncovered_area_mm2=row['area_mm2']))
    for index,row in enumerate(actual_self['rows']):
        warnings.append(dict(region_id=f'actual_support_{index}',xyz=[*row['representative_xy'],row['z_mm']],classification='SUPPORT_UNSUPPORTED',uncovered_area_mm2=row['area_mm2'],source='ACTUAL_EXTRUSION_LAYER_GROWTH'))
    for row in bc['rows']:
        if not row['covered']:warnings.append(dict(region_id=row['region_id'],xyz=[*row['representative_xy'],row['z_mm']],classification='UNKNOWN',detail='Unresolved box overhang',uncovered_area_mm2=row['uncovered_area_mm2']))
    native=result.get('warning_audit',{}).get('critical',[])
    # Never clear a native critical warning solely because it has no location.
    unresolved=sum(r['classification']!='INTENTIONAL_AIR_GAP' for r in warnings)
    if native and not warnings:unresolved+=len(native);warnings.append(dict(classification='UNKNOWN',detail=str(native)))
    table(out/'warning_localization.csv',warnings)
    # Coverage is exact polygon evidence; near-vessel burden is a separate
    # descriptive projection, not a percentage improvement claim against V2.
    burden=float(unary_union([p for _,p in supported]).area)
    contact=contact_burden(kept,trimesh.load_mesh(folder/'vascular_with_ports.stl'),profile)
    result=dict(slicer_success=True,auto_support_off=str(settings.get('enable_support')) in ['0','False','false'],native_support_toolpaths=int((paths['kind']==1).sum()),
        actual_vascular_coverage={k:v for k,v in vc.items() if k!='rows'},actual_box_coverage={k:v for k,v in bc.items() if k!='rows'},
        all_features_preserved=features['passed'],feature_preservation=features,actual_self_buildability=actual_self,unresolved_warning_count=unresolved,native_warnings=native,
        actual_near_vessel_projection_mm2=burden,contact_area_mm2=contact['contact_area_mm2'],contact_metric=contact,estimated_print_time_seconds=result.get('estimated_print_time_seconds'),
        actual_profile={k:settings.get(k) for k in ['layer_height','line_width','nozzle_diameter','enable_support','filament_type']},
        archive_crc_passed=crc is None,production_geometry=result.get('production_geometry'),gcode_model_segments=len(segments))
    save(out/'bambu_validation.json',result)
    return result


def slice_candidate(record,profile,discovery):
    folder=Path(record['folder']);engine=discover_engine('Bambu',discovery,folder)
    profiles=engine_profiles(discovery,'Bambu',{'enable_support':'0'})
    source=folder/'mold_with_custom_support_v3.stl'
    r=sliced(engine,profiles,source,source,folder/'bambu',900)
    print('V3 BAMBU FINISHED',record['candidate_id'],r.get('slicer_success'),flush=True)
    key=hashlib.sha256(json.dumps(dict(source=digest(source),support=digest(folder/'custom_support_v3.stl'),ribs=digest(folder/'rib_geometry_cache.pkl'),
        archive=digest(r['archive']) if r.get('slicer_success') else None,profile=profile,code=digest(Path(__file__)),coverage_policy=POLICY),sort_keys=True).encode()).hexdigest()
    cached=folder/'bambu_validation.json'
    if cached.exists():
        previous=json.loads(cached.read_text())
        if previous.get('cache_fingerprint')==key:
            print('V3 ACTUAL AUDIT CACHE',record['candidate_id'],flush=True);return previous
    metrics=audit(record,r,profile)
    if r.get('slicer_success'):shutil.copy2(r['archive'],folder/'X1C_project.3mf')
    support=trimesh.load_mesh(folder/'custom_support_v3.stl');box=trimesh.load_mesh(folder/'box.stl')
    path=folder/'support_with_box_foundation.stl';union([support,box]).export(path)
    solo=sliced(engine,profiles,path,path,folder/'bambu_support_only',900)
    metrics['support_only_slicer_success']=solo.get('slicer_success',False)
    metrics['cache_fingerprint']=key
    save(folder/'bambu_validation.json',metrics)
    print('V3 ACTUAL COVERAGE',record['candidate_id'],metrics.get('actual_vascular_coverage',{}).get('percent'),flush=True)
    return metrics


def thin_wall_trial(out,profile,discovery):
    out.mkdir(parents=True,exist_ok=True);ribs=[]
    for lines in [1,2,3]:
        xy=np.c_[np.linspace(80,110,61),np.full(61,80+lines*12.)]
        rib=curtain(xy,np.full(len(xy),12.),0,profile,lines,(lines+2)*profile['line_width'],5.,3.)
        rib['rib_id']=f'{lines}_LINE';ribs.append(rib);rib['mesh'].export(out/f'rib_{lines}_line.stl')
    path=out/'thin_wall_preservation_trial.stl';union([r['mesh'] for r in ribs]).export(path)
    engine=discover_engine('Bambu',discovery,out);profiles=engine_profiles(discovery,'Bambu',{'enable_support':'0'})
    result=sliced(engine,profiles,path,path,out/'bambu',900)
    checks=[]
    if result.get('slicer_success'):
        paths=parse_actual_paths(result['archive']);segments=paths['segments'][paths['kind']==0]
        footprints,_=extrusion_sections(segments,trimesh.load_mesh(path),profile)
        checks=preserved_features(ribs,footprints,profile)['rows']
    perline={i:bool(checks) and all(r['preserved'] for r in checks if r['feature_id'].startswith(f'{i}_LINE_')) for i in [1,2,3]}
    selected=next((i for i in [1,2,3] if perline[i]),None)
    result=dict(slicer_success=result.get('slicer_success',False),line_preservation=perline,selected_minimum_preserved_lines=selected,
        one_line_upgrade_required=not perline[1],checks=checks,physical_strength_validated=False,physical_status='NO PHYSICAL COUPON AVAILABLE')
    save(out/'thin_wall_result.json',result);table(out/'feature_preservation.csv',checks)
    return result
