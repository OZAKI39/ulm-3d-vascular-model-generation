"""Semantic actual-Gcode coverage and fixed-support release gates."""
from __future__ import annotations
import json,pickle,copy,shutil,zipfile,io,re
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from .support_inputs import save,table,transform
from .support_v3_validation import slice_candidate,coverage,extrusion_sections
from .support_v3_geometry import union
from .support_v3_1_orientation import vertical_clear
from .support_v3_1_release import plate_release
from .support_v3_removal import run as fragment_removal


def semantic_fields(folder,inputs):
    with (folder/'frozen_detection_cache.pkl').open('rb') as f:det,field,bf,intervals,branches=pickle.load(f)
    xyz=np.vstack([b['xyz'] for b in branches]);flag=np.concatenate([np.full(len(b['xyz']),b['is_port']) for b in branches]);tree=cKDTree(xyz)
    labels=[]
    for row in field['rows']:
        idx=tree.query([*row['representative_xy'],row['z_mm']])[1];labels.append('PORT' if flag[idx] else 'VASCULAR')
    record=json.loads((folder/'candidate.json').read_text());restore=np.asarray(record['orientation']['casting_restore_transform'])
    box=trimesh.load_mesh(inputs['box']);bottom=box.triangles_center[:,2]<=inputs['inner_bounds'][0][2]+1e-4
    bp=np.array([[*r['representative_xy'],r['z_mm']] for r in bf['rows']]).reshape(-1,3)
    if len(bp):
        _,_,face=trimesh.proximity.closest_point(box,transform(bp,restore));blabels=np.where(bottom[face],'BOX_BOTTOM','BOX_WALLS').tolist()
    else:blabels=[]
    result={}
    for name,base,ls in [('VASCULAR',field,labels),('PORT',field,labels),('BOX_BOTTOM',bf,blabels),('BOX_WALLS',bf,blabels)]:
        result[name]=dict(rows=[{**r,'component':name} for r,l in zip(base['rows'],ls) if l==name],polygons=[p for p,l in zip(base['polygons'],ls) if l==name])
    return result,field,bf,branches


def semantic_audit(record,inputs,profile):
    folder=Path(record['folder']);fields,field,bf,branches=semantic_fields(folder,inputs)
    with (folder/'bambu/actual_support_footprints.pkl').open('rb') as f:fp=pickle.load(f)
    from .support_v3_1_residual import area_tolerance
    metrics={};rows=[];geometry={}
    for name,fi in fields.items():
        c,hit,miss=coverage(fi,fp,profile)
        for rr in c['rows']:
            tol=area_tolerance(profile['line_width'],max(abs(x) for x in rr['bounds']))
            rr['area_tolerance_mm2']=tol;rr['covered']=rr['uncovered_area_mm2']<=tol
        c['all_required_covered']=all(rr['covered'] for rr in c['rows']);c['uncovered_region_count']=sum(not rr['covered'] for rr in c['rows'])
        c['tolerance_derivation']='min(1e-5 mm² legacy ceiling, 4 * float32 coordinate ULP * actual line width)'
        metrics[name]={k:v for k,v in c.items() if k!='rows'};rows+=c['rows'];geometry[name]=(hit,miss)
    table(folder/'actual_coverage.csv',rows);table(folder/'vascular_residuals.csv',[r for r in rows if not r['covered']]);
    table(folder/'support_self_residuals.csv',record['slice']['actual_self_buildability']['rows'])
    with (folder/'semantic_coverage_geometry.pkl').open('wb') as f:pickle.dump(geometry,f)
    record['semantic_coverage']=metrics
    # Recompute direct opportunity over frozen required layer polygons, not
    # over the earlier face proxy. Weight by each polygon's area.
    required=[r for r in field['rows'] if r['classification']=='SUPPORT_REQUIRED']
    points=np.array([[*r['representative_xy'],r['z_mm']] for r in required]).reshape(-1,3)
    weights=np.array([r['area_mm2'] for r in required]);model=trimesh.load_mesh(folder/'mold_without_support.stl')
    clear=vertical_clear(model,points,profile['roof_gap_mm'])
    record['direct_buildplate_support_fraction']=float(weights[clear].sum()/weights.sum()) if weights.sum() else 1.
    table(folder/'direct_buildplate_support.csv',[dict(region_id=r['region_id'],area_mm2=r['area_mm2'],direct_ray_clear=bool(ok),scope='REPRESENTATIVE_RAY_RANKING_ONLY') for r,ok in zip(required,clear)])
    from .support_qc import parse_actual_paths
    paths=parse_actual_paths(folder/'X1C_project.3mf');segments=paths['segments'][paths['kind']==0]
    from .support_inputs import numeric
    with zipfile.ZipFile(folder/'X1C_project.3mf') as archive:settings=json.loads(archive.read('Metadata/project_settings.config'))
    first_width=numeric(settings.get('initial_layer_line_width',profile['line_width']),profile['nozzle_diameter'])
    first={};firstpolys={}
    for name,file in [('BOX','box.stl'),('CUSTOM_SUPPORT','custom_support_v3.stl')]:
        mesh=trimesh.load_mesh(folder/file);foot,_=extrusion_sections(segments[abs(segments[:,:,2].mean(1)-profile['first_layer_height'])<.001],mesh,{**profile,'initial_layer_line_width':first_width})
        poly=union_polys(foot.values());firstpolys[name]=poly;first[name]=dict(actual_contact_area_mm2=float(poly.area),islands=len(poly.geoms) if hasattr(poly,'geoms') else int(not poly.is_empty))
    brim=actual_brim(folder/'X1C_project.3mf',profile,first_width)
    first['brim_connection']={name:dict(brim_present=not brim.is_empty,minimum_gap_mm=float(poly.distance(brim)) if not brim.is_empty and not poly.is_empty else None,
        footprint_connected=bool(not brim.is_empty and not poly.is_empty and poly.distance(brim)<=4*float(np.spacing(np.float32(256))))) for name,poly in firstpolys.items()}
    first['actual_brim_area_mm2']=float(brim.area);first['adhesion_strength_validated']=False
    record['actual_first_layer']=first
    # RS[n] is the frozen geometry detector's support-required polygon, while
    # CS[n] comes from the previous ACTUAL support/approved-box toolpaths.
    # Keep this distinct from current-actual versus previous-actual growth.
    from .support_v2_layers import analyse
    from .support_v3_1_residual import difference
    support=trimesh.load_mesh(folder/'custom_support_v3.stl');boxmesh=trimesh.load_mesh(folder/'box.stl')
    sf=analyse(support,profile,'support',boxmesh,bridge_limit=0.) if record['self_buildability']['unsupported_region_count'] else dict(rows=[],polygons=[])
    geometry_rows=[]
    if any(r['classification']=='SUPPORT_REQUIRED' for r in sf['rows']):
        boxfp,_=extrusion_sections(segments,boxmesh,{**profile,'initial_layer_line_width':first_width})
        for rr,poly in zip(sf['rows'],sf['polygons']):
            if rr['classification']!='SUPPORT_REQUIRED':continue
            # Frozen sections are layer centres, actual G-code Z is its top.
            previous_top=profile['first_layer_height']+(rr['layer_index']-1)*profile['layer_height']
            actual=union_polys([p for z,p in fp.items() if abs(z-previous_top)<.001]+[p for z,p in boxfp.items() if abs(z-previous_top)<.001])
            bad=difference(poly,actual.buffer(min(profile['line_width']/2,profile['layer_height'])))
            geometry_rows.append(dict(**rr,uncovered_area_mm2=float(bad.area),uncovered_wkt=bad.wkt,actual_previous_layer_top_mm=previous_top,
                covered=bad.area<=area_tolerance(profile['line_width'],max(abs(x) for x in poly.bounds))))
    table(folder/'support_self_layer_closure.csv',geometry_rows)
    record['geometry_support_self_actual_coverage']=dict(required_regions=len(geometry_rows),unsupported_region_count=sum(not x['covered'] for x in geometry_rows),rows=geometry_rows,
        method='RS[n] frozen support-required polygons minus CS[n] previous actual support and box extrusion with frozen growth reach')
    save(folder/'semantic_validation.json',dict(coverage=metrics,first_layer=first,direct_buildplate_support_fraction=record['direct_buildplate_support_fraction'],semantic_method='Nearest frozen SWC branch for vessel/port; nearest casting box face for bottom/walls'))
    return record


def actual_brim(path,profile,width):
    """Keep original coordinates/extrusion; filter only original FEATURE roles.

    Reuse the established G-code interpreter for offsets/arcs/extrusion modes.
    The in-memory view never replaces or modifies the real sliced archive.
    """
    from .bambu_casting_mold import parse_paths
    from shapely.geometry import LineString
    stream=io.BytesIO()
    with zipfile.ZipFile(path) as original,zipfile.ZipFile(stream,'w') as filtered:
        for name in original.namelist():
            if name=='Metadata/project_settings.config':filtered.writestr(name,original.read(name))
            elif name.endswith('.gcode'):
                lines=original.read(name).decode(errors='replace').splitlines()
                lines=[line if not line.startswith('; FEATURE:') or 'brim' in line.lower() else '; FEATURE:' for line in lines]
                filtered.writestr(name,'\n'.join(lines))
    stream.seek(0);paths=parse_paths(stream);segs=paths['segments'];segs=segs[abs(segs[:,:,2].mean(1)-profile['first_layer_height'])<.001]
    return union_polys(LineString(s[:,:2]).buffer(width/2,cap_style=2) for s in segs)


def union_polys(items):
    from shapely.ops import unary_union
    return unary_union(list(items))


def solid_components(mesh):
    """Group negative cavity shells with their containing positive solid.

    Trimesh.split returns surface-connected shells, not material components;
    treating a negative internal void as a free support fragment is wrong.
    All original triangles/coordinates are preserved in this partition.
    """
    shells=list(mesh.split(only_watertight=False));positive=[p for p in shells if p.volume>0]
    groups=[[p] for p in positive]
    for cavity in [p for p in shells if p.volume<=0]:
        point=cavity.vertices.mean(0)
        owners=[i for i,p in enumerate(positive) if np.all(point>=p.bounds[0]-1e-6) and np.all(point<=p.bounds[1]+1e-6) and p.contains([point])[0]]
        if not owners:raise ValueError('CAVITY_SHELL_OWNER_UNRESOLVED')
        owner=min(owners,key=lambda i:positive[i].volume);groups[owner].append(cavity)
    return [trimesh.util.concatenate(g) for g in groups]


def removal_audit(record,inputs):
    folder=Path(record['folder'])
    with (folder/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
    from .support_custom_qc import intersection_volume
    from .support_v2_mesh import exact_stl
    # Partition ACTUAL exported connected components, not pre-union ribs.
    # A ceiling rib fused to a plate rib is fixed with it; pretending these
    # are independent pieces would incorrectly certify release.
    actual=trimesh.load_mesh(folder/'custom_support_v3.stl')
    components=solid_components(actual)
    plate_parts=[p for p in components if p.bounds[0,2]<.001]
    other_parts=[p for p in components if p.bounds[0,2]>=.001]
    fixed=trimesh.util.concatenate(plate_parts) if plate_parts else trimesh.Trimesh()
    nonplate=trimesh.util.concatenate(other_parts) if other_parts else trimesh.Trimesh()
    other=[]
    for r in ribs:
        if len(nonplate.faces) and intersection_volume(r['mesh'],nonplate)>1e-5:
            local=copy.deepcopy(r)
            local['cuts']=[c for c in local['cuts'] if intersection_volume(c['intended'],nonplate)>1e-7]
            other.append(local)
    fixed.export(folder/'plate_rooted_support.stl')
    result=plate_release(trimesh.load_mesh(folder/'mold_without_support.stl'),fixed)
    result.update(plate_rooted_components=len(plate_parts),nonplate_components=len(other_parts),nonplate_ribs=len(other),partition='ACTUAL_EXPORTED_CONNECTED_COMPONENTS');save(folder/'plate_release.json',result);record['plate_release']=result
    if not other:
        removal=dict(trapped_pieces=0,tool_access_failed=0,tool_access_passed=True,all_internal_exit_top=True,fragments=[],scope='ALL_ROOTS_BUILD_PLATE; separate plate release gate is mandatory')
    else:
        dest=folder/'nonplate_removal';dest.mkdir(exist_ok=True)
        with (dest/'rib_geometry_cache.pkl').open('wb') as f:pickle.dump(other,f)
        # Concatenating already disjoint exported components preserves their
        # coordinates exactly and avoids introducing new Boolean/STL seams.
        exact_stl(nonplate).export(dest/'custom_support_v3.stl')
        sub={**record,'folder':str(dest)};removal=fragment_removal(sub,inputs)
        removal['scope']='Nonplate ribs after separation from fixed plate supports; meaningful only if plate-release gate passes'
    record['removal']=removal;save(folder/'support_removal_plan.json',removal)
    return record


def gates(record,phase_a=True):
    s=record.get('slice',{});rem=record.get('removal',{});cov=record.get('semantic_coverage',{})
    checks={
        'ORIENTATION_LABEL_MATRIX_MISMATCH':record['orientation'].get('label_matrix_matches',True),
        'BUILD_VOLUME_EXCEEDED':record.get('build_volume_fit',False),
        'SOURCE_PROTECTION_NOT_VERIFIED':record.get('source_unchanged',False),
        'NON_RIGID_TRANSFORM':record.get('scale_factor')==1.,
        'ACTUAL_SUPPORT_SELF_UNSUPPORTED':s.get('actual_self_buildability',{}).get('unsupported_region_count')==0,
        'STL_SUPPORT_SELF_UNSUPPORTED':record.get('self_buildability',{}).get('unsupported_region_count')==0,
        'GEOMETRY_SUPPORT_SELF_ACTUAL_GAP':record.get('geometry_support_self_actual_coverage',{}).get('unsupported_region_count',0)==0,
        'ILLEGAL_INTERSECTIONS':record.get('vascular_intersection_mm3',1.)<1e-5 and record.get('exported_geometry_audit',{}).get('illegal_box_intersection_mm3',1.)<1e-4,
        'EXPORTED_SUPPORT_GEOMETRY_AUDIT_FAILED':record.get('exported_geometry_audit',{}).get('passed',False),
        'AUTO_SUPPORT_NOT_OFF':s.get('auto_support_off') and s.get('native_support_toolpaths')==0,
        'CRITICAL_FEATURE_DROPPED':s.get('all_features_preserved',False),
        'UNRESOLVED_CRITICAL_WARNINGS':s.get('unresolved_warning_count',-1)==0,
        'REAL_BAMBU_SLICE_FAILED':s.get('slicer_success',False) and s.get('support_only_slicer_success',False),
        'TRAPPED_FRAGMENT':rem.get('trapped_pieces')==0,
        'CUT_POINT_TOOL_ACCESS_FAILED':rem.get('tool_access_failed')==0,
        'CUT_SECTION_CROSSES_UNWEAKENED_SUPPORT':record.get('cuts_valid',False),
        'ASSEMBLY_STL_PRECISION_FAILURE':not record.get('assembly_export_error')}
    for name in ('VASCULAR','PORT','BOX_BOTTOM','BOX_WALLS'):checks['ACTUAL_'+name+'_COVERAGE_INCOMPLETE']=cov.get(name,{}).get('all_required_covered',False)
    if phase_a:checks['PLATE_RELEASE_FAILED']=record.get('plate_release',{}).get('passed',False)
    return [name for name,value in checks.items() if not value]
