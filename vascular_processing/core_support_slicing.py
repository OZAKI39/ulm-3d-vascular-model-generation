"""Native Bambu slicing and strict actual-extrusion coverage for s1-11."""
from __future__ import annotations
import copy,json,pickle,shutil,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import shapely
from shapely.geometry import Point,GeometryCollection
from shapely.ops import unary_union
from .support_inputs import save,table,digest
from .all_abs_print_plate import write_project as write_objects,NS
from .support_slicer_adapter import slice_local,production_geometry
from .support_rfus_travel import events
from .support_detection import contiguous_runs,section_polygon


def clip_extrusion_lines(lines,region):
    """Exact clipping with prepared predicates for fully in/out paths.

    Avoid repeating a large polygon overlay for every already-contained
    extrusion. Only boundary-crossing lines need GEOS intersection; no
    tolerance or geometric approximation is introduced by this shortcut.
    """
    shapely.prepare(region)
    contained=shapely.covers(region,lines);hits=shapely.intersects(region,lines)
    result=np.asarray(lines,dtype=object).copy();result[~hits]=GeometryCollection()
    mixed=hits&~contained
    if mixed.any():result[mixed]=shapely.intersection(lines[mixed],region)
    return result


def actual_sections(parsed,mesh,profile,progress_label=None):
    """Round bead footprints from actual positive extrusion and variable width.

    Individual G1 subdivisions are not flat-ended physical beads: a continuous
    turn retains the circular nozzle sweep at each joining endpoint.
    """
    keep=np.array([not any(k in t.lower() for k in ['brim','skirt','support','wipe','prime','custom']) for t in parsed['features']])
    seg=parsed['extrusion'][keep];width=parsed['extrusion_widths'][keep];zs=np.round(seg[:,:,2].mean(1),5);levels=np.unique(zs);h=profile['layer_height'];first=profile['first_layer_height'];result={}
    if progress_label:print('ACTUAL_SECTIONS_BEGIN',progress_label,'layers',len(levels),'faces',len(mesh.faces),flush=True)
    paths=mesh.section_multiplane([0,0,0],[0,0,1],levels-np.where(levels<=first+.001,first,h)/2)
    if progress_label:print('ACTUAL_SECTIONS_PLANES_READY',progress_label,flush=True)
    for index,(z,path) in enumerate(zip(levels,paths)):
        if progress_label and index%25==0:print('ACTUAL_SECTIONS_LAYER',progress_label,index,float(z),flush=True)
        poly=section_polygon(path)
        if poly.is_empty:continue
        ix=np.flatnonzero(zs==z);lines=shapely.linestrings(seg[ix,:,:2]);clipped=clip_extrusion_lines(lines,poly.buffer(.001));valid=~shapely.is_empty(clipped)
        if valid.any():result[float(z)]=shapely.union_all(shapely.buffer(clipped[valid],width[ix][valid]/2,cap_style='round',quad_segs=16))
    return result


def coordinate_tolerance(point,profile):
    # XY positions in this native archive are emitted on a 0.001 mm grid.
    # Half a coordinate cell per axis plus float32 STL precision, recorded in
    # the audit; this is unrelated to a desired coverage percentage.
    return 4*float(np.spacing(np.float32(max(1.,*np.abs(point)))))+profile.get('gcode_xy_half_cell_mm',0.)*np.sqrt(2)


def write_project(path,vascular,support,profiles):
    """One assembly with two normal parts preserves the deliberate Z gap."""
    ET.register_namespace('',NS)
    records=[dict(name='FROZEN_VASCULAR_CORE',support_enabled=False),dict(name='EXPLICIT_REMOVABLE_SUPPORT',support_enabled=False)]
    write_objects(path,[vascular,support],records,profiles)
    with zipfile.ZipFile(path) as z:data={n:z.read(n) for n in z.namelist()}
    root=ET.fromstring(data['3D/3dmodel.model']);resources=root.find(f'{{{NS}}}resources');build=root.find(f'{{{NS}}}build')
    for item in list(build):build.remove(item)
    group=ET.SubElement(resources,f'{{{NS}}}object',{'id':'10','type':'model','name':'BG001 vascular core and removable support'})
    components=ET.SubElement(group,f'{{{NS}}}components')
    for oid in ['1','3']:ET.SubElement(components,f'{{{NS}}}component',{'objectid':oid})
    ET.SubElement(build,f'{{{NS}}}item',{'objectid':'10','printable':'1'})
    cfg=ET.fromstring(data['Metadata/model_settings.config']);old=cfg.findall('object');obj=ET.Element('object',{'id':'10'})
    for k,v in dict(name='BG001 vascular core + custom support',extruder='1',enable_support='0').items():ET.SubElement(obj,'metadata',{'key':k,'value':v})
    for o in old:
        for part in o.findall('part'):obj.append(part)
        cfg.remove(o)
    cfg.insert(0,obj);plate=cfg.find('plate')
    for mi in plate.findall('model_instance'):plate.remove(mi)
    inst=ET.SubElement(plate,'model_instance')
    for k,v in dict(object_id='10',instance_id='0',identify_id='1').items():ET.SubElement(inst,'metadata',{'key':k,'value':v})
    data['3D/3dmodel.model']=ET.tostring(root,encoding='utf-8',xml_declaration=True);data['Metadata/model_settings.config']=ET.tostring(cfg,encoding='utf-8',xml_declaration=True)
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
        for n,b in data.items():z.writestr(n,b)


def actual_underbody(branches,footprints,profile):
    zs=np.array(sorted(footprints));h=profile['layer_height'];gap=profile['roof_gap_mm'];rows=[];cache={};maxarc=0.
    for b in branches:
        flags=[]
        for i,low in enumerate(b['lower']):
            valid=bool(np.isfinite(low).all());ulp=coordinate_tolerance(low,profile) if valid else 0
            eligible=zs[(zs>=low[2]-gap-h-ulp)&(zs<=low[2]-gap+ulp)] if valid else np.array([])
            key=tuple(eligible)
            if key not in cache:cache[key]=unary_union([footprints[float(z)] for z in eligible])
            fp=cache[key];distance=float(fp.distance(Point(low[:2]))) if valid and not fp.is_empty else None
            ok=distance is not None and distance<=ulp
            higher=zs[(zs>low[2]-gap+ulp)&(zs<=low[2]+ulp)] if valid else np.array([])
            too_close=[float(z) for z in higher if footprints[float(z)].distance(Point(low[:2]))<=ulp]
            ok=ok and not too_close;flags.append(ok)
            rows.append(dict(branch_id=b['branch_id'],is_port=bool(b['is_port']),sample_index=i,arc_mm=float(b['arc'][i]),lower_surface=low.tolist() if valid else None,supported=ok,actual_custom_support=ok,eligible_layer_tops=eligible.tolist(),too_close_support_layer_tops=too_close,distance_mm=distance,float_tolerance_mm=ulp,gap_mm=gap,layer_quantization_mm=h))
        for a,z in contiguous_runs(~np.asarray(flags)):
            arc=b['arc'];lo=arc[a] if a==0 else (arc[a-1]+arc[a])/2;hi=arc[z] if z==len(arc)-1 else (arc[z]+arc[z+1])/2;maxarc=max(maxarc,float(hi-lo))
    ports=[r for r in rows if r['is_port']];anatomy=[r for r in rows if not r['is_port']]
    return dict(passed=bool(rows) and all(r['supported'] for r in rows),centerline_underbody_percent=100*sum(r['supported'] for r in rows)/max(1,len(rows)),port_underbody_percent=100*sum(r['supported'] for r in ports)/max(1,len(ports)),anatomical_underbody_percent=100*sum(r['supported'] for r in anatomy)/max(1,len(anatomy)),sample_count=len(rows),branch_count=len(branches)-sum(b['is_port'] for b in branches),port_count=sum(b['is_port'] for b in branches),unsupported_samples=sum(not r['supported'] for r in rows),max_uncovered_arc_mm=maxarc,maximum_sampling_step_mm=max(np.diff(b['arc']).max() for b in branches),no_vascular_self_support_credit=True,rows=rows)


def actual_self(footprints,profile):
    h=profile['layer_height'];reach=min(profile['line_width']/2,h);rows=[];allarea=0.
    for z,poly in sorted(footprints.items()):
        if z<=profile['first_layer_height']+.001:continue
        prev=unary_union([p for zz,p in footprints.items() if abs(zz-(z-h))<1e-5]);bad=poly.difference(prev.buffer(reach))
        tol=coordinate_tolerance(poly.bounds,profile)
        for p in shapely.get_parts(bad):
            if p.is_empty:continue
            # Only perimeter x float32 error is disregarded, never percentage.
            area_tol=p.length*tol+np.pi*tol**2
            if p.area>area_tol:rows.append(dict(z_mm=z,area_mm2=float(p.area),area_tolerance_mm2=area_tol,wkt=p.wkt));allarea+=p.area
    return dict(passed=not rows,self_unsupported_regions=len(rows),unsupported_area_mm2=float(allarea),rows=rows,method='Actual custom support footprints minus previous custom support footprints buffered by min(line_width/2, layer_height); no vascular foundation')


def layer_coverage(vascular,footprints,profile):
    h=profile['layer_height'];first=profile['first_layer_height'];gap=profile['roof_gap_mm'];tops=np.arange(first,vascular.bounds[1,2]+h,h);cuts=vascular.section_multiplane([0,0,0],[0,0,1],tops-np.where(tops<=first+.001,first,h)/2);prev=GeometryCollection();zs=np.array(sorted(footprints));rows=[];area=missing=0.
    # Only exposed layer tops are roofs. A wall's lower interior layers must
    # not receive support credit when its actual top is too close to a vessel.
    roof_footprints={}
    for q,fp in footprints.items():
        above=unary_union([p for zz,p in footprints.items() if abs(zz-q-h)<1e-5])
        roof_footprints[q]=fp.difference(above)
    triangles=vascular.triangles[vascular.face_normals[:,2]<-1e-8]
    zlo=triangles[:,:,2].min(1);zhi=triangles[:,:,2].max(1)
    def projected_band(lo,hi):
        from shapely.geometry import Polygon
        polygons=[]
        for triangle in triangles[(zlo<=hi)&(zhi>=lo)]:
            pts=list(triangle)
            for level,above in [(lo,True),(hi,False)]:
                old=pts;pts=[]
                if not old:break
                for a,b in zip(old,old[1:]+old[:1]):
                    ai=a[2]>=level if above else a[2]<=level;bi=b[2]>=level if above else b[2]<=level
                    if ai:pts.append(a)
                    if ai!=bi:pts.append(a+(b-a)*(level-a[2])/(b[2]-a[2]))
            if len(pts)>=3:
                p=Polygon(np.asarray(pts)[:,:2])
                if p.area>1e-14:polygons.append(p)
        return unary_union(polygons)
    for z,path in zip(tops,cuts):
        poly=section_polygon(path);required=poly.difference(prev.buffer(min(h,profile['line_width']/2)));prev=poly
        if z<=first+.001 or required.is_empty:continue
        # A newly occupied XY region enters the frozen surface between the
        # previous and current slicing planes, NOT necessarily at z-h/2.
        # Clip actual downward triangles to each exact eligible surface-Z band
        # and project them; no arbitrary expanded gap window is credited.
        tol=coordinate_tolerance(poly.bounds,profile);eligible=zs[(zs+gap+h>=z-1.5*h-tol)&(zs+gap<=z-h/2+tol)]
        pieces=[]
        for q in eligible:
            lo=max(q+gap-tol,z-1.5*h-tol);hi=min(q+gap+h+tol,z-h/2+tol)
            pieces.append(roof_footprints[float(q)].intersection(projected_band(lo,hi)))
        actual=unary_union(pieces);bad=required.difference(actual);numerical=required.length*tol+np.pi*tol**2
        passed=bad.area<=numerical;rows.append(dict(z_mm=float(z),required_mm2=float(required.area),missing_mm2=float(bad.area),tolerance_mm2=numerical,passed=passed,required_wkt=required.wkt,uncovered_wkt=bad.wkt,eligible_layer_tops=eligible.tolist()));area+=required.area;missing+=bad.area
    return dict(passed=bool(rows) and all(r['passed'] for r in rows),percent=max(0.,100*(1-missing/area)) if area else 100.,required_area_sum_mm2=float(area),uncovered_area_sum_mm2=float(missing),uncovered_region_count=sum(not r['passed'] for r in rows),rows=rows,bridge_credit=False,method='Entire required polygons intersect actual extrusion footprints AND exact frozen downward-triangle gap bands between consecutive slice planes')


def vascular_features(branches,footprints,profile):
    zs=np.array(sorted(footprints));rows=[];h=profile['layer_height']
    for b in branches:
        for i,xyz in enumerate(b['xyz']):
            eligible=zs[np.abs(zs-h/2-xyz[2])<=h/2+1e-5];fp=unary_union([footprints[float(z)] for z in eligible]);d=float(fp.distance(Point(xyz[:2]))) if not fp.is_empty else None
            # A centreline is not an infill requirement; at least one real
            # vascular extrusion must lie inside the local physical radius.
            ok=d is not None and d<=b['radius'][i]+2e-5;rows.append(dict(branch_id=b['branch_id'],sample_index=i,preserved=ok,distance_mm=d))
    return dict(passed=bool(rows) and all(r['preserved'] for r in rows),missing_samples=sum(not r['preserved'] for r in rows),rows=rows)


def evaluate(record,actual,engine,reuse_native=False):
    folder=Path(record['folder']);out=folder/'bambu';out.mkdir(exist_ok=True)
    with (folder/'geometry.pkl').open('rb') as f:g=pickle.load(f)
    path=out/'input.3mf'
    if not reuse_native:write_project(path,g['vascular'],g['support'],actual)
    elif not path.exists():raise FileNotFoundError('NATIVE_RECOVERY_REQUIRES_ORIGINAL_INPUT_PROJECT')
    signature=dict(project=digest(path),profiles=actual,binary=engine['binary_sha256']);cache=out/'slice_inputs.json';reuse=cache.exists() and json.loads(cache.read_text())==signature and (out/'slice_result.json').exists();rec=json.loads((out/'slice_result.json').read_text()) if reuse else {}
    if reuse_native and (not reuse or not rec.get('archive') or not Path(rec['archive']).exists()):raise ValueError('NATIVE_RECOVERY_INPUTS_OR_ARCHIVE_MISMATCH_NO_RESLICE')
    if not rec.get('archive') or not Path(rec['archive']).exists():
        # The shared adapter derives its Windows staging name from out.name.
        # A candidate-specific directory prevents simultaneous fixtures or
        # candidates from overwriting each other's native slicing inputs.
        native_out=out/('s1_11_'+folder.name)
        rec=slice_local(engine,actual,path,native_out)
        for item in native_out.iterdir():
            if item.is_file():shutil.copy2(item,out/item.name)
        if rec.get('archive'):rec['archive']=str(out/'sliced.3mf')
        save(out/'slice_result.json',rec)
    save(cache,signature);print('SLICE',folder.name,rec.get('slicer_success'),rec.get('warning_audit'),flush=True)
    if not rec.get('slicer_success'):save(folder/'evaluation.json',dict(passed=False,slice=rec));return dict(passed=False,slice=rec)
    import re
    with zipfile.ZipFile(rec['archive']) as z:
        code=z.read(next(n for n in z.namelist() if n.endswith('.gcode'))).decode()
    motion='\n'.join(line for line in code.splitlines() if re.match(r'^G[0123]\s',line))
    decimals=[len(m.group(1)) for m in re.finditer(r'(?:^|\s)[XY][-+]?\d*\.(\d+)',motion)]
    digits=max(decimals) if decimals else 0
    g['profile']['gcode_xy_half_cell_mm']=.5*10.**(-digits)
    save(folder/'gcode_precision.json',dict(maximum_xy_decimal_places=digits,half_coordinate_cell_mm=g['profile']['gcode_xy_half_cell_mm'],footprint_cap='ROUND_NOZZLE_SWEEP'))
    ev=events(rec['archive']);fp=actual_sections(ev,g['support'],g['profile'],'support');vf=actual_sections(ev,g['vascular'],g['profile'],'vascular')
    with (folder/'actual_footprints.pkl').open('wb') as f:pickle.dump(dict(support=fp,vascular=vf),f)
    cov=actual_underbody(g['branches'],fp,g['profile']);selfcheck=actual_self(fp,g['profile']);layer=layer_coverage(g['vascular'],fp,g['profile']);feat=vascular_features(g['branches'],vf,g['profile'])
    for name,obj in [('coverage',cov),('support_self',selfcheck),('layer_coverage',layer),('vascular_features',feat)]:save(folder/(name+'.json'),obj)
    with zipfile.ZipFile(rec['archive']) as z:
        cfg=json.loads(z.read('Metadata/project_settings.config'));mc=ET.fromstring(z.read('Metadata/model_settings.config'));parts=mc.findall('.//part');crc=z.testzip() is None
    geom=production_geometry(rec['archive'],folder/'vascular_core_plus_support_review.stl');off=str(cfg.get('enable_support'))=='0' and not any('support' in tag.lower() for tag in ev['features'])
    # A floating warning is explained only when all corresponding targets AND
    # layer regions are covered. Unknown warning classes remain unresolved.
    warnings=rec['warning_audit']['critical'];unresolved=[w for w in warnings if not ('floating' in w.lower() and cov['passed'] and layer['passed'] and selfcheck['passed'])]
    localization=dict(native_warnings=rec['warning_audit']['all'],unresolved_critical=unresolved,intentional_gap_explanation_accepted=bool(warnings) and not unresolved,uncovered_samples=[dict(branch_id=r['branch_id'],sample_index=r['sample_index'],xyz=r['lower_surface']) for r in cov['rows'] if not r['supported']],uncovered_layer_regions=[dict(z_mm=r['z_mm'],representative_xy=list(shapely.from_wkt(r['uncovered_wkt']).representative_point().coords[0]),area_mm2=r['missing_mm2']) for r in layer['rows'] if not r['passed']],self_unsupported_regions=[dict(z_mm=r['z_mm'],representative_xy=list(shapely.from_wkt(r['wkt']).representative_point().coords[0]),area_mm2=r['area_mm2']) for r in selfcheck['rows']],method='Bambu emits an object-level floating warning; its unresolved regions are localized by independent actual extrusion and frozen-surface audits, not assumed from the warning text')
    save(folder/'warning_localization.json',localization)
    result=dict(passed=cov['passed'] and selfcheck['passed'] and layer['passed'] and feat['passed'] and off and crc and geom['passed'] and not unresolved,slice=rec,auto_support_off=off,independent_part_count=len(parts),archive_geometry=geom,unresolved_warnings=unresolved,coverage={k:v for k,v in cov.items() if k!='rows'},support_self={k:v for k,v in selfcheck.items() if k!='rows'},layer_coverage={k:v for k,v in layer.items() if k!='rows'},feature_preservation={k:v for k,v in feat.items() if k!='rows'})
    save(folder/'evaluation.json',result);print('ACTUAL',folder.name,'coverage',cov['centerline_underbody_percent'],'maxarc',cov['max_uncovered_arc_mm'],'layer',layer['percent'],'selfbad',selfcheck['self_unsupported_regions'],'features',feat['missing_samples'],flush=True);return result
