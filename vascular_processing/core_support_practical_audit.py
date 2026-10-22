"""Practical s1-11 acceptance from actual extrusion, separate from gap quality.

Profile-scale macro tests are geometric screening rules for a first trial,
not calibrated ABS strength predictions. Historical strict audits are retained.
"""
from pathlib import Path
from collections import Counter
import json,pickle
import numpy as np
import shapely
from shapely.geometry import Point,LineString,GeometryCollection
from shapely.ops import nearest_points,unary_union
from .core_support_slicing import coordinate_tolerance
from .core_support_geometry import bed_roots
from .support_inputs import save,table


def structural_samples(branches,footprints,profile,old_rows=(),archive=''):
    """Use the highest actual support top; a close roof is still a roof.

    The maximum useful gap stays at nominal gap + one actual layer. We do
    not credit distant walls or numerical exceptions without local evidence.
    """
    levels=np.array(sorted(footprints));nominal=profile['roof_gap_mm'];maximum=nominal+profile['layer_height']
    for fp in footprints.values():shapely.prepare(fp)
    old={(r['branch_id'],r['sample_index']):r for r in old_rows};rows=[];maxarc=0.
    for b in branches:
        local=[]
        for i,q in enumerate(b['lower']):
            q=np.asarray(q);tol=coordinate_tolerance(q,profile);pt=Point(q[:2]);best=None;witness=None
            usable=levels[(levels<=q[2]+tol)&(levels>=q[2]-maximum-tol)]
            for z in usable[::-1]:
                fp=footprints[float(z)];dist=0. if fp.covers(pt) else float(fp.distance(pt))
                if best is None or dist<best[0]:best=(dist,float(z))
                if dist<=tol:witness=(dist,float(z));break
            previous=old.get((b['branch_id'],i),{})
            old_reason='PASS' if previous.get('supported',False) else 'SUPPORT_TOO_CLOSE' if previous.get('too_close_support_layer_tops') else 'OLD_STRICT_GAP_OR_FOOTPRINT_FAILURE'
            supported=witness is not None;new='ACCEPTABLE';gap=None;far_witness=None
            if supported:
                best=witness;gap=float(q[2]-best[1]);new='SUPPORT_TOO_CLOSE' if gap<nominal-tol else 'ACCEPTABLE'
            else:
                new='TRUE_MISSING_SUPPORT'
                # A farther wall is evidence of TOO_FAR, never support credit.
                farther=levels[levels<q[2]-maximum-tol]
                for z in farther[::-1]:
                    if footprints[float(z)].distance(pt)<=tol:
                        new='SUPPORT_TOO_FAR';gap=float(q[2]-z);far_witness=dict(layer_top_mm=float(z),gap_mm=gap,lateral_distance_mm=float(footprints[float(z)].distance(pt)));break
                if best is not None and gap is None:gap=float(q[2]-best[1])
            path=None
            if best is not None:
                near=pt if best[0]==0 else nearest_points(pt,footprints[best[1]])[1]
                path=dict(archive=str(archive),role='ACTUAL_CUSTOM_SUPPORT_EXTRUSION_FOOTPRINT',layer_top_mm=best[1],nearest_footprint_xy=list(near.coords[0]),lateral_distance_mm=best[0])
            row=dict(sample_id=f"{b['branch_id']}:{i}",branch=b['branch_id'],sample_index=i,port_or_vessel='PORT' if b['is_port'] else 'VESSEL',is_port=b['is_port'],arc_mm=float(b['arc'][i]),xyz=q.tolist(),old_reason=old_reason,old_pass=bool(previous.get('supported',False)),new_class=new,gap_mm=gap,nearest_support_path=path,structurally_supported=supported,true_missing_support=not supported,numerical_exception=False,coordinate_tolerance_mm=tol,action_required='LOCAL_MISSING_SUPPORT_REPAIR' if not supported else 'REMOVAL_AND_SURFACE_REVIEW' if new=='SUPPORT_TOO_CLOSE' else 'NONE')
            local.append(row)
            row['farther_support_witness']=far_witness
        # An isolated boundary point may use two independently rounded bead
        # endpoints, but both adjacent samples AND their actual connected
        # footprint must remain continuous. Never exempt a multi-point gap.
        for i in range(1,len(local)-1):
            r=local[i]
            if r['structurally_supported'] or not local[i-1]['structurally_supported'] or not local[i+1]['structurally_supported']:continue
            path=r['nearest_support_path']
            if path is None or path['lateral_distance_mm']>2*r['coordinate_tolerance_mm']:continue
            a=local[i-1]['nearest_support_path'];c=local[i+1]['nearest_support_path']
            if a['layer_top_mm']!=c['layer_top_mm'] or path['layer_top_mm']!=a['layer_top_mm']:continue
            line=LineString([local[i-1]['xyz'][:2],r['xyz'][:2],local[i+1]['xyz'][:2]])
            fp=footprints[path['layer_top_mm']]
            connected=any(part.buffer(2*r['coordinate_tolerance_mm']).covers(line) for part in shapely.get_parts(fp))
            if connected:r.update(new_class='NUMERICAL_BOUNDARY_ONLY',structurally_supported=True,true_missing_support=False,numerical_exception=True,action_required='VISUAL_REVIEW_NUMERICAL_EXCEPTION',exception_evidence=dict(neighbors=[local[i-1]['sample_id'],local[i+1]['sample_id']],same_actual_connected_footprint=True,bound_mm=2*r['coordinate_tolerance_mm']))
        bad=[r['true_missing_support'] for r in local]
        start=None
        for i,flag in enumerate([*bad,False]):
            if flag and start is None:start=i
            if not flag and start is not None:
                left=max(0,start-1);right=min(len(local)-1,i)
                maxarc=max(maxarc,float(b['arc'][right]-b['arc'][left]));start=None
        rows.extend(local)
    counts=Counter(r['new_class'] for r in rows);missing=sum(r['true_missing_support'] for r in rows)
    return dict(passed=missing==0,true_missing_support_count=missing,max_true_unsupported_arc_mm=maxarc,structural_support_percent=100*sum(r['structurally_supported'] for r in rows)/len(rows),gap_quality_counts=dict(TOO_CLOSE=counts['SUPPORT_TOO_CLOSE'],ACCEPTABLE=counts['ACCEPTABLE']+counts['NUMERICAL_BOUNDARY_ONLY'],TOO_FAR=counts['SUPPORT_TOO_FAR']),class_counts=dict(counts),numerical_exceptions=counts['NUMERICAL_BOUNDARY_ONLY'],sample_count=len(rows),maximum_usable_vertical_gap_mm=maximum,rows=rows,method='Highest actual custom-support extrusion below each frozen lower-envelope target; close gaps are separate removal/surface risks; no distant-wall credit')


def macro_printability(footprints,profile,mesh=None,strict_rows=()):
    """Material connectivity, floating starts, and profile-scale overhang depth.

    A deficit is micro only while attached to a grounded extrusion component,
    shallower than one bead width, and not propagating into a growing deep
    overhang on successive layers. Perimeter polygon area alone is no gate.
    """
    h=profile['layer_height'];w=profile['support_line_width'];reach=min(h,w/2);rows=[];floating=[];prev=GeometryCollection();previous_bad=[];previous_z=None
    strict_by_z={}
    for r in strict_rows:strict_by_z.setdefault(round(r['z_mm'],5),[]).append(r)
    for z,fp in sorted(footprints.items()):
        if previous_z is not None and z-previous_z>h+1e-5:prev=GeometryCollection();previous_bad=[]
        first=z<=profile['first_layer_height']+1e-5;tol=coordinate_tolerance(fp.bounds,profile)
        foundation=prev.buffer(reach);shapely.prepare(foundation)
        if not first:
            for j,part in enumerate(shapely.get_parts(fp)):
                if part.is_empty:continue
                grounded=foundation.intersects(part)
                if not grounded and part.area>part.length*tol+np.pi*tol**2:
                    item=dict(region_id=f'FLOAT_{z:.5f}_{j}',z_mm=z,status='MACRO_FAIL',reason='FLOATING_EXTRUSION_ISLAND_OR_FEATURE_STARTS_IN_AIR',area_mm2=float(part.area),wkt=part.wkt,maximum_depth_mm=None,propagation_layers=1)
                    floating.append(item);rows.append(item)
        current_bad=[]
        for j,item in enumerate(strict_by_z.get(round(z,5),[])):
            p=shapely.from_wkt(item['wkt']);xy=shapely.get_coordinates(p);depth=float(shapely.distance(shapely.points(xy),prev).max()) if not prev.is_empty else float('inf')
            old=[r for r in previous_bad if p.intersects(r['polygon'].buffer(reach))]
            chain=1+max((r['chain'] for r in old),default=0)
            growing=any(p.area>r['area']+w*h and depth>r['depth']+tol for r in old)
            deep=depth>w+tol;propagating=chain>=3 and growing and depth>reach+w/2
            island=any(p.intersects(shapely.from_wkt(r['wkt'])) for r in floating if abs(r['z_mm']-z)<1e-5)
            macro=island or deep or propagating
            row=dict(region_id=f'DEF_{z:.5f}_{j}',z_mm=z,status='MACRO_FAIL' if macro else 'MICRO_DEFICIT',reason='DETACHED_MATERIAL' if island else 'DEEP_UNSUPPORTED_SUPPORT_SPAN' if deep else 'PROPAGATING_STRUCTURAL_DEFICIT' if propagating else 'ATTACHED_SHALLOW_NONPROPAGATING_PERIMETER_DEFICIT',area_mm2=float(p.area),maximum_depth_mm=depth if np.isfinite(depth) else None,propagation_layers=chain,growing=growing,wkt=p.wkt)
            rows.append(row);current_bad.append(dict(polygon=p,chain=chain,depth=depth,area=p.area))
        previous_bad=current_bad;prev=fp;previous_z=z
    roots=bed_roots(mesh) if mesh is not None else []
    for root in roots:
        if not root['bed_rooted']:rows.append(dict(region_id=f"MESH_{root['component_id']}",status='MACRO_FAIL',reason='DETACHED_SUPPORT_GEOMETRY',z_mm=root['minimum_z_mm'],area_mm2=None,wkt=None,maximum_depth_mm=None,propagation_layers=None))
    if not rows:rows=[dict(region_id='ALL',status='PASS',reason='NO_DEFICIT',z_mm=None,area_mm2=0.,wkt=None,maximum_depth_mm=0.,propagation_layers=0)]
    fails=sum(r['status']=='MACRO_FAIL' for r in rows)
    return dict(passed=fails==0,macro_fail_count=fails,floating_support_islands=len(floating)+sum(not r['bed_rooted'] for r in roots),micro_deficit_count=sum(r['status']=='MICRO_DEFICIT' for r in rows),bed_roots=roots,rows=rows,criteria=dict(maximum_shallow_deficit_depth_mm=w,previous_layer_reach_mm=reach,propagating_layers=3,growth_area_mm2=w*h),physical_stability_validated=False)


def vascular_layer_starts(vascular_footprints,support_footprints,profile):
    """Locate new actual vascular extrusion islands without a prior anchor.

    Existing vascular layer material may anchor its own next layer, but is
    never credited as custom support in the separate lower-envelope audit.
    Dense branch sampling independently catches long partially rooted spans.
    """
    h=profile['layer_height'];reach=min(h,profile['line_width']/2);maximum=profile['roof_gap_mm']+h;levels=np.array(sorted(support_footprints));rows=[];prev=GeometryCollection();last=None
    for z,fp in sorted(vascular_footprints.items()):
        if last is not None and z-last>h+1e-5:prev=GeometryCollection()
        foundation=prev.buffer(reach);shapely.prepare(foundation)
        for i,part in enumerate(shapely.get_parts(fp)):
            if foundation.intersects(part):continue
            tol=coordinate_tolerance(part.bounds,profile)
            if part.area<=part.length*tol+np.pi*tol**2:continue
            eligible=levels[(levels>=z-h-maximum-tol)&(levels<=z+tol)]
            roofs=unary_union([support_footprints[float(q)] for q in eligible]);overlap_area=float(part.intersection(roofs.buffer(tol)).area)
            rows.append(dict(region_id=f'V_START_{z:.5f}_{i}',z_mm=float(z),area_mm2=float(part.area),actual_support_overlap_mm2=overlap_area,passed=overlap_area>0,eligible_support_tops=eligible.tolist(),wkt=part.wkt))
        prev=fp;last=z
    return dict(passed=all(r['passed'] for r in rows),floating_vascular_region_count=sum(not r['passed'] for r in rows),rows=rows,method='Actual vascular layer connected-island starts, anchored by preceding vascular material or actual nearby custom support; dense lower-envelope test remains independent')


def warning_classification(warnings,structural,macro,floating_vascular=0):
    rows=[]
    for warning in warnings:
        text=str(warning)
        if 'floating' in text.lower():
            critical=bool(structural['true_missing_support_count'] or macro['macro_fail_count'] or floating_vascular)
            reason='LOCALIZED_REAL_MISSING_OR_MACRO_FAILURE' if critical else 'INTENTIONAL_GAP_AND_ATTACHED_MICRO_DEFICITS_WITH_ACTUAL_EXTRUSION_EVIDENCE'
        else:critical=True;reason='UNCLASSIFIED_NATIVE_CRITICAL_WARNING'
        rows.append(dict(warning=text,status='CRITICAL_UNRESOLVED_WARNING' if critical else 'NONCRITICAL_WARNING',reason=reason))
    return dict(critical_unresolved=sum(r['status']=='CRITICAL_UNRESOLVED_WARNING' for r in rows),rows=rows)


def audit_candidate(folder,out=None):
    folder=Path(folder);out=Path(out or folder/'practical');out.mkdir(parents=True,exist_ok=True)
    g=pickle.load((folder/'geometry.pkl').open('rb'));fp=pickle.load((folder/'actual_footprints.pkl').open('rb'));e=json.loads((folder/'evaluation.json').read_text());old=json.loads((folder/'coverage.json').read_text());sc=json.loads((folder/'support_self.json').read_text())
    profile=dict(g['profile']);precision=folder/'gcode_precision.json'
    if precision.exists():profile['gcode_xy_half_cell_mm']=json.loads(precision.read_text())['half_coordinate_cell_mm']
    structural=structural_samples(g['branches'],fp['support'],profile,old['rows'],e['slice']['archive'])
    macro=macro_printability(fp['support'],profile,g['support'],sc['rows'])
    # Dense lower-envelope sampling covers every source branch/port; missing
    # native vascular features remain a separate hard gate, not gap quality.
    floating=vascular_layer_starts(fp['vascular'],fp['support'],profile)
    warning=warning_classification(e['slice']['warning_audit']['critical'],structural,macro,floating['floating_vascular_region_count'])
    result=dict(candidate=folder.name,structural={k:v for k,v in structural.items() if k!='rows'},macro={k:v for k,v in macro.items() if k!='rows'},warnings=warning,floating_vascular_regions_unresolved=floating['floating_vascular_region_count'],vascular_features_preserved=e['feature_preservation']['passed'],native_slice_pass=e['slice']['slicer_success'],auto_support_off=e['auto_support_off'])
    result['support_vascular_intersection_mm3']=json.loads((folder/'geometry.json').read_text())['intersection_mm3']
    save(out/'floating_vascular_regions.json',floating)
    save(out/'final_true_missing_support.json',structural);save(out/'support_self_macro_audit.json',macro);save(out/'practical_audit.json',result)
    table(out/'support_failure_reclassification.csv',[r for r in structural['rows'] if not r['old_pass'] or r['true_missing_support'] or r['numerical_exception']])
    table(out/'support_self_macro_audit.csv',macro['rows'])
    print('PRACTICAL_AUDIT',folder.name,'missing',structural['true_missing_support_count'],'close',structural['gap_quality_counts']['TOO_CLOSE'],'macro',macro['macro_fail_count'],flush=True)
    return result,structural,macro
