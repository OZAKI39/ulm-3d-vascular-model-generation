"""Measured roof-gap candidates and conservative bridge-limit sensitivity."""
import copy,json
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import save,table
from .support_v2_layers import analyse,coverage,bridge_evidence
from .support_v2_primitives import roof_pad
from .support_v2_planner import manager
from .support_v2_mesh import support_stl
from .support_v2_inputs import public

def unbraced_trunk_metrics(primitives):
    """A break neck is not a lateral brace: retain the full root-trunk height."""
    groups={}
    for p in primitives:
        if 'root_id' in p and p['kind'] in ['VERTICAL_PILLAR','BREAKAWAY_SEGMENT_NECK','BUTTRESS']:
            groups.setdefault(p['root_id'],[]).append(p)
    rows=[]
    for identity,members in sorted(groups.items()):
        height=max(p['bbox_mm'][1][2] for p in members)-min(p['bbox_mm'][0][2] for p in members)
        width=min(p['width_mm'] for p in members if 'width_mm' in p)
        rows.append(dict(root_id=identity,unbraced_height_mm=height,minimum_section_width_mm=width,conservative_slenderness=height/width))
    return dict(rows=rows,max_unbraced_height_mm=max([r['unbraced_height_mm'] for r in rows],default=0.),max_slenderness=max([r['conservative_slenderness'] for r in rows],default=0.),method='Continuous root trunk height / narrowest designed neck; segmentation alone supplies no lateral bracing; not a physical strength prediction')

def run(record,profile):
    folder=Path(record['folder']);out=folder/'sensitivity';out.mkdir(exist_ok=True)
    v=trimesh.load_mesh(folder/'vascular_with_ports.stl');b=trimesh.load_mesh(folder/'box.stl');vman=manager(v);bman=manager(b)
    fields=[analyse(v,profile,'vascular',b),analyse(b,profile,'box')];roofs=json.loads((folder/'roofs.json').read_text());gap_rows=[];bridge_rows=[]
    for gap in [.1,.15,.2,.25]:
        variants=[];shapes=[]
        for r in roofs:
            q=copy.deepcopy(r);dz=gap-r['gap_mm'] if r['geometry_class']!='box' else 0.
            top=np.r_[r['xyz'][:2],r['top_z_mm']-dz]
            shape,footprint=roof_pad(top,r['length_mm'],r['width_mm'],r['thickness_mm'],0,r['angle_rad'])
            measured=0. if vman.in_collision_single(shape) else float(vman.min_distance_single(shape))
            q.update(gap_mm=gap if r['geometry_class']!='box' else 0.,top_z_mm=float(top[2]),measured_min_gap_mm=measured,_footprint=footprint)
            q['valid']=not vman.in_collision_single(shape) and (r['geometry_class']=='box' or (measured>=gap-.002 and not bman.in_collision_single(shape)))
            if q['valid']:variants.append(q);shapes.append(shape)
        row=dict(gap_mm=gap,valid_roofs=len(variants),invalid_roofs=len(roofs)-len(variants),status='PHYSICAL_GAP_UNVALIDATED',scope='CONTACT_GEOMETRY_ONLY_NOT_COMPLETE_PRODUCTION_SCAFFOLD')
        cov=coverage(fields[0],variants,profile);row.update(vascular_polygon_coverage_percent=cov['percent'],uncovered_regions=cov['uncovered_region_count']);gap_rows.append(row)
        if shapes:
            raw=trimesh.boolean.union(shapes,engine='manifold')
            try:support_stl(raw).export(out/f'roof_gap_{gap:.2f}.stl');row['contact_mesh_closed']=True
            except ValueError as exc:raw.export(out/f'roof_gap_{gap:.2f}.stl');row.update(contact_mesh_closed=False,export_error=str(exc))
        save(out/f'roof_gap_{gap:.2f}.json',public(variants))
    for limit in [2.,4.,6.,8.]:
        for f in fields:
            bridge_count=0;required_area=0.
            for row,poly in zip(f['rows'],f['polygons']):
                previous=f['layer_polygons'][max(0,row['layer_index']-1)]
                evidence=bridge_evidence(poly,previous,limit,profile['line_width']);bridge_count+=int(evidence['bridgeable'])
                if not evidence['bridgeable']:required_area+=poly.area
            bridge_rows.append(dict(geometry_class=f['rows'][0]['geometry_class'] if f['rows'] else 'box',bridge_limit_mm=limit,bridgeable_regions=bridge_count,required_area_sum_mm2=required_area,physical_validation=False,status='BRIDGE_LIMIT_UNCALIBRATED' if profile['bridge_calibration_status']!='PHYSICAL_COUPON_CALIBRATED' else 'SENSITIVITY_NOT_CALIBRATION'))
    table(out/'gap_sensitivity.csv',gap_rows);table(out/'bridge_sensitivity.csv',bridge_rows)
    return dict(gaps=gap_rows,bridge_limits=bridge_rows)

def add_metrics(record):
    folder=Path(record['folder']);graph=json.loads((folder/'support_graph_v2.json').read_text());roofs=json.loads((folder/'roofs.json').read_text());support=trimesh.load_mesh(folder/'custom_support.stl')
    record['near_vessel_roof_area_mm2']=sum(r['length_mm']*r['width_mm'] for r in roofs if r['geometry_class']!='box')
    record['support_component_count']=len(support.split())
    angles=[np.degrees(np.arctan2(np.linalg.norm(np.array(e['end'])[:2]-np.array(e['start'])[:2]),abs(e['end'][2]-e['start'][2]))) for e in graph['edges']]
    record['mechanical_proxy']['max_branch_inclination_degrees']=max(angles,default=0.)
    # Exact segment direction from the generated primitive, not its AABB.
    import pickle
    with (folder/'scaffold_local_cache.pkl').open('rb') as f:_,_,pieces,_=pickle.load(f)
    for p,m in zip(graph['primitives'],pieces):
        if p['kind'] in ['SHORT_BRIDGE','LONGITUDINAL_RAIL']:
            _,vectors=np.linalg.eigh(np.cov(m.vertices[:,:2].T));p['bridge_direction_xy']=vectors[:,-1].tolist()
            p['span_mm']=max(0.,p['span_mm']);p['bridge_state']='VALIDATED_BRIDGE_CANDIDATE' if p['span_mm']>0 else 'SELF_SUPPORTED_RAIL'
    graph['bridge_count']=sum(p.get('span_mm',0)>0 for p in graph['primitives'] if p['kind'] in ['SHORT_BRIDGE','LONGITUDINAL_RAIL'])
    record['bridge_count']=graph['bridge_count'];record['mechanical_proxy']['bridge_count']=graph['bridge_count']
    trunks=unbraced_trunk_metrics(graph['primitives']);table(folder/'unbraced_trunk_metrics.csv',trunks['rows'])
    record['mechanical_proxy']['unbraced_trunk_metrics']=trunks
    record['mechanical_proxy']['max_slenderness']=max(record['mechanical_proxy']['max_slenderness'],trunks['max_slenderness'])
    record['mechanical_proxy']['max_pillar_height_mm']=max(record['mechanical_proxy']['max_pillar_height_mm'],trunks['max_unbraced_height_mm'])
    component_rows=[];centers=np.array([p.bounds.mean(0) for p in pieces])
    for identity,part in enumerate(support.split()):
        possible=np.flatnonzero(np.all(centers>=part.bounds[0]-.001,axis=1)&np.all(centers<=part.bounds[1]+.001,axis=1))
        inside=part.contains(centers[possible]);assigned=[graph['primitives'][j] for j in possible[inside]]
        pillars=[p for p in assigned if p['kind']=='VERTICAL_PILLAR'];branches=[p for p in assigned if p['kind']=='SELF_SUPPORTING_BRANCH'];spans=[p.get('span_mm',0.) for p in assigned if p['kind'] in ['SHORT_BRIDGE','LONGITUDINAL_RAIL']]
        slants=[np.degrees(np.arctan2(np.linalg.norm(np.array(p['end'])[:2]-np.array(p['start'])[:2]),abs(p['end'][2]-p['start'][2]))) for p in branches if 'start' in p]
        pads=[p for p in assigned if p['kind']=='ROOF_PAD']
        component_rows.append(dict(component_id=identity,volume_mm3=float(part.volume),max_pillar_height_mm=max([p['bbox_mm'][1][2]-p['bbox_mm'][0][2] for p in pillars],default=0),min_width_mm=min([p['width_mm'] for p in assigned if 'width_mm' in p],default=None),max_slenderness=max([(p['bbox_mm'][1][2]-p['bbox_mm'][0][2])/p['width_mm'] for p in pillars],default=0),max_branch_inclination_degrees=max(slants,default=0),bridge_count=len(spans),bridge_max_span_mm=max(spans,default=0),supported_load_center=np.mean([np.mean(p['bbox_mm'],axis=0) for p in pads],axis=0).tolist() if pads else None,primitive_membership_method='Representative interior primitive centers inside closed union component; boundary-center ambiguity retained',physical_strength_validated=False))
        unbraced=unbraced_trunk_metrics(assigned);component_rows[-1]['max_unbraced_trunk_height_mm']=unbraced['max_unbraced_height_mm'];component_rows[-1]['conservative_trunk_slenderness']=unbraced['max_slenderness']
        component_rows[-1]['rail_and_bridge_member_count']=len(spans);component_rows[-1]['bridge_count']=sum(span>0 for span in spans)
    table(folder/'mechanical_components.csv',component_rows)
    table(folder/'support_primitives.csv',graph['primitives'])
    save(folder/'support_graph_v2.json',graph);save(folder/'candidate.json',record)
    return record
