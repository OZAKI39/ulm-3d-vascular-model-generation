"""GEOS slice differences, branch-aware grains and bounded closure controls."""
from __future__ import annotations
import numpy as np
import shapely
from shapely.geometry import GeometryCollection
from shapely.ops import unary_union
from scipy.spatial import cKDTree


def area_tolerance(width,coordinate_magnitude):
    # Two independently rounded XY endpoints: four float32 coordinate ULPs
    # times one extrusion width. Never relax the historical 1e-5 mm² ceiling.
    ulp=float(np.spacing(np.float32(max(1.,coordinate_magnitude))))
    return min(1e-5,4*ulp*width)


def difference(required,actual):
    return shapely.difference(shapely.make_valid(required),shapely.make_valid(actual))


def quantized_width(width,line_width):
    return float(np.ceil((width-1e-10)/line_width)*line_width)


def cluster(rows,branches,profile,self_support=False):
    points=np.vstack([b['xyz'] for b in branches]);tree=cKDTree(points)
    names=np.concatenate([np.full(len(b['xyz']),b['branch_id'],dtype=object) for b in branches])
    arcs=np.concatenate([b['arc'] for b in branches]);tangent=np.vstack([b['tangent'] for b in branches])
    pieces=[];w=profile['line_width'];h=profile['layer_height']
    for index,row in enumerate(rows):
        poly=shapely.from_wkt(row['polygon_wkt'] if self_support else row['uncovered_wkt'])
        for part in shapely.get_parts(shapely.make_valid(poly)):
            if part.area<=area_tolerance(w,max(abs(x) for x in part.bounds)):continue
            xy=np.array(part.representative_point().coords[0]);z=row['z_mm'];_,j=tree.query([*xy,z])
            pieces.append(dict(source_index=index,component='SUPPORT_SELF' if self_support else row['component'],
                               z=z,layer=int(round((z-profile['first_layer_height'])/h)),poly=part,branch=str(names[j]),
                               arc=float(arcs[j]),tangent=tangent[j],point=[*xy,z]))
    parent=list(range(len(pieces)))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    order=sorted(range(len(pieces)),key=lambda i:pieces[i]['layer'])
    for pos,i in enumerate(order):
        a=pieces[i]
        for j in order[pos+1:]:
            b=pieces[j]
            if b['layer']-a['layer']>2:break
            if a['component']!=b['component'] or a['branch']!=b['branch']:continue
            if abs(a['arc']-b['arc'])>max(4*w,1.) or abs(np.dot(a['tangent'],b['tangent']))<np.cos(np.radians(30)):continue
            if a['poly'].distance(b['poly'])<=w+1e-9:parent[root(j)]=root(i)
    groups={}
    for i,p in enumerate(pieces):groups.setdefault(root(i),[]).append(p)
    grains=[]
    for group in groups.values():
        group.sort(key=lambda p:(p['z'],p['arc']));xy=unary_union([p['poly'] for p in group]);zs=[p['z'] for p in group]
        grains.append(dict(grain_id=f"{'S' if self_support else 'V'}G{len(grains):04}",component=group[0]['component'],
            layer_start=min(p['layer'] for p in group),layer_end=max(p['layer'] for p in group),
            bbox=[list(xy.bounds[:2])+[min(zs)],list(xy.bounds[2:])+[max(zs)]],area_mm2=sum(p['poly'].area for p in group),
            branch=group[0]['branch'],arc_interval_mm=[min(p['arc'] for p in group),max(p['arc'] for p in group)],
            centroid_path=[p['point'] for p in group],source_rows=sorted({p['source_index'] for p in group}),footprint_wkt=xy.wkt,
            required_repair_type='TOP_DOWN_THIN_WEB' if self_support else 'EXTEND_EXISTING_ROOF_RIDGE',
            tangent=np.mean([p['tangent'] for p in group],axis=0).tolist()))
    return grains


def top_down_layers(footprint,start_z,actual,profile,root_z=0.):
    """Propagate the same local grain, stopping only at printable material/root."""
    z=float(start_z);h=profile['layer_height'];rows=[];levels=np.array(sorted(actual))
    while z>root_z+h/2:
        previous=unary_union([actual[float(q)] for q in levels[abs(levels-(z-h))<h/2]])
        missing=difference(footprint,previous.buffer(min(profile['line_width']/2,h)))
        rows.append(dict(z_mm=z,missing_area_mm2=float(missing.area)))
        if missing.area<=area_tolerance(profile['line_width'],max(abs(x) for x in footprint.bounds)):
            return dict(base_z_mm=max(root_z,z-h),connected_to='ACTUAL_PREVIOUS_MATERIAL',layers=rows)
        z=round(z-h,6)
    return dict(base_z_mm=root_z,connected_to='LEGAL_ROOT_REQUIRED',layers=rows)


def iteration_stop(areas,all_covered=False,max_repairs=5):
    if all_covered:return 'ALL_COVERAGE_PASSED'
    if len(areas)>=3:
        gains=[(a-b)/max(a,1e-15) for a,b in zip(areas[-3:-1],areas[-2:])]
        if all(g<.01 for g in gains):return 'COVERAGE_REPAIR_STAGNATED'
    if len(areas)-1>=max_repairs:return 'REPAIR_ITERATION_LIMIT'
    return None
