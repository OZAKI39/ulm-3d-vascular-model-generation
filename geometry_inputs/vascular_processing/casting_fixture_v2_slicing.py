"""Actual local Bambu evidence for wall fins and flat separate retainers."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .bambu_casting_mold import parse_paths,sample_segments
from .casting_fixture_profiles import slice_part
from .casting_fixture_v2_geometry import port_transform
from .support_inputs import save,table
from .support_slicer_adapter import production_geometry


def blade_queries(l,g,layer=.1,first=.2):
    """Interior fin centre-lines on every occupied layer; 0.25 mm along web.

    These are mandatory blade/rib features, deliberately away from walls.
    Cradle points above the analytic lower arc are not material queries.
    """
    result=[];t=g['lines']*g['line_width']
    for p in l['ports']:
        R=p['radius_mm']+p['cradle_clearance_mm'];inner=R+.55;w=p['wall_distance_from_trim_mm'];z0=-p['trim_origin'][2];retz=p['radius_mm']+p['anti_float_gap_mm'];rett=3*g['line_width'];top=l['high'][2]-p['trim_origin'][2]
        zs=np.arange(first, l['high'][2]+.01,layer);zcentres=zs-layer/2-p['trim_origin'][2]
        features={};cradles=g['local'][p['port_id']]['cradle_ys'];yfront=min(cradles)-t/2
        def line_points(name,xy,zmax):
            pts=[]
            for x,y in xy:
                ok=(zcentres<zmax(x,y)-.035)&(zcentres>z0+.025)
                pts.extend(np.c_[np.full(ok.sum(),x),np.full(ok.sum(),y),zs[ok]-p['trim_origin'][2]])
            features[name]=np.asarray(pts).reshape(-1,3)
        ys=np.linspace(yfront+.12,w-.4,max(2,int((w-yfront)/.25)))
        for sign in [-1,1]:line_points('radial_web_'+str(sign),[(sign*(inner+t/2),y) for y in ys],lambda x,y:retz+rett+1.8)
        for j,y in enumerate(cradles):
            xs=np.linspace(-inner,inner,max(3,int(2*inner/.25)))
            line_points('cradle_'+str(j),[(x,y) for x in xs],lambda x,y:-np.sqrt(max(0,R*R-x*x)))
        ro=g['local'][p['port_id']]['rib_outer_mm']
        line_points('alignment_rib',[(x,t/2) for x in np.linspace(inner+t+.2,ro-.05,5)],lambda x,y:top-.5)
        line_points('rib_wall_connection',[(ro-t/2,y) for y in np.linspace(t,w-.4,18)],lambda x,y:top-.5)
        m=port_transform(p)
        for name,ps in features.items():result.append((p['port_id'],name,ps@m[:3,:3].T+m[:3,3]))
    return result


def validate_blades(record,l,g,out):
    out=Path(out)
    if not record.get('slicer_success'):return dict(passed=False,reason='NO_ACTUAL_SLICE')
    paths=parse_paths(record['archive']);segments=paths['segments'][paths['kind']==0]
    m=np.array(json.loads((out/'print_transform.json').read_text())['transform_4x4'])
    # A common STEP/STL placement error must fail before any path comparison.
    archive_geometry=production_geometry(record['archive'],out/'print_input.stl')
    points=sample_segments(segments,.10);tree=cKDTree(points);rows=[];worst=[]
    tolerance=g['line_width']*.75 # 0.165 mm: includes half-line + sample chord.
    settings=record['actual_project_settings']
    for port,name,qs in blade_queries(l,g,float(settings['layer_height']),float(settings['initial_layer_print_height'])):
        qs=qs@m[:3,:3].T+m[:3,3];d=tree.query(qs)[0];good=d<=tolerance
        layers=np.unique(np.round(qs[:,2],4));missing_layers=sum(not good[np.isclose(qs[:,2],z)].any() for z in layers)
        rows.append(dict(port_id=port,feature=name,query_count=len(qs),all_points_preserved=bool(good.all()),missing_query_count=int((~good).sum()),max_path_distance_mm=float(d.max()),layer_count=len(layers),missing_layers=int(missing_layers),tolerance_mm=tolerance))
        if not good.all():worst.extend([dict(port_id=port,feature=name,point=x.tolist(),distance_mm=float(dd)) for x,dd in zip(qs[~good][:20],d[~good][:20])])
    r=dict(passed=all(x['all_points_preserved'] for x in rows) and archive_geometry['passed'],rows=rows,missing_samples=worst,archive_geometry=archive_geometry,model_extrusion_segment_count=len(segments),auto_support_extrusion_segments=int((paths['kind']==1).sum()),method='Every occupied layer: mandatory blade/rib centre-line samples spaced <=0.25 mm, nearest actual positive-extrusion model G-code; 0.10 mm path sampling. Finite-resolution preservation, not physical print certification.')
    save(out/'blade_preservation.json',r);table(out/'blade_preservation.csv',rows);return r
