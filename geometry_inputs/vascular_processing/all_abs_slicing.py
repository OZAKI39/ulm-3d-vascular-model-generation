"""Actual isolated ABS Bambu slices; no printer delivery or profile mutation."""
from __future__ import annotations
import json,zipfile
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .support_inputs import save,digest,table,numeric,resolve_profile
from .support_slicer_adapter import engine_profiles,discover_engine,slice_local,production_geometry
from .casting_fixture_profiles import bed_mesh
from .bambu_casting_mold import parse_paths,sample_segments
from .casting_fixture_v2_geometry import port_transform
from .casting_fixture_geometry import shape_mesh


def profiles(inputs,out):
    values,d=resolve_profile(inputs,out/'profile_discovery');actual=engine_profiles(d,'Bambu',{'enable_support':'0'});reference=engine_profiles(d,'Bambu');engine=discover_engine('Bambu',d,out/'engine')
    if actual['filament'].get('filament_type')!=['ABS']:raise ValueError('ACTUAL_ABS_PROFILE_REQUIRED')
    save(out/'actual_ABS_profiles.json',dict(values=values,support_off=actual,vascular_reference=reference,source=d));return values,actual,reference,engine


def slice_abs(part,engine,actual,out,support_off=True):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);bed,m=bed_mesh(part);path=out/'print_input.stl';bed.export(path);save(out/'print_transform.json',dict(transform_4x4=m,scale_factor=1.))
    signature=dict(source=digest(path),profiles=actual,binary_hash=engine.get('binary_sha256'));cache=out/'slice_inputs.json';result=out/'slice_result.json'
    reuse=cache.exists() and result.exists() and json.loads(cache.read_text())==signature
    rec=json.loads(result.read_text()) if reuse else {}
    if not rec.get('archive') or not Path(rec['archive']).exists():rec={}
    save(cache,signature)
    if not rec:rec=slice_local(engine,actual,path,out)
    if rec.get('archive'):
        with zipfile.ZipFile(rec['archive']) as z:
            cfg=json.loads(z.read('Metadata/project_settings.config'));crc=z.testzip() is None
        paths=parse_paths(rec['archive']);native=int((paths['kind']==1).sum());off=str(cfg.get('enable_support')) in ('0','False','false');geom=production_geometry(rec['archive'],path)
        rec.update(archive_crc_ok=crc,native_support_extrusion_segments=native,auto_support_off=off,archive_geometry=geom,actual_profile={k:cfg.get(k) for k in ['printer_model','printer_settings_id','nozzle_diameter','filament_type','filament_settings_id','print_settings_id','line_width','layer_height','initial_layer_print_height','enable_support']})
        ok=bool(rec.get('slicer_success') and crc and geom['passed'] and cfg.get('filament_type')==['ABS'] and not rec['warning_audit']['critical'])
        rec['slice_passed']=ok and (not support_off or (off and native==0))
    else:rec['slice_passed']=False
    rec.update(support_off_required=support_off,physical_print_validated=False,no_printer_delivery=True);save(result,rec);return rec


def feature_preservation(l,g,mold_record,ret_record,packing,out):
    rows=[];detail=[]
    for kind,record,folder in [('MOLD',mold_record,out/'outer_mold'),('RETAINERS',ret_record,out/'retainers')]:
        if not record.get('slice_passed'):return dict(passed=False,reason='ACTUAL_ABS_SLICE_MISSING_OR_FAILED')
        paths=parse_paths(record['archive']);segments=paths['segments'][paths['kind']==0];tree=cKDTree(sample_segments(segments,.08));m=np.asarray(json.loads((folder/'print_transform.json').read_text())['transform_4x4']);settings=record['actual_profile'];layer=numeric(settings['layer_height']);first=numeric(settings['initial_layer_print_height']);t=g['blade_thickness_mm'];queries=[]
        if kind=='MOLD':
            for p in l['ports']:
                R=p['radius_mm']+p['cradle_clearance_mm'];inner=R+.55;w=p['wall_distance_from_trim_mm'];k=p['port_id'];retz=p['radius_mm']+p['anti_float_gap_mm'];top=retz+3*g['line_width_mm']+.7+1.1
                levels=np.arange(first,l['high'][2]+layer,layer);zc=levels-layer/2-p['trim_origin'][2];zf=l['low'][2]-p['trim_origin'][2]
                def vertical(name,xy,zmax,zmin=zf+.03):
                    q=[]
                    for x,y in xy:
                        mask=(zc>zmin)&(zc<zmax-.04);q.extend(np.c_[np.full(mask.sum(),x),np.full(mask.sum(),y),levels[mask]-p['trim_origin'][2]])
                    if q:
                        mm=port_transform(p);qs=np.asarray(q)@mm[:3,:3].T+mm[:3,3];queries.append((k,name,qs))
                ys=np.linspace(min(g['local_parts'][k]['cradle_ys'])-.2,w-.35,20)
                for sign in [-1,1]:vertical('wall_fin_'+str(sign),[(sign*(inner+t/2),y) for y in ys],top)
                vertical('solvent_web',[(0,y) for y in np.linspace(2.3,w-.4,15)],-R-.08)
                rec=next(r for r in g['access_records'] if r['port_id']==k);contactz=min(np.asarray(rec['contact_witnesses'])[:,2])-p['trim_origin'][2]
                vertical('conformal_access_land',[(0,c) for c in [2.,2.2,2.4]],contactz,zmin=-R-.12+.03)
                # Shallow bottom guide: evaluate both 0.2 mm relief layers,
                # outside fixture blades and entirely on the trim side.
                vertical('transferred_trim_guide',[(x,t/2) for x in np.linspace(-3.8,-2.8,5)],zf+.2,zmin=zf+.02)
        else:
            for pk,p in zip(packing,l['ports']):
                k=p['port_id'];part=shape_mesh(g['local_parts'][k]['retainer']);R=np.asarray(pk['local_rotation']);shift=np.asarray(pk['translation']);part.apply_transform(np.block([[R,shift[:,None]],[np.zeros((1,3)),np.ones((1,1))]]))
                # Stem centre-line across its full height in assembly, sampled
                # every 0.25 mm, repeated on every occupied flat-print layer.
                w=p['wall_distance_from_trim_mm'];retz=p['radius_mm']+p['anti_float_gap_mm'];ztop=l['high'][2]-p['trim_origin'][2]+.8
                line=np.array([[0,g['local_parts'][k]['stem_center_y_mm'],z] for z in np.arange(retz+.7+.2,ztop-.1,.25)]);pts=line@R.T+shift;levels=np.arange(first,part.bounds[1,2],layer);q=[]
                for z in levels:
                    if z-layer/2>part.bounds[0,2]+.04 and z-layer/2<part.bounds[1,2]-.04:
                        a=pts.copy();a[:,2]=z;q.extend(a)
                queries.append((k,'top_exposed_retainer_stem',np.asarray(q)))
        for k,name,qs in queries:
            qs=qs@m[:3,:3].T+m[:3,3];dist=tree.query(qs)[0];tol=g['line_width_mm']*.80;good=dist<=tol
            row=dict(part=kind,port_id=k,feature=name,query_count=len(qs),missing_query_count=int((~good).sum()),max_toolpath_distance_mm=float(dist.max()),tolerance_mm=tol,passed=bool(good.all()));rows.append(row)
            if not good.all():detail.append(dict(**row,missing_points=qs[~good][:12].tolist(),missing_distances=dist[~good][:12].tolist()))
    r=dict(passed=bool(rows) and all(r['passed'] for r in rows),rows=rows,failures=detail,method='Mandatory thin-feature centre-line samples on every occupied layer versus actual positive model extrusion; finite-resolution geometric preservation, no assertion of actual zero-gap solvent contact.',physical_contact_preservation_calibrated=False);save(out/'feature_preservation.json',r);table(out/'feature_preservation.csv',rows);return r
