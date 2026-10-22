"""Local Bambu slicing and audits of real model/support extrusion toolpaths."""
import json
import re
from pathlib import Path
import zipfile
import numpy as np
from scipy.spatial import cKDTree
from . import bambu_slice_adapter as adapter
from .sacrificial_fixture_review import save_json
from .sacrificial_fixture import write_csv


def discover(cfg,out):return adapter.discover_current(cfg,out)


def extruder_coordinate_frame(archive,allow_missing_settings=False):
    """Restore model/plate coordinates using the actual archived machine offset.

    Bambu GCode::point_to_gcode subtracts extruder_offset.  Its inverse adds
    it; estimating this offset by fitting paths to the STL would hide errors.
    This adapter intentionally rejects ambiguous multi-extruder configurations.
    """
    name='Metadata/project_settings.config'
    settings=json.loads(archive.read(name)) if name in archive.namelist() else {}
    entries=settings.get('extruder_offset')
    if entries is None:
        if not allow_missing_settings:raise ValueError('GCODE_EXTRUDER_OFFSET_UNRESOLVED')
        return np.zeros(3),dict(source='EXPLICIT_ZERO_OFFSET_TEST_FIXTURE',offset_added_mm=[0,0,0])
    if isinstance(entries,str):entries=[entries]
    try:
        offsets=np.array([[float(v) for v in entry.split('x')] for entry in entries])
    except (TypeError,ValueError,AttributeError) as exc:
        raise ValueError('GCODE_EXTRUDER_OFFSET_UNRESOLVED') from exc
    if offsets.ndim!=2 or offsets.shape[1]!=2 or not len(offsets) or not np.isfinite(offsets).all():
        raise ValueError('GCODE_EXTRUDER_OFFSET_UNRESOLVED')
    if not np.allclose(offsets,offsets[0],rtol=0,atol=1e-12):
        raise ValueError('MULTI_EXTRUDER_COORDINATE_MAPPING_UNSUPPORTED')
    offset=np.r_[offsets[0],0.]
    return offset,dict(source=name,archived_extruder_offset=entries,offset_added_mm=offset.tolist(),
        output_frame='MODEL_BUILD_PLATE_COORDINATES',raw_gcode_modified=False,
        reference='https://github.com/bambulab/BambuStudio/blob/master/src/libslic3r/GCode.cpp',
        rule='GCode::point_to_gcode subtracts extruder_offset; restore it before comparing paths to the placed model.')


def parse_paths(archive_path,allow_missing_settings=False):
    segments=[];kinds=[];extrusions=[];tags=set();arcs=0
    with zipfile.ZipFile(archive_path) as z:
        offset,coordinate_frame=extruder_coordinate_frame(z,allow_missing_settings)
        for name in z.namelist():
            if not name.endswith('.gcode'):continue
            position=np.zeros(3);previous_e=0.;relative_e=True;absolute=True;feature='';plane='G17'
            for raw in z.read(name).decode(errors='replace').splitlines():
                if raw.startswith('; FEATURE:'):feature=raw.split(':',1)[1].strip();tags.add(feature)
                line=raw.split(';',1)[0].strip()
                if line=='M83':relative_e=True
                elif line=='M82':relative_e=False
                elif line=='G90':absolute=True
                elif line=='G91':absolute=False
                elif line in ('G17','G18','G19'):plane=line
                match=re.match(r'^(G0|G1|G2|G3|G92)\s+(.*)',line)
                if not match:continue
                values={k:float(v) for k,v in re.findall(r'([XYZEIJR])([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)',match.group(2))}
                code=match[1]
                if code=='G92':
                    previous_e=values.get('E',previous_e)
                    for k,axis in enumerate('XYZ'):
                        if axis in values:position[k]=values[axis]
                    continue
                new=position.copy()
                for k,axis in enumerate('XYZ'):
                    if axis in values:new[k]=values[axis]+(0 if absolute else position[k])
                extrusion=values.get('E',0.) if relative_e else values.get('E',previous_e)-previous_e
                if 'E' in values:previous_e=previous_e+values['E'] if relative_e else values['E']
                if extrusion>0 and feature:
                    if code in ('G2','G3'):
                        if plane!='G17':raise ValueError('UNSUPPORTED_GCODE_ARC_PLANE')
                        points=adapter.arc_polyline(position,new,values,code=='G2');arcs+=1
                    else:points=np.array([position,new])
                    pieces=np.stack([points[:-1],points[1:]],axis=1);length=np.linalg.norm(np.diff(pieces,axis=1)[:,0],axis=1)
                    if length.sum()>1e-9:
                        segments.extend(pieces);extrusions.extend(extrusion*length/length.sum())
                        kind=1 if 'support' in feature.lower() else 2 if any(x in feature.lower() for x in ('skirt','brim','prime','wipe','custom')) else 0
                        kinds.extend([kind]*len(pieces))
                position=new
    return dict(segments=(np.array(segments).reshape(-1,2,3)+offset).astype(np.float32),kind=np.array(kinds,dtype=np.uint8),
        extrusion_filament_mm=np.array(extrusions),arc_command_count=arcs,feature_tags=sorted(tags),coordinate_frame=coordinate_frame)


def sample_segments(segments,step):
    chunks=[]
    for batch in np.array_split(segments,max(1,int(np.ceil(len(segments)/50000)))):
        if not len(batch):continue
        lengths=np.linalg.norm(batch[:,1]-batch[:,0],axis=1)
        counts=np.maximum(1,np.ceil(lengths/step).astype(int))
        # Group equal segment sample counts to avoid a Python loop per point.
        for count in np.unique(counts):
            selected=batch[counts==count];u=np.linspace(0,1,count+1)
            chunks.append((selected[:,0,None,:]*(1-u[None,:,None])+selected[:,1,None,:]*u[None,:,None]).reshape(-1,3))
    return np.vstack(chunks) if chunks else np.empty((0,3))


def feature_preservation(paths,mesh,partition,box,transform,cfg):
    model=paths['segments'][paths['kind']==0]
    points=sample_segments(model,cfg['feature_qc']['toolpath_sample_step_mm']);rows=[]
    if not len(points):return dict(passed=False,status='SLICER_DROPPED_VASCULAR_FEATURE',rows=[],reason='No model extrusion paths')
    tree=cKDTree(points);centers=mesh.triangles_center[partition['face_ids']];labels=partition['labels']
    rotation=np.array(transform)[:3,:3];translation=np.array(transform)[:3,3];step=cfg['feature_qc']['spatial_bin_mm'];tolerance=cfg['feature_qc']['maximum_surface_distance_mm']
    for label in sorted(set(labels)):
        p=centers[labels==label];distance=tree.query(p@rotation.T+translation)[0]
        bins=np.floor((p-p.min(0))/step).astype(int);unique,inverse=np.unique(bins,axis=0,return_inverse=True)
        nearest=np.full(len(unique),np.inf);np.minimum.at(nearest,inverse,distance)
        covered=nearest<=tolerance;fraction=float(covered.mean())
        rows.append(dict(feature=str(label),type='VASCULAR_BRANCH_OR_PORT',tested_bins=len(unique),covered_bins=int(covered.sum()),
            coverage_fraction=fraction,maximum_bin_nearest_path_distance_mm=float(nearest.max()),
            passed=fraction>=cfg['feature_qc']['minimum_bin_coverage_fraction'],
            method='Actual non-support positive-extrusion G-code within 3 mm spatial bins; finite-resolution feature-presence test'))
    lo,hi=box['inner']
    for wall in ('-X','+X','-Y','+Y','BOTTOM'):
        if wall=='BOTTOM':
            queries=np.array([[x,y,lo[2]] for x in np.linspace(lo[0]+3,hi[0]-3,5) for y in np.linspace(lo[1]+3,hi[1]-3,5)])
        else:
            axis='XY'.index(wall[-1]);side=0 if wall[0]=='-' else 1;other=1-axis;queries=[]
            for a in np.linspace(lo[other]+3,hi[other]-3,5):
                for z in np.linspace(lo[2]+3,hi[2]-.4,5):
                    p=np.zeros(3);p[axis]=(lo if side==0 else hi)[axis];p[other]=a;p[2]=z;queries.append(p)
            queries=np.array(queries)
        distances=tree.query(queries@rotation.T+translation)[0]
        rows.append(dict(feature=wall,type='BOX_WALL',tested_bins=len(queries),covered_bins=int((distances<=tolerance).sum()),
            coverage_fraction=float((distances<=tolerance).mean()),maximum_bin_nearest_path_distance_mm=float(distances.max()),
            passed=bool(np.all(distances<=tolerance)),method='5 x 5 inner-surface samples checked against actual model extrusion'))
    passed=all(r['passed'] for r in rows)
    return dict(passed=passed,status='SAMPLED_GCODE_FEATURES_PRESENT' if passed else 'SLICER_DROPPED_VASCULAR_FEATURE',rows=rows,
        spatial_resolution_mm=step,maximum_surface_distance_mm=tolerance,full_bead_volume_reconstruction_performed=False)


def support_quantity(paths,transform,box,discovery):
    inverse=np.linalg.inv(transform);middle=paths['segments'].mean(axis=1).astype(float)@inverse[:3,:3].T+inverse[:3,3]
    inside=np.all((middle>=box['inner'][0])&(middle<=box['inner'][1]),axis=1)
    support=paths['kind']==1;internal=support&inside
    filament=discovery['flattened']['filament']
    def value(key):
        entry=filament.get(key);entry=entry[0] if isinstance(entry,list) and entry else entry
        try:return float(entry)
        except (TypeError,ValueError):return None
    diameter=value('filament_diameter');density=value('filament_density')
    length=float(paths['extrusion_filament_mm'][internal].sum())
    volume=length*np.pi*(diameter/2)**2 if diameter else None
    return dict(actual_support_segment_count=int(support.sum()),internal_support_segment_count=int(internal.sum()),
        total_support_extrusion_filament_mm=float(paths['extrusion_filament_mm'][support].sum()),
        internal_support_extrusion_filament_mm=length,estimated_internal_support_volume_mm3=volume,
        estimated_internal_support_mass_g=volume/1000*density if volume is not None and density else None,
        filament_diameter_mm=diameter,filament_density_g_cm3=density,
        internal_support_quantity_method='Actual positive E extrusion assigned by segment midpoint inside the casting cavity; uses active filament diameter/density.',
        support_geometry_available=False,support_geometry_status='BAMBU_SUPPORT_GEOMETRY_UNAVAILABLE')


def slice_one(discovery,cfg,row,out,logger,mesh,partition,box):
    result=adapter.slice_top_five(discovery,cfg,[row],out,logger)
    result['support_geometry_status']='BAMBU_SUPPORT_GEOMETRY_UNAVAILABLE'
    for record in result.get('results',[]):
        if record.get('archive'):audit_slice(record,cfg,row,mesh,partition,box)
    save_json(Path(out)/'slice_results.json',result)
    return result


def native_warning_lines(record):
    messages=[]
    def visit(value):
        if isinstance(value,dict):
            for key,item in value.items():
                if 'warning' in key.lower() and isinstance(item,str):messages.extend(item.splitlines())
                elif isinstance(item,(dict,list)):visit(item)
        elif isinstance(value,list):
            for item in value:visit(item)
    visit(record.get('native_cli_result',{}))
    messages.extend((record.get('stdout','')+'\n'+record.get('stderr','')).splitlines())
    return sorted(set(line for line in messages if re.search(r'unsupported|floating|thin.feature',line,re.I)))


def audit_slice(record,cfg,row,mesh,partition,box):
    """Recheck an unchanged actual 3MF; no slicing or printer action."""
    folder=Path(record['archive']).parent
    paths=parse_paths(record['archive'])
    np.savez_compressed(folder/'actual_extrusion_paths.npz',segments=paths['segments'],kind=paths['kind'],extrusion_filament_mm=paths['extrusion_filament_mm'])
    record['actual_paths_file']=str(folder/'actual_extrusion_paths.npz')
    record['gcode_coordinate_frame']=paths['coordinate_frame']
    record['feature_preservation']=feature_preservation(paths,mesh,partition,box,row['transform_4x4'],cfg)
    with zipfile.ZipFile(record['archive']) as archive:
        archived_settings=json.loads(archive.read('Metadata/project_settings.config'))
    record['support_quantity']=support_quantity(paths,row['transform_4x4'],box,dict(flattened=dict(filament=archived_settings)))
    record['gcode_arc_count']=paths['arc_command_count'];record['gcode_feature_tags']=paths['feature_tags']
    # Keep genuine slicer warnings and native results; do not infer physical
    # printability merely from a successful file export.
    record['unsupported_floating_thin_warning_lines']=native_warning_lines(record)
    with zipfile.ZipFile(record['archive']) as z:
        if 'Metadata/plate_1.png' in z.namelist():
            path=folder/'native_bambu_plate_preview.png';path.write_bytes(z.read('Metadata/plate_1.png'));record['native_preview']=str(path)
        for name in z.namelist():
            if name.endswith('.gcode'):
                path=folder/Path(name).name;path.write_bytes(z.read(name));record['gcode_file']=str(path)
    save_json(folder/'casting_slice_audit.json',record)
    return record


def write_slice_comparison(records,path):
    rows=[]
    for r in records:
        q=r.get('support_quantity',{})
        rows.append(dict(candidate_id=r['orientation_id'],support_mode=r['support_mode'],status=r['status'],
            print_time_seconds=r.get('estimated_print_time_seconds'),total_filament_g=r.get('filament_used_g'),
            internal_support_estimated_g=q.get('estimated_internal_support_mass_g'),
            internal_support_estimated_mm3=q.get('estimated_internal_support_volume_mm3'),
            total_support_filament_mm=q.get('total_support_extrusion_filament_mm'),
            sampled_features_preserved=r.get('feature_preservation',{}).get('passed'),
            model_geometry_preserved=r.get('geometry_qc',{}).get('rotation_scale_position_preserved'),
            actual_support_type=r.get('actual_project_settings',{}).get('support_type'),
            warning_lines=r.get('unsupported_floating_thin_warning_lines',[])))
    # Geometry-only runs have no slice records; still emit a valid empty table.
    fields=['candidate_id','support_mode','status','print_time_seconds','total_filament_g',
        'internal_support_estimated_g','internal_support_estimated_mm3','total_support_filament_mm',
        'sampled_features_preserved','model_geometry_preserved','actual_support_type','warning_lines']
    write_csv(path,rows,fields=fields)
