"""Native Rhino 7 3dm files, with one closed mesh and layer per printed part.

The assembly and separated plate are distinct documents. No duplicated
reference surfaces, casting substrate, or imagined slicer supports are added.
The source core's analytic CAD remains available as STEP beside the meshes.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import rhino3dm as rhino
import trimesh

from . import y_channel_casting as shared
from .bottomless_channel_fixture import casting_height_qc, part_height_limit


def role_color(name):
    if 'core' in name: return 'lumen_core', (196, 81, 66, 255)
    if 'box' in name: return 'bottomless_frame', (132, 169, 187, 255)
    if 'retainer' in name: return 'dovetail_retainer', (212, 166, 76, 255)
    return 'cross_key', (63, 128, 101, 255)


def write_model(path, parts, purpose, cfg):
    model=rhino.File3dm()
    model.ApplicationName='vascular_printing'
    model.ApplicationDetails='P1S 0.4 mm; closed mesh components; rhino3dm '+rhino.__version__
    model.Settings.ModelUnitSystem=rhino.UnitSystem.Millimeters
    model.Settings.ModelAbsoluteTolerance=.001
    model.Strings['purpose']=purpose
    model.Strings['geometry_type']='Closed mesh objects; analytic lumen core in separate STEP'
    model.Strings['height_limits']=json.dumps(cfg['height_limits'],sort_keys=True)
    model.Strings['casting_reference']='External substrate top Z0; substrate is not a printed component'
    model.Strings['temporary_print_supports']='Not included; generate them when manually slicing the editable Bambu project'
    for name, source in parts.items():
        if not source.is_volume or len(source.split())!=1:
            raise ValueError('RHINO_SOURCE_NOT_SINGLE_CLOSED_SOLID: '+name)
        role,color=role_color(name)
        layer=rhino.Layer();layer.Name=name;layer.Color=color
        layer_index=model.Layers.Add(layer)
        attributes=rhino.ObjectAttributes();attributes.Name=name;attributes.LayerIndex=layer_index
        attributes.ColorSource=rhino.ObjectColorSource.ColorFromLayer
        attributes.SetUserString('units','mm');attributes.SetUserString('role',role)
        attributes.SetUserString('purpose',purpose)
        if cfg.get('part_materials'):
            attributes.SetUserString('print_material',cfg['part_materials']['core' if 'core' in name else 'fixture'])
        mesh=rhino.Mesh();mesh.Vertices.UseDoublePrecisionVertices=True
        for xyz in source.vertices: mesh.Vertices.Add(*map(float,xyz))
        for face in source.faces: mesh.Faces.AddFace(*map(int,face))
        mesh.Normals.ComputeNormals();mesh.Compact()
        if not mesh.IsValid or not mesh.IsClosed:
            raise ValueError('RHINO_INVALID_MESH: '+name+' '+str(mesh.IsValidWithLog))
        model.Objects.AddMesh(mesh,attributes)
    if not model.Write(str(path),7): raise ValueError('RHINO_WRITE_FAILED: '+str(path))


def read_and_check(path, expected, cfg, purpose):
    """Read the actual 3dm again and compare every mesh with its source STL."""
    model=rhino.File3dm.Read(str(path))
    if model is None: raise ValueError('RHINO_READ_FAILED: '+str(path))
    if model.Settings.ModelUnitSystem!=rhino.UnitSystem.Millimeters:
        raise ValueError('RHINO_MM_UNITS_REQUIRED')
    names=[o.Attributes.Name for o in model.Objects]
    if len(names)!=len(expected) or set(names)!=set(expected):
        raise ValueError('RHINO_PART_SET_MISMATCH')
    if model.Strings['purpose']!=purpose: raise ValueError('RHINO_PURPOSE_MISMATCH')
    recovered={};rows=[]
    for obj in model.Objects:
        name=obj.Attributes.Name;geometry=obj.Geometry;source=expected[name]
        if not isinstance(geometry,rhino.Mesh) or not geometry.IsValid or not geometry.IsClosed:
            raise ValueError('RHINO_ROUNDTRIP_NOT_CLOSED_MESH: '+name)
        vertices=np.array([[v.X,v.Y,v.Z] for v in geometry.Vertices])
        faces=np.array([geometry.Faces[i] for i in range(geometry.Faces.Count)])
        if not np.array_equal(faces[:,2],faces[:,3]): raise ValueError('RHINO_EXPECTED_TRIANGLES')
        if vertices.shape!=source.vertices.shape or not np.array_equal(faces[:,:3],source.faces):
            raise ValueError('RHINO_TOPOLOGY_CHANGED: '+name)
        error=float(abs(vertices-source.vertices).max())
        if error>2e-5: raise ValueError('RHINO_COORDINATES_CHANGED: '+name)
        actual=trimesh.Trimesh(vertices,faces[:,:3],process=False)
        if not actual.is_volume or len(actual.split())!=1: raise ValueError('RHINO_INVALID_VOLUME: '+name)
        if model.Layers[obj.Attributes.LayerIndex].Name!=name:
            raise ValueError('RHINO_COMPONENT_LAYER_MISMATCH: '+name)
        recovered[name]=actual
        if cfg.get('part_materials') and obj.Attributes.GetUserString('print_material')!=cfg['part_materials']['core' if 'core' in name else 'fixture']:
            raise ValueError('RHINO_MATERIAL_ASSIGNMENT_MISMATCH')
        rows.append(dict(name=name,vertices=len(vertices),triangles=len(faces),closed=True,
            layer=name,bounds_mm=actual.bounds,vertex_max_error_mm=error))
    if purpose=='casting_assembly':
        height=casting_height_qc(cfg,recovered)
    else:
        # This is the separated CAD plate, before adding slicer raft/support.
        height=dict(passed=all(abs(m.bounds[0,2])<2e-5 and m.bounds[1,2]<=
                    part_height_limit(cfg,n)+2e-5
                    for n,m in recovered.items()),reference='Print bed Z0; raft not included')
    if not height['passed']: raise ValueError('RHINO_HEIGHT_LIMIT_FAILED')
    return dict(passed=True,units='mm',archive_version=model.ArchiveVersion,
        purpose=purpose,part_count=len(rows),parts=rows,height_qc=height,
        rhino_desktop_open_tested=False,reader='rhino3dm '+rhino.__version__)


def export_models(cfg,parts,summary,out):
    if summary['status']!='GEOMETRY_PASSED': raise ValueError('RHINO_REQUIRES_PASSED_GEOMETRY')
    out=Path(out);kind=summary['variant'];records=[]
    # Use saved production STL meshes, making the exported geometry auditable
    # against exactly the same surfaces that enter the native slicer.
    for purpose,folder in [('casting_assembly','assembly'),('print_layout','print_parts')]:
        sources={name:out/folder/(name+'.stl') for name in parts}
        meshes={name:trimesh.load_mesh(path) for name,path in sources.items()}
        path=out/(kind+'_'+cfg.get('rhino_export_tag','H15')+'_'+purpose+'.3dm')
        write_model(path,meshes,purpose,cfg)
        record=read_and_check(path,meshes,cfg,purpose)
        record.update(file=path.name,sha256=shared.digest(path),
            sources_sha256={str(p.relative_to(out)):shared.digest(p) for p in sources.values()})
        records.append(record)
    result=dict(passed=True,geometry_type='Closed meshes, not NURBS/Brep',
        exported_files=records,fixture_height_qc=summary['manufacturing_height_qc'],
        actual_dish_inner_diameter_mm=None,physical_dish_fit_validated=False)
    shared.save(out/'rhino_qc.json',result)
    print('RHINO_3DM_VERIFIED',kind,[r['file'] for r in records],flush=True)
    return result


def verify_exports(out):
    """Reject stale/corrupt Rhino files instead of silently packaging them."""
    out=Path(out);result=json.loads((out/'rhino_qc.json').read_text())
    cfg=json.loads((out/'effective_config.json').read_text())
    if not result['passed'] or len(result['exported_files'])!=2:
        raise ValueError('RHINO_QC_PASS_REQUIRED')
    for record in result['exported_files']:
        path=out/record['file']
        if shared.digest(path)!=record['sha256']: raise ValueError('RHINO_FILE_HASH_MISMATCH')
        meshes={}
        for relative,expected_hash in record['sources_sha256'].items():
            source=out/relative
            if shared.digest(source)!=expected_hash: raise ValueError('RHINO_SOURCE_STALE')
            meshes[source.stem]=trimesh.load_mesh(source)
        read_and_check(path,meshes,cfg,record['purpose'])
    return result
