"""Print-pose/profile revision of a frozen two-port core. No slice/toolpath API.

Only a rigid transform and editable Bambu project settings change. The original
assembly, all fixtures and the other two channel designs are read-only inputs.
"""
from __future__ import annotations
import copy
import hashlib
import json
from pathlib import Path
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import trimesh
import yaml
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from .all_abs_print_plate import NS, read_objects
from .channel_manual_delivery import validate_project, PROJECT_MEMBERS
from .straight_channel_casting import load_config, design
from .y_channel_casting import ROOT, save, digest

NAME = 'Straight_lumen_core'


def geometry_audit(before, after, old_pose, new_pose, cfg):
    """Separate build height from casting height; preserve every source triangle."""
    back = after.copy()
    back.apply_transform(old_pose @ np.linalg.inv(new_pose))
    error = float(abs(back.vertices-before.vertices).max())
    if not np.array_equal(before.faces, after.faces) or error > 1e-7:
        raise ValueError('CORE_GEOMETRY_CHANGED')
    if not after.is_volume or not np.isclose(before.volume, after.volume, rtol=1e-8):
        raise ValueError('CORE_VOLUME_CHANGED')
    d = design(cfg)
    axis = new_pose[:3, :3] @ np.array([1., 0., 0.])
    # The unchanged straight cylindrical segment is parallel to build Z, so
    # its radial vectors have no Z component. This is its complete Z extent.
    endpoints = trimesh.transform_points(np.array([
        [d['bend_x'], 0., d['center_z']],
        [d['outlet_bend_x'], 0., d['center_z']]]), new_pose)
    if not np.allclose(axis, [0., 0., 1.], atol=1e-9):
        raise ValueError('STRAIGHT_AXIS_MUST_BE_VERTICAL')
    clearance = float(endpoints[:, 2].min())
    if clearance < 3. or abs(after.bounds[0, 2]) > 1e-8:
        raise ValueError('UNSAFE_CORE_BED_PLACEMENT')
    if np.any(after.bounds[0, :2] < 20) or np.any(after.bounds[1, :2] > 236):
        raise ValueError('CORE_OUTSIDE_RESERVED_BED')
    if after.bounds[1, 2] > 256:
        raise ValueError('CORE_EXCEEDS_P1S_BUILD_HEIGHT')
    return dict(passed=True,rigid_transform_only=True,maximum_inverse_vertex_error_mm=error,
        triangle_count=len(after.faces),source_volume_mm3=before.volume,
        print_volume_mm3=after.volume,diameter_mm=cfg['diameter_mm'],
        main_straight_length_mm=d['straight_centerline_length_mm'],
        main_straight_axis_in_print=axis,main_straight_bed_clearance_mm=clearance,
        main_straight_z_range_mm=sorted(endpoints[:, 2].tolist()),
        print_bounds_mm=after.bounds,print_height_mm=float(after.extents[2]),
        casting_height_limits_unchanged=True,
        casting_to_print_transform=new_pose,print_to_casting_transform=np.linalg.inv(new_pose),
        main_segment_contacts_bed=False,
        support_geometry_generated=False,toolpaths_generated=False,
        physical_surface_finish_validated=False,physical_support_removal_validated=False)


def write_refined_project(source, destination, mesh, overrides, *, plate_name=None):
    with zipfile.ZipFile(source) as z:
        if set(z.namelist()) != PROJECT_MEMBERS or z.testzip() is not None:
            raise ValueError('EXPECTED_UNSLICED_SOURCE_PROJECT')
        members = {n:z.read(n) for n in z.namelist()}
    root = ET.fromstring(members['3D/3dmodel.model'])
    meshes = root.findall(f'{{{NS}}}resources/{{{NS}}}object/{{{NS}}}mesh')
    if len(meshes) != 1:
        raise ValueError('EXPECTED_EXACTLY_ONE_CORE_MESH')
    vertices = meshes[0].findall(f'{{{NS}}}vertices/{{{NS}}}vertex')
    faces = np.array([[int(t.get(k)) for k in ('v1','v2','v3')]
                      for t in meshes[0].findall(f'{{{NS}}}triangles/{{{NS}}}triangle')])
    if len(vertices) != len(mesh.vertices) or not np.array_equal(faces, mesh.faces):
        raise ValueError('SOURCE_MESH_STORAGE_UNEXPECTED')
    for row, xyz in zip(vertices, mesh.vertices):
        for k,v in zip('xyz',xyz):row.set(k,format(float(v),'.12g'))
    ET.register_namespace('',NS)
    members['3D/3dmodel.model'] = ET.tostring(root,encoding='utf-8',xml_declaration=True)
    settings = json.loads(members['Metadata/project_settings.config'])
    before = copy.deepcopy(settings)
    settings.update(overrides)
    settings.pop('adaptive_layer_height',None)  # Not a Bambu 2.7 serialized option.
    settings['inherits_group'] = ['0.08mm High Quality @BBL X1C','','']
    # PresetCollection::load_external_preset replaces every unlisted setting
    # with the installed system preset, EVEN when the project contains a value.
    # This is a complete frozen snapshot: pin every present option in each of
    # the process/filament/printer groups (the loader intersects with its keys).
    # GUI save/read-back verifies that values survive, not just that it opens.
    pinned=';'.join(sorted(k for k in settings if k!='different_settings_to_system'))
    settings['different_settings_to_system']=[pinned,pinned,pinned]
    # Object overrides take precedence over global settings. Update the old
    # normal/snug/two-raft-layer settings here too, otherwise the GUI keeps them.
    config = ET.fromstring(members['Metadata/model_settings.config'])
    objects = config.findall('object')
    if len(objects) != 1:raise ValueError('EXPECTED_ONE_CORE_OBJECT')
    for key in ['enable_support','support_type','support_style','raft_layers']:
        node = next((m for m in objects[0].findall('metadata') if m.get('key')==key),None)
        if node is None:node=ET.SubElement(objects[0],'metadata',key=key)
        node.set('value',str(settings[key]))
    for node in config.findall('.//plate/metadata'):
        if node.get('key')=='plater_name':
            node.set('value',plate_name or 'Straight 2Top - ABS Fine Vertical 0.08 - UNSLICED')
    members['Metadata/model_settings.config'] = ET.tostring(config,encoding='utf-8',xml_declaration=True)
    members['Metadata/project_settings.config'] = json.dumps(settings).encode()
    destination.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(destination,'w',zipfile.ZIP_DEFLATED) as z:
        for name,data in members.items():z.writestr(name,data)
    return {k:dict(before=before.get(k),after=v) for k,v in settings.items() if before.get(k)!=v}


def verify_gui_roundtrip(project, saved, quality):
    """Validate GUI-effective settings via an unsliced Save As, not CLI --info."""
    original,_,_,_=read_objects(project)
    recovered,settings,plates,crc=read_objects(saved)
    if not crc or plates!=1 or len(recovered)!=1:
        raise ValueError('GUI_ROUNDTRIP_OBJECT_MISMATCH')
    meta=recovered[0]['metadata']
    # Object support overrides are effective even if the GUI normalizes the
    # global preset's support type/style. Check both scopes, never silently
    # substitute a changed global value for an explicit object setting.
    effective=dict(settings)
    for k in ['enable_support','support_type','support_style','raft_layers']:
        if meta.get(k)!=quality['settings'][k]:raise ValueError('GUI_OBJECT_OVERRIDE_RESET: '+k)
        effective[k]=meta[k]
    differences={k:dict(expected=v,actual=effective.get(k)) for k,v in quality['settings'].items()
                 if k!='print_settings_id' and effective.get(k)!=v}
    if differences:raise ValueError('GUI_SETTINGS_RESET: '+json.dumps(differences))
    error=float(abs(original[0]['mesh'].bounds-recovered[0]['mesh'].bounds).max())
    if error>2e-5:raise ValueError('GUI_PRINT_POSE_CHANGED')
    with zipfile.ZipFile(saved) as z:
        if any(n.endswith(('.gcode','.bgcode')) for n in z.namelist()):
            raise ValueError('GUI_ROUNDTRIP_HAS_TOOLPATHS')
    return dict(passed=True,application='Bambu Studio 02.07.01.57 Linux GUI',
        method='Open project, Save As without slicing, read and compare effective settings',
        project_sha256=digest(project),gui_saved_sha256=digest(saved),
        settings_verified={k:effective[k] for k in quality['settings'] if k!='print_settings_id'},
        global_support_settings={k:settings.get(k) for k in ['enable_support','support_type','support_style','raft_layers']},
        effective_object_support_settings={k:meta[k] for k in ['enable_support','support_type','support_style','raft_layers']},
        object_overrides_verified=True,maximum_bounds_difference_mm=error,
        slicing_performed=False,physical_print_tested=False)


def preview(before, after, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    fig, axes = plt.subplots(1,2,figsize=(11,6))
    for ax,m,title in zip(axes,[before,after],['Previous: horizontal / 0.20 mm','Revised: upright main tube / 0.08 mm']):
        v=m.vertices.copy();v[:,0]-=m.bounds[:,0].mean()
        # Front projection, with actual mesh contours. No imagined support mesh.
        ax.add_collection(PolyCollection(v[m.faces][:,:,[0,2]],facecolors='#c45142',edgecolors='none'))
        ax.axhline(0,color='#46546a',linewidth=2,label='Build plate')
        ax.set(xlim=(-20,20),ylim=(-2,36),xlabel='Print X (mm)',ylabel='Print Z (mm)',title=title)
        ax.set_aspect('equal');ax.legend(loc='upper left');ax.grid(alpha=.15)
    fig.suptitle('Same 1.60 mm core; casting geometry unchanged\nSupport settings only - supports are NOT generated')
    fig.tight_layout();fig.savefig(path,dpi=180);plt.close(fig)


def package(out):
    """Package this core revision only. Never bundle the earlier core by accident."""
    out=Path(out);summary=json.loads((out/'quality_summary.json').read_text())
    project=out/summary['project']['file']
    if digest(project)!=summary['project']['sha256']:raise ValueError('STALE_QUALITY_SUMMARY')
    files=[project,out/'Straight_lumen_core_fine_vertical.stl',out/'quality_summary.json',
           out/'effective_quality_config.yaml',out/'orientation_comparison.png',out/'使用说明.md']
    if (out/'bambu_gui_verified.png').exists():
        files += [out/'bambu_gui_verified.png',out/'gui_verification.json']
    hashes={p.name:digest(p) for p in files}
    dest=ROOT/'deliverables/Straight_2Top_P1S04_ABS_fine_vertical_unsliced.zip'
    dest.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED) as z:
        for p in files:z.write(p,p.name)
        z.writestr('SHA256.json',json.dumps(hashes,indent=2,ensure_ascii=False))
    with zipfile.ZipFile(dest) as z:
        if z.testzip() is not None:raise ValueError('QUALITY_ZIP_CRC_FAILURE')
        for n,h in hashes.items():
            if hashlib.sha256(z.read(n)).hexdigest()!=h:raise ValueError('QUALITY_ZIP_HASH_FAILURE')
    save(out/'delivery_receipt.json',dict(zip=str(dest.relative_to(ROOT)),sha256=digest(dest),
         files_sha256=hashes,sliced=False,generated_support_meshes=0))
    return dest


def run(config_path):
    config_path=Path(config_path);quality=yaml.safe_load(config_path.read_text())
    cfg=load_config(ROOT/quality['source_config']);source_out=cfg['output']
    if cfg.get('channel_topology')!='two_vertical_ports':raise ValueError('ONLY_TWO_TOP_PORT_CORE')
    source=source_out/'Straight_2Top_P1S_04_core_ABS_project.3mf'
    assembly=source_out/'assembly/Straight_lumen_core.stl'
    for p,key in [(source,'source_project_sha256'),(assembly,'source_assembly_sha256')]:
        if digest(p)!=quality[key]:raise ValueError('FROZEN_SOURCE_HASH_MISMATCH: '+str(p))
    protected={p:digest(p) for p in [*source_out.glob('*.3mf'),*source_out.glob('*.3dm'),
        *source_out.glob('assembly/*'),source_out/'effective_config.json',source_out/'print_layout.json'] if p.is_file()}
    objects,settings,plates,crc=read_objects(source)
    if not crc or plates!=1 or len(objects)!=1 or objects[0]['metadata']['name']!=NAME:
        raise ValueError('EXPECTED_ORIGINAL_SINGLE_CORE')
    before=objects[0]['mesh'];after=before.copy()
    old_record=next(r for r in json.loads((source_out/'print_layout.json').read_text())['records'] if r['name']==NAME)
    old_pose=np.array(old_record['casting_to_plate_transform']);original=before.copy()
    original.apply_transform(np.linalg.inv(old_pose))
    saved=trimesh.load_mesh(assembly)
    if max(cKDTree(saved.vertices).query(original.vertices)[0].max(),
           cKDTree(original.vertices).query(saved.vertices)[0].max())>2e-5:
        raise ValueError('PROJECT_DOES_NOT_MATCH_FROZEN_ASSEMBLY')
    new_pose=np.eye(4);new_pose[:3,:3]=Rotation.from_euler('y',quality['rotation_y_degrees'],degrees=True).as_matrix()
    after=original.copy();after.apply_transform(new_pose)
    translation=np.r_[quality['plate_center_mm'],0.]-np.r_[after.bounds.mean(0)[:2],after.bounds[0,2]]
    new_pose[:3,3]=translation;after.apply_translation(translation)
    audit=geometry_audit(before,after,old_pose,new_pose,cfg)
    out=(ROOT/quality['output']).resolve()
    if out==source_out or source_out in out.parents:raise ValueError('USE_SEPARATE_QUALITY_OUTPUT')
    out.mkdir(parents=True,exist_ok=True)
    project=out/quality['project_filename']
    differences=write_refined_project(source,project,after,quality['settings'])
    qc=validate_project(project,{NAME:after},'ABS',cfg)
    after.export(out/'Straight_lumen_core_fine_vertical.stl')
    preview(before,after,out/'orientation_comparison.png')
    (out/'effective_quality_config.yaml').write_text(config_path.read_text())
    (out/'使用说明.md').write_text((ROOT/'docs/Straight_2Top血管芯精细打印.md').read_text())
    summary=dict(status='UNSLICED_PROCESS_CANDIDATE_PHYSICAL_TEST_REQUIRED',project=qc,
        geometry=audit,settings_changes=differences,
        source_project=str(source.relative_to(ROOT)),source_project_sha256=digest(source),
        source_assembly_sha256=digest(assembly),
        original_files_unchanged=all(digest(p)==h for p,h in protected.items()),
        support_settings_enabled=True,slicer_executed=False,generated_support_meshes=0,
        physical_print_tested=False,gui_open_tested=False)
    if not summary['original_files_unchanged']:raise ValueError('SOURCE_FILES_MODIFIED')
    # Old GUI evidence must never be reused for newly generated bytes.
    for n in ['bambu_gui_verified.png','gui_verification.json']:
        (out/n).unlink(missing_ok=True)
    save(out/'quality_summary.json',summary)
    dest=package(out)
    print('UNSLICED_CORE_QUALITY',project,'ZIP',dest,flush=True)
    return summary
