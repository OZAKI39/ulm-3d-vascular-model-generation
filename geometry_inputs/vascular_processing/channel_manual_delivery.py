"""Editable single-material projects and model bundles. Never run a slicer."""
from pathlib import Path
import copy
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import trimesh

from . import y_channel_casting as y, straight_channel_casting as s
from .all_abs_print_plate import read_objects

PROJECT_MEMBERS={'[Content_Types].xml','_rels/.rels','3D/3dmodel.model',
                 'Metadata/model_settings.config','Metadata/project_settings.config'}

# These labels are also used as filenames by Bambu Studio. XML escaping does
# not make a forbidden filename character valid after the project is loaded.
FORBIDDEN_NAME_CHARS=re.compile(r'[<>:"/\\|?*\x00-\x1f\x7f]')


def validate_bambu_names(filename,model_xml,settings_xml):
    model=ET.fromstring(model_xml);settings=ET.fromstring(settings_xml)
    names=[('project filename',Path(filename).name)]
    for obj in model.iter():
        if obj.tag.rsplit('}',1)[-1]=='object' and 'name' in obj.attrib:
            names.append(('model object',obj.get('name')))
        if obj.tag.rsplit('}',1)[-1]=='metadata' and obj.get('name')=='Title':
            names.append(('project title',obj.text or ''))
    for meta in settings.iter('metadata'):
        if meta.get('key') in ('name','plater_name'):
            names.append((meta.get('key'),meta.get('value','')))
    for field,value in names:
        if not value.strip() or FORBIDDEN_NAME_CHARS.search(value):
            raise ValueError(f'BAMBU_INVALID_NAME: {field}={value!r}')
    return dict(passed=True,checked_name_count=len(names),
        plate_names=[v for k,v in names if k=='plater_name'])


def validate_gui_filaments(settings):
    """Match the GUI spool/variant checks that CLI --info does not exercise."""
    colors=settings.get('filament_colour')
    if (not isinstance(colors,list) or len(colors)!=1
            or not isinstance(colors[0],str) or not re.fullmatch(r'#[0-9a-fA-F]{6}',colors[0])):
        raise ValueError('BAMBU_GUI_FILAMENT_COLOUR_REQUIRED')
    variants=settings.get('filament_extruder_variant') or ['Direct Drive Standard']
    indices=settings.get('filament_self_index')
    if (not isinstance(variants,list) or not isinstance(indices,list)
            or len(indices)!=len(variants) or any(str(i)!='1' for i in indices)):
        raise ValueError('BAMBU_GUI_FILAMENT_VARIANT_MAPPING_INVALID')
    for field in ['filament_ids','filament_settings_id','filament_type']:
        value=settings.get(field)
        if not isinstance(value,list) or len(value)!=1 or not value[0]:
            raise ValueError('BAMBU_GUI_SINGLE_FILAMENT_METADATA_REQUIRED: '+field)
    return dict(passed=True,spool_count=1,extruder_variant_count=len(variants),
                variant_spool_indices=indices,filament_colour=colors)


def validate_project(path,expected,material,cfg):
    with zipfile.ZipFile(path) as z:
        if set(z.namelist())!=PROJECT_MEMBERS or z.testzip() is not None:
            raise ValueError('UNSLICED_PROJECT_HAS_UNEXPECTED_MEMBERS')
        name_qc=validate_bambu_names(path.name,z.read('3D/3dmodel.model'),
                                    z.read('Metadata/model_settings.config'))
    objects,settings,nplates,crc=read_objects(path)
    filament_qc=validate_gui_filaments(settings)
    if (not crc or nplates!=1 or len(objects)!=len(expected)
            or {o['metadata']['name'] for o in objects}!=set(expected)):
        raise ValueError('MANUAL_PROJECT_OBJECTS_MISMATCH')
    if (settings.get('printer_model')!='Bambu Lab P1S' or settings.get('nozzle_diameter')!=['0.4']
            or settings.get('filament_type')!=[material]):
        raise ValueError('MANUAL_PROJECT_MATERIAL_OR_MACHINE_MISMATCH')
    if not y.slicer_scale_qc(cfg,settings)['passed']:raise ValueError('MANUAL_PROJECT_UNEXPECTED_SCALING')
    rows=[]
    for obj in objects:
        name=obj['metadata']['name'];mesh=obj['mesh'];source=expected[name]
        if mesh.vertices.shape!=source.vertices.shape or mesh.faces.shape!=source.faces.shape:
            raise ValueError('MANUAL_PROJECT_GEOMETRY_MISMATCH')
        error=float(abs(mesh.vertices-source.vertices).max())
        if error>1e-5 or not np.array_equal(mesh.faces,source.faces) or not mesh.is_volume:
            raise ValueError('MANUAL_PROJECT_GEOMETRY_MISMATCH')
        actual_material='ABS' if 'core' in name else 'PLA'
        if material!=actual_material or obj['metadata']['extruder']!='1':
            raise ValueError('MANUAL_PROJECT_PART_MATERIAL_MISMATCH')
        rows.append(dict(name=name,material=material,mesh_error_mm=error,bounds_mm=mesh.bounds,
            suggested_support_enabled=obj['metadata']['enable_support'],
            suggested_raft_layers=obj['metadata']['raft_layers']))
    return dict(passed=True,file=path.name,sha256=y.digest(path),units='mm',plate_count=1,
        material=material,parts=rows,printer='Bambu Lab P1S',nozzle_mm=.4,zip_crc_passed=True,
        has_slice_or_toolpath_members=False,settings_are_editable=True,
        name_qc=name_qc,gui_filament_qc=filament_qc,gui_open_tested=False,slicer_executed=False)


def export_delivery(cfg,parts,summary,out,kind):
    if cfg.get('delivery_mode')!='manual_slicing' or cfg.get('part_materials')!={'core':'ABS','fixture':'PLA'}:
        raise ValueError('MANUAL_ABS_CORE_PLA_FIXTURE_POLICY_REQUIRED')
    if summary['status']!='GEOMETRY_PASSED' or not all(summary['gates'].values()):
        raise ValueError('MANUAL_PROJECT_REQUIRES_VALID_GEOMETRY')
    # Do not hide stale sliced outputs inside a new user-slicing delivery.
    if (out/'bambu').exists() or (out/'print_qc.json').exists() or any(out.glob('*sliced*.3mf')):
        raise ValueError('ARCHIVE_OLD_SLICED_OUTPUTS_BEFORE_MANUAL_EXPORT')
    module=y if kind=='Y' else s
    meshes,records=module.arrange(parts,out)
    prefix='Y' if kind=='Y' else ('Straight_2Top' if s.two_top_ports(cfg) else 'Straight_3Port')
    projects=[];assignments=[]
    for role,material in [('core','ABS'),('fixture','PLA')]:
        group=[(mesh,row) for mesh,row in zip(meshes,records) if ('core' in row['name'])==(role=='core')]
        gm,gr=map(list,zip(*group));c=copy.deepcopy(cfg)
        if material=='PLA':c['printer']['filament']=cfg['fixture_filament']
        profiles,engine=y.get_profiles(c,out/'profiles'/material,material=material,discover=False)
        assert engine is None
        path=out/(prefix+'_P1S_04_'+role+'_'+material+'_project.3mf')
        module.write_project(path,gm,gr,profiles,c)
        with zipfile.ZipFile(path) as z:data={n:z.read(n) for n in z.namelist()}
        model=ET.fromstring(data['Metadata/model_settings.config'])
        for meta in model.findall('.//plate/metadata'):
            if meta.get('key')=='plater_name':meta.set('value',prefix+' - '+role+' '+material+' - UNSLICED')
        data['Metadata/model_settings.config']=ET.tostring(model,encoding='utf-8',xml_declaration=True)
        with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
            for n,b in data.items():z.writestr(n,b)
        projects.append(validate_project(path,dict(zip([r['name'] for r in gr],gm)),material,cfg))
        assignments.extend(dict(part=r['name'],material=material,project=path.name) for r in gr)
    # Write saved configuration first, because Rhino readback checks it.
    y.save(out/'effective_config.json',cfg)
    if cfg.get('rhino_export'):
        from .rhino_channel_export import export_models
        export_models(cfg,parts,summary,out)
    result=dict(status='READY_FOR_USER_SLICING',slicer_executed=False,toolpaths_generated=False,
        printer_upload=False,printer_start=False,part_materials=assignments,projects=projects,
        geometry_status=summary['status'],existing_cad_geometry_preserved=True,
        actual_print_time_seconds=None,material_mass_g=None,
        material_policy='Only the lumen core is ABS; the complete frame including magnets, all retainers and keys are PLA',
        process_policy='Two separate single-material projects; editable starting presets; supports generated by user slicing')
    y.save(out/'manual_slicing_summary.json',result)
    package(cfg,out,result,kind)
    print('MANUAL_DELIVERY',kind,result['status'],[p['file'] for p in projects],flush=True)
    return result


def package(cfg,out,result,kind):
    if result['status']!='READY_FOR_USER_SLICING' or result['slicer_executed']:
        raise ValueError('MANUAL_DELIVERY_STATUS_MISMATCH')
    names=[r['part'] for r in result['part_materials']]
    files=[out/p['file'] for p in result['projects']]
    files += [out/'assembly'/(n+'.stl') for n in names]
    files += [out/'print_parts'/(n+'.stl') for n in names]
    files += sorted((out/'assembly').glob('*.step'))+sorted(out.glob('labelled_*_mm.vtp'))
    files += [out/n for n in ['geometry_qc.json','magnet_slots_qc.json','manual_slicing_summary.json',
                              'print_layout.json','effective_config.json','rhino_qc.json']]
    files += sorted(out.glob('*.3dm'))+sorted((out/'reference').iterdir())
    files += [p for p in (out/'QC').glob('*.png') if not p.name.startswith('native_')]
    files += [out/'profiles'/m/'profiles.json' for m in ['ABS','PLA']]
    from .rhino_channel_export import verify_exports
    verify_exports(out)
    members={str(p.relative_to(out)):p for p in files}
    members['使用说明.md']=y.ROOT/'docs/ABS血管芯_PLA夹具_自行切片使用说明.md'
    hashes={n:y.digest(p) for n,p in members.items()}
    dest=y.ROOT/'deliverables'/(( 'Y_channel_' if kind=='Y' else 'Straight_')+out.name+'.zip')
    dest.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED) as z:
        for n,p in members.items():z.write(p,kind+'_channel/'+n)
        z.writestr(kind+'_channel/SHA256.json',json.dumps(hashes,indent=2,ensure_ascii=False))
    with zipfile.ZipFile(dest) as z:
        if z.testzip() is not None:raise ValueError('MANUAL_ZIP_CRC_FAILED')
        for n,h in hashes.items():
            if hashlib.sha256(z.read(kind+'_channel/'+n)).hexdigest()!=h:raise ValueError('MANUAL_ZIP_HASH_FAILED')
        if any('sliced.3mf' in n.lower() or n.lower().endswith(('.gcode','.bgcode')) for n in z.namelist()):
            raise ValueError('MANUAL_ZIP_CONTAINS_TOOLPATH')
    y.save(out/'delivery_receipt.json',dict(zip=dest,zip_sha256=y.digest(dest),zip_bytes=dest.stat().st_size,
        members_sha256=hashes,crc_and_hashes_passed=True,status='READY_FOR_USER_SLICING',
        sliced_artifacts_included=False))
