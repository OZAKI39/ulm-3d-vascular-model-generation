"""Discover local Bambu presets, flatten official inheritance and slice locally.

No printer/device/network commands are issued. Active GUI presets are read only.
"""
from __future__ import annotations

import json
import io
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import zipfile
import xml.etree.ElementTree as ET

import numpy as np

from .topbrain_qc import sha256

OFFICIAL_CLI = 'https://github.com/bambulab/BambuStudio/wiki/Command-Line-Usage'
A1_MATERIAL_SOURCE = 'https://bambulab.com/en/a1/tech-specs'


def json_object(path):
    # BambuStudio.conf can have trailing data after its first complete JSON object.
    return json.JSONDecoder().raw_decode(Path(path).read_text(encoding='utf-8-sig').lstrip())[0]


def flatten(path, seen=()):
    path=Path(path)
    if path in seen:
        raise ValueError('Cyclic Bambu preset inheritance')
    data=json_object(path);chain=[]
    parent=data.get('inherits')
    if parent:
        base,chain=flatten(path.parent/(parent+'.json'),(*seen,path))
        base.update(data);data=base
    data.pop('inherits',None)
    return data,[*chain,dict(path=str(path),sha256=sha256(path))]


def discover(config):
    binaries=[];dirs=[]
    for name in ['bambu-studio','BambuStudio','bambu-studio.exe']:
        found=shutil.which(name)
        if found:binaries.append(Path(found))
    for base in [Path.home()/'.config/BambuStudio',Path.home()/'.local/share/BambuStudio']:
        if base.is_dir():dirs.append(base)
    for mount in sorted(Path('/mnt').glob('[a-z]')):
        for base in [mount/'Program Files/Bambu Studio',mount/'Program Files (x86)/Bambu Studio']:
            if (base/'bambu-studio.exe').is_file():binaries.append(base/'bambu-studio.exe')
        users=mount/'Users'
        if users.is_dir():
            for user in users.iterdir():
                base=user/'AppData/Roaming/BambuStudio'
                try:
                    if base.is_dir():dirs.append(base)
                except PermissionError:
                    continue
    for base in [Path('/opt'),Path.home()/'Applications',Path.home()/'.local/bin']:
        if base.is_dir():
            binaries += [p for p in base.glob('*Bambu*') if p.is_file() and os.access(p,os.X_OK)]
    report=dict(binary_candidates=[str(p) for p in dict.fromkeys(binaries)],config_directories=[str(p) for p in dirs],
        model='GENERIC_BAMBU_256_'+str(config['printer']['nozzle_diameter_mm']).replace('.','P'),status='BAMBU_PROFILE_NOT_CONFIRMED',binary=None,version=None,
        nozzle_diameter_mm=config['printer']['nozzle_diameter_mm'],filament=config['printer']['filament'],
        build_volume_mm=config['printer']['build_volume_mm'],warnings=[],official_cli_reference=OFFICIAL_CLI)
    if binaries:
        report['binary']=str(binaries[0])
        try:
            result=subprocess.run([str(binaries[0]),'--help'],cwd=binaries[0].parent,
                capture_output=True,timeout=config['slicer']['help_timeout_seconds'])
            report['help_probe']=dict(command=[str(binaries[0]),'--help'],exit_code=result.returncode,
                stdout=result.stdout.decode(errors='replace'),stderr=result.stderr.decode(errors='replace'))
        except (OSError,subprocess.TimeoutExpired) as exc:
            report['help_probe']=dict(error=str(exc))
        if binaries[0].suffix.lower()=='.exe':
            ps=Path('/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe')
            if ps.is_file():
                win=windows_path(binaries[0]).replace("'","''")
                result=subprocess.run([str(ps),'-NoProfile','-NonInteractive','-Command',
                    f"(Get-Item '{win}').VersionInfo.FileVersion"],capture_output=True,timeout=20)
                report['version']=result.stdout.decode(errors='replace').strip() or None
    active=None;root=None
    for directory in dirs:
        path=directory/'BambuStudio.conf'
        if not path.is_file():continue
        try:
            presets=json_object(path).get('presets',{})
        except (OSError,ValueError):continue
        if presets.get('machine'):
            active={k:presets.get(k) for k in ['machine','process','filaments']}
            root=directory;report['active_preset_source']=str(path);break
    if not active:return report
    report['active_gui_presets']=active
    system=root/'system/BBL'
    machine_path=system/'machine'/(active['machine']+'.json')
    if not machine_path.is_file():return report
    machine,mc=flatten(machine_path)
    nozzle=float(machine['nozzle_diameter'][0]);model=machine['printer_model']
    # Explicit alternative nozzle profiles keep the real model but select a
    # matching preset; this does not assert a hardware nozzle has been changed.
    requested=float(config['printer']['nozzle_diameter_mm'])
    if abs(nozzle-requested)>1e-8:
        alternate=system/'machine'/f'{model} {requested:g} nozzle.json'
        if not alternate.is_file():return report
        machine_path=alternate;machine,mc=flatten(machine_path)
    process_name=machine.get('default_print_profile') if abs(nozzle-requested)>1e-8 else active['process']
    process_path=system/'process'/(process_name+'.json')
    filament_path=system/'filament'/f'Bambu {config["printer"]["filament"]} @BBL {model.removeprefix("Bambu Lab ")} {requested:g} nozzle.json'
    if not filament_path.is_file():
        fallback=system/'filament'/f'Bambu {config["printer"]["filament"]} @BBL {model.removeprefix("Bambu Lab ")}.json'
        filament_path=fallback
    if not process_path.is_file() or not filament_path.is_file():
        report['warnings'].append('Matching official process/filament preset unavailable')
        return report
    process,pc=flatten(process_path);filament,fc=flatten(filament_path)
    area=np.array([[float(x) for x in point.split('x')] for point in machine['printable_area']])
    volume=[*np.ptp(area,axis=0).tolist(),float(machine['printable_height'])]
    report.update(model=model,status='BAMBU_PROFILE_CONFIRMED',nozzle_diameter_mm=requested,
        active_hardware_nozzle_matches=requested==nozzle,build_volume_mm=volume,
        machine_profile=str(machine_path),process_profile=str(process_path),filament_profile=str(filament_path),
        profile_sources=mc+pc+fc,flattened=dict(machine=machine,process=process,filament=filament),
        active_gui_presets_modified=False,material_selection='Explicit task ABS; GUI active filament retained unchanged')
    if 'A1' in model and config['printer']['filament']=='ABS':
        report['warnings'].append('ABS is not recommended by the manufacturer for this printer profile.')
        report['material_warning_source']=A1_MATERIAL_SOURCE
    return report


def windows_path(path):
    return subprocess.check_output(['wslpath','-w',str(Path(path).resolve())],text=True).strip()


def prepare_slicer_profiles(discovery,config,directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    profiles=discovery['flattened']
    process=dict(profiles['process'])
    process['enable_support']='1' if config['slicer']['enable_support'] else '0'
    profiles={**profiles,'process':process}
    paths={}
    for key,data in profiles.items():
        path=directory/(key+'.json');path.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n');paths[key]=path
    return paths


def parse_sliced_archive(path):
    with zipfile.ZipFile(path) as archive:
        names=archive.namelist();gcodes=[n for n in names if n.endswith('.gcode')]
        metadata={};filaments=[];warnings_found=[]
        if 'Metadata/slice_info.config' in names:
            tree=ET.fromstring(archive.read('Metadata/slice_info.config'))
            for element in tree.iter():
                key=element.attrib.get('key')
                if key:
                    metadata[key]=element.attrib.get('value')
                    if 'warning' in key.lower():warnings_found.append(dict(element.attrib))
                if element.tag=='filament':filaments.append(dict(element.attrib))
        headers={};feature_names=set();extrusion_lines=0
        for name in gcodes:
            text=archive.read(name).decode(errors='replace')
            for line in text.splitlines():
                if line.startswith(';') and '=' in line:
                    key,value=line[1:].split('=',1)
                    if any(token in key.lower() for token in ['estimated','filament used','support']):headers[key.strip()]=value.strip()
                if line.startswith('; FEATURE:'):feature_names.add(line.split(':',1)[1].strip())
                if re.match(r'G[01]\s.*\bE[-+\d.]',line):extrusion_lines+=1
        prediction=metadata.get('prediction')
        seconds=float(prediction) if prediction and re.fullmatch(r'[\d.]+',prediction) else None
        # Return actual exposed fields. Missing support quantities stay null.
        support_settings={k:v for k,v in headers.items() if 'support' in k.lower()}
        support=dict(support_used=metadata.get('support_used'),support_type=headers.get('support_type'),
            filament_usage_g=None,quantity_source='Separate support quantity not exposed by this archive',
            feature_tags=sorted(x for x in feature_names if 'support' in x.lower()))
        project=json.loads(archive.read('Metadata/project_settings.config')) if 'Metadata/project_settings.config' in names else {}
        first_layer=metadata.get('first_layer_time')
        if first_layer is not None and seconds is not None:
            try:
                if not 0<=float(first_layer)<=seconds:
                    warnings_found.append(dict(code='UNRELIABLE_FIRST_LAYER_TIME',value=first_layer,used_for_selection=False))
            except ValueError:pass
        usage=sum(float(f['used_g']) for f in filaments if 'used_g' in f) if any('used_g' in f for f in filaments) else None
        return dict(slice_success=bool(gcodes and extrusion_lines),gcode_files=gcodes,extrusion_lines=extrusion_lines,
            estimated_print_time_seconds=seconds,filament_used_g=usage,filament_records=filaments,
            support_information=support,support_settings=support_settings,
            actual_project_settings={k:project.get(k) for k in ['printer_model','printer_settings_id','print_settings_id',
                'nozzle_diameter','filament_type','filament_settings_id','printable_height','printable_area','enable_support','support_type']},
            actual_metadata=metadata,actual_gcode_headers=headers,warnings=warnings_found)


def sliced_geometry_qc(path,stl):
    """Apply actual 3MF component/build matrices and compare all vertices to STL."""
    import pyvista as pv
    from scipy.spatial import cKDTree
    def matrix(value):
        result=np.eye(4)
        if value:
            array=np.fromstring(value,sep=' ').reshape(4,3)
            result[:3,:3]=array[:3].T;result[:3,3]=array[3]
        return result
    namespace='{http://schemas.microsoft.com/3dmanufacturing/core/2015/02}'
    production='{http://schemas.microsoft.com/3dmanufacturing/production/2015/06}'
    with zipfile.ZipFile(path) as archive:
        cache={}
        def objects(filename):
            if filename in cache:return cache[filename]
            result={};current=None
            for event,element in ET.iterparse(io.BytesIO(archive.read(filename)),events=('start','end')):
                tag=element.tag.rsplit('}',1)[-1]
                if event=='start' and tag=='object':
                    current={'points':[],'components':[]};result[element.attrib['id']]=current
                elif event=='end':
                    if tag=='vertex':current['points'].append([float(element.attrib[k]) for k in ('x','y','z')])
                    elif tag=='component':current['components'].append(dict(element.attrib))
                    if tag in ['vertex','triangle']:element.clear()
            cache[filename]=result;return result
        def collect(filename,identity,parent):
            obj=objects(filename)[identity];chunks=[]
            if obj['points']:
                chunks.append(np.asarray(obj['points'])@parent[:3,:3].T+parent[:3,3])
            for child in obj['components']:
                linked=child.get(production+'path',filename).lstrip('/')
                chunks.extend(collect(linked,child['objectid'],parent@matrix(child.get('transform'))))
            return chunks
        model=ET.fromstring(archive.read('3D/3dmodel.model'));chunks=[]
        for item in model.find(namespace+'build'):
            if item.attrib.get('printable','1')=='1':
                chunks.extend(collect('3D/3dmodel.model',item.attrib['objectid'],matrix(item.attrib.get('transform'))))
        points=np.vstack(chunks)
    expected=pv.read(stl).points
    forward=cKDTree(expected).query(points,workers=1)[0]
    backward=cKDTree(points).query(expected,workers=1)[0]
    distance=float(max(forward.max(),backward.max()))
    return dict(rotation_scale_position_preserved=distance<2e-4,max_bidirectional_vertex_distance_mm=distance,
        sliced_bbox_extents_mm=np.ptp(points,axis=0).tolist(),input_bbox_extents_mm=np.ptp(expected,axis=0).tolist(),
        tolerance_mm=2e-4,input_stl_sha256=sha256(stl),sliced_3mf_sha256=sha256(path))


def slice_candidates(discovery,config,candidates,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    if not discovery.get('binary') or discovery['status']!='BAMBU_PROFILE_CONFIRMED':
        return dict(status='BAMBU_SLICE_VALIDATION_SKIPPED',reason='CLI or confirmed printer/process/filament profile unavailable',results=[])
    binary=Path(discovery['binary']);is_windows=binary.suffix.lower()=='.exe'
    if is_windows:
        conf=Path(discovery['active_preset_source'])
        staging=conf.parents[2]/'Local/Temp'/('bg001_manufacturing_'+str(time.time_ns()))
    else:staging=output/'staging'
    staging.mkdir(parents=True)
    paths=prepare_slicer_profiles(discovery,config,staging)
    # Also save only the three official flattened presets in the project, never
    # the GUI config (which may contain account/device information).
    local_profiles=output/'profiles';local_profiles.mkdir(exist_ok=True)
    for p in paths.values():shutil.copy2(p,local_profiles/p.name)
    argpath=windows_path if is_windows else lambda p:str(p.resolve())
    results=[]
    for item in candidates:
        rank=item['rank'];stem=f'orientation_candidate_{rank:02d}'
        model=staging/(stem+'.stl');shutil.copy2(item['stl'],model)
        destination=staging/(stem+'.3mf')
        command=[str(binary),'--debug','2','--arrange','0',
            '--load-settings',argpath(paths['machine'])+';'+argpath(paths['process']),
            '--load-filaments',argpath(paths['filament']),'--curr-bed-type',config['slicer']['bed_type'],
            '--slice','0','--export-3mf',argpath(destination),argpath(model)]
        assert '--orient' not in command and '--scale' not in command
        print('BAMBU SLICE START',stem,flush=True)
        record=dict(rank=rank,orientation_id=item['candidate_id'],command=command,
                    rotation_overridden=False,status='BAMBU_SLICE_VALIDATION_SKIPPED')
        try:
            start=time.monotonic()
            result=subprocess.run(command,cwd=binary.parent,capture_output=True,timeout=config['slicer']['timeout_seconds'])
            record.update(exit_code=result.returncode,elapsed_seconds=time.monotonic()-start,
                stdout=result.stdout.decode(errors='replace'),stderr=result.stderr.decode(errors='replace'))
            if destination.is_file():
                local=output/destination.name;shutil.copy2(destination,local)
                record['archive']=str(local);record.update(parse_sliced_archive(local))
                record['status']='BAMBU_SLICED' if result.returncode==0 and record['slice_success'] else 'BAMBU_SLICE_FAILED'
                record['geometry_qc']=sliced_geometry_qc(local,item['stl'])
                settings=record['actual_project_settings']
                record['printer_filament_profile_matches']=(settings['printer_model']==discovery['model']
                    and config['printer']['filament'] in (settings['filament_type'] or [])
                    and float(settings['nozzle_diameter'][0])==discovery['nozzle_diameter_mm'])
                if not record['geometry_qc']['rotation_scale_position_preserved'] or record['actual_metadata'].get('outside')=='true':
                    record['status']='BAMBU_SLICE_FAILED'
                if not record['printer_filament_profile_matches']:
                    record['status']='BAMBU_PROFILE_MISMATCH'
            else:record['reason']='CLI produced no 3MF; headless slicing not confirmed'
        except (OSError,subprocess.TimeoutExpired,ValueError,zipfile.BadZipFile) as exc:
            record['reason']=str(exc)
        results.append(record)
        (output/(stem+'_slice.json')).write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')
        print('BAMBU SLICE END',stem,record['status'],flush=True)
    successes=[r for r in results if r['status']=='BAMBU_SLICED']
    return dict(status='BAMBU_SLICED' if successes else 'BAMBU_SLICE_VALIDATION_SKIPPED',
                staging=str(staging),results=results,successful_count=len(successes))
