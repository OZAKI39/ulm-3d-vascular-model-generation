"""External engines only. Audited local CLI and 3MF volume roles; no G-code writing."""
from __future__ import annotations
import copy,json,re,shutil,subprocess,time,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from .support_inputs import ROOT,save,digest
from . import bambu_manufacturing as native

ALLOW={'--help','--version','--debug','--arrange','--load-settings','--load-filaments','--curr-bed-type','--slice','--export-3mf'}

def local_command(command):
    options={s.split('=')[0] for s in command if s.startswith('--')}
    if not options<=ALLOW:raise ValueError('NONLOCAL_OR_UNAUDITED_SLICER_COMMAND')
    return command

def cli_options(text):return sorted(set(re.findall(r'--[a-z][a-z0-9-]+',text)))

def discover_engine(name,discovery,out):
    paths=[]
    for command in (['orca-slicer','OrcaSlicer'] if name=='Orca' else ['bambu-studio']):
        if shutil.which(command):paths.append(Path(shutil.which(command)))
    if name=='Orca':
        for mount in Path('/mnt').glob('[a-z]'):
            paths.extend(p for p in [mount/'Program Files/OrcaSlicer/orca-slicer.exe',mount/'Program Files (x86)/OrcaSlicer/orca-slicer.exe'] if p.exists())
        paths.extend(sorted(Path('/mnt/c/Users').glob('*/AppData/Local/Temp/vascular_support_orca_*/orca-slicer.exe')))
        local=ROOT/'tools/external/OrcaSlicer/run-orca.sh'
        if local.exists():paths.append(local)
    elif discovery.get('binary') and Path(discovery['binary']).is_file():paths.append(Path(discovery['binary']))
    result=dict(engine=name,status='ENGINE_UNAVAILABLE',binary=None)
    if paths:
        binary=paths[-1];command=local_command([str(binary),'--help'])
        try:
            p=subprocess.run(command,cwd=binary.parent,capture_output=True,timeout=30)
            text=p.stdout.decode(errors='replace')+'\n'+p.stderr.decode(errors='replace')
            version=None
            if name=='Orca':
                v=re.search(r'OrcaSlicer-([\d.]+)',text)
                if v:version=v[1]
            elif binary.suffix!='.exe':
                v=re.search(r'BambuStudio-([\d.]+)',text)
                if v:version=v[1]
            if binary.suffix=='.exe':
                ps=Path('/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe');win=native.windows_path(binary)
                v=subprocess.run([str(ps),'-NoProfile','-NonInteractive','-Command',f"(Get-Item '{win}').VersionInfo.FileVersion"],capture_output=True,timeout=25)
                version=v.stdout.decode(errors='replace').strip()
            result.update(status='ENGINE_FOUND',binary=str(binary),binary_sha256=digest(binary),help_command=command,help_returncode=p.returncode,help_output=text,options_from_help=cli_options(text),version=version,
                cli_evidence='Actual help probe; Windows GUI build may redirect to CONOUT$. CLI flags also read in pinned PrintConfig.cpp and verified by real slice.')
            if discovery.get('native_binary'):
                # A headless wrapper is the launcher; cache identity must
                # follow the actual slicer ELF, not just an unchanged script.
                native_binary=Path(discovery['native_binary'])
                result.update(launcher_sha256=digest(binary),native_binary=str(native_binary),
                              binary_sha256=digest(native_binary))
        except (OSError,subprocess.TimeoutExpired) as exc:result.update(status='ENGINE_PROBE_FAILED',error=str(exc))
    save(Path(out)/(name.lower()+'_engine.json'),result);return result

def engine_profiles(discovery,engine,updates=None):
    profiles=copy.deepcopy(discovery['flattened'])
    if engine=='Orca' and (ROOT/'tools/external/OrcaSlicer/squashfs-root/resources/profiles/BBL').is_dir():
        # Translate against the installed engine's actual presets. New Bambu
        # multi-nozzle vectors are not a valid Orca single-nozzle schema.
        base=ROOT/'tools/external/OrcaSlicer/squashfs-root/resources/profiles/BBL'
        translated={};names={'machine':'Bambu Lab X1 Carbon 0.2 nozzle','process':'0.10mm Standard @BBL X1C 0.2 nozzle','filament':'Bambu ABS @BBL X1C 0.2 nozzle'}
        for kind,name in names.items():
            native_base,_=native.flatten(base/kind/(name+'.json'),search_dirs=[base/kind],resolve_includes=True)
            for key,original in list(native_base.items()):
                value=profiles[kind].get(key)
                if value is None or value==[] or key in ('name','type','from','instantiation','version','setting_id') or 'gcode' in key:continue
                if isinstance(original,list) and not isinstance(value,list):value=[value]
                if not isinstance(original,list) and isinstance(value,list):value=value[0]
                native_base[key]=value
            translated[kind]=native_base
        profiles=translated
    p=profiles['process'];p['enable_support']='1';p['support_type']='tree(auto)';p['support_style']='tree_organic'
    for kind,settings in profiles.items():
        settings['type']=kind;settings['instantiation']='true'
    if engine=='Orca':
        p['support_style']='organic'
        for settings in profiles.values():
            for key in [k for k,v in settings.items() if v==[]]:del settings[key]
        # Orca's traditional single-extruder process options are scalar. Bambu
        # 2.7 stores speed pairs for different nozzle systems; slot 0 is X1C.
        for k,v in list(p.items()):
            if isinstance(v,list) and k not in ('compatible_printers','compatible_printers_condition') and v:p[k]=v[0]
        p['tree_support_branch_angle_organic']=p.get('tree_support_branch_angle','40')
        p['tree_support_branch_diameter_organic']=p.get('tree_support_branch_diameter','2.5')
    p.update(updates or {})
    return profiles

def write_project(path,model,enforcer,blocker,profiles):
    """3MF core components + BBS part subtype, verified against pinned importer.

    Each mesh is in the SAME bed frame; modifiers never enter build as separate
    printable objects. No scaling/rearrangement or G-code is generated here.
    """
    ns='http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
    modelxml=ET.Element('model',{'unit':'millimeter','xmlns':ns})
    ET.SubElement(modelxml,'metadata',{'name':'Application'}).text='OrcaSlicer-2.4.2'
    ET.SubElement(modelxml,'metadata',{'name':'BambuStudio:3mfVersion'}).text='1'
    resources=ET.SubElement(modelxml,'resources');config=ET.Element('config');obj=ET.SubElement(config,'object',{'id':'4'})
    ET.SubElement(obj,'metadata',{'key':'name','value':'BG001 frozen 2 mm mold + support modifiers'})
    ET.SubElement(obj,'metadata',{'key':'extruder','value':'1'})
    for identity,mesh,role in [(1,model,'normal_part'),(2,enforcer,'support_enforcer'),(3,blocker,'support_blocker')]:
        o=ET.SubElement(resources,'object',{'id':str(identity),'type':'model'});m=ET.SubElement(o,'mesh');vs=ET.SubElement(m,'vertices');fs=ET.SubElement(m,'triangles')
        for v in mesh.vertices:ET.SubElement(vs,'vertex',dict(zip('xyz',[format(float(x),'.10g') for x in v])))
        for face in mesh.faces:ET.SubElement(fs,'triangle',dict(zip(['v1','v2','v3'],map(str,face))))
        part=ET.SubElement(obj,'part',{'id':str(identity),'subtype':role});ET.SubElement(part,'metadata',{'key':'name','value':role})
        ET.SubElement(part,'metadata',{'key':'matrix','value':'1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1'})
        ET.SubElement(part,'mesh_stat',{'face_count':str(len(mesh.faces)),'edges_fixed':'0','degenerate_facets':'0','facets_removed':'0','facets_reversed':'0','backwards_edges':'0'})
    group=ET.SubElement(resources,'object',{'id':'4','type':'model'});components=ET.SubElement(group,'components')
    for i in [1,2,3]:ET.SubElement(components,'component',{'objectid':str(i)})
    ET.SubElement(ET.SubElement(modelxml,'build'),'item',{'objectid':'4','printable':'1'})
    settings={}
    for k in ['machine','filament','process']:settings.update(profiles[k])
    settings.update(curr_bed_type='Textured PEI Plate',printer_settings_id=profiles['machine']['name'],print_settings_id=profiles['process']['name'],filament_settings_id=[profiles['filament']['name']])
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml','<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>')
        z.writestr('_rels/.rels','<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>')
        z.writestr('3D/3dmodel.model',ET.tostring(modelxml,encoding='utf-8',xml_declaration=True));z.writestr('Metadata/model_settings.config',ET.tostring(config,encoding='utf-8',xml_declaration=True));z.writestr('Metadata/project_settings.config',json.dumps(settings))
    return project_roles(path)

def project_roles(path):
    with zipfile.ZipFile(path) as z:
        cfg=ET.fromstring(z.read('Metadata/model_settings.config'))
        roles=[p.attrib.get('subtype') for p in cfg.findall('.//part')]
    return roles

def production_geometry(path,expected):
    """Ignore non-printable modifier meshes when checking immutable geometry."""
    import trimesh
    from scipy.spatial import cKDTree
    ns='{http://schemas.microsoft.com/3dmanufacturing/core/2015/02}';prod='{http://schemas.microsoft.com/3dmanufacturing/production/2015/06}'
    def matrix(text):
        m=np.eye(4)
        if text:
            a=np.fromstring(text,sep=' ').reshape(4,3);m[:3,:3]=a[:3].T;m[:3,3]=a[3]
        return m
    with zipfile.ZipFile(path) as z:
        config=ET.fromstring(z.read('Metadata/model_settings.config'));parts={p.attrib['id']:p.attrib.get('subtype') for p in config.findall('.//part')};cache={}
        def objects(file):
            if file not in cache:cache[file]={o.attrib['id']:o for o in ET.fromstring(z.read(file)).findall(ns+'resources/'+ns+'object')}
            return cache[file]
        def visit(file,identity,m):
            o=objects(file)[identity];mesh=o.find(ns+'mesh');points=[]
            if mesh is not None and parts.get(identity,'normal_part')=='normal_part':
                v=np.array([[float(p.attrib[k]) for k in 'xyz'] for p in mesh.findall(ns+'vertices/'+ns+'vertex')]);points.append(v@m[:3,:3].T+m[:3,3])
            for c in o.findall(ns+'components/'+ns+'component'):
                points.extend(visit(c.attrib.get(prod+'path',file).lstrip('/'),c.attrib['objectid'],m@matrix(c.attrib.get('transform'))))
            return points
        root=ET.fromstring(z.read('3D/3dmodel.model'));allpoints=[]
        for item in root.findall(ns+'build/'+ns+'item'):allpoints.extend(visit('3D/3dmodel.model',item.attrib['objectid'],matrix(item.attrib.get('transform'))))
        points=np.vstack(allpoints);target=trimesh.load_mesh(expected).vertices
        distance=float(max(cKDTree(points).query(target)[0].max(),cKDTree(target).query(points)[0].max()))
        settings=json.loads(z.read('Metadata/project_settings.config'))
    return dict(passed=distance<2e-4,max_vertex_distance_mm=distance,expected_sha256=digest(expected),roles=list(parts.values()),
        actual_settings={k:settings.get(k) for k in ['layer_height','initial_layer_print_height','line_width','nozzle_diameter','support_top_z_distance','support_object_xy_distance','filament_type','nozzle_temperature','textured_plate_temp','support_type','support_style']})

def warnings(record):
    native_result=record.get('native_result',{});items=[]
    def visit(v):
        if isinstance(v,dict):
            for k,val in v.items():
                if any(word in k.lower() for word in ['warning','error']) and val and val!='Success.':items.append(str(val))
                elif isinstance(val,(dict,list)):visit(val)
        elif isinstance(v,list):
            for x in v:visit(x)
    visit(native_result)
    items.extend(str(x) for x in record.get('warnings',[]))
    critical=[x for x in items if any(w in x.lower() for w in ['floating','unsupported','empty layer','overhang'])]
    return dict(all=items,critical=critical)

def slice_local(engine,profiles,source,out,timeout=900):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if engine['status']!='ENGINE_FOUND':
        r=dict(status='ENGINE_UNAVAILABLE',slicer_success=False,engine=engine['engine']);save(out/'slice_result.json',r);return r
    binary=Path(engine['binary']);win=binary.suffix=='.exe'
    staging=binary.parent/('bg001_support_'+out.name) if win else out/'staging';staging.mkdir(exist_ok=True)
    argpath=native.windows_path if win else str
    paths={}
    for k,v in profiles.items():paths[k]=staging/(k+'.json');save(paths[k],v);save(out/'profiles'/(k+'.json'),v)
    model=staging/('input'+Path(source).suffix);shutil.copy2(source,model);dest=staging/'sliced.3mf'
    if dest.exists():dest.unlink()
    for p in [staging/'result.json']:
        if p.exists():p.unlink()
    command=local_command([str(binary),'--debug','2','--arrange','0','--load-settings',argpath(paths['machine'])+';'+argpath(paths['process']),
        '--load-filaments',argpath(paths['filament']),'--curr-bed-type','Textured PEI Plate','--slice','0','--export-3mf',argpath(dest),argpath(model)])
    save(out/'command.json',dict(command=command,local_only=True,input_sha256=digest(source)))
    start=time.monotonic();record=dict(engine=engine['engine'],command=command,status='SLICE_FAILED',slicer_success=False,printer_upload=False,printer_start=False)
    try:
        run=subprocess.run(command,cwd=staging,capture_output=True,timeout=timeout)
        record.update(exit_code=run.returncode,elapsed_seconds=time.monotonic()-start)
        (out/'stdout.log').write_bytes(run.stdout);(out/'stderr.log').write_bytes(run.stderr)
        if (staging/'result.json').exists():record['native_result']=json.loads((staging/'result.json').read_text())
        if dest.exists():
            archive=out/'sliced.3mf';shutil.copy2(dest,archive);record.update(native.parse_sliced_archive(archive));record['archive']=str(archive)
            record['slicer_success']=bool(record['slice_success'] and run.returncode==0);record['status']='SLICED' if record['slicer_success'] else 'SLICE_FAILED'
            record['roles']=project_roles(archive)
            with zipfile.ZipFile(archive) as z:
                for name in z.namelist():
                    if name.endswith('.gcode'): (out/Path(name).name).write_bytes(z.read(name))
                    elif name=='Metadata/plate_1.png':(out/'native_preview.png').write_bytes(z.read(name))
    except (OSError,subprocess.TimeoutExpired,ValueError) as exc:record['error']=str(exc)
    record['warning_audit']=warnings(record);save(out/'slice_result.json',record);return record
