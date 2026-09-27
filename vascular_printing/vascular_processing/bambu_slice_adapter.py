"""Read active Bambu profiles and run local slicing only; never send a job."""
import copy
import json
from pathlib import Path
import re
import shutil
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
from . import bambu_manufacturing as existing
from .sacrificial_fixture_review import save_json

FORBIDDEN={'--send','--send-to-printer','--upload','--print','--start-print','--host','--access-code','--device','--device-id'}
LOCAL_OPTIONS={'--debug','--arrange','--load-settings','--load-filaments','--curr-bed-type','--slice','--export-3mf'}


def check_local_command(command):
    options={arg.split('=')[0] for arg in command if arg.startswith('--')}
    if options & FORBIDDEN or not options <= LOCAL_OPTIONS:
        raise ValueError('NONLOCAL_OR_UNAUDITED_SLICER_COMMAND')
    return True


def discover_current(cfg,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    config=dict(printer=dict(nozzle_diameter_mm=.4,filament=cfg['printer']['filament'],build_volume_mm=cfg['printer']['provisional_build_volume_mm']),
        slicer=dict(help_timeout_seconds=cfg['bambu']['help_timeout_seconds']))
    report=existing.discover(config)
    # .4 above is only a discovery seed, never an unconfirmed hardware choice.
    if report.get('active_gui_presets') and not report.get('active_hardware_nozzle_matches',False):
        root=Path(report['active_preset_source']).parent
        path=root/'system/BBL/machine'/(report['active_gui_presets']['machine']+'.json')
        if path.is_file():
            actual,_=existing.flatten(path);nozzles=actual.get('nozzle_diameter',[])
            if len(nozzles)==1:
                config['printer']['nozzle_diameter_mm']=float(nozzles[0]);report=existing.discover(config)
    if not report.get('binary'):report['status']='BAMBU_STUDIO_NOT_FOUND'
    elif report.get('status')!='BAMBU_PROFILE_CONFIRMED' or not report.get('active_hardware_nozzle_matches'):
        report['status']='PROFILE_UNRESOLVED'
    # Record what help ACTUALLY returned; Windows GUI launchers may be silent.
    help_record=report.get('help_probe',{})
    help_text=help_record.get('stdout','')+'\n'+help_record.get('stderr','')
    report['available_cli_options_from_help']=sorted(set(re.findall(r'--[a-z][a-z0-9-]+',help_text)))
    report['help_output_available']=bool(help_text.strip())
    report['cli_options_evidence']='Actual --help probe; locally executed slice flags independently verified in each result.'
    report['official_cli_reference']=existing.OFFICIAL_CLI
    save_json(out/'cli_probe.json',help_record)
    (out/'cli_help.txt').write_text(help_text if help_text.strip() else '本机 Windows 启动器执行 --help 返回 0，但未返回标准输出；没有伪造帮助文本。\n',encoding='utf-8')
    save_json(out/'discovery.json',{k:v for k,v in report.items() if k!='flattened'})
    return report


def arc_polyline(start,end,values,clockwise,absolute_center=False,chord_error_mm=.02):
    """G17 XY I/J arc; approximation is for evidence rendering, never slicing.

    Marlin G2/G3 reference: https://marlinfw.org/docs/gcode/G002-G003.html
    The actual Bambu output uses I/J centers, millimetres and the XY plane.
    """
    if 'R' in values or not ('I' in values or 'J' in values):
        raise ValueError('UNSUPPORTED_SUPPORT_ARC_CENTER: expected audited I/J format')
    center=np.array([values.get('I',0.),values.get('J',0.)])+(0 if absolute_center else start[:2])
    a=start[:2]-center;b=end[:2]-center;r0=np.linalg.norm(a);r1=np.linalg.norm(b)
    if min(r0,r1)<=0:raise ValueError('INVALID_SUPPORT_ARC_RADIUS')
    angle0=np.arctan2(a[1],a[0]);angle1=np.arctan2(b[1],b[0])
    sweep=(angle1-angle0)%(2*np.pi)
    if clockwise:sweep=-((angle0-angle1)%(2*np.pi))
    if abs(sweep)<1e-10:sweep=-2*np.pi if clockwise else 2*np.pi
    max_step=2*np.arccos(np.clip(1-chord_error_mm/max(r0,r1),-1,1))
    count=max(1,int(np.ceil(abs(sweep)/max(max_step,1e-5))))
    u=np.linspace(0,1,count+1);angle=angle0+u*sweep;r=r0+(r1-r0)*u
    points=np.column_stack([center[0]+r*np.cos(angle),center[1]+r*np.sin(angle),start[2]+u*(end[2]-start[2])])
    points[0]=start;points[-1]=end
    return points


def support_paths_from_archive(path):
    """Actual G-code extrusion segments, NOT a reconstructed support solid."""
    segments=[];tags=set();gcodes=[];support_arcs=0;support_lines=0;support_e=0.;mesh_objects=[];part_types=[]
    with zipfile.ZipFile(path) as archive:
        names=archive.namelist()
        for name in names:
            if name.endswith('.model'):
                tree=ET.fromstring(archive.read(name))
                for obj in tree.findall('.//{*}object'):
                    mesh=obj.find('{*}mesh')
                    if mesh is not None:mesh_objects.append(dict(file=name,id=obj.attrib.get('id'),name=obj.attrib.get('name'),
                        triangle_count=len(mesh.findall('.//{*}triangle'))))
            elif name=='Metadata/model_settings.config':
                tree=ET.fromstring(archive.read(name));part_types=[p.attrib.get('subtype') for p in tree.findall('.//part')]
        for name in names:
            if not name.endswith('.gcode'):continue
            gcodes.append(name);position=np.zeros(3);e_previous=0.;relative_e=True;absolute_xyz=True;feature='';plane='G17';absolute_center=False
            for raw in archive.read(name).decode(errors='replace').splitlines():
                if raw.startswith('; FEATURE:'):feature=raw.split(':',1)[1].strip();tags.add(feature)
                line=raw.split(';',1)[0].strip()
                if line=='M83':relative_e=True
                elif line=='M82':relative_e=False
                elif line=='G90':absolute_xyz=True
                elif line=='G91':absolute_xyz=False
                elif line=='G90.1':absolute_center=True
                elif line=='G91.1':absolute_center=False
                elif line in ('G17','G18','G19'):plane=line
                match=re.match(r'^(G0|G1|G2|G3|G92)\s+(.*)',line)
                if not match:continue
                values={k:float(v) for k,v in re.findall(r'([XYZEIJR])([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)',match.group(2))}
                if match.group(1)=='G92':
                    if 'E' in values:e_previous=values['E']
                    for k,axis in enumerate('XYZ'):
                        if axis in values:position[k]=values[axis]
                    continue
                new=position.copy()
                for k,axis in enumerate('XYZ'):
                    if axis in values:new[k]=values[axis]+(0 if absolute_xyz else position[k])
                extrusion=values.get('E',0.) if relative_e else values.get('E',e_previous)-e_previous
                if 'E' in values:e_previous=values['E'] if not relative_e else e_previous+values['E']
                if extrusion>0 and 'support' in feature.lower():
                    if match.group(1) in ('G2','G3'):
                        if plane!='G17':raise ValueError('UNSUPPORTED_SUPPORT_ARC_PLANE: '+plane)
                        points=arc_polyline(position,new,values,match.group(1)=='G2',absolute_center)
                        segments.extend(np.stack([points[:-1],points[1:]],axis=1));support_arcs+=1;support_e+=extrusion
                    elif np.linalg.norm(new-position)>1e-8:
                        segments.append(np.vstack([position,new]));support_lines+=1;support_e+=extrusion
                position=new
    segments=np.array(segments).reshape(-1,2,3)
    return dict(archive_members=names,gcode_files=gcodes,support_feature_tags=sorted(t for t in tags if 'support' in t.lower()),
        support_segment_count=len(segments),support_toolpath_length_mm=float(np.linalg.norm(np.diff(segments,axis=1)[:,0],axis=1).sum()) if len(segments) else 0.,
        support_linear_move_count=support_lines,support_arc_move_count=support_arcs,
        support_extrusion_filament_mm=support_e,arc_preview_chord_error_mm=.02,
        mesh_objects=mesh_objects,model_part_types=part_types,
        support_geometry_available=False,support_geometry_status='BAMBU_SUPPORT_GEOMETRY_NOT_AVAILABLE',
        support_evidence='G-code paths available; no separable support solid in the exported archive'),segments


def slice_top_five(discovery,cfg,top,out,logger):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if discovery['status']!='BAMBU_PROFILE_CONFIRMED':
        return dict(status=discovery['status'],results=[],successful_count=0,
            support_geometry_available=False,support_removal_status='GEOMETRIC_SUPPORT_ACCESS_ONLY')
    if len(top)>5:raise ValueError('More than Top 5 orientations submitted')
    config=dict(printer=dict(filament=cfg['printer']['filament'],nozzle_diameter_mm=discovery['nozzle_diameter_mm']),
        slicer=dict(enable_support=False,bed_type=cfg['bambu']['bed_type'],timeout_seconds=cfg['bambu']['timeout_seconds']))
    results=[]
    for item in top:
        for enabled in (False,True):
            if enabled and item['downward_overhang_area_mm2']<=cfg['bambu']['support_on_overhang_area_mm2'] and any(
                r['orientation_id']==item['candidate_id'] and r['status']=='BAMBU_SLICED' for r in results):continue
            mode='ON' if enabled else 'OFF';config['slicer']['enable_support']=enabled
            folder=out/f"orientation_{item['rank']:02}_support_{mode}"
            logger.info('Bambu local slice rank=%s support=%s',item['rank'],mode)
            # Use the existing audited local command builder, with an explicit
            # allowlist check BEFORE its subprocess boundary.
            original=existing.subprocess.run
            def guarded(command,*args,**kwargs):
                if isinstance(command,(list,tuple)) and command and Path(command[0]).name.lower().startswith('bambu-studio'):
                    check_local_command(command)
                    destination=command[command.index('--export-3mf')+1]
                    if str(command[0]).lower().endswith('.exe'):
                        local=existing.subprocess.check_output(['wslpath','-u',destination],text=True).strip()
                        kwargs['cwd']=Path(local).parent
                return original(command,*args,**kwargs)
            existing.subprocess.run=guarded
            try:response=existing.slice_candidates(discovery,config,[item],folder)
            finally:existing.subprocess.run=original
            native_result=Path(response.get('staging',out))/'result.json'
            if native_result.is_file():shutil.copy2(native_result,folder/'native_result.json')
            for record in response['results']:
                if native_result.is_file():record['native_cli_result']=json.loads(native_result.read_text())
                record['support_mode']=mode
                if record.get('archive'):
                    audit,segments=support_paths_from_archive(record['archive']);record['support_archive_audit']=audit
                    np.save(folder/'support_toolpath_segments.npy',segments)
                    record['support_paths_file']=str(folder/'support_toolpath_segments.npy')
                    save_json(folder/'support_archive_audit.json',audit)
                results.append(record)
            save_json(out/'slice_attempts.json',results)
    on=[r for r in results if r['status']=='BAMBU_SLICED' and r['support_mode']=='ON']
    off=[r for r in results if r['status']=='BAMBU_SLICED' and r['support_mode']=='OFF']
    # Frame geometry ranks the Top 5. Prefer the first successful supported
    # orientation because overhangs are explicitly present in this geometry.
    selected=min(on,key=lambda r:r['rank']) if on else None
    return dict(status='BAMBU_SLICED' if selected else 'BAMBU_VALIDATION_PENDING',results=results,
        successful_count=len(on)+len(off),support_off_success_count=len(off),support_on_success_count=len(on),selected=selected,
        support_geometry_available=False,support_geometry_status='BAMBU_SUPPORT_GEOMETRY_NOT_AVAILABLE',
        support_removal_status='GEOMETRIC_SUPPORT_ACCESS_ONLY',printer_job_sent=False)
