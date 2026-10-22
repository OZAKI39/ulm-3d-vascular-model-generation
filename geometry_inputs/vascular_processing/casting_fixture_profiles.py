"""Installed X1C presets; isolated local slicing only, no GUI/profile mutation."""
import copy,json,zipfile
from pathlib import Path
import numpy as np
import trimesh
from . import bambu_manufacturing as b
from .support_inputs import save,digest,numeric,resolve_profile
from .support_slicer_adapter import discover_engine,engine_profiles,slice_local,warnings

def profiles(inputs,cfg,out):
    _,abs_discovery=resolve_profile(inputs,Path(out)/'ABS_reference')
    root=Path(abs_discovery['active_preset_source']).parent
    names=dict(machine=cfg['printer'],process=cfg['PLA_process'],filaments=[cfg['PLA_filament']])
    try:
        p=b.active_profiles(root,names,'PLA');p['binary']=abs_discovery['binary'];p['material_selection']='INSTALLED_PLA_TEMPLATE_SELECTED_FOR_BCF_NOT_GUI_ACTIVE'
        actual=engine_profiles(p,'Bambu',{'enable_support':'0'})
        save(Path(out)/'PLA_discovery.json',p)
    except (FileNotFoundError,ValueError,KeyError) as e:
        save(Path(out)/'PLA_discovery.json',dict(status='PLA_SLICE_VALIDATION_PENDING',error=str(e)));return abs_discovery,None,None
    engine=discover_engine('Bambu',p,Path(out))
    return abs_discovery,actual,engine

def bed_mesh(mesh,rotation=None,center=(128,128)):
    result=mesh.copy();m=np.eye(4)
    if rotation is not None:m[:3,:3]=rotation
    points=result.vertices@m[:3,:3].T;lo=points.min(0);hi=points.max(0);m[:3,3]=[center[0]-(lo[0]+hi[0])/2,center[1]-(lo[1]+hi[1])/2,-lo[2]];result.apply_transform(m)
    return result,m

def slice_part(mesh,engine,actual,out,rotation=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);bed,m=bed_mesh(mesh,rotation);path=out/'print_input.stl';bed.export(path);save(out/'print_transform.json',dict(transform_4x4=m,scale_factor=1.,input_mesh_hash=digest(path)))
    signature=dict(source=digest(path),profiles=actual,binary_hash=engine.get('binary_sha256'));save(out/'slice_inputs.json',signature)
    record=slice_local(engine,actual,path,out)
    if record.get('archive'):
        with zipfile.ZipFile(record['archive']) as z:
            cfg=json.loads(z.read('Metadata/project_settings.config'));native=0
            for name in z.namelist():
                if name.endswith('.gcode'):
                    feature=''
                    for ln in z.read(name).decode(errors='replace').splitlines():
                        if ln.startswith('; FEATURE:'):feature=ln.split(':',1)[1].lower()
                        if ln.startswith(('G0 ','G1 ','G2 ','G3 ')) and ' E' in ln and 'support' in feature:native+=1
            record['archive_crc_ok']=z.testzip() is None
        record['actual_profile']={k:cfg.get(k) for k in ['printer_model','printer_settings_id','nozzle_diameter','filament_type','filament_settings_id','print_settings_id','line_width','layer_height','enable_support']}
        record['native_support_moves']=native;record['auto_support_off']=str(cfg.get('enable_support')) in ('0','False','false')
        record['support_free_slice_passed']=bool(record['slicer_success'] and record['auto_support_off'] and native==0 and not record['warning_audit']['critical'] and cfg.get('filament_type')==['PLA'])
    else:record['support_free_slice_passed']=False
    record['physical_print_validated']=False;save(out/'slice_result.json',record);return record
