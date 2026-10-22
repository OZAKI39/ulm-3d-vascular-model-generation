"""Fresh Linux Bambu slice of the actual two-part P1S vascular assembly."""
import sys,json,zipfile
from pathlib import Path
root=Path('/workspace/vascular_printing');sys.path.insert(0,str(root))
from vascular_processing.core_support_p1s import default_output,profiles
from vascular_processing.support_slicer_adapter import slice_local,production_geometry
name=sys.argv[1] if len(sys.argv)>1 else 'native_slice'
if Path(name).name!=name:raise ValueError('Output name must be a basename')
out=root/'migration/20261009'/name;out.mkdir(parents=True,exist_ok=True)
bundle=profiles(out/'profiles')
archived=default_output();summary=json.loads((archived/'regeneration_summary.json').read_text());selected=archived/'selected'
r=slice_local(bundle['engine'],bundle['actual'],selected/'BG001_RMCA_core_support_P1S_04_ABS_project.3mf',out/'bambu')
assert r['slicer_success'],r
with zipfile.ZipFile(r['archive']) as z:
 assert z.testzip() is None
 settings=json.loads(z.read('Metadata/project_settings.config'))
 assert settings['printer_model']=='Bambu Lab P1S'
 assert settings['nozzle_diameter']==['0.4'] and settings['filament_type']==['ABS']
 assert str(settings['enable_support'])=='0'
 assert len([n for n in z.namelist() if n.endswith('.gcode')])==1
import trimesh
combined=trimesh.util.concatenate([trimesh.load_mesh(selected/name) for name in ['BG001_RMCA_vascular_core_P1S_04_ABS.stl','BG001_RMCA_custom_support_P1S_04_ABS.stl']])
combined.export(out/'expected_assembly.stl')
check=production_geometry(r['archive'],out/'expected_assembly.stl');assert check['passed']
result=dict(status='NATIVE_LINUX_SLICE_VERIFIED',engine=bundle['engine'],slice=r,geometry=check,
    no_printer_upload_or_start=not r['printer_upload'] and not r['printer_start'],
    unchanged_manufacturing_acceptance=summary['status'],
    note='Environment validation only; Linux slice does not replace historical deliverables or waive support/removal failures.')
(out/'native_smoke_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ['engine','slice','geometry']},indent=2))
