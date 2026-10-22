"""Clean only this deployment's redundant downloads and summarize real checks."""
from pathlib import Path
import json,hashlib,subprocess,sys,xml.etree.ElementTree as ET,shutil
root=Path('/workspace/vascular_printing');out=root/'migration/20261009'
manifest=json.loads((out/'project_files.json').read_text());removed=[]
# Initial resumable transfer started with a wider list. Remove only those
# explicitly narrowed temporary/cache classes, never unknown user files.
for base in [root/'references/image',root/'outputs']:
 if not base.exists():continue
 for p in sorted(base.rglob('*'),key=lambda x:len(x.parts),reverse=True):
  if not p.is_file() or p.is_symlink():continue
  rel=str(p.relative_to(root))
  temporary=rel.startswith('references/image/') or p.name in {'support_toolpath_segments.npy','model_toolpath_segments.npy','travel_segments.npy','extrusion_segments.npy'} or any(x.startswith('s1_11_') for x in p.parts) or 'staging' in p.parts
  if temporary and rel not in manifest:
   removed.append(dict(path=rel,bytes=p.stat().st_size,reason='Explicitly omitted regenerated cache/duplicate'));p.unlink()
installer=root/'tools/bambu-studio-2.7.1.57/BambuStudio.AppImage'
expected='85b053853f1a238775cd7ae2d4d19544e35260592428d2276ca27dcbd5f8eb53'
if installer.exists():
 assert hashlib.sha256(installer.read_bytes()).hexdigest()==expected
 removed.append(dict(path=str(installer.relative_to(root)),bytes=installer.stat().st_size,reason='Verified installer; extracted native runtime retained'));installer.unlink()
env=dict(__import__('os').environ,UV_CACHE_DIR=str(root/'tools/runtime/uv-cache'))
clean=subprocess.run(['uv','cache','clean'],env=env,capture_output=True,text=True,check=True)
(out/'cache_cleanup.log').write_text(clean.stdout+clean.stderr)
(out/'removed_deployment_caches.json').write_text(json.dumps(removed,indent=2))
subprocess.run([str(root/'.venv/bin/python'),'-m','pip','check'],check=True,stdout=(out/'pip_check.txt').open('w'))
subprocess.run([str(root/'.venv/bin/python'),'-m','pip','freeze'],check=True,stdout=(out/'pip-freeze.txt').open('w'))
tests=ET.parse(out/'pytest.xml').getroot().find('testsuite')
failures=[x.attrib['name'] for x in tests.findall('testcase') if x.find('failure') is not None or x.find('error') is not None]
assert failures==['test_practical_delivery_integrity_and_no_false_production_freeze'],failures
gui=list((out/'gui_smoke').rglob('*ui_compatibility.json'));assert len(gui)==1
gui_result=json.loads(gui[0].read_text());assert gui_result['gui_smoke_passed'],gui_result
native=json.loads((out/'native_slice_linux_verified/native_smoke_summary.json').read_text())
assert native['slice']['slicer_success'] and native['geometry']['passed']
native_path=Path(native['engine']['native_binary']);assert hashlib.sha256(native_path.read_bytes()).hexdigest()==native['engine']['binary_sha256']
verification=json.loads((out/'transfer_verification.json').read_text());assert verification['status']=='VERIFIED'
env_check=json.loads((out/'environment_smoke.json').read_text());assert env_check['python']['isolated'] and env_check['boolean_and_collision']
topbrain=list((out/'topbrain_smoke').rglob('case_manifest.json'));assert len(topbrain)==1
case=json.loads(topbrain[0].read_text());assert case['status']=='PASS',case['status']
report=dict(status='MIGRATED_AND_ENVIRONMENT_VERIFIED_WITH_PREEXISTING_MANUFACTURING_FAILURES',client_date='2026-10-09',
 host='95.3.33.46',ssh_port=15283,project=str(root),environment=str(root/'.venv'),python_version='3.13.11',
 file_verification=verification,environment_smoke=env_check,
 tests=dict(total=int(tests.attrib['tests']),passed=int(tests.attrib['tests'])-len(failures),failed=len(failures),preexisting_failures=failures,new_failures=0),
 gui_smoke=dict(passed=True,report=str(gui[0]),static_view_count=6),topbrain_smoke=dict(passed=True,case_id='001',modality='mr',case_manifest=str(topbrain[0])),
 native_slice=dict(passed=True,version=native['engine']['version'],native_binary=str(native_path),native_binary_sha256=native['engine']['binary_sha256'],
                  report=str(out/'native_slice_linux_verified/native_smoke_summary.json'),geometry_preserved=True,printer='Bambu Lab P1S',nozzle_mm=.4,material='ABS'),
 manufacturing_status='NEEDS_ADJUSTMENT',workspace_is_volume=False,no_printer_upload_or_start=True,
 original_project_deleted=False,original_full_backup_untouched=True,base_server_venv_used=False,
 physical_project_data_bytes=verification['source_bytes'],removed_cache_files=len(removed),removed_cache_bytes=sum(x['bytes'] for x in removed))
(out/'migration_result.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in report.items() if k not in ['environment_smoke']},ensure_ascii=False,indent=2))
