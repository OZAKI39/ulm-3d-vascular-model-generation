"""Close and hash this bounded run; validate deliverables without invoking a solver."""
from pathlib import Path
import sys,json,hashlib,shutil,time,re,py_compile
from html.parser import HTMLParser
import numpy as np
from PIL import Image
base=Path(sys.argv[1]).resolve();figs=Path(sys.argv[2]).resolve();root=Path.cwd();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();J=lambda p:json.loads(p.read_text())
old=root/'verification/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z';assert sha(old/'FILE_MANIFEST.json')=='609fd7997a857bf30d72fd1b1ce0fbc09b9c539424e9a3d41fbe018bacf14323'
prior=J(old/'FILE_MANIFEST.json');subset={k:v for k,v in prior['files'].items() if k.startswith(('M2/','M5/'))};assert all(sha(old/k)==v['sha256'] for k,v in subset.items())
# Match the formal solver source and the version embedded in the launch archive.
launch=J(base/'provenance/LAUNCH_SOURCE_SHA256.json');assert sha(root/'pd_clot/lammps_fragmentation/macro_lab_v1.py')==launch['pd_clot/lammps_fragmentation/macro_lab_v1.py'];assert sha(root/'analysis/shear_boundary_notch_v2/run.py')==launch['analysis/shear_boundary_notch_v2/run.py']
postdir=base/'provenance/postprocess_source';postdir.mkdir(exist_ok=True);sources={}
for p in sorted((root/'analysis/shear_boundary_notch_v2').glob('*.py')):
 py_compile.compile(str(p),doraise=True);shutil.copyfile(p,postdir/p.name);sources[str(p.relative_to(root))]=sha(p)
sources['pd_clot/lammps_fragmentation/macro_lab_v1.py']=sha(root/'pd_clot/lammps_fragmentation/macro_lab_v1.py');(base/'provenance/FINAL_SOURCE_IDENTITY.json').write_text(json.dumps(dict(sha256=sources,formal_sources='LAUNCH_SOURCE_SHA256.json',scope='Formal macro/run unchanged since launch; later audit metadata/figures/report-only code snapshots stored here.'),indent=2))
# Keep report input capture consistent with final metadata.
inputs=J(base/'REPORT_INPUTS.json');inputs['shear']=J(base/'SHEAR_BOUNDARY_AUDIT.json');(base/'REPORT_INPUTS.json').write_text(json.dumps(inputs,indent=2))
execution=J(base/'CAMPAIGN_EXECUTION.json');assert len(execution)==2 and all(v['returncode']==0 for v in execution);assert {v['case'] for v in execution}=={'M2P','M5B'}
status=J(base/'FINAL_STATUS.json');assert status['new_nontrivial_LAMMPS_cases']==2 and status['FRACTURE_WORK_IDENTIFIABILITY']['Gc_J_m2'] is None
receipt=J(base/'provenance/WSL_DOCUMENT_UPDATE.json');assert sha(base/'provenance/clot_rupture.before.md')==receipt['before_sha256'];assert sha(base/'provenance/clot_rupture.after.md')==receipt['after_sha256'];assert receipt['state']=='ACTUALLY_UPDATED_AND_READ_BACK'
# Decode, rather than only checking extensions, of exact state and figure-source archives.
npz_count=0
for p in base.rglob('*.npz'):
 with np.load(p,allow_pickle=False) as z:
  for k in z.files:z[k]
 npz_count+=1
pngs={}
for p in figs.glob('*.png'):
 with Image.open(p) as im:pngs[p.name]=im.size;assert min(im.size)>=900;im.verify()
assert len(pngs)==9
class Links(HTMLParser):
 def __init__(self):super().__init__();self.urls=[]
 def handle_starttag(self,tag,attrs):
  for k,v in attrs:
   if k in ['href','src']:self.urls.append(v)
missing=[];link_count=0
for f in [base/'OPEN_RESULTS.html',figs/'OPEN_RESULTS.html']:
 p=Links();p.feed(f.read_text())
 for x in p.urls:
  if x.startswith(('http:','https:','#')):continue
  link_count+=1
  if x!='FILE_MANIFEST.json' and not (f.parent/x).is_file():missing.append([str(f),x])
assert not missing,missing
for f in base.rglob('*.json'):J(f)
check=dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),status='DELIVERY_CHECKS_PASSED_WITH_SCIENTIFIC_LIMITATIONS_RETAINED',new_solver_cases=2,zero_step_preflights=2,old_manifest_matches_preexisting_WSL_copy=True,old_M2_M5_artifact_hashes_verified=len(subset),launch_solver_source_unchanged=True,python_compile=True,npz_archives_decoded=npz_count,png_dimensions=pngs,html_local_links=link_count,missing_links=missing,WSL_document_receipt_verified=True,formal_execution=execution,viewed_figures=['A1_shear_definitions_and_profiles','B1_notch_vs_grip_first_break','B3_notch_and_grip_projection','C2_actual_fracture_states','D2_crossing_connectivity_limits'],scientific_gates='FINAL_STATUS.json; independent of execution/delivery checks')
(base/'DELIVERY_VALIDATION.json').write_text(json.dumps(check,indent=2))
# SHA manifest of accepted artifacts. Exclude self and transfer receipts to avoid recursive hashes.
excluded={'FILE_MANIFEST.json','DELIVERY_VALIDATION.json','COMPACT_DOWNLOAD.json','LOCAL_FINAL_VALIDATION.json'};files={}
for directory in [base,figs]:
 for p in sorted(directory.rglob('*')):
  if not p.is_file() or p.name in excluded or '__pycache__' in p.parts:continue
  files[str(p.relative_to(root))]=dict(bytes=p.stat().st_size,sha256=sha(p),WSL_selection='SERVER_ONLY_FULL_NATIVE_DATA' if p.name in ['history.csv','body.data','native_energy_density.dat'] else 'COMPACT_REVIEW')
manifest=dict(server_project_root=str(root),run_id=base.name,utc=check['utc'],scope='Closed task verification and figure artifacts; excludes manifest itself, validation/export receipts and Python cache',files=files,total_bytes=sum(v['bytes'] for v in files.values()),compact_bytes=sum(v['bytes'] for v in files.values() if v['WSL_selection']=='COMPACT_REVIEW'))
(base/'FILE_MANIFEST.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(dict(checks=check,files=len(files),total_bytes=manifest['total_bytes'],compact_bytes=manifest['compact_bytes']),indent=2))
