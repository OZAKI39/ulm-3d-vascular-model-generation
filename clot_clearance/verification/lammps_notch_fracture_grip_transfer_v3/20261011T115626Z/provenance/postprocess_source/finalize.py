"""Close the bounded task; verify immutability, provenance, outputs and links."""
from pathlib import Path
import sys,json,hashlib,shutil,subprocess,py_compile,datetime
from html.parser import HTMLParser
import numpy as np
from PIL import Image
from design import js,save

def main(base,figs):
 root=Path(__file__).resolve().parents[2];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();launch=js(base/'provenance/LAUNCH_SOURCE_IDENTITY.json')['sha256'];assert all(sha(root/k)==v for k,v in launch.items()),'Formal launch source/design changed'
 before=js(base/'provenance/REUSED_IDENTITY.json');allowed='pd_clot/lammps_fragmentation/macro_lab_v1.py';changed=[k for k,v in before.items() if sha(root/k)!=v];assert changed==[allowed]
 old=root/'verification/lammps_shear_boundary_and_notch_localization_v2/20261011T001536Z';manifest=js(old/'FILE_MANIFEST.json');subset={k:v for k,v in manifest['files'].items() if '/M5B/' in k};assert len(subset)>0;assert all(sha(root/k)==v['sha256'] for k,v in subset.items())
 archive=base/'provenance/postprocess_source';archive.mkdir(exist_ok=True);sources={}
 for p in sorted((root/'analysis/notch_grip_transfer_v3').glob('*.py')):
  py_compile.compile(str(p),doraise=True);shutil.copyfile(p,archive/p.name);sources[str(p.relative_to(root))]=sha(p)
 sources[allowed]=sha(root/allowed);save(base/'provenance/FINAL_SOURCE_IDENTITY.json',dict(sha256=sources,formal_launch_source_and_design_unchanged=True,launch_record='LAUNCH_SOURCE_IDENTITY.json',old_identity_matches=len(before)-1,old_M5B_manifest_files_verified=len(subset)))
 receipt=js(base/'provenance/WSL_DOCUMENT_UPDATE.json');assert receipt['state']=='ACTUALLY_UPDATED_AND_READ_BACK';assert sha(base/'provenance/clot_rupture.before.md')==receipt['before_sha256'];assert sha(base/'provenance/clot_rupture.after.md')==receipt['after_sha256']
 status=js(base/'FINAL_STATUS.json');assert status['new_formal_cases']==1 and status['bug_retests']==0 and status['MATERIAL_GC_STATUS']=='GC_UNAVAILABLE';assert status['IMPLEMENTATION_VALIDITY']=='PASS';assert js(base/'RUN_RESOURCE.json')['exit_code']==0;assert not list((base/'M5C').glob('*ERROR*'))
 status['WSL_document']=receipt;status['raw_server_root']=str(base);status['figures_server_root']=str(figs);status['delivery_manifest']='FILE_MANIFEST.json';save(base/'FINAL_STATUS.json',status)
 git=lambda *args:subprocess.check_output(['git',*args],cwd=root,text=True).strip();sup=subprocess.run(['supervisorctl','status','y_clot_grip_v3_20261011T115626Z'],capture_output=True,text=True);assert 'EXITED' in sup.stdout;assert git('rev-parse','HEAD')==js(base/'provenance/SERVER_START.json')['head'];save(base/'provenance/SERVER_END.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),git_head=git('rev-parse','HEAD'),git_branch=git('branch','--show-current'),git_status=git('status','--short'),own_supervisor=sup.stdout.strip(),all_preexisting_identities_unchanged_except_lab_adapter=True,checked=len(before)-1,free_disk_bytes=shutil.disk_usage(root).free))
 npz=0
 for p in base.rglob('*.npz'):
  with np.load(p,allow_pickle=False) as z:
   for k in z.files:z[k]
  npz+=1
 png={}
 for p in figs.glob('*.png'):
  with Image.open(p) as im:png[p.name]=list(im.size);assert min(im.size)>=900;im.verify()
 assert len(png)==6
 required=['REVIEW_LAMMPS_NOTCH_FRACTURE_GRIP_TRANSFER_V3.md','FINAL_STATUS.json','M5B_KINEMATIC_ROOT_CAUSE.md','M5C_DESIGN_FREEZE.json','GRIP_KINEMATIC_AUDIT.json','CRACK_PROPAGATION_AUDIT.json','CRACK_AREA_ATTEMPT.json','ENERGY_LEDGER_AUDIT.json','M5B_M5C_COMPARISON.csv','IMPLEMENTATION_CHECKS.json','OPEN_RESULTS.html'];assert all((base/f).is_file() for f in required)
 class Links(HTMLParser):
  def __init__(self):super().__init__();self.urls=[]
  def handle_starttag(self,tag,attrs):
   self.urls.extend(v for k,v in attrs if k in ['href','src'])
 missing=[];links=0
 for p in [base/'OPEN_RESULTS.html',figs/'OPEN_RESULTS.html']:
  parser=Links();parser.feed(p.read_text())
  for x in parser.urls:
   if x.startswith(('http:','https:','#')):continue
   links+=1
   if x!='FILE_MANIFEST.json' and not (p.parent/x).is_file():missing.append([str(p),x])
 assert not missing,missing
 for p in base.rglob('*.json'):js(p)
 check=dict(status='DELIVERY_VALIDATED_SCIENTIFIC_FAILURE_RETAINED',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),formal_solver_cases=1,zero_step_cases=1,bug_retests=0,old_immutable_identities_verified=len(before)-1,old_M5B_manifest_entries_verified=len(subset),formal_launch_source_and_design_unchanged=True,NPZ_archives_decoded=npz,PNG_dimensions=png,html_links_verified=links,missing_links=missing,WSL_document_backup_diff_and_hash_receipt_verified=True,scientific_gates='Separate FINAL_STATUS.json; delivery success does not upgrade failed notch/path gates.')
 save(base/'DELIVERY_VALIDATION.json',check)
 excluded={'FILE_MANIFEST.json','DELIVERY_VALIDATION.json','COMPACT_DOWNLOAD.json','LOCAL_FINAL_VALIDATION.json','EXPLICIT_SYNC.json','finalize.py.stdout.log','finalize.py.stderr.log'};files={}
 for directory in [base,figs]:
  for p in sorted(directory.rglob('*')):
   if not p.is_file() or p.name in excluded or '__pycache__' in p.parts:continue
   heavy=p.name in ['body.data','history.csv','native_energy_density.dat','boundary_dynamics_audit.csv'] or p.parent.name=='M5C' and p.name.startswith('frame_')
   files[str(p.relative_to(root))]=dict(bytes=p.stat().st_size,sha256=sha(p),WSL_selection='SERVER_ONLY_FULL_NATIVE_DATA' if heavy else 'COMPACT_REVIEW')
 result=dict(server_root=str(root),run_id=base.name,files=files,total_bytes=sum(v['bytes'] for v in files.values()),compact_bytes=sum(v['bytes'] for v in files.values() if v['WSL_selection']=='COMPACT_REVIEW'),scope='Only this V3 raw/analysis/provenance/figures; manifest, transfer-validation receipts and finalizer stdout/stderr excluded to avoid self-changing hashes.',excluded_names=sorted(excluded))
 save(base/'FILE_MANIFEST.json',result);print(json.dumps(dict(checks=check,files=len(files),total_bytes=result['total_bytes'],compact_bytes=result['compact_bytes']),indent=2))

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
