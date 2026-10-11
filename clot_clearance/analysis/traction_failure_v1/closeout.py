"""Read-only baseline verification and manifests for this isolated run."""
from analysis.traction_failure_v1.common import *
import subprocess,sys,importlib.metadata as md,difflib,platform
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def load(p):return json.loads(Path(p).read_text())
def entry(p):return dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=digest(p))
def verify():
 baseline=load(COMMON/'provenance/PROTECTED_START.json');changed=[]
 for name,expected in baseline.items():
  p=Path(name);actual=digest(p) if p.is_file() else None
  if actual!=expected:changed.append(dict(path=name,before=expected,after=actual))
 after={x.metadata['Name']:x.version for x in md.distributions()};before=load(COMMON/'provenance/ENVIRONMENT_BEFORE.json')
 head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
 envsame=before==after;headsame=head=='cc5ee77149f55e58d715bb5a4660c6467231a826'
 save(COMMON/'provenance/ENVIRONMENT_AFTER.json',after)
 save(COMMON/'PROTECTION_AUDIT.json',dict(checked_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),registered_files=len(baseline),changed=changed,all_registered_sha_unchanged=not changed,environment_unchanged=envsame,git_HEAD=head,git_HEAD_unchanged=headsame,baseline_manifest_sha256=digest(COMMON/'provenance/PROTECTED_START.json'),production_edits=False,git_mutations=False))
 assert not changed and envsame and headsame
 st=load(A/'FINAL_STATUS.json');st['states']['ORIGINAL_LPS_BASELINE_PROTECTED']=f'PASS_SHA_{len(baseline)}_REGISTERED_FILES';save(A/'FINAL_STATUS.json',st)
 co=load(COMMON/'FINAL_STATUS.json');co['A']=st['states'];co['protection']='PASS';save(COMMON/'FINAL_STATUS.json',co)
 inputs={ROOT/k for k in load(B/'INPUT_SHA.json')}
 inputs.update([M5C/'geometry.npz',PREV/'refined/SV_MESH/mesh-complete.mesh.vtu',PREV/'refined/CFD_RESULT.json',ROOT/'analysis/material_macro_audit_v1/energy.cpp',ROOT/'analysis/material_macro_audit_v1/energy.so',ROOT/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so'])
 for n in ['M1','M2']:inputs.add(ROOT/'verification/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z'/n/'history.csv')
 for name in baseline:
  p=Path(name)
  if p.name in ['bond_bpm_spring.cpp','pair_bpm_spring.cpp','pair_y_clot_wall.cpp'] or ('src/' in name and 'bpm' in p.name):inputs.add(p)
 save(COMMON/'INPUT_MANIFEST.json',[entry(p) for p in sorted(inputs)])
 code=[p for p in (ROOT/'analysis/traction_failure_v1').rglob('*') if p.is_file() and '__pycache__' not in p.parts]
 sources=[p for p in (A/'sources').rglob('*') if p.is_file() and '__pycache__' not in p.parts]
 save(COMMON/'SOURCE_MANIFEST.json',dict(new_code=[entry(p) for p in sorted(code)],literature_and_pinned_source=[entry(p) for p in sorted(sources)],source_input_audit='INPUT_MANIFEST.json; complete preexisting protection ledger in provenance/PROTECTED_START.json',python=platform.python_version()))
 patch=[]
 for p in sorted(code):
  if p.suffix in ['.py','.json','.txt']:
   patch.extend(difflib.unified_diff([],p.read_text().splitlines(True),fromfile='/dev/null',tofile=str(p.relative_to(ROOT))))
 (COMMON/'CODE_ADDITIONS.diff').write_text(''.join(patch))
 (COMMON/'provenance/GIT_STATUS_AFTER.txt').write_text(subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True))
 jobs={}
 for p in COMMON.glob('*_LAUNCH.json'):
  name=load(p)['name'];r=subprocess.run(['supervisorctl','status',name],capture_output=True,text=True);jobs[name]=r.stdout.strip()
 assert all('RUNNING' not in v for v in jobs.values())
 save(COMMON/'OWN_JOBS_FINAL.json',dict(jobs=jobs,automatic_restart=False,original_figures_failure='Matplotlib Axes.set did not accept linthresh; plotting-only fix in figures_fixed. Numerical analyses not rerun.',other_jobs_modified=False))
 print('protected',len(baseline),'files; environment and HEAD unchanged')
def manifest():
 assert load(COMMON/'PROTECTION_AUDIT.json')['all_registered_sha_unchanged']
 assert load(COMMON/'document_edit/AUDIT.json')['readback_match']
 excluded={'SERVER_ARTIFACT_MANIFEST.json','COMPACT_MANIFEST.json','LOCAL_RETRIEVAL_AUDIT.json'}
 files=[]
 for folder in [A,B,COMMON,FA,FB,ROOT/'analysis/traction_failure_v1']:
  for p in folder.rglob('*'):
   if p.is_file() and '__pycache__' not in p.parts and p.name not in excluded:files.append(entry(p))
 files.sort(key=lambda d:d['path']);base=dict(run=RUN,files=files,count=len(files),total_bytes=sum(x['bytes'] for x in files),scope='Only current-run code, three verification trees and seven PNG/SVG groups. No historical raw CFD or production dependencies. Self-referential manifests/retrieval receipts excluded from their own inventory.')
 save(COMMON/'SERVER_ARTIFACT_MANIFEST.json',base);save(COMMON/'COMPACT_MANIFEST.json',base)
 print('manifest',len(files),'files',base['total_bytes'],'bytes')
if __name__=='__main__':
 {'verify':verify,'manifest':manifest}[sys.argv[1]]()
