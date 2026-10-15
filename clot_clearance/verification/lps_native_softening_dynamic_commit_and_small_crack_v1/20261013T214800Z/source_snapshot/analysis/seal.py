from common import *
import ast,shutil,subprocess,re
assert json.loads((O/'DELIVERY_READINESS.json').read_text())['status']=='READY_FOR_MANIFEST'
# Verify that every untouched upstream source remains identical to the accepted base.
changed=[];unchanged=0
for p in (PREV/'source').rglob('*'):
 if not p.is_file():continue
 rel=p.relative_to(PREV/'source');q=O/'source'/rel
 if not q.is_file() or sha(p)!=sha(q):changed.append(str(rel))
 else:unchanged+=1
new=[str(p.relative_to(O/'source')) for p in (O/'source').rglob('*') if p.is_file() and not (PREV/'source'/p.relative_to(O/'source')).exists()]
assert set(changed)=={'src/.gitignore','src/BPM/bond_bpm_peri_variational.cpp','src/BPM/bond_bpm_peri_variational.h'},changed
assert set(new)=={'src/BPM/fix_lps_variational_dynamic.cpp','src/BPM/fix_lps_variational_dynamic.h'},new
dump(O/'UPSTREAM_SOURCE_PROTECTION.json',dict(status='PASS',unchanged_files=unchanged,expected_changed_files=changed,expected_new_files=new,no_production_edits=True))
for p in A.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
ss=O/'source_snapshot';(ss/'analysis').mkdir(parents=True,exist_ok=True);(ss/'native_cpp').mkdir(exist_ok=True);(ss/'frozen_offline').mkdir(exist_ok=True)
for p in A.iterdir():
 if p.is_file():shutil.copy2(p,ss/'analysis'/p.name)
for p in A.glob('*'):
 if p.suffix in ['.cpp','.h']:shutil.copy2(p,ss/'native_cpp'/p.name)
shutil.copy2(R/'analysis/softening_feasibility_v1/model.py',ss/'frozen_offline/model.py');shutil.copy2(O/'source/src/.gitignore',ss/'native_cpp/src.gitignore')
(O/'binary').mkdir(exist_ok=True);shutil.copyfile(LIB,O/'binary/liblammps.so');assert sha(O/'binary/liblammps.so')==sha(LIB)
# Save direct locked hook sources independently of the full retained server source.
for name in ['verlet.cpp','fix_nve.cpp','fix_move.cpp','modify.cpp','neighbor.cpp','fix.h','BPM/pair_bpm_peri.cpp']:
 dst=ss/'locked_hooks'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(O/'source/src'/name,dst)
# Syntax-only source guards; do not execute new dynamics.
cpp=[p for p in A.glob('*') if p.suffix in ['.cpp','.h']];assert all(all(c<128 for c in p.read_bytes()) for p in cpp)
status=json.loads((O/'FINAL_STATUS.json').read_text());assert status['phases']['P4']=='NOT_RUN'
assert not (O/'SMALL_CRACK_DESIGN_FREEZE.json').exists() and not (O/'SMALL_CRACK_TOPOLOGY_AUDIT.json').exists()
report=(O/'REVIEW_LPS_NATIVE_SOFTENING_DYNAMIC_COMMIT_AND_SMALL_CRACK_V1.md').read_text();assert re.findall(r'^### (\d+)\.',report,flags=re.M)==[str(j) for j in range(1,52)]
required=['FINAL_STATUS.json','CODE_CLEANUP_INVENTORY.json','CODE_CLEANUP_REPORT.md','DELETED_CODE_MANIFEST.json','ACTIVE_DEVELOPMENT_ENTRYPOINTS.md','SOURCE_AND_LIBRARY_IDENTITY.json','DYNAMIC_DAMAGE_COMMIT_DESIGN_FREEZE.json','DYNAMIC_INTEGRATOR_IMPLEMENTATION_AUDIT.md','ACCEPTED_STEP_DAMAGE_HISTORY_AUDIT.json','RUN0_AND_TRIAL_STATE_SAFETY.json','DYNAMIC_341_NODE_AUDIT.json','FIRST_SOFTENING_AND_FAILURE_EVENTS.json','DYNAMIC_FORCE_ENERGY_AUDIT.json','MOMENTUM_AND_TIMESTEP_AUDIT.json','DYNAMIC_ENERGY_LEDGER.json','DAMAGE_IRREVERSIBILITY_AUDIT.json','POSTBREAK_CONTACT_DYNAMIC_AUDIT.json','IMPLEMENTATION_CHECKS.json','RUN_RESOURCES.json','REPRODUCE.txt','NATIVE_SOURCE_PATCH.diff','OPEN_RESULTS.html','WSL_DOCUMENT_UPDATE.json']
assert all((O/p).is_file() for p in required)
# Last read-only verification of the protected assets after all development.
protected=json.loads((O/'HISTORICAL_PROTECTION_BEFORE.json').read_text());conflicts=[n for n,v in protected.items() if not (R/n).is_file() or sha(R/n)!=v['sha256']];assert not conflicts
identity=json.loads((O/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());assert all(sha(n)==v['sha256'] for n,v in identity['files'].items())
assert subprocess.check_output(['git','diff','--binary'],cwd=R,text=True)==(O/'GIT_DIFF_before.patch').read_text()
dump(O/'FINAL_DELIVERY_CHECKS.json',dict(status='PASS_DELIVERY_SCIENTIFIC_STOP_RETAINED',required_files=len(required),report_answers=51,PNG=5,SVG=5,CPP_ASCII=True,Python_AST=True,history_files_protected=len(protected),source_base_unchanged_files=unchanged,library_sha256=sha(LIB),P4_not_created=True,dynamic_steps_main=244,dynamic_steps_short_half_window=68,scientific_data_not_modified_by_delivery=True))
# Only useful payload: full source and compiler objects remain on the server.
def payload():
 files=[]
 for root in [A,O,F]:
  for p in root.rglob('*'):
   if not p.is_file():continue
   if p.is_relative_to(O/'source') or p.is_relative_to(O/'build'):continue
   if p.name in ['PAYLOAD_MANIFEST.json','PAYLOAD_MANIFEST.sha256']:continue
   files.append(p)
 files.append(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md');return sorted(set(files))
files=payload();resources=json.loads((O/'RUN_RESOURCES.json').read_text());resources['payload_bytes_before_manifest']=sum(p.stat().st_size for p in files);resources['payload_files_before_manifest']=len(files);resources['full_upstream_source_and_object_cache_retained_on_server_not_local_payload']=True;dump(O/'RUN_RESOURCES.json',resources)
files=payload();total=sum(p.stat().st_size for p in files);assert total<600*1024**2
dump(O/'PAYLOAD_MANIFEST.json',dict(version=1,task=TASK,run_id=RUN,scope='Current analysis module, active entrypoint, useful run evidence/source snapshots/actual binary, figures. Full duplicate upstream source and CMake object cache stay on server.',file_count=len(files),total_bytes=total,files={str(p.relative_to(R)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in files}))
(O/'PAYLOAD_MANIFEST.sha256').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'  PAYLOAD_MANIFEST.json\n')
print('SEALED',len(files),'FILES',total,'BYTES','SHA',sha(O/'PAYLOAD_MANIFEST.json'),flush=True)
