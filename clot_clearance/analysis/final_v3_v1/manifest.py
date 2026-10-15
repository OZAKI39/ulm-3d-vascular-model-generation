from common import *
import json,re,subprocess,html,time
from PIL import Image
from urllib.parse import unquote,urlsplit
required=['REVIEW_LPS_FINAL_NATIVE_BUGFIX_FULL_V3_REGRESSION_V1.md','FINAL_STATUS.json','DESIGN_FREEZE.json','SOURCE_AND_LIBRARY_IDENTITY.json','PATCH_PROVENANCE.json','V3_INPUT_IDENTITY_AUDIT.json','WVOLUME_RUNTIME_IDENTITY_AUDIT.json','MPI_SORT_AND_MIGRATION_EVENTS.json','CROSS_RANK_BREAK_AUDIT.json','NATIVE_FORCE_AND_TOPOLOGY_AUDIT.json','EARLY_FLOW_COUPLING_AUDIT.json','V3_OLD_NEW_COMPARISON.csv','FRAGMENT_ATTACHMENT_AWARE_V2.json','IMPLEMENTATION_CHECKS.json','OPEN_RESULTS.html','REPRODUCE.txt','DOCUMENT_UPDATE_AUDIT.json']
assert all((O/n).is_file() for n in required);assert not (O/'REGRESSION_FAILURE_DIAGNOSTIC.json').exists()
status=json.loads((O/'FINAL_STATUS.json').read_text());assert status['WSL_DOCUMENT_UPDATE']=='PASS_ACTUAL_WSL_UPDATED' and not any('PENDING' in v for v in status.values())
run=json.loads((O/'RUN_COMMAND.json').read_text());core=['common.py','driver.py','launch.py','runtime_audit.py','observer.cpp'];coreproof={}
for name in core:
 p=R/'analysis/final_v3_v1'/name;actual=sha(p);old=run['code_sha256'][str(p)];assert actual==old;coreproof[name]={'sha256':actual,'matches_executed_source':True}
assert sha(O/'observer.so')==run['observer_sha256'];assert sha(B/'liblammps.so')==EXPECTED
protected=json.loads((O/'PROTECTED_START.json').read_text());bad=[p for p,q in protected.items() if not Path(p).is_file() or sha(p)!=q['sha256']];assert not bad
pngs={}
for p in F.glob('*.png'):
 with Image.open(p) as im:im.verify()
 with Image.open(p) as im:pngs[p.name]={'width':im.width,'height':im.height};assert im.width>=3000 and im.height>=2000
assert len(pngs)==4
links=[]
for path in [O/'OPEN_RESULTS.html',F/'OPEN_RESULTS.html']:
 for target in re.findall(r'(?:href|src)="([^"]+)"',path.read_text()):
  if urlsplit(target).scheme or target.startswith('#'):continue
  p=(path.parent/unquote(target)).resolve();assert p.is_file(),p;links.append(str(p))
end={'status':'PASS','required_deliverables':required,'executed_core_code_identity':coreproof,'observer_binary_sha256':sha(O/'observer.so'),'final_native_library_sha256':sha(B/'liblammps.so'),'protected_files':len(protected),'changed_protected_files':bad,'PNG_resolution':pngs,'static_index_links_verified':len(links),'production_promoted':False,'formal_runs':1,'native_driver_invocations':3}
dump(O/'FINAL_IMPLEMENTATION_AND_DELIVERY_CHECKS.json',end)
code="""import json,xmlrpc.client,supervisor.xmlrpc
r=xmlrpc.client.ServerProxy('http://127.0.0.1',transport=supervisor.xmlrpc.SupervisorTransport(None,None,'unix:///var/run/supervisor.sock'))
print(json.dumps(r.supervisor.getProcessInfo('clot_final_v3_20261011T211309Z')))
""";res=subprocess.run(['/usr/bin/python3','-c',code],capture_output=True,text=True,check=True);job=json.loads(res.stdout);assert job['statename']=='EXITED' and job['exitstatus']==0;dump(O/'MANAGED_FINAL_STATUS.json',job)
# Archive the exact executed core and current postprocessing source without extra libraries.
arch=O/'reproduction_code';arch.mkdir(exist_ok=True)
import shutil
for p in (R/'analysis/final_v3_v1').glob('*'):
 if p.is_file() and p.suffix in ['.py','.cpp','.json']:shutil.copy2(p,arch/p.name)
excluded={'SERVER_OUTPUT_SHA256.json','DELIVERY_VERIFICATION.json','LOCAL_DELIVERY_SHA256.json'}
manifest={}
for root in [O,F]:
 for p in sorted(root.rglob('*')):
  if not p.is_file() or p.name in excluded or '__pycache__' in p.parts:continue
  manifest[str(p.relative_to(R))]={'sha256':sha(p),'bytes':p.stat().st_size}
dump(O/'SERVER_OUTPUT_SHA256.json',{'files':manifest,'file_count':len(manifest),'total_bytes':sum(q['bytes'] for q in manifest.values()),'exclusions':sorted(excluded),'sha256_scope':'complete new verification/results payload including raw native trajectory/body/restart, excluding recursive manifest/transfer audit'});print(json.dumps({'status':'PASS','files':len(manifest),'bytes':sum(q['bytes'] for q in manifest.values()),'source_matches_executed':True,'history_unchanged':True,'managed_job':job['statename']}))
