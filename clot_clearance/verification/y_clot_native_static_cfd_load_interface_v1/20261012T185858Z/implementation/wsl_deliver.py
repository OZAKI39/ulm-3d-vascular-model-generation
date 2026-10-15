from pathlib import Path
import json,hashlib,subprocess,shlex,sys,re
from remote import C,L,OUT,REL,SSH
O=OUT;F=L/'results'/C['task']/C['run'];P=L/'analysis/static_cfd_v1';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
# Local required-output and relative HTML link check, no additional solve.
required=['REVIEW_Y_CLOT_NATIVE_STATIC_CFD_LOAD_INTERFACE_V1.md','FINAL_STATUS.json','STATIC_LOAD_DESIGN_FREEZE.json','SOURCE_AND_LIBRARY_IDENTITY.json','CFD_LOAD_SOURCE_AUDIT.json','NATIVE_FORCE_UNIT_AUDIT.json','ATOM_ID_LOAD_APPLICATION_AUDIT.json','UNLOADED_REFERENCE_AUDIT.json','STATIC_EQUILIBRIUM_AUDIT.json','ANCHOR_REACTION_AUDIT.json','FORCE_MOMENT_BALANCE.json','BULK_INTERFACE_STRETCH_AUDIT.json','INDEPENDENT_LPS_ENERGY_AUDIT.json','CLOT_STATIC_DISPLACEMENTS.npz','CLOT_STATIC_REACTIONS.npz','IMPLEMENTATION_CHECKS.json','OPEN_RESULTS.html']
assert all((O/p).is_file() for p in required)
for page in [O/'OPEN_RESULTS.html',F/'OPEN_RESULTS.html']:
 for target in re.findall(r'(?:href|src)="([^"]+)"',page.read_text()):assert (page.parent/target).resolve().exists(),target
st=json.loads((O/'FINAL_STATUS.json').read_text());assert len(st['independent_statuses'])==19;assert all(x['status'].startswith('PASS') for x in st['independent_statuses'].values());assert len(re.findall(r'^## \d+\.',(O/'REVIEW_Y_CLOT_NATIVE_STATIC_CFD_LOAD_INTERFACE_V1.md').read_text(),re.M))==33
(O/'DELIVERABLE_CHECKS.json').write_text(json.dumps({'required_files':required,'required_files_present':True,'HTML_relative_links_valid':True,'four_PNG_four_SVG':len(list(F.glob('Fig*.png')))==4 and len(list(F.glob('Fig*.svg')))==4,'question_count':33,'independent_status_count':19,'actual_document_SHA_matches':sha(Path('/home/lzy/projects/clot_rupture.md'))==json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text())['after_sha256']},indent=2))
shutil=__import__('shutil');shutil.copy2(Path(__file__),O/'implementation/wsl_deliver.py')
# The receipt cannot hash itself. All substantive data and code are hashed.
exclude={'PAYLOAD_MANIFEST.json','WSL_DELIVERY_MANIFEST.json','DELIVERY_VERIFICATION.json'};files={str(p.relative_to(L)):{'sha256':sha(p),'bytes':p.stat().st_size} for base in [O,F,P] for p in sorted(base.rglob('*')) if p.is_file() and p.name not in exclude and '__pycache__' not in p.parts};total=sum(z['bytes'] for z in files.values());manifest={'run':C['run'],'payload_files':len(files),'payload_bytes':total,'root_server':C['server'],'root_local':str(L),'excluded_self_referential_files':sorted(exclude),'files':files};(O/'PAYLOAD_MANIFEST.json').write_text(json.dumps(manifest,indent=2));(O/'WSL_DELIVERY_MANIFEST.json').write_text(json.dumps({'local_root':str(L),'server_root':C['server'],'payload_files':len(files),'payload_bytes':total,'verification':str(O),'figures':str(F),'code':str(P),'scope':'Essential source snapshots, mapping/geometry/surface inputs, full precision native/static audits, logs, figures, report, document backup/diff. No entire LAMMPS build, old trajectories or new CFD solve.','canonical_inventory':'PAYLOAD_MANIFEST.json','canonical_inventory_sha256':sha(O/'PAYLOAD_MANIFEST.json')},indent=2))
for base in [O,F,P]:subprocess.run(['rsync','-a','--exclude=__pycache__','-e',shlex.join(SSH[:-1]),str(base)+'/',SSH[-1]+':'+C['server']+'/'+str(base.relative_to(L))+'/'],check=True)
script="""from pathlib import Path
import hashlib,json,subprocess
R=Path(ROOT);O=R/REL;m=json.loads((O/'PAYLOAD_MANIFEST.json').read_text());bad=[]
for rel,v in m['files'].items():
 p=R/rel
 if not p.is_file() or p.stat().st_size!=v['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=v['sha256']:bad.append(rel)
code=\"\"\"import json,xmlrpc.client,supervisor.xmlrpc
r=xmlrpc.client.ServerProxy('http://127.0.0.1',transport=supervisor.xmlrpc.SupervisorTransport(None,None,'unix:///var/run/supervisor.sock'))
print(json.dumps([r.supervisor.getProcessInfo(n) for n in NAMES]))
\"\"\".replace('NAMES',repr(['clot_static_cfd_'+RUN+'_formal','clot_static_cfd_'+RUN+'_repair1']))
q=subprocess.run(['/usr/bin/python3','-c',code],capture_output=True,text=True);assert q.returncode==0,q.stderr;jobs=json.loads(q.stdout);assert all(x['statename']=='EXITED' for x in jobs);assert jobs[1]['exitstatus']==0
receipt={'status':'PASS_SERVER_AND_LOCAL_PAYLOAD_SHA','payload_files':m['payload_files'],'payload_bytes':m['payload_bytes'],'mismatches':bad,'manifest_sha256':hashlib.sha256((O/'PAYLOAD_MANIFEST.json').read_bytes()).hexdigest(),'managed_jobs':[{'name':x['name'],'statename':x['statename'],'exitstatus':x['exitstatus']} for x in jobs],'failed_formal_job_reason':'Preserved initial sorting assertion, completed from saved native states in repair1; no second static solve','no_active_task_solver':True}
(O/'DELIVERY_VERIFICATION.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2));assert not bad
""".replace('ROOT',repr(C['server'])).replace('REL',repr(REL)).replace('RUN',repr(C['run']))
remote_cmd='source /workspace/y_clot/env.sh\npython - <<\'VERIFY\'\n'+script+'\nVERIFY\n'
subprocess.run(SSH+['bash -lc '+shlex.quote(remote_cmd)],check=True)
subprocess.run(['rsync','-a','-e',shlex.join(SSH[:-1]),SSH[-1]+':'+C['server']+'/'+REL+'/DELIVERY_VERIFICATION.json',str(O)+'/'],check=True)
receipt=json.loads((O/'DELIVERY_VERIFICATION.json').read_text());assert receipt['manifest_sha256']==sha(O/'PAYLOAD_MANIFEST.json');assert all((L/k).stat().st_size==v['bytes'] and sha(L/k)==v['sha256'] for k,v in files.items());print('LOCAL_READBACK_AND_ALL_PAYLOAD_SHA_PASS',len(files),total)
