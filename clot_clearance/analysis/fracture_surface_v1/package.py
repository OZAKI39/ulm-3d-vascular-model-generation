"""Final readback, packaging and reproducibility inventory. No numerical re-run."""
from common import *
from datetime import datetime,timezone
import re,shutil,subprocess,xml.etree.ElementTree as ET,resource,time
from html.parser import HTMLParser

start=time.perf_counter();CODE=R/'analysis/fracture_surface_v1'
required=['REVIEW_LPS_3D_FRACTURE_SURFACE_ENERGY_IDENTIFICATION_V1.md','FINAL_STATUS.json','DESIGN_FREEZE.json','THEORY_FORMULA_AUDIT.json','PAPER_EQUATION_AUDIT.md','SOURCE_LPS_CONSISTENCY_AUDIT.json','PRISTINE_REFERENCE_TOPOLOGY_AUDIT.json','CROSS_PLANE_REFERENCE_TOPOLOGY_AUDIT.json','LAMBDA_AND_ENERGY_DENSITY.npz','DISCRETE_SURFACE_ENERGY_RESULTS.json','INDEPENDENT_ENERGY_FORM_CROSSCHECK.json','GEOMETRIC_SENSITIVITY.json','MODEL_VS_EXPERIMENT_SCOPE.md','IMPLEMENTATION_CHECKS.json','REPRODUCE.txt','OPEN_RESULTS.html','WSL_DOCUMENT_UPDATE.json']
assert all((O/n).is_file() for n in required)
report=(O/required[0]).read_text(); assert [int(x) for x in re.findall(r'^### (\d+)\.',report,re.M)]==list(range(1,29))
status=json.loads((O/'FINAL_STATUS.json').read_text()); assert len(status['statuses'])==18 and status['GC_EFFECTIVE'] is None and not status['GC_EXPERIMENTALLY_CALIBRATED']
for p in O.rglob('*.json'):json.loads(p.read_text())
for p in F.glob('*.svg'):ET.parse(p)
png=[]
for p in F.glob('*.png'):
    b=p.read_bytes();assert b[:8]==b'\x89PNG\r\n\x1a\n';width=int.from_bytes(b[16:20],'big');height=int.from_bytes(b[20:24],'big'); assert min(width,height)>1000;png.append({'file':p.name,'width':width,'height':height})
class Parser(HTMLParser):
    def __init__(self):super().__init__();self.targets=[]
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ['href','src'] and v and not v.startswith(('https:','http:','#')):self.targets.append(v)
linkresults=[]
for p in [O/'OPEN_RESULTS.html',F/'OPEN_RESULTS.html']:
    parser=Parser();parser.feed(p.read_text())
    for v in parser.targets:assert (p.parent/v).resolve().is_file(),(p,v)
    linkresults.append({'page':str(p),'all_local_targets_exist':True,'targets':parser.targets})
# Small reproducible evaluator artifact, not a new build or a replacement of the native library.
shutil.copy2(R/'analysis/dual_track_v1/energy_y.so',O/'provenance/energy_y.so')
sourcefiles=[{'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(CODE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts]
dump(O/'TASK_SOURCE_MANIFEST.json',{'files':sourcefiles,'native_CXX_changes':False,'existing_energy_evaluator_recompiled':False})
old=json.loads((O/'HISTORICAL_PROTECTION_BEFORE.json').read_text());changed=[p for p,v in old.items() if not Path(p).is_file() or sha(p)!=v['sha256']];assert not changed
dump(O/'FINAL_PROTECTION_READBACK.json',{'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),'files_checked':len(old),'changed_or_missing':changed,'native_library_sha256':sha(B/'liblammps.so'),'production_library_sha256':sha(PROD),'env_sha256':sha('/workspace/y_clot/env.sh'),'document_is_intentional_WSL_update_outside_server_historical_set':True})
stages={k:json.loads((O/k).read_text()) for k in ['INTEGRATION_RESOURCES.json','PREPARATION_RESOURCES.json','VISUALIZATION_RESOURCES.json','CGROUP_AUDIT.json']}
vis=stages['VISUALIZATION_RESOURCES.json'];renders=[]
while vis:
    renders.append({'wall_s':vis['wall_s'],'peak_RSS_MiB':vis['peak_RSS_MiB']});vis=vis.get('previous_render')
stages['render_summary']={'passes':len(renders),'total_wall_s':sum(r['wall_s'] for r in renders),'peak_recorded_render_RSS_MiB':max(r['peak_RSS_MiB'] for r in renders)}
stages['task_resource_scope']='Per-process recorded peaks, not summed RAM or a continuously sampled all-session peak. Reading/analysis/editing/SSH wall time is separate from the 0.068s formal integration.'
stages['elapsed_since_run_id_s']=(datetime.now(timezone.utc)-datetime.strptime(C['run'],'%Y%m%dT%H%M%SZ').replace(tzinfo=timezone.utc)).total_seconds()
stages['native_instances_this_task']=0;stages['solver_jobs_started']=0;stages['GPU_used']=False;dump(O/'RESOURCE_USAGE.json',stages)
procs=subprocess.run(['ps','-eo','pid,ppid,etimes,pcpu,pmem,comm'],capture_output=True,text=True,check=True)
dump(O/'TASK_CLOSURE.json',{'closed_utc':datetime.now(timezone.utc).isoformat(),'status':'OFFLINE_WORK_COMPLETE_NO_SOLVER_JOB_TO_STOP','managed_jobs_created':0,'background_jobs_created':0,'formal_plane_integrals':2,'new_native_runs':0,'next_trial_not_started':True,'process_snapshot':procs.stdout,'other_services_untouched':True})
dump(O/'PACKAGE_CONTENT_CHECK.json',{'status':'PASS','required_files':required,'report_questions':28,'status_items':18,'figure_groups':4,'PNG_dimensions':png,'SVG_XML_valid':True,'HTML_link_checks':linkresults,'package_check_wall_s':time.perf_counter()-start,'package_process_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024})
# Exclude post-transfer proof/self hashes to avoid recursive manifests; record exclusions explicitly.
exclude={'ARTIFACT_MANIFEST.json','WSL_DELIVERY_VERIFICATION.json','SERVER_FINAL_PACKAGE_CHECK.json','FINAL_READBACK.json','package.log'}
files=[]
for base in [CODE,O,F]:
    for p in sorted(base.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and p.name not in exclude:files.append({'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)})
manifest={'format':'SHA256-relative-project-path-v1','created_utc':datetime.now(timezone.utc).isoformat(),'server_root':str(R),'WSL_root':C['local'],'files':files,'file_count':len(files),'total_bytes':sum(x['bytes'] for x in files),'excludes':sorted(exclude)+['__pycache__'],'scope':'All new-task source, full-precision input/result, figure, review, source proof, logs and WSL document backup artifacts; no historical bulk libraries or unrelated project outputs.'}
dump(O/'ARTIFACT_MANIFEST.json',manifest)
print(json.dumps({'status':'PACKAGE_READY','files':len(files),'bytes':manifest['total_bytes'],'manifest_sha256':sha(O/'ARTIFACT_MANIFEST.json')}))
