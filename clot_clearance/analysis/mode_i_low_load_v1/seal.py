from common import *
import time
assert json.loads((O/'FINAL_STATUS.json').read_text())['stopped']
excluded={'ARTIFACT_MANIFEST.json','WSL_DELIVERY_VERIFICATION.json','SERVER_DELIVERY_RECEIPT.json'}
def files():return sorted(p for base in [O,F,Path(__file__).parent] for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in excluded)
res=json.loads((O/'RUN_RESOURCES.json').read_text());res['visualization_render_passes']=2;res['first_render_resources']=json.loads((O/'FIRST_FIGURE_PASS_RESOURCES.json').read_text());res['visualization_revision_scope']='Readable arrows, true circular front-region outline, subplot spacing and annotation placement only. No state, criteria or result changed.'
res['delivery_scope']='All 989 full-precision accepted coordinate/residual states plus endpoint/native files, code, logs, five PNG/SVG sets and documentation audit. No unrelated projects, binaries or raw historical trajectories copied.';dump(O/'RUN_RESOURCES.json',res)
sources=json.loads((O/'SOURCE_INPUT_SHA256.json').read_text())
for p in Path(__file__).parent.glob('*.py'):sources[str(p)]=sha(p)
dump(O/'SOURCE_INPUT_SHA256.json',sources)
required=['REVIEW_LPS_MODE_I_LOW_LOAD_ADMISSION_AND_SOLVER_RECOVERY_V1.md','FINAL_STATUS.json','SOURCE_AND_LIBRARY_IDENTITY.json','SPECIMEN_IDENTITY_AUDIT.json','SOLVER_FAILURE_ROOT_CAUSE_AUDIT.md','STATIC_SOLVER_RECOVERY_PLAN.json','LOW_LOAD_DESIGN_FREEZE.json','LOW_LOAD_FORCE_APPLICATION_AUDIT.json','STATIC_SOLVER_TRACE.csv','STATIC_EQUILIBRIUM_AUDIT.json','NATIVE_FORCE_GRADIENT_CROSSCHECK.json','RIGID_MODE_AND_GAUGE_AUDIT.json','NOTCH_VS_LOADING_END_STRETCH.json','LOW_LOAD_COMPLIANCE_AND_ENERGY.json','LOADING_TIMESCALE_AUDIT.md','RUN_RESOURCES.json','IMPLEMENTATION_CHECKS.json','REPRODUCE.txt','OPEN_RESULTS.html','WSL_DOCUMENT_UPDATE.json']
assert all((O/k).is_file() for k in required)
entries={str(p.relative_to(R)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in files()}
manifest=dict(task=C['task'],run_id=C['run'],created_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),root=str(R),files=entries,count=len(entries),bytes=sum(v['bytes'] for v in entries.values()),excluded_control_records=sorted(excluded),required_deliverables_present=True,closure='Payload sealed after science, figures, history SHA and WSL doc update; transfer/receipt records cite this manifest SHA without self-hash recursion.')
dump(O/'ARTIFACT_MANIFEST.json',manifest);print('SEALED',manifest['count'],manifest['bytes'],sha(O/'ARTIFACT_MANIFEST.json'),flush=True)
