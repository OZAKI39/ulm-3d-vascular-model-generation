from common import *
import datetime
assert json.loads((O/'DELIVERY_VALIDATION.json').read_text())['status']=='PASS'
assert (F/'OPEN_RESULTS.html').exists()
freeze=json.loads((O/'SOURCE_FREEZE.json').read_text())
for p,h in freeze['implementation_sha256'].items():assert sha(R/p)==h,p
for p,h in freeze['input_sha256'].items():assert sha(O/p)==h,p
dump(O/'POSTPROCESS_SOURCE_FREEZE.json',dict(files={str(p.relative_to(R)):sha(p) for p in A.glob('*.py')},runtime_freeze_sha256=sha(O/'SOURCE_FREEZE.json')))
resources=json.loads((O/'RUN_RESOURCES.json').read_text())
resources['figure_render_calls']=json.loads((O/'FIGURE_RESOURCES.json').read_text())
resources['figure_detail_fix']=json.loads((O/'FIGURE_DETAIL_FIX_RESOURCES.json').read_text())
resources['postprocess_driver']=json.loads((O/'POSTPROCESS_RESOURCES.json').read_text())
dump(O/'RUN_RESOURCES.json',resources)
excludes={'DELIVERY_MANIFEST.json','LOCAL_READBACK.json'}
files={}
for root in [A,O,F]:
    for p in sorted(root.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.name in excludes or p.name.endswith('.tmp'):continue
        item=dict(bytes=p.stat().st_size,sha256=sha(p))
        if p.is_symlink():item['symlink_target']=os.readlink(p)
        files[str(p.relative_to(R))]=item
dump(O/'DELIVERY_MANIFEST.json',dict(task=TASK,run=RUN,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=files,total_bytes=sum(q['bytes'] for q in files.values()),file_count=len(files),native_library='Frozen server dependency, not copied into this output bundle'))
print('SEALED',len(files),sum(q['bytes'] for q in files.values()),sha(O/'DELIVERY_MANIFEST.json'))
