"""Seal only this round's deliverables, no large historical tree scan."""
from common import *
import shutil,subprocess
assert (O/'FINAL_STATUS.json').exists()
protected=json.loads((O/'PROTECTED_FILES_BEFORE.json').read_text())['files']
for p,h in protected.items():
 if p!=str(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md'):assert sha(p)==h,(p,'unexpected historical mutation')
assert subprocess.check_output(['git','diff','--binary'],cwd=R)==(O/'tracked_diff.before.patch').read_bytes()
for p in A.glob('*.py'):shutil.copy2(p,O/'source_snapshot'/p.name)
files={}
paths=[R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md']+sorted(A.glob('*.py'))+[p for base in [O,F] for p in base.rglob('*') if p.is_file() and p.name not in ['PAYLOAD_MANIFEST.json','PAYLOAD_MANIFEST.sha256','SERVER_DELIVERY_CHECK.json']]
for p in sorted(set(paths)):
 if '__pycache__' not in str(p):files[str(p.relative_to(R))]=dict(bytes=p.stat().st_size,sha256=sha(p))
manifest=dict(version=1,task=TASK,run_id=RUN,scope='Current new analysis/audits/FP64 evidence/figures/document backup and mutable active entrypoint. Reuses prior verified library instead of duplicating source/build caches.',file_count=len(files),total_bytes=sum(q['bytes'] for q in files.values()),files=files,library_dependency=dict(server_path=str(LIB),sha256=sha(LIB),local_prior_binary_relative='verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/binary/liblammps.so'))
assert manifest['total_bytes']<600*1024**2
dump(O/'PAYLOAD_MANIFEST.json',manifest);(O/'PAYLOAD_MANIFEST.sha256').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'\n')
print('SEALED',manifest['file_count'],manifest['total_bytes'],sha(O/'PAYLOAD_MANIFEST.json'),flush=True)
