"""Finish existing authorized campaign only. No CFD or extra cohort creation."""
from pathlib import Path
import hashlib,json,subprocess,time
HERE=Path(__file__).resolve().parents[1]
while not (HERE/'data/final_summary.json').exists():
    status=subprocess.run(['supervisorctl','-c',str(HERE/'supervisord.conf'),'status','production'],
                          text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).stdout
    if any(word in status for word in ['FATAL','STOPPED','EXITED']):
        (HERE/'data/DELIVERY_BLOCKED.json').write_text(json.dumps(dict(reason='Production exited before final summary',supervisor_status=status),indent=2))
        raise RuntimeError(status)
    time.sleep(5)
subprocess.run(['/root/particle8_2_runs/env/bin/python',str(HERE/'scripts/render_results.py')],check=True)
files={}
for folder in ['data','figures','animations','scripts']:
    for path in sorted((HERE/folder).rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts and path.suffix!='.so':
            with path.open('rb') as stream:files[str(path.relative_to(HERE))]=hashlib.file_digest(stream,'sha256').hexdigest()
for name in ['MICROBUBBLE_RESULTS_ZH.md','OPEN_RESULTS.html']:
    with (HERE/name).open('rb') as stream:files[name]=hashlib.file_digest(stream,'sha256').hexdigest()
(HERE/'DELIVERY_COMPLETE.json').write_text(json.dumps(dict(files=files,count=1500,completed_unix_s=time.time()),indent=2)+'\n')
print('DELIVERY_COMPLETE',flush=True)
