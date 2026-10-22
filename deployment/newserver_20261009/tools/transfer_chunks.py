import json, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

root=Path(__file__).resolve().parents[1]
ssh='ssh -i /home/lzy/.ssh/vast_step3b_ed25519 -p 15283 -o IdentitiesOnly=yes -o BatchMode=yes -o ConnectTimeout=15 -o IPQoS=none'
def transfer(path):
    start=time.monotonic()
    args=['rsync','-a','--no-owner','--no-group','--compress','--partial','--partial-dir=.rsync-partial','--stats',f'--files-from={path}','-e',ssh,str(root/'stage')+'/', 'root@95.3.33.46:/workspace/brava_cfd/']
    with (root/'logs'/f'transfer_chunk_{path.stem}.log').open('w') as log:
        p=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT)
    row=dict(chunk=path.stem,returncode=p.returncode,seconds=time.monotonic()-start)
    print(json.dumps(row),flush=True)
    return row
results=[]
with ThreadPoolExecutor(max_workers=2) as pool:
    for job in as_completed([pool.submit(transfer,p) for p in sorted((root/'inventory/chunks').glob('*.txt'))]):
        results.append(job.result())
        (root/'inventory/CHUNK_TRANSFER_RESULTS.json').write_text(json.dumps(results,indent=2)+'\n')
raise SystemExit(any(r['returncode'] for r in results))
