from pathlib import Path
import csv,json,hashlib
S=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out=[]
for name in ['LOW','MEDIUM','LONG_TRANSPORT']:
 a=S/'runs'/f'{name}_mpi1';b=S/'runs'/f'{name}_mpi4';files=[]
 for p in sorted(a.glob('*.csv')):
  q=b/p.name
  if p.name in ['TRAJECTORIES.csv','SOLVER_HISTORY.csv']:
   with p.open() as f:left=list(csv.DictReader(f))
   with q.open() as f:right=list(csv.DictReader(f))
   excluded={'owner_rank'} if p.name=='TRAJECTORIES.csv' else {'crossrank_pairs'}
   assert len(left)==len(right)
   assert [{k:v for k,v in r.items() if k not in excluded} for r in left]==[{k:v for k,v in r.items() if k not in excluded} for r in right],p.name
   files.append(dict(name=p.name,rows=len(left),comparison='all physical fields exact',excluded=sorted(excluded)))
  else:
   assert sha(p)==sha(q),(name,p.name);files.append(dict(name=p.name,sha256=sha(p),comparison='byte identical'))
 sa=json.loads((a/'RUN_STATE.json').read_text());sb=json.loads((b/'RUN_STATE.json').read_text());assert sa.pop('mpi_ranks')==1 and sb.pop('mpi_ranks')==4;assert sa==sb
 ra=json.loads((a/'EXECUTION_RECEIPT.json').read_text());rb=json.loads((b/'EXECUTION_RECEIPT.json').read_text());assert ra['returncode']==rb['returncode']==0;assert all(ra[k]==rb[k] for k in ['config_sha256','binary_sha256','library_sha256'])
 out.append(dict(case=name,status='PASS',time=sa['time_s'],reason=sa['reason'],files=files,terminal_state_except_rank_count='EXACT',same_binary_and_config=True))
(S/'validation/MPI_DETERMINISM.json').write_text(json.dumps(dict(status='PASS',cases=out),indent=2)+'\n');print('MPI_DETERMINISM_PASS',[(r['case'],r['time']) for r in out])
