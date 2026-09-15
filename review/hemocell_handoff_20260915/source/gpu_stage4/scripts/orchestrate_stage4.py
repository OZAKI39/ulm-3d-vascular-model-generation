#!/usr/bin/env python3
"""Bounded, no-retry Stage4 schedule. Supervisor writes durable stdout locally."""
from pathlib import Path
import subprocess,sys,json,csv,hashlib,gzip,os,time,statistics,shutil,fcntl,datetime
sys.dont_write_bytecode=True
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'scripts'))
from evaluate_case import compare
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def j(p):return json.loads(p.read_text())
def compress(D):
 assert D.is_relative_to(R/'runs') and j(D/'RUN_TERMINAL.json')['status']=='PASS'
 store=R/'blobs';store.mkdir(exist_ok=True);index=[]
 for p in sorted((D/'diagnostics/field_samples').glob('*.bin')):
  rawhash=sha(p);target=store/(rawhash+'.bin.gz')
  if not target.exists():
   with p.open('rb') as a,target.open('wb') as f:
    with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0,compresslevel=6) as b:shutil.copyfileobj(a,b,2**20)
  with gzip.open(target,'rb') as f:
   h=hashlib.sha256()
   for b in iter(lambda:f.read(2**20),b''):h.update(b)
  assert h.hexdigest()==rawhash
  dest=Path(str(p)+'.gz');os.link(target,dest)
  index.append(dict(path=str(p.relative_to(D)),raw_sha256=rawhash,raw_bytes=p.stat().st_size,gzip_sha256=sha(dest),gzip_bytes=dest.stat().st_size))
  p.unlink()
 save(D/'provenance/LOSSLESS_STORAGE.json',dict(status='PASS',after_recorded_end_to_end=True,codec='gzip',entries=index))
 (D/'SHA256SUMS_RAW_PRE_STORAGE').write_text((D/'SHA256SUMS').read_text())
 (D/'SHA256SUMS').write_text(''.join(sha(p)+'  '+str(p.relative_to(D))+'\n' for p in sorted(D.rglob('*')) if p.is_file() and p.name!='SHA256SUMS'))
 print('LOSSLESS_ARCHIVE',D.name,'freeMiB',shutil.disk_usage(R).free/2**20,flush=True)
def run(name):
 print('SCHEDULE_CASE',name,flush=True)
 p=subprocess.run(['/usr/bin/python3','-B','-u',str(R/'scripts/run_case.py'),name],cwd=R)
 assert p.returncode==0,dict(case=name,returncode=p.returncode)
 D=R/'runs'/name;assert j(D/'RUN_TERMINAL.json')['status']=='PASS';compress(D);return D
def check(a,b,label,points=None):
 d=compare(a,b,R/'verification'/label,points)
 save(R/'verification'/(label+'.json'),{k:v for k,v in d.items() if k!='maxima'})
 print('NUMERICAL',label,d['status'],'byte_identical',d['BYTE_IDENTICAL'],'checks',d['comparisons'],flush=True)
 assert d['status']=='PASS',dict(numerical_failure=label)
 return d
def rate(name):return j(R/'runs'/name/'RUN_TERMINAL.json')['solver_timing']['steps_per_second']
def phase(n,reps,threshold,label):
 rows=[]
 for rep in range(1,reps+1):
  order=['baseline','batch'] if rep%2 else ['batch','baseline']
  for v in order:run(f'{v}_{n}_r{rep}')
  check(R/'runs'/f'baseline_{n}_r{rep}',R/'runs'/f'batch_{n}_r{rep}',f'pair_{n}_r{rep}')
  for v in ['baseline','batch']:
   t=j(R/'runs'/f'{v}_{n}_r{rep}'/'RUN_TERMINAL.json')
   rows.append(dict(variant=v,steps=n,repetition=rep,steps_per_second=rate(f'{v}_{n}_r{rep}'),end_to_end_seconds=t['end_to_end_seconds'],initialization_seconds=t['initialization_seconds'],status='PASS'))
 base=statistics.median(x['steps_per_second'] for x in rows if x['variant']=='baseline')
 batch=statistics.median(x['steps_per_second'] for x in rows if x['variant']=='batch')
 gate=dict(status='PASS' if batch/base-1>=threshold else 'REJECT_LOW_BENEFIT',baseline=base,batch=batch,speedup=batch/base,improvement=batch/base-1,threshold=threshold,rows=rows)
 save(R/'verification'/(label+'.json'),gate)
 out={200:'STAGE4_SHORT_BENCHMARK.csv',1000:'STAGE4_1000_BENCHMARK.csv',5000:'STAGE4_5000_BENCHMARK.csv'}[n]
 with (R/out).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 print('PERFORMANCE_GATE',label,json.dumps({k:v for k,v in gate.items() if k!='rows'}),flush=True)
 return gate
def terminal(status,reason,**extra):
 save(R/'provenance/SCHEDULE_TERMINAL.json',dict(status=status,reason=reason,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),automatic_retries=0,**extra))
 print('SCHEDULE_TERMINAL',status,reason,flush=True)
def main():
 lock=(R/'schedule.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 assert not (R/'provenance/SCHEDULE_STARTED.json').exists(),'Never repeat schedule'
 save(R/'provenance/SCHEDULE_STARTED.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),pid=os.getpid(),automatic_retries=0))
 a=run('base_correctness')
 check(R.parent/'20260914_stage3_sustained/runs/gpu_1000_r1',a,'stage3_baseline_equivalence',[0,100,200])
 b=run('batch_correctness');correct=check(a,b,'initial_strong_correctness')
 save(R/'verification/INITIAL_CORRECTNESS_GATE.json',dict(status='PASS',BYTE_IDENTICAL=correct['BYTE_IDENTICAL'],checkpoints=correct['checkpoints']))
 cpu=run('cpu_sanity');check(cpu,b,'cpu_sanity_comparison',[0,100,200])
 short=phase(200,3,.05,'SHORT_GATE')
 prof=run('batch_profile');check(R/'runs/baseline_200_r1',prof,'profile_correctness')
 result=subprocess.run(['/usr/bin/python3','-B',str(R/'scripts/analyze_post_profile.py')],cwd=R)
 assert result.returncode==0,'Profile analysis failed; evidence preserved'
 if short['status']!='PASS':terminal('REJECTED','REJECT_LOW_BENEFIT_SHORT');return
 mem=j(R/'POST_OPT_MEMORY_MIGRATION.json');pre=j(R.parent/'20260914_stage3_sustained/POST_OPT_MEMORY_MIGRATION.json')
 raw=mem['MANAGED_MIGRATION_BYTES_PER_STEP']/pre['MANAGED_MIGRATION_BYTES_PER_STEP']
 ordinary=mem['ordinary_steps_101_199']['MANAGED_MIGRATION_BYTES_PER_STEP']/pre['ordinary_steps_101_199']['MANAGED_MIGRATION_BYTES_PER_STEP']
 save(R/'verification/MIGRATION_GATE.json',dict(status='FAIL' if raw>1.25 and ordinary>1.25 else 'PASS',raw_ratio=raw,ordinary_ratio=ordinary))
 if raw>1.25 and ordinary>1.25:terminal('REJECTED','MAJOR_MIGRATION_REGRESSION');return
 middle=phase(1000,2,.10,'HORIZON_1000_GATE')
 if middle['status']!='PASS':terminal('REJECTED','REJECT_LOW_BENEFIT_1000');return
 long=phase(5000,2,.10,'HORIZON_5000_GATE')
 rates=[rate(f'batch_5000_r{i}') for i in [1,2]]
 spread=abs(rates[0]-rates[1])/statistics.mean(rates)
 drop=1-long['batch']/middle['batch']
 stable=spread<=.05 and drop<=.05
 save(R/'verification/STABILITY_GATE.json',dict(status='PASS' if stable else 'FAIL',repeat_relative_spread=spread,relative_drop_vs_1000=drop,threshold=.05))
 if long['status']!='PASS' or not stable:terminal('REJECTED','5000_BENEFIT_OR_STABILITY_GATE');return
 terminal('ACCEPTED','ALL_FROZEN_GATES_PASS',benefit_class='ACCEPT_MAJOR' if long['improvement']>=.4 else 'ACCEPT_GOOD' if long['improvement']>=.2 else 'ACCEPT_MINOR')
if __name__=='__main__':
 try:main()
 except Exception as e:
  terminal('FAILED_STOPPED',repr(e));raise
