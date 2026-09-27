#!/usr/bin/env python3
"""Frozen Stage2 tolerances at Stage3 horizons. Every scalar/safety row is tested."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import csv,json,math,struct,sys,hashlib,gzip
import numpy as np
R=Path(__file__).resolve().parents[1]
EXPECTED='5e78c61974be7ca2a52b3add288858330804c8d74d47ae909fced717d1f02456'
def raw_sha(p,step):
 fp=p/'diagnostics/field_samples'/f'fields_{step}.bin'
 with (fp.open('rb') if fp.exists() else gzip.open(str(fp)+'.gz','rb')) as f:
  h=hashlib.sha256()
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
  return h.hexdigest()
def compare(cpu,gpu,out,checkpoints_override=None):
 contract=R/'verification/GPU_NUMERICAL_COMPARISON_CONTRACT.json'
 assert hashlib.sha256(contract.read_bytes()).hexdigest()==EXPECTED
 c=json.loads(contract.read_text());out.mkdir(parents=True,exist_ok=False)
 cr=json.loads((cpu/'RUN_TERMINAL.json').read_text());gr=json.loads((gpu/'RUN_TERMINAL.json').read_text())
 N=gr['total_steps'];milestones=[x for x in ([0,1,10,100,200] if gr.get('strong_checkpoints') else [0,100,200,1000,5000]) if x<=N]
 if checkpoints_override is not None:milestones=checkpoints_override
 hist=lambda p:list(csv.DictReader((p/'diagnostics/flow_history.csv').open()))
 a,b=hist(cpu),hist(gpu)
 assert len(a)==cr['total_steps']+1 and cr['total_steps']>=N
 assert [int(x['iteration']) for x in a]==list(range(len(a)))
 maxima={};checkpoints=[];failed=[];total=0
 fields=['quantity','step','cpu_value','gpu_value','error','tolerance','status','detail']
 stream=gzip.open(out/'ALL_NUMERICAL_CHECKS.csv.gz','wt',newline='')
 w=csv.DictWriter(stream,fieldnames=fields);w.writeheader()
 def add(q,step,error,tol,detail='',cv='',gv=''):
  nonlocal total
  ok=math.isfinite(float(error)) and error<=tol
  row=dict(quantity=q,step=step,cpu_value=cv,gpu_value=gv,error=error,tolerance=tol,status='PASS' if ok else 'FAIL',detail=detail)
  w.writerow(row);total+=1
  if q not in maxima or error>maxima[q]['error']:maxima[q]=dict(row)
  if step in milestones:checkpoints.append(dict(row))
  if not ok:failed.append(dict(row))
 add('required_history_length',N,0 if len(b)==N+1 and [int(x['iteration']) for x in b]==list(range(N+1)) else float('inf'),0)
 params=(cpu/'contracts/solver_parameters.txt').read_text().splitlines()[2].split()
 scale=float(params[6])/float(params[7]);qt=c['physical_inputs']['Qtarget_m3_s']
 nominal=json.loads((cpu/'diagnostics/inlet_command.json').read_text())['nominal_velocity_LU']*scale
 for i,(x,y) in enumerate(zip(a,b)):
  for key in ['rho_min','rho_max','rho_mean']:
   add(key,i,abs(float(x[key])-float(y[key])),c['fields']['rho_min_max_abs_difference_max'],cv=x[key],gv=y[key])
  for key in ['total_mass','control_volume_mass']:
   add(key,i,abs(float(x[key])-float(y[key]))/abs(float(x[key])),c['mass']['relative_difference_max'],cv=x[key],gv=y[key])
  for key in x:
   if key.startswith(('Qin_','Qout01_','Qout02_','Qout03_')):
    add(key,i,abs(float(x[key])-float(y[key]))/qt,c['flux']['abs_Q_gpu_minus_cpu_over_Qtarget_max'],'normalized by frozen Qtarget',x[key],y[key])
 dt=np.dtype([('id','<u8'),('rho','<f8'),('u','<f8',(3,))])
 def field(p,step):
  fp=p/'diagnostics/field_samples'/f'fields_{step}.bin'
  with (fp.open('rb') if fp.exists() else gzip.open(str(fp)+'.gz','rb')) as f:
   n=struct.unpack('<Q',f.read(8))[0];v=np.frombuffer(f.read(),dtype=dt)
  assert n==len(v)==182694 and np.all(v['id'][1:]>v['id'][:-1])
  assert np.isfinite(v['rho']).all() and np.isfinite(v['u']).all()
  return v
 for step in milestones:
  if not ((gpu/'diagnostics/field_samples'/f'fields_{step}.bin').exists() or (gpu/'diagnostics/field_samples'/f'fields_{step}.bin.gz').exists()):
   add('required_field_present',step,float('inf'),0);continue
  x,y=field(cpu,step),field(gpu,step);assert np.array_equal(x['id'],y['id'])
  add('rho_normalized_L2',step,np.linalg.norm(y['rho']-x['rho'])/np.linalg.norm(x['rho']),c['fields']['rho_normalized_l2_max'])
  den=max(np.linalg.norm(x['u']*scale),math.sqrt(len(x))*nominal*1e-12)
  add('velocity_normalized_L2',step,np.linalg.norm((y['u']-x['u'])*scale)/den,c['fields']['velocity_normalized_l2_max'])
  xm=np.linalg.norm(x['u'],axis=1).max()*scale;ym=np.linalg.norm(y['u'],axis=1).max()*scale
  add('velocity_max_normalized_difference',step,abs(ym-xm)/max(xm,nominal*1e-12),c['fields']['velocity_max_normalized_abs_difference_max'],cv=xm,gv=ym)
 for name,p,receipt in [('cpu',cpu,cr),('gpu',gpu,gr)]:
  s=json.loads((p/'diagnostics/solver_status.json').read_text())
  n=receipt['total_steps']
  add(name+'_runtime_safety',N,0 if receipt['status']=='PASS' and s.get('runtime_safety')=='PASS' and s.get('timesteps')==n and s.get('safety_checks')==n+1 else float('inf'),0)
  add(name+'_MPI_ownership_reduction',N,0 if s.get('mpi_runtime_correctness')=='PASS' and s.get('quadrature_ownership_checks',0)>=n//100+1 and receipt.get('binding_status')=='PASS' else float('inf'),0)
  for fn in ['runtime_geometry_check.json','quadrature_ownership_check.json','control_volume_ownership_check.json']:
   item=json.loads((p/'diagnostics'/fn).read_text())
   add(name+'_'+fn,N,0 if item.get('status')=='PASS' else float('inf'),0)
  safety=list(csv.DictReader((p/'diagnostics/safety_counts.csv').open()))
  add(name+'_safety_evidence_rows',N,0 if len(safety)==n+1 and [int(x['iteration']) for x in safety]==list(range(n+1)) else float('inf'),0)
  for row in safety[:N+1]:
   for key in ['fluid_nan','fluid_inf','ghost_nan','ghost_inf']:add(name+'_'+key,int(row['iteration']),int(row[key]),0)
 stream.close()
 by={(x['quantity'],x['step']):x for x in checkpoints}
 growth=[]
 for row in checkpoints:
  if row['step']<5000 or not row['tolerance']:continue
  ref=by.get((row['quantity'],1000))
  if ref and row['error']>max(4*ref['error'],0.01*row['tolerance']):
   growth.append(dict(quantity=row['quantity'],step=row['step'],error=row['error'],step1000_error=ref['error'],tolerance=row['tolerance']))
 trend='FAIL' if failed else ('GROWING' if growth else 'STABLE')
 for path,rows in [('CHECKPOINT_CORRECTNESS.csv',checkpoints),('MAXIMUM_ERRORS.csv',list(maxima.values()))]:
  with (out/path).open('w',newline='') as f:
   wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(rows)
 result=dict(BYTE_IDENTICAL=all(raw_sha(cpu,s)==raw_sha(gpu,s) for s in milestones),status='PASS' if not failed else 'FAIL',comparisons=total,failed_count=len(failed),failed_examples=failed[:20],contract_sha256=EXPECTED,
 CPU=str(cpu),GPU=str(gpu),horizon=N,thresholds_changed=False,checkpoints=milestones,CPU_history_records=len(a),GPU_history_records=len(b),
 NUMERICAL_DIVERGENCE_TREND=trend,growing_quantities=growth,velocity_scale_m_s_per_LU=scale,nominal_inlet_speed_m_s=nominal,
 maxima=maxima,comparison_definition='All Stage2 field/scalar/safety tests extended to covered Stage3 horizons; only diagnostic trend rule added before results.')
 (out/'correctness_summary.json').write_text(json.dumps(result,indent=2)+'\n')
 return result
