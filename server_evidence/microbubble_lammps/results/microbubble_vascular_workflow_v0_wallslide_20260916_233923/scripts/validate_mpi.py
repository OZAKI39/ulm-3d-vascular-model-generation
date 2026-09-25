from pathlib import Path
import csv,json,hashlib
S=Path(__file__).resolve().parents[1];results=[]
for name in ['LOW','MEDIUM','LONG_TRANSPORT']:
 a=S/'runs'/f'{name}_release_mpi1';b=S/'runs'/f'{name}_release_mpi4';checks={}
 for file in ['injection_events.csv','flux_timeseries.csv','INTEGRATION_STAGES.csv','PORT_EVENTS.csv','NEAR_WALL_EVENTS.csv','PAIR_HYDRODYNAMIC_EVENTS.csv','RETRY_HISTORY.csv','wall_constraint_timeseries.csv','GEOMETRIC_PROJECTION_EVENTS.csv']:
  checks[file]=a.joinpath(file).read_bytes()==b.joinpath(file).read_bytes()
 def rows(p):
  with p.open() as f:return [r for r in csv.DictReader(f) if r.get('disclaimer')=='NOT EXPERIMENTAL CONCENTRATION' and None not in r and all(v is not None for v in r.values())]
 r1=rows(a/'TRAJECTORIES.csv');r4=rows(b/'TRAJECTORIES.csv');checks['trajectory_complete_rows_same']=len(r1)==len(r4);maxerr=0.;ownership=0
 for x,y in zip(r1,r4):
  ownership+=x.pop('owner_rank')!=y.pop('owner_rank')
  for k in ['x_m','y_m','z_m']:maxerr=max(maxerr,abs(float(x[k])-float(y[k])))
  assert x==y,('STOP_MPI_DIVERGENCE',name,x,y)
 checks['trajectory_all_other_fields_bitwise_equal']=True
 if name!='LONG_TRANSPORT':
  checks['PENDING_FINAL.csv']=a.joinpath('PENDING_FINAL.csv').read_bytes()==b.joinpath('PENDING_FINAL.csv').read_bytes()
  s1=json.loads((a/'RUN_STATE.json').read_text());s4=json.loads((b/'RUN_STATE.json').read_text());s1.pop('mpi_ranks');s4.pop('mpi_ranks');checks['terminal_state_equal_except_rank_count']=s1==s4
 else:
  e=json.loads((S/'validation/LONG_TRANSPORT_ABORT_EVIDENCE.json').read_text())['cases'];checks['same_stop_reason']=e[0]['stop_reason']==e[1]['stop_reason'];checks['same_last_flushed_accounting']=all(e[0][k]==e[1][k] for k in ['last_flushed_accounting_step','last_flushed_accounting_time_s','last_flushed_active','last_flushed_admitted','last_flushed_source_drawn','last_flushed_pending','outlet_exits'])
 assert all(checks.values()),('STOP_MPI_DIVERGENCE',checks)
 results.append(dict(case=name,status='PASS',scope='persisted outputs and identical STOP only; missing terminal state and fatal trial unverified' if name=='LONG_TRANSPORT' else 'complete terminal run',checks=checks,complete_trajectory_rows=len(r1),max_position_difference_m=maxerr,owner_rank_label_differences=ownership))
(S/'validation/MPI_CONSISTENCY.json').write_text(json.dumps(dict(status='PASS',fatal_trial_state_verified=False,cases=results),indent=2)+'\n');print(json.dumps(results,indent=2))
