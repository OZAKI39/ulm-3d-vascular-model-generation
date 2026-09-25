"""Read-only audit of aborted outputs. Never fabricates a terminal engine state."""
from pathlib import Path
import csv,json,collections,hashlib
S=Path(__file__).resolve().parents[1]
def complete(p):
 with p.open() as f:
  d=csv.DictReader(f);header=d.fieldnames;allrows=list(d)
 rows=[r for r in allrows if r.get('disclaimer')=='NOT EXPERIMENTAL CONCENTRATION' and None not in r and all(v is not None for v in r.values())]
 return header,rows,len(allrows)-len(rows)
results=[]
for rank in [1,4]:
 D=S/'runs'/f'LONG_TRANSPORT_release_mpi{rank}';assert 'STOP_WALL_NORMAL_AMBIGUITY' in (D/'engine.log').read_text();assert not (D/'RUN_STATE.json').exists()
 files={};data={}
 for p in sorted(D.glob('*.csv')):
  h,r,bad=complete(p);data[p.name]=(h,r);files[p.name]=dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),complete_rows=len(r),incomplete_rows=bad,last_step=max((int(x['step']) for x in r if 'step' in x),default=None))
 hist=data['SOLVER_HISTORY.csv'][1];flux=data['flux_timeseries.csv'][1];hs={int(r['step']):r for r in hist};st=data['INTEGRATION_STAGES.csv'][1];wall=data['wall_constraint_timeseries.csv'][1];tr=data['TRAJECTORIES.csv'][1]
 groups=[]
 for r,keyadd in [(st,1),(wall,1),(tr,0)]:
  counts=collections.defaultdict(set)
  for x in r:counts[(int(x['step'])+keyadd,int(x.get('stage',0)))].add(x['particle_id'])
  complete_steps={k for (k,stage),ids in counts.items() if k in hs and len(ids)==int(hs[k]['active_before_removal']) and (not keyadd or (k,1-stage) in counts and len(counts[(k,1-stage)])==len(ids))}
  groups.append(complete_steps)
 common=set(hs).intersection(*groups);cut=max(common);assert set(range(1,cut+1)).intersection(hs).issuperset({k for k in hs if k<=cut})
 last=flux[-1];prefix=next(r for r in flux if int(r['step'])==cut)
 out=dict(status='ABORTED_AT_REQUIRED_STOP',stop_reason='STOP_WALL_NORMAL_AMBIGUITY',mpi_ranks=rank,last_flushed_accounting_step=int(last['step']),last_flushed_accounting_time_s=float(last['time']),last_flushed_active=int(last['active_count']),last_flushed_admitted=int(last['admitted_cumulative']),last_flushed_source_drawn=int(last['N_source_drawn']),last_flushed_pending=int(last['pending_count']),outlet_exits=sum(int(last[f'outlet_{i}_cumulative']) for i in range(3)),common_diagnostic_prefix_step=cut,common_diagnostic_prefix_time_s=float(prefix['time']),terminal_state_written=False,fatal_trial_unverified=True,files=files,normal_ambiguity_classification='Detector trigger only; physical nonuniqueness versus numerical tie not resolved',particle_at_fault='not recorded by fatal exception')
 results.append(out)
 if rank==1:
  V=S/'validation/abort_prefix_mpi1';V.mkdir(exist_ok=True)
  for name,(h,r) in data.items():
   if not h:(V/name).write_bytes(b'');continue
   if name in ['INTEGRATION_STAGES.csv','wall_constraint_timeseries.csv','GEOMETRIC_PROJECTION_EVENTS.csv','injection_events.csv','RETRY_HISTORY.csv','PAIR_HYDRODYNAMIC_EVENTS.csv']:rr=[x for x in r if int(x['step'])<cut]
   else:rr=[x for x in r if int(x.get('step',0))<=cut]
   with (V/name).open('w') as f:w=csv.DictWriter(f,fieldnames=h);w.writeheader();w.writerows(rr)
  # Empty pending list is supported by the native flushed accounting row.
  assert int(prefix['pending_count'])==0
  (V/'PENDING_FINAL.csv').write_text('particle_id,radius_m,born_s,source_u,disclaimer\n')
  for name in ['case.cfg','EXECUTION_RECEIPT.json']:(V/name).write_bytes((D/name).read_bytes())
  state=dict(reason='ABORT_PREFIX_ONLY_STOP_WALL_NORMAL_AMBIGUITY',time_s=float(prefix['time']),accepted_steps=cut,committed_LAMMPS_insertions=int(prefix['admitted_cumulative']),committed_LAMMPS_removals=0,active=int(prefix['active_count']),pending=int(prefix['pending_count']))
  (V/'DERIVED_PREFIX_STATE.json').write_text(json.dumps(state,indent=2)+'\n')
  (V/'AUDIT_VIEW_README.md').write_text('This is a derived readback view of the common complete CSV prefix, not a simulation run or recovered terminal state. Original aborted files remain unchanged under runs/LONG_TRANSPORT_release_mpi1. Near-wall event buffering and the fatal trial cannot be validated. Pending list is empty per native accounting.\n')
(S/'validation/LONG_TRANSPORT_ABORT_EVIDENCE.json').write_text(json.dumps(dict(status='BLOCKED',cases=results),indent=2)+'\n')
print(json.dumps([{k:v for k,v in x.items() if k!='files'} for x in results],indent=2))
