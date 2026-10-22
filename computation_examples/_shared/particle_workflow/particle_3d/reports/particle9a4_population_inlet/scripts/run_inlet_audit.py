"""Prepare diagnostics/contract, generate physical source ledger, or replay workers."""
from pathlib import Path
import argparse,hashlib,json,math,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.population_inlet_p9a4 import load_new_environment,make_contract,make_source,generate,REPORT_RELATIVE
from particle_3d.continuous_infusion import canonical_bytes,PopulationLedger,sha256
from particle_3d.particle82a_geometry import lower_gap,maximum_handoff_radius,InletGeometryAudit
R=ROOT/REPORT_RELATIVE
CONTRACT=ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json'


def write(path,value):
 with Path(path).open('xb') as f:f.write(canonical_bytes(value))

def wilson(k,n):
 z=1.959963984540054;p=k/n;den=1+z*z/n;center=(p+z*z/(2*n))/den;half=z/den*math.sqrt(p*(1-p)/n+z*z/(4*n*n))
 return [max(0.,center-half),min(1.,center+half)]

def prepare(env):
 if CONTRACT.exists():raise FileExistsError('Immutable contract already exists')
 t=time.perf_counter();count=20000;seed=2026092595
 xyz,ids=env.sampler.sample(np.random.default_rng(seed),count)
 distances=np.array([env.wall.nearest_center_triangle(x)[1] for x in xyz]);owners=np.array([env.field.locate(x)[0] for x in xyz]);assert min(owners)>=0
 diameter_grid=np.array([.5,.8,1.,1.2,1.6,2.,2.6,3.,3.5,4.])
 qacc=[];checked=0
 for d in diameter_grid:
  a=d*.5e-6;accepted=distances-a-lower_gap(a)>=-env.wall.roundoff_m
  # Independently validate this diagnostic shortcut against original check.
  for j in range(0,count,250):
   e=dict(species='MB',particle_id=j+1,radius_m=a,q=[1,0,0,0]);particle,status,_=env.checker.check(e,xyz[j],{})
   assert (particle is not None)==bool(accepted[j]),(j,d,status);checked+=1
  n=int(accepted.sum());qacc.append(dict(diameter_um=float(d),accepted=n,samples=count,fraction=n/count,ci95=wilson(n,count),Q_acc_m3_s=n/count*env.sampler.Q_m3_s))
 assert np.all(np.diff([x['fraction'] for x in qacc])<=0)
 maxD=2*maximum_handoff_radius(distances+env.wall.roundoff_m)
 weight=np.asarray(env.distribution.cdf(maxD));p=float(weight.mean());se=float(weight.std(ddof=1)/np.sqrt(count));ci=[max(0,p-1.96*se),min(1,p+1.96*se)]
 provisional=make_contract(ROOT,env);rate=provisional['lambda_source_s_inv']
 estimate=dict(s_inv=rate*p,ci95_s_inv=[rate*x for x in ci],acceptance_probability=p,ci95_probability=ci,
  method='INDEPENDENT_FULL_FLUX_MONTE_CARLO_WITH_EXACT_CONDITIONAL_CDF_SIZE_INTEGRATION',seed=seed,samples=count,authoritative_check_parity_cases=checked)
 contract=make_contract(ROOT,env,lambda_enter_estimate=estimate);write(CONTRACT,contract)
 geometry=InletGeometryAudit(env)
 write(R/'data/qacc_diagnostic.json',dict(rows=qacc,rate_estimate=estimate,wall_roundoff_m=env.wall.roundoff_m,
  common_positions=True,monotone=True,source_generation_uses_this_diagnostic=False,runtime_s=time.perf_counter()-t))
 np.savez_compressed(R/'data/qacc_positions.npz',position_m=xyz,wall_distance_m=distances,owner=owners,triangle_id=ids)
 np.savez_compressed(R/'data/inlet_geometry.npz',triangles_m=env.sampler.triangles,q_m_s=env.sampler.q,origin=geometry.origin,basis=geometry.basis)
 write(R/'data/open_inlet_geometry.json',geometry.semantics())
 print('PREPARED',json.dumps(estimate),flush=True)


def summary(ledger,contract,runtime):
 rows=ledger.rows;accepted=np.array([r['particle_id'] is not None for r in rows]);d=np.array([r['diameter_um'] for r in rows]);time_s=np.array([r['proposal_time_s'] for r in rows]);dt=np.array([r['arrival_delta_t_s'] for r in rows]);clearance=np.array([r['clearance_m'] for r in rows])[accepted]
 def stats(a):return dict(mean=float(np.mean(a)),p10=float(np.quantile(a,.1)),p50=float(np.quantile(a,.5)),p90=float(np.quantile(a,.9)))
 from collections import Counter
 reasons=dict(Counter(r['rejection_reason'] for r in rows if r['particle_id'] is None))
 windows=np.arange(0,time_s[-1],20/ledger.identity['lambda_source_s_inv']);counts=np.histogram(time_s,bins=windows)[0];accepted_counts=np.histogram(time_s[accepted],bins=windows)[0]
 return dict(N_source_proposals=len(rows),N_accepted=int(accepted.sum()),N_rejected=int((~accepted).sum()),acceptance_fraction=float(accepted.mean()),rejection_fraction=float((~accepted).mean()),
  physical_source_time_span_s=float(time_s[-1]),lambda_source_s_inv=ledger.identity['lambda_source_s_inv'],lambda_enter_measured=float(accepted.sum()/time_s[-1]),lambda_enter_expected=contract['lambda_enter']['estimate']['s_inv'],
  acceptance_ci95=wilson(int(accepted.sum()),len(rows)),source_diameter_um=stats(d),enter_diameter_um=stats(d[accepted]),accepted_clearance_m=stats(clearance),
  source_mean_interarrival_s=float(dt.mean()),source_variance_interarrival_s2=float(dt.var()),accepted_mean_interarrival_s=float(np.diff(np.r_[0,time_s[accepted]]).mean()),
  fano_source_audit=float(counts.var()/counts.mean()),fano_accepted_audit=float(accepted_counts.var()/accepted_counts.mean()),fano_window_expected_source_count=20,fano_is_strict_ci_gate=False,
  rejection_reasons=reasons,runtime_s=runtime,workers=None,all_proposals_recorded=True,trajectory_executed=False)


def main():
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','generate','replay']);p.add_argument('--workers',type=int,default=6);p.add_argument('--count',type=int,default=100000);p.add_argument('--label',default='inlet100k');a=p.parse_args()
 env=load_new_environment(ROOT)
 if a.mode=='prepare':prepare(env);return
 contract=json.loads(CONTRACT.read_text());source=make_source(env,contract);t=time.perf_counter()
 ledger=generate(source,a.count,workers=a.workers,progress=lambda n:print('PROPOSALS',n,flush=True) if n%4096==0 else None)
 if a.mode=='generate':
  folder=R/'data'/a.label;ledger.checkpoint(folder);write(folder/'accepted_births.json',dict(identity=ledger.identity,events=ledger.birth_events()))
  write(folder/'smoke30_births.json',dict(identity=ledger.identity,events=ledger.birth_events(30)))
  report=summary(ledger,contract,time.perf_counter()-t);report['workers']=a.workers;write(folder/'summary.json',report);print(json.dumps(report),flush=True)
 else:
  original=R/'data/inlet100k/proposal_ledger.jsonl';digest=hashlib.sha256(b''.join(canonical_bytes(r) for r in ledger.rows)).hexdigest();assert digest==sha256(original)
  births=json.loads((R/'data/inlet100k/accepted_births.json').read_text())['events'];assert canonical_bytes(births)==canonical_bytes(ledger.birth_events())
  write(R/f'data/worker_replay_{a.workers}.json',dict(workers=a.workers,count=a.count,proposal_ledger_sha256=digest,accepted_births_byte_identical=True,all_proposals_byte_identical=True,runtime_s=time.perf_counter()-t))
  print('ALL_SOURCE_AND_ACCEPTED_EVENTS_BITWISE_IDENTICAL',a.workers,flush=True)
if __name__=='__main__':main()
