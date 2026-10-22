"""Numerical receipt for the same fixed-seed analytic fixture as permanent tests."""
from pathlib import Path
import importlib.util,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
p=ROOT/'particle_3d/tests/particle9a4_population_inlet/conftest.py'
spec=importlib.util.spec_from_file_location('p9a4_benchmark_fixture',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
start=time.perf_counter();ledger=m.poiseuille.__wrapped__();rows=ledger.rows
d=np.array([r['diameter_m'] for r in rows]);rad=np.array([np.linalg.norm(r['position_m'][:2])/2e-6 for r in rows]);keep=np.array([r['particle_id'] is not None for r in rows])
sizes=np.array([.6e-6,2.4e-6]);x=(2e-6-sizes/2-2e-9)/2e-6;fraction=2*x*x-x**4;p=np.array([.4,.6])@fraction
source_errors=[abs(float(np.mean(rad<=r))-(2*r*r-r**4)) for r in [.2,.4,.6,.8]]
by_size=[]
for size,xc,target in zip(sizes,x,fraction):
 selected=d==size;accepted=selected&keep
 radial_errors=[abs(float(np.mean(rad[accepted]<=xc*t))-(2*(xc*t)**2-(xc*t)**4)/target) for t in [.25,.5,.75]]
 by_size.append(dict(diameter_um=size*1e6,source_count=int(selected.sum()),accepted_count=int(accepted.sum()),
  fraction_observed=float(keep[selected].mean()),fraction_expected=float(target),max_conditional_radial_cdf_error=max(radial_errors)))
times=np.array([r['proposal_time_s'] for r in rows]);bins=np.arange(0,times[-1],4.);counts=np.histogram(times,bins)[0];accepted_counts=np.histogram(times[keep],bins)[0]
result=dict(N=len(rows),accepted=ledger.accepted,lambda_source=5.,lambda_enter_expected=5*p,lambda_enter_observed=ledger.accepted/ledger.time_s,
 max_source_radial_cdf_error=max(source_errors),small_size_entering_probability_expected=.4*fraction[0]/p,
 small_size_entering_probability_observed=float((d[keep]==sizes[0]).mean()),by_size=by_size,
 source_fano=float(counts.var()/counts.mean()),entering_fano=float(accepted_counts.var()/accepted_counts.mean()),
 tolerances=dict(source_radial_cdf_absolute=.015,conditional_radial_cdf_absolute=.025,size_marginal_absolute=.018,entering_rate_relative=.035),
 fano_role='DESCRIPTIVE_NOT_STRICT_GATE',runtime_s=time.perf_counter()-start)
assert max(source_errors)<.015 and all(r['max_conditional_radial_cdf_error']<.025 for r in by_size)
assert abs(result['small_size_entering_probability_observed']-result['small_size_entering_probability_expected'])<.018
assert abs(result['lambda_enter_observed']/result['lambda_enter_expected']-1)<.035
out=ROOT/'particle_3d/reports/particle9a4_population_inlet/data/synthetic_poiseuille_audit.json'
with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
print(json.dumps(result,indent=2))
