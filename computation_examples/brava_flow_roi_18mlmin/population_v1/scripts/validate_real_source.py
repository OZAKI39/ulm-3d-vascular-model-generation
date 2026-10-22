"""Unconditioned source statistics on the real FEM inlet (not accepted-only)."""
from pathlib import Path
import argparse,csv,json,sys
import numpy as np
from scipy.stats import kstest,chisquare
ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);args=ap.parse_args()
sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.source import PerfusionSource
from particle_3d.population_v1.provenance import write_json
root=Path(__file__).resolve().parents[1];env=load_environment(json.loads((root/'data/environment.json').read_text()))
s=PerfusionSource(env,concentration_m3=4e9,seed=2026092992);n=30000
rows=[s.marks(i+1) for i in range(n)];dt=np.array([r['arrival_delta_t_s'] for r in rows]);diameter=np.array([r['diameter_m'] for r in rows])
xyz=np.array([r['position_m'] for r in rows]);triangles=np.array([r['inlet_triangle_id'] for r in rows])
weights=env.sampler.weights/env.sampler.Q_m3_s;counts=np.bincount(triangles,minlength=len(weights))
# Predeclared ten cumulative-flux groups avoid invalid low-expected-count chi-square bins.
groups=np.minimum((np.cumsum(weights)*10).astype(int),9)
expected=np.bincount(groups,weights=n*weights,minlength=10);observed=np.bincount(groups,weights=counts,minlength=10)
poisson=kstest(dt,'expon',args=(0,1/s.rate));size=kstest(diameter,env.distribution.cdf);flux=chisquare(observed,expected)
error=(xyz[0]+(xyz-xyz[0]).mean(0))-env.sampler.expectation()
centroid_budget=6*xyz.std(0)/np.sqrt(n)+64*np.finfo(float).eps*np.abs(xyz).max()
ok=abs((1/dt.mean())/s.rate-1)<6/np.sqrt(n) and min(poisson.pvalue,size.pvalue,flux.pvalue)>1e-5 and np.all(abs(error)<centroid_budget)
record=dict(PASS=bool(ok),phase=4,flow_sha256=env.flow_sha256,n=n,source_identity=s.identity,
    empirical_arrival_rate_s_inv=float(1/dt.mean()),poisson_ks_pvalue=float(poisson.pvalue),
    size_ks_pvalue=float(size.pvalue),flux_chi_square_pvalue=float(flux.pvalue),
    flux_group_expected=expected,flux_group_observed=observed,centroid_error_m=error,
    source_draws_unconditioned_on_admission=True,flux_audit=env.flux_audit,
    test_role='DIAGNOSTIC_SOURCE_MARKS_SEPARATE_SEED_NOT_PRODUCTION_EVENTS')
write_json(root/'validation/phase4_real_inlet_statistics.json',record)
with (root/'validation/real_inlet_triangle_sampling.csv').open('w',newline='') as stream:
    writer=csv.writer(stream);writer.writerow(['inlet_triangle_id','theoretical_flux_probability','expected_count','observed_count','flux_group'])
    for i in range(len(weights)):writer.writerow([i,weights[i],n*weights[i],counts[i],groups[i]])
print((root/'validation/phase4_real_inlet_statistics.json').read_text(),flush=True)
if not ok:raise SystemExit(1)
