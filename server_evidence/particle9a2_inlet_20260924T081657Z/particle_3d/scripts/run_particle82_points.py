#!/usr/bin/env python3
"""Remote-only point-tracer basin diagnostic, with raw per-seed path storage."""
from pathlib import Path
import argparse,json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy.stats import binomtest
from particle_3d.particle82_provenance import atomic_json, require_remote, sha256
from particle_3d.particle82_tracers import init_environment,run


def main(args):
    provenance=json.loads(Path(args.host_provenance).read_text())
    require_remote(provenance,provenance['hostname'])
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    start=time.time();env=init_environment()
    if args.positions:
        positions=np.load(args.positions)['positions']
        source=dict(role='EXPLICIT_DIAGNOSTIC_POSITIONS_NO_RESAMPLING',sha256=sha256(args.positions))
    else:
        rng=np.random.default_rng(2026092182)
        positions,faces=env.sampler.sample(rng,args.count)
        source=dict(seed=2026092182,role='P7_EXACT_POSITIVE_P1_FLUX_WEIGHTED_POINT_TRACER_PROPOSALS',rng_final=rng.bit_generator.state)
        np.savez_compressed(out/'birth_positions.npz',positions=positions,faces=faces)
    config=dict(step_m=args.step_m,error=args.error,horizon_m=2e-3)
    shards=run(positions,out/'shards',provenance,workers=args.workers,shard_size=args.shard_size,config=config)
    rows=sorted([r for s in shards for r in s['rows']],key=lambda r:r['tracer_id'])
    counts={o:sum(r['outlet']==o for r in rows) for o in ['OUTLET_01','OUTLET_02','OUTLET_03']}
    unresolved=len(rows)-sum(counts.values());fraction={}
    for outlet,n in counts.items():
        p=env.audit[outlet]['positive_Q_m3_s']/env.sampler.Q_m3_s
        test=binomtest(n,len(rows),p)
        ci=test.proportion_ci(confidence_level=1-.05/3,method='exact')
        fraction[outlet]=dict(count=n,observed_fraction=n/len(rows),frozen_flux_fraction=p,
            binomial_two_sided_p=float(test.pvalue),familywise_alpha=.05,per_outlet_alpha=.05/3,
            simultaneous_bonferroni_exact_interval=[float(ci.low),float(ci.high)],
            sampling_compatibility=bool(test.pvalue>=.05/3))
    summary=dict(role='DIAGNOSTIC_ONLY_ZERO_RADIUS_NO_P65_NO_SONOVUE',count=len(rows),outlet_counts=counts,
        unresolved=unresolved,three_outlets_observed=all(counts.values()),
        all_resolved=unresolved==0,fraction_vs_frozen=fraction,source=source,config=config,
        numerical_fraction_convergence='REQUIRES_REFINED_POINT_TRACER_COMPARISON',
        wall_seconds=time.time()-start,compute_host=provenance['hostname'],workers=args.workers,
        git_commit=provenance['source_git_commit'],shards=[dict(path=f'shards/tracer_{s["first_id"]:07d}_{s["first_id"]+s["count"]-1:07d}.json',sha256=sha256(out/'shards'/f'tracer_{s["first_id"]:07d}_{s["first_id"]+s["count"]-1:07d}.json')) for s in shards])
    atomic_json(out/'POINT_TRACER_SUMMARY.json',summary)
    print(json.dumps({k:summary[k] for k in ['count','outlet_counts','unresolved','wall_seconds','fraction_vs_frozen']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',required=True);p.add_argument('--host-provenance',required=True)
    p.add_argument('--count',type=int,default=100000);p.add_argument('--workers',type=int,default=8)
    p.add_argument('--shard-size',type=int,default=512);p.add_argument('--step-m',type=float,default=.2e-6)
    p.add_argument('--error',type=float,default=1e-11);p.add_argument('--positions')
    main(p.parse_args())
