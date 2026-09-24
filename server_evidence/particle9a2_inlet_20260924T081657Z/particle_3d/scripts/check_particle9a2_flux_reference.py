#!/usr/bin/env python3
"""Independent original-sampler rejection reference, used only after births.

This audit fixes each diameter in advance; future trajectories do not feed back.
It also checks real-vessel acceleration against a direct, unaccelerated sampler.
"""
from concurrent.futures import ProcessPoolExecutor
import time,json
import numpy as np
from scipy.stats import ks_2samp
from run_particle9a2 import setup,R,D,ROLES
from particle_3d.particle81_simulation import environment,dump
from particle_3d.injection_method_c import sample_flux_weighted_feasible_position
from particle_3d.particle82_point_native import NativePointTracer

SOURCE=None;POINT=None

def probe(spec):
    k,diameter_um=spec;start=time.time();s=SOURCE;d=diameter_um*1e-6;r=d/2
    e=dict(species='MB',particle_id=k+1,diameter_m=d,radius_m=r,q=[1.,0.,0.,0.],diameter_draw_id=[2026092495,k,920,0])
    original=s.sampler;rng=np.random.default_rng([2026092495,k,923]);centers=[];attempts=0;counts={v:0 for v in ROLES}
    while len(centers)<500:
        p,tri=original.sample(rng,n=512)
        for x in p:
            attempts+=1
            # Same exact sphere predicate used by legacy B, followed by full checker.
            distance=s.wall.nearest_center_triangle(x)[1]
            if distance-r-2e-9 < -s.wall.roundoff_m:continue
            particle,status,_=s.checker.check(e,x,{})
            if particle is None:continue
            centers.append(x.tolist())
            trace=POINT.trace(x);counts[trace['outlet'] or 'NO_EXIT']+=1
            if len(centers)==500:break
        if attempts>2000000:raise RuntimeError('AUDIT_ORIGINAL_PROPOSAL_GUARD')
    proposal,mapping,bounds=s.tree.feasible_proposal(r);accelerated=[]
    for pid in range(1,501):
        e['particle_id']=pid
        p,_,_=sample_flux_weighted_feasible_position(e,sampler=proposal,mapping=mapping,checker=s.checker,seed=2026092496+k,guard=bounds['position_guard'],bounds=bounds)
        accelerated.append(p.tolist())
    orig=np.array(centers);acc=np.array(accelerated)
    # Fixed projected axes from inlet SVD; no basin labels in comparison.
    allp=original.triangles.reshape(-1,3);origin=allp.mean(axis=0);basis=np.linalg.svd(allp-origin,full_matrices=False)[2][:2]
    a=(orig-origin)@basis.T;b=(acc-origin)@basis.T
    ks=[float(ks_2samp(a[:,j],b[:,j]).statistic) for j in range(2)]
    return dict(diameter_um=diameter_um,original_attempts=attempts,accepted_count=500,
        original_feasible_flux_fraction_estimate=500/attempts,conditional_point_counts=counts,
        original_positions_m=centers,accelerated_positions_m=accelerated,
        acceleration_bounds=bounds,position_coordinate_two_sample_KS=ks,
        empirical_position_agreement=bool(max(ks)<.14),
        diagnostic_threshold='.14: conservative simultaneous empirical CDF tolerance; no outlet target',
        wall_seconds=time.time()-start)

def main():
    global SOURCE,POINT
    SOURCE=setup(2026092495);POINT=NativePointTracer(environment());start=time.time()
    with ProcessPoolExecutor(max_workers=6) as pool:rows=list(pool.map(probe,enumerate([.8,1.,1.2,1.6,2.,2.6])))
    dump(D/'independent_flux_reference.json',dict(rows=rows,workers=6,wall_seconds=time.time()-start,
        all_position_agreement=all(r['empirical_position_agreement'] for r in rows),
        role='AUDIT_ONLY_FIXED_DIAMETER_UNACCELERATED_ORIGINAL_INLET_FLUX_SAMPLER'))
    print(json.dumps([{k:v for k,v in r.items() if 'positions' not in k} for r in rows],indent=2),flush=True)
if __name__=='__main__':main()
