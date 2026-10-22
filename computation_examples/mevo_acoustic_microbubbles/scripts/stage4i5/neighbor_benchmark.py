"""Zero-miss and independent broad-phase runtime stress tests with real gold states."""
import time
import h5py
import numpy as np

from mevo_mb.bubbles.contact import contact_candidates
from mevo_mb.bubbles.neighbor_cache import ContactNeighborCache
from scripts.stage4i5.common import GOLD,new_run,finish


def main():
    out=new_run("p1_vtk_contact",dict(test="neighbor_cache_stress",seed=20261005,skin_m=7.12e-6))
    rng=np.random.default_rng(20261005)
    rows=[]
    for kind in ["sparse","dense","wall_near","high_velocity","variable_radius","gold_B2","gold_B3"]:
        samples=[]
        r=rng.uniform(1e-6,3e-6,256)
        if kind.startswith("gold"):
            with h5py.File(GOLD/"cases"/kind[-2:]/"particles/trajectories.h5") as h:
                for i in range(len(h["time_s"])):
                    mask=h["active"][i]
                    samples.append((h["positions_m"][i][mask],h["radii_m"][:][mask]))
        else:
            x=rng.normal(size=(256,3))*(2e-5 if kind=="sparse" else 6e-6)
            if kind=="wall_near":
                x[:,2]=r
            for i in range(400):
                x+=rng.normal(size=x.shape)*(1e-5 if kind=="high_velocity" else 5e-8)
                if kind=="variable_radius" and i%11==0:
                    r=r*rng.uniform(.95,1.05,len(r))
                samples.append((x.copy(),r.copy()))
        cache=ContactNeighborCache(7.12e-6,1e-11)
        missing=extra=0
        begin=time.perf_counter()
        original=[contact_candidates(x,r,1e-11,"ckdtree") for x,r in samples]
        reference_s=time.perf_counter()-begin
        begin=time.perf_counter()
        candidate=[cache.candidates(x,r,1e-11) for x,r in samples]
        candidate_s=time.perf_counter()-begin
        for a,b in zip(original,candidate):
            aa={tuple(p) for p in a}
            bb={tuple(p) for p in b}
            missing+=len(aa-bb)
            extra+=len(bb-aa)
        rows.append(dict(kind=kind,queries=len(samples),missed_true_contacts=missing,extra_narrow_phase_contacts=extra,
            sorted_pairs_identical=all(np.array_equal(a,b) for a,b in zip(original,candidate)),
            original_s=reference_s,candidate_s=candidate_s,speedup=reference_s/candidate_s,
            rebuilds=cache.rebuilds,reuses=cache.reuses,maximum_broad_phase_candidates=cache.maximum_candidates))
    finish(out,dict(status="PASS" if all(r["sorted_pairs_identical"] for r in rows) else "FAIL",rows=rows,
        scope="Independent search-only cost, not overall trajectory speedup. Real-state cadence often forces rebuild; actual projection-loop reuse is measured by the S1/S2 replay."))


if __name__=="__main__":
    main()
