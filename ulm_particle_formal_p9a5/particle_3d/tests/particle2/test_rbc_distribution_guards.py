from math import erf,sqrt,log
import numpy as np
from particle_3d.rbc_distribution import classify_candidate


def cdf_normal(x,mean,sd):
    return np.array([.5*(1+erf(float((v-mean)/(sd*sqrt(2))))) for v in x])


def ks_distance(x,cdf):
    n=len(x);return max(np.max(np.arange(1,n+1)/n-cdf),np.max(cdf-np.arange(n)/n))


def test_paired_rng_raw_draws_source_order_and_counts(population,contract):
    ledger=population.candidates;meta=population.metadata
    rng=np.random.default_rng(meta["seed"])
    reference=np.concatenate([rng.normal([6.79,47.9],[.93,7.0892],size=(4096,2))
                              for _ in range(meta["generated_candidate_count"]//4096)])[:len(ledger)]
    np.testing.assert_array_equal(reference[:,0],ledger["D_raw_um"])
    np.testing.assert_array_equal(reference[:,1],ledger["V_raw_fL"])
    mask=ledger["status"]=="ACCEPT"
    np.testing.assert_array_equal(population.samples["candidate_id"],ledger["candidate_id"][mask])
    np.testing.assert_array_equal(population.samples["D_um"],ledger["D_raw_um"][mask])
    np.testing.assert_array_equal(population.samples["V_fL"],ledger["V_raw_fL"][mask])
    assert meta["candidate_count"]==meta["N"]+meta["guard_rejection_count"]+meta["shape_rejection_count"]
    assert meta["acceptance_rate"]==meta["N"]/len(ledger)
    assert ledger[-1]["status"]=="ACCEPT" and mask.sum()==100000
    for code,key in [("D_GUARD","diameter_guard_rejections"),("V_GUARD","volume_guard_rejections"),("SHAPE","shape_rejection_count")]:
        assert (ledger["status"]==code).sum()==meta[key]


def test_latent_and_guard_pass_distributions_with_independent_erf_cdf(population):
    ledger=population.candidates
    # Pre-shape population includes shape rejections. DKW alpha=1e-6 fixed.
    guarded=ledger[np.isin(ledger["status"],["ACCEPT","SHAPE"])]
    lo,hi=.5*(1+erf(-3/sqrt(2))),.5*(1+erf(3/sqrt(2)))
    for name,mean,sd in [("D_raw_um",6.79,.93),("V_raw_fL",47.9,7.0892)]:
        raw=np.sort(ledger[name]);passed=np.sort(guarded[name])
        assert ks_distance(raw,cdf_normal(raw,mean,sd))<=sqrt(log(2/1e-6)/(2*len(raw)))
        cdf=(cdf_normal(passed,mean,sd)-lo)/(hi-lo)
        assert ks_distance(passed,cdf)<=sqrt(log(2/1e-6)/(2*len(passed)))


def test_guard_order_is_disjoint_and_inclusive(contract):
    assert classify_candidate(3.,100.,contract)[0]=="D_GUARD"
    assert classify_candidate(6.,100.,contract)[0]=="V_GUARD"
    assert classify_candidate(4.,26.6324,contract)[0]=="ACCEPT"
    assert classify_candidate(9.58,69.1676,contract)[0]=="ACCEPT"
    for d,v in [(np.nan,40),(6,np.inf)]:
        assert classify_candidate(d,v,contract)[0] in ["D_GUARD","V_GUARD"]
