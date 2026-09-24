import json
import pytest
from scipy.stats import binomtest

def test_exact_sampling_statistics_and_unresolved_denominator(point_root,evidence_root):
    s=json.loads((point_root/'POINT_TRACER_SUMMARY.json').read_text())
    assert sum(s['outlet_counts'].values())+s['unresolved']==s['count']
    for outlet,r in s['fraction_vs_frozen'].items():
        test=binomtest(r['count'],s['count'],r['frozen_flux_fraction']);ci=test.proportion_ci(1-.05/3,method='exact')
        assert r['observed_fraction']==r['count']/s['count']
        assert r['binomial_two_sided_p']==pytest.approx(test.pvalue)
        assert r['simultaneous_bonferroni_exact_interval']==pytest.approx([ci.low,ci.high])
        assert r['sampling_compatibility']==(test.pvalue>=.05/3)
    comparison=json.loads((evidence_root/'POINT_REFINEMENT_AUDIT.json').read_text())
    assert comparison['identical_seed_positions'] and comparison['count']==s['count']
    if comparison['refined_summary']['unresolved'] or not comparison['sampling_only_flux_match']:
        assert comparison['flux_sanity_status']=='NOT_ESTABLISHED'
