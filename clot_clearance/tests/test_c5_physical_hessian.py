import inspect
import numpy as np
from scripts import c5_stability_audit as a
from scripts import c5_precise_energy as precise


def rows():return a.read(a.OUT/'ALL_STATES.json')['states']


def test_physical_gradient_is_negative_internal_force():
    d=a.old.ca.read_npz(a.OUT/'states/S1/STATE.npz');p=a.problem(float(d['strain']),d['lower']);z=p.pack(d['x'],d['s']);_,g=p.objective_gradient(z)
    expected=-p.state(*p.unpack(z)).force[p.free].ravel()*p.length_scale/p.energy_scale
    np.testing.assert_array_equal(g[:p.nx],expected)


def test_Hxx_symmetry_and_three_scale_convergence():
    for r in rows():
        assert r['Hxx_gate'] and r['Hxx_FD_relative_spread']<.02
        assert len(r['assembly'])==3
        assert max(x['Hxx_raw_asymmetry'] for x in r['assembly'])<1e-5


def test_ten_mechanical_curvatures_three_amplitudes():
    for r in rows():
        o=a.read(a.Path(r['cache'])/'MECHANICAL_ORACLE.json');assert o['passed'] and len(o['directions'])==10
        assert all(len(d['rows'])==3 and d['passed'] for d in o['directions'])
        assert o['max_relative_error']<.02


def test_Hss_analytic_matches_physical_gradient_FD():
    for label in ['S1','S5','S7']:
        r=a.read(a.OUT/f'states/{label}/STABILITY.json')
        assert r['Hss_gate'] and max(x['Hss_relative_error'] for x in r['assembly'])<1e-6


def test_cross_blocks_reciprocal():
    for r in rows():assert max(x['cross_relative_error'] for x in r['assembly'])<1e-5


def test_extended_energy_reuses_original_code_objects():
    from pd_clot.pfpd_candidate import baqp
    copied=precise.rebind(baqp.evaluate,dict(baqp.evaluate.__globals__))
    assert copied.__code__ is baqp.evaluate.__code__
    d=a.old.ca.read_npz(a.OUT/'states/S7/STATE.npz');p=a.problem(float(d['strain']),d['lower']);z=p.pack(d['x'],d['s']);x,s=p.unpack(z)
    assert abs(precise.energy(p,x,s)/p.objective_gradient(z)[0]-1)<1e-12


def test_spectral_switch_detector_preserves_raw_and_resolved_flags():
    before=np.array([[-.01,1e-30,.02]]);after=np.array([[-.01,-1e-30,.02]])
    r=a.spectral_switch(before,after);assert r['raw_count']==1 and r['resolved_count']==0
    r=a.spectral_switch(before,np.array([[-.01,-1e-5,.02]]));assert r['resolved_count']==1


def test_switching_probes_are_not_used_as_centered_evidence():
    count=0
    for r in rows():
        d=a.read(a.OUT/f'states/{r["label"]}/MODE_ENERGY.json')
        for q in d['rows']:
            if any(p.get('spectral',{}).get('raw_count',0) for p in q['points']):
                count+=1;assert not q['centered_used'] and q['validation_scheme']=='one_sided_physical_energy'
    assert count>0


def test_zero_spectral_surfaces_cancel_on_fixed_intact_phase_bonds():
    for r in rows():
        q=r['spectral'];assert q['near_zero_bonds']==68
        assert q['smooth_restricted_energy'] and q['exposed_switches']==0
        assert q['near_zero_in_critical_phase']==0 and q['near_zero_with_unequal_split_weight']==0
