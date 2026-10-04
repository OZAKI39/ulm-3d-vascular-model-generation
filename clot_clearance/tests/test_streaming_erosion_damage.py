"""Analytic slow-time checks, deliberately independent of mechanics."""
import numpy as np
import pytest
from pd_clot.regularization.erosion import StreamingErosionDamage, tensile_strain, SettlingMonitor


def config(**kw):
    return dict(C_t_s_inv=.01,m_t=1.,epsilon_ref=1e-4,threshold_ratio=0.,
                D_break=1.,tension_only=True,spatial_feedback=False,demonstration_parameter=True,**kw)


def test_zero_driver_any_exposure_preserves_initial_damage():
    law=StreamingErosionDamage(3,config(),[0.,.4,1.])
    for dt in [1e-12,.01,1e9]:law.advance(np.zeros(3),dt)
    np.testing.assert_array_equal(law.D,[0.,.4,1.])


def test_constant_driver_matches_independent_analytic_solution():
    c=config();c.update(m_t=1.5,threshold_ratio=.25)
    eps=np.array([0.,.5e-4,1e-4,2e-4]);initial=np.array([0.,.1,.2,.3])
    law=StreamingErosionDamage(4,c,initial)
    for _ in range(10):law.advance(eps,.01)
    expected=initial+.1*.01*np.maximum(eps/1e-4-.25,0.)**1.5
    np.testing.assert_allclose(law.D,expected,rtol=2e-15,atol=2e-16)


def test_compression_does_not_damage():
    law=StreamingErosionDamage(4,config())
    law.advance(tensile_strain([.01,.8,.999,1.]),100.)
    np.testing.assert_array_equal(law.D,np.zeros(4))


def test_irreversible_capped_and_loss_of_connectivity_only_at_one():
    law=StreamingErosionDamage(3,config(),[.9,.99,1.])
    assert law.active.tolist()==[True,True,False]
    ids,old,increment=law.advance(np.array([1e-4,1.,1.]),1.)
    np.testing.assert_array_equal(ids,[1]);assert np.all(increment>=0)
    assert law.D[0]==pytest.approx(.91);assert law.active.tolist()==[True,False,False]
    np.testing.assert_array_equal(law.g,1.-law.D)
    previous=law.D.copy();law.advance(np.zeros(3),1e6)
    np.testing.assert_array_equal(law.D,previous);np.testing.assert_array_equal(law.g[1:],[0.,0.])


@pytest.mark.parametrize('value',[.5,.999,1.01])
def test_break_threshold_must_equal_one(value):
    c=config();c['D_break']=value
    with pytest.raises(ValueError,match='D_break'):StreamingErosionDamage(1,c)


@pytest.mark.parametrize('dt',[.02,.01,.005])
def test_constant_driver_time_step_consistency_without_mechanical_clock(dt):
    eps=np.array([0.,1e-5,1e-4,2e-4]);law=StreamingErosionDamage(4,config())
    for _ in range(round(.1/dt)):law.advance(eps,dt)
    np.testing.assert_allclose(law.D,.1*.01*eps/1e-4,rtol=4e-15,atol=1e-18)


@pytest.mark.parametrize('key,value',[('epsilon_ref',0.),('m_t',0.),('C_t_s_inv',-1.),
                                    ('spatial_feedback',True),('tension_only',False),('demonstration_parameter',False)])
def test_invalid_or_unimplemented_parameters_rejected(key,value):
    c=config();c[key]=value
    with pytest.raises(ValueError):StreamingErosionDamage(1,c)


def test_settling_monitor_diagnoses_bond_changes_hidden_by_scalar_max():
    monitor=SettlingMonitor([True,True,False])
    monitor.add(0.,0.,[1e-4,0.,99.],1e-10)
    monitor.add(.002,.002,[0.,1e-4,0.],1e-11)
    summary=monitor.summary(.05)
    assert summary['maximum_driver_relative_drift']==0
    assert summary['mean_driver_relative_drift']==0
    assert summary['driver_relative_drift']==1
    assert summary['settling_status']=='MECHANICAL_STEADY_STATE_NOT_ESTABLISHED'
    assert summary['settling_final_kinetic_energy_J']==1e-11


def test_constant_driver_settling_passes_without_creating_damage():
    monitor=SettlingMonitor([True,True])
    for t in [0.,.001,.002]:monitor.add(t,t,[1e-4,2e-4],1e-20)
    assert monitor.summary(.05)['settling_status']=='DRIVER_SETTLING_CHECK_PASSED'


def test_energy_increase_is_preserved_signed_and_flagged_not_dissipation():
    from pd_clot.regularization.erosion_energy import ErosionEnergyLedger
    ledger=ErosionEnergyLedger(1.)
    allocation=ledger.damage_drop(1.,1.2,np.array([1.,3.]))
    np.testing.assert_allclose(allocation,[-.05,-.15])
    row=ledger.row(1.2,0.,0.)
    assert abs(row['numerical_energy_residual_J'])<1e-15
    assert row['damage_energy_increase_count']==1
    assert row['damage_energy_status']=='ENERGY_DISSIPATION_NOT_ESTABLISHED'
    assert row['damage_stored_energy_drop_J']==pytest.approx(-.2)


def test_positive_energy_drop_is_accounted_without_coupon_energy():
    from pd_clot.regularization.erosion_energy import ErosionEnergyLedger
    ledger=ErosionEnergyLedger(1.)
    np.testing.assert_allclose(ledger.damage_drop(1.,.7,np.array([1.,2.])),[.1,.2])
    row=ledger.row(.7,0.,0.)
    assert abs(row['numerical_energy_residual_J'])<1e-15
    assert row['damage_energy_increase_count']==0
