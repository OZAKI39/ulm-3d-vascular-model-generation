import copy,json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.regularization.equilibrium import make_fixture
from pd_clot.regularization.damage_energy_audit import audit_state
from pd_clot.regularization.fast_mechanics import evaluate_supported
from pd_clot.fragment_mechanics import prepare_supported_shape
from pd_clot.regularization.energy import potential

ROOT=Path(__file__).resolve().parents[1]


def setup():
    c=json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text());c['clot']['cells']=[4,3,3]
    cl=make_fixture(c);x=cl.X@np.diag([1.01,.9975,.9975]);return c,cl,x


def test_total_and_particle_energy_components():
    c,cl,x=setup();a=audit_state(cl,x,c['material'],c['safety'],g=np.ones(len(cl.pairs)))
    assert a.global_energy['U_total_J']==a.global_energy['U_constitutive_J']+a.global_energy['U_stabilization_J']
    for name in ['constitutive','stabilization']:
        assert np.isclose(a.particle[name+'_energy_J'].sum(),a.global_energy['U_'+name+'_J'],rtol=1e-14,atol=1e-30)


@pytest.mark.parametrize('pattern',['none','uniform','single','localized'])
def test_undeformed_damage_cannot_create_resolved_energy(pattern):
    c,cl,_=setup();D=np.zeros(len(cl.pairs))
    if pattern=='uniform':D[:]=.1
    elif pattern=='single':D[0]=.1
    elif pattern=='localized':D[np.any(cl.pairs==5,axis=1)]=.1
    a=audit_state(cl,cl.X,c['material'],c['safety'],D=D)
    assert abs(a.global_energy['U_total_J'])<1e-20
    assert np.all(a.particle['support_rank']==3)


def test_uniform_full_rank_scaling_preserves_F_and_scales_energy():
    c,cl,x=setup();baseline=audit_state(cl,x,c['material'],c['safety'],g=np.ones(len(cl.pairs)))
    energies=[]
    for factor in [1.,.99,.95,.9,.8]:
        a=audit_state(cl,x,c['material'],c['safety'],g=np.full(len(cl.pairs),factor))
        np.testing.assert_allclose(a.particle['F'],baseline.particle['F'],rtol=0,atol=1e-12)
        np.testing.assert_allclose(a.global_energy['U_total_J'],factor*baseline.global_energy['U_total_J'],rtol=1e-8,atol=1e-25)
        energies.append(a.global_energy['U_total_J'])
    assert np.all(np.diff(energies)<0)


def test_nonaffine_uniform_scaling_includes_stabilization():
    c,cl,x=setup();x[15,0]+=.002*cl.spacing
    a=audit_state(cl,x,c['material'],c['safety'],g=np.ones(len(cl.pairs)))
    b=audit_state(cl,x,c['material'],c['safety'],g=np.full(len(cl.pairs),.8))
    assert a.global_energy['U_stabilization_J']>0
    for key in a.global_energy:np.testing.assert_allclose(b.global_energy[key],.8*a.global_energy[key],rtol=1e-8,atol=1e-25)


def test_audit_inputs_and_production_results_are_unchanged():
    c,cl,x=setup();g=np.linspace(.9,1,len(cl.pairs));D=1-g;oldc=copy.deepcopy(c)
    oldx=x.copy();oldg=g.copy();oldD=D.copy();oldX=cl.X.copy()
    shape=prepare_supported_shape(cl,g,c['safety']);before=evaluate_supported(cl,x,g,c['material'],c['safety'],shape)
    a=audit_state(cl,x,c['material'],c['safety'],g=g);audit_state(cl,x,c['material'],c['safety'],D=D)
    after=evaluate_supported(cl,x,g,c['material'],c['safety'],shape)
    for u,v in zip(before,after):np.testing.assert_array_equal(u,v)
    np.testing.assert_array_equal(a.force_N,before[0]);assert a.global_energy['U_total_J']==sum(potential(cl,before))
    for u,v in [(x,oldx),(g,oldg),(D,oldD),(cl.X,oldX)]:np.testing.assert_array_equal(u,v)
    assert c==oldc


def test_api_rejects_ambiguous_or_invalid_integrity():
    c,cl,x=setup()
    with pytest.raises(ValueError,match='exactly one'):audit_state(cl,x,c['material'],c['safety'],g=np.ones(len(cl.pairs)),D=np.zeros(len(cl.pairs)))
    with pytest.raises(ValueError,match='Integrity'):audit_state(cl,x,c['material'],c['safety'],g=np.full(len(cl.pairs),1.1))
