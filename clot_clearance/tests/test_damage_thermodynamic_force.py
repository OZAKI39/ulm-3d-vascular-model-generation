import json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.regularization.equilibrium import make_fixture
from pd_clot.regularization.damage_energy_audit import audit_state,thermodynamic_force

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('initial_damage',[0.,.2])
def test_FD_converges_to_affine_analytic_damage_force(initial_damage):
    c=json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text());c['clot']['cells']=[4,3,3]
    cl=make_fixture(c);x=cl.X@np.diag([1.01,.9975,.9975]);D=np.full(len(cl.pairs),initial_damage)
    a=audit_state(cl,x,c['material'],c['safety'],D=D);ids=np.arange(0,len(D),55)
    f=thermodynamic_force(cl,x,D,ids,c['material'],c['safety'],[1e-6,3e-6,1e-5,3e-5,1e-4])
    i,j=cl.pairs.T
    denominator=np.bincount(i,weights=cl.volume[j],minlength=len(x))+np.bincount(j,weights=cl.volume[i],minlength=len(x))
    W=a.particle['W_undegraded_Pa'];analytic=cl.volume[i]*cl.volume[j]*(W[i]/denominator[i]+W[j]/denominator[j])
    assert f['stable'].all() and f['sign_resolved'].all()
    np.testing.assert_allclose(f['estimate_J'],analytic[ids],rtol=1e-4,atol=1e-20)
    assert np.all(f['Y_J']>0)
    assert np.all(f['scheme']==('forward' if initial_damage==0 else 'central'))


def test_finite_difference_does_not_advance_damage_or_move_particles():
    c=json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text());c['clot']['cells']=[4,3,3]
    cl=make_fixture(c);x=cl.X@np.diag([1.01,1.,1.]);D=np.zeros(len(cl.pairs));xx=x.copy();dd=D.copy()
    thermodynamic_force(cl,x,D,[0,1],c['material'],c['safety'],[1e-6,3e-6,1e-5])
    np.testing.assert_array_equal(x,xx);np.testing.assert_array_equal(D,dd)
