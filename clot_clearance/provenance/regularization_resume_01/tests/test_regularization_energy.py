import json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.geometry import make_cloud
from pd_clot.fragment_mechanics import prepare_supported_shape, evaluate_supported
from pd_clot.regularization.energy import EnergyLedger, potential, kinetic, bond_quadrature_volume


def fixture():
    c=json.loads((Path(__file__).parents[1]/'configs/fragmentation/pilot_001.json').read_text())
    c['clot']['cells']=[4,4,4]
    return c,make_cloud(c['clot'])


def test_energy_partition_and_homogeneous_softening():
    c,g=fixture();w=bond_quadrature_volume(g)
    assert np.all(w>0)
    np.testing.assert_allclose(w.sum(),g.volume.sum(),rtol=1e-14)
    x=g.X@np.diag([1.01,.999,1.002]);one=np.ones(len(g.pairs))
    def energy(b):
        return sum(potential(g,evaluate_supported(g,x,b,c['material'],c['safety'],prepare_supported_shape(g,b,c['safety']))))
    a,b=energy(one),energy(.7*one)
    np.testing.assert_allclose(b,.7*a,rtol=1e-9,atol=1e-20)
    ledger=EnergyLedger(a);allocation=ledger.damage_drop(a,b,w)
    np.testing.assert_allclose(allocation.sum(),a-b)
    assert abs(ledger.row(b,0.,0.)['numerical_energy_residual_J'])<1e-20


def test_energy_increase_is_not_hidden():
    with pytest.raises(FloatingPointError):EnergyLedger(1.).damage_drop(1.,1.1,np.ones(3))


def test_exact_damping_kinetic_ledger():
    m=np.array([1.,2.]);v=np.array([[1.,0,0],[0,2.,0]])
    a=kinetic(m,v);b=kinetic(m,v*np.exp(-.1))
    ledger=EnergyLedger(a);ledger.damping_dissipation_J=a-b
    assert abs(ledger.row(0.,0.,b)['numerical_energy_residual_J'])<1e-14
