from dataclasses import FrozenInstanceError,replace
import numpy as np
import pytest
from particle_3d.rbc import RBCGeometry


def test_geometry_is_immutable_and_consistent(geometries):
    for g in geometries:
        assert 0<g.c_m<g.a_m==g.b_m
        assert 0<g.r<1 and -1<g.jeffery_lambda<0
        with pytest.raises(FrozenInstanceError):g.a_m=1.
        with pytest.raises(FrozenInstanceError):g.provenance.seed=2
        with pytest.raises(ValueError):replace(g,a_m=2*g.a_m)
        with pytest.raises(ValueError):replace(g,volume_m3=2*g.volume_m3)
        with pytest.raises(ValueError):replace(g,provenance=None)


@pytest.mark.parametrize("a,c",[(0,1e-6),(-1e-6,1e-6),(1e-6,2e-6),(1e-6,1e-6),(np.nan,1e-6),(2e-6,np.inf)])
def test_invalid_geometry_rejected(a,c):
    with pytest.raises(ValueError):RBCGeometry.synthetic_only(a,c)


def test_artificial_geometry_explicitly_synthetic_only():
    g=RBCGeometry.synthetic_only(2e-6,1e-6)
    assert g.provenance.role==g.distribution_id=="SYNTHETIC_ONLY"
