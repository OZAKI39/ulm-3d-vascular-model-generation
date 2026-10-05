"""Independent ordered arithmetic and full constraint regression for an opt-in backend."""
import shutil

import numpy as np
import pytest

from mevo_mb.bubbles.contact import project_contacts
from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep
from mevo_mb.coupling.geometry import Cylinder


@pytest.fixture(scope="module")
def sweep(tmp_path_factory):
    if shutil.which("cc") is None:
        pytest.skip("Optional compiler unavailable; Python backend unaffected")
    return OrderedContactSweep(tmp_path_factory.mktemp("ordered-contact"))


def test_ordered_scalar_updates_match_numpy(sweep):
    rng=np.random.default_rng(20261005)
    # Overlapping, separated, variable-radius and order-dependent connected pairs.
    for scale in [1e-9,1e-6,1e-3,1.]:
        for _ in range(100):
            x=rng.normal(size=(32,3))*scale
            radii=rng.uniform(.2,2.,32)*scale
            pairs=np.array([(i,j) for i in range(31) for j in range(i+1,32)],dtype=np.int64)
            expected=x.copy()
            for i,j in pairs:
                delta=expected[j]-expected[i]
                d=np.linalg.norm(delta)
                amount=max(0.,radii[i]+radii[j]-d)
                fraction=radii[j]/(radii[i]+radii[j])
                expected[i]-=fraction*amount*delta/d
                expected[j]+=(1-fraction)*amount*delta/d
            sweep.sweep(x,radii,pairs)
            np.testing.assert_array_equal(x,expected)


def test_full_wall_and_contact_projection_matches(sweep):
    geometry=Cylinder(radius_m=10e-6,length_m=100e-6)
    radii=np.array([1e-6,1.2e-6,1.4e-6])
    positions=np.array([[8e-6,0,1e-5],[8e-6,0,1.1e-5],[8e-6,0,1.2e-5]])
    a=project_contacts(positions,radii,geometry,1e-11,1000,"ckdtree")
    b=sweep.project(positions,radii,geometry,1e-11,1000,"ckdtree")
    for x,y in zip(a[:2],b[:2]):
        np.testing.assert_array_equal(x,y)
    assert a[2:]==b[2:]


def test_invalid_pair_and_coincident_centers_fail(sweep):
    x=np.zeros((2,3))
    radii=np.ones(2)
    with pytest.raises(ValueError,match="out of bounds"):
        sweep.sweep(x,radii,np.array([[0,2]]))
    with pytest.raises(RuntimeError,match="Coincident"):
        sweep.sweep(x,radii,np.array([[0,1]]))
