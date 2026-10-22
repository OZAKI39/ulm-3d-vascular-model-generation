import h5py
import numpy as np
import pytest
from scripts.straight_tube_exp1.f10z_scaling.metrics import series, validate, edge_events
from scripts.straight_tube_exp1.f10z_scaling.common import BASE, HISTORICAL, R
from scripts.straight_tube_exp1.f10z.state import series as legacy

@pytest.mark.parametrize('n', [32, 64, 128])
def test_tensor_normalization_translation_and_permutation(n):
    rng = np.random.default_rng(1)
    x = rng.normal(size=(n, 3))*10e-6
    a = np.array([x, x + np.array([1e-6, 2e-6, -3e-6])])
    m = series(np.array([0., .001]), a)
    np.testing.assert_allclose(m['rg_m'], np.sqrt(np.sum((x-x.mean(0))**2)/n), rtol=1e-14)
    np.testing.assert_allclose(m['centroid_m'][1]-m['centroid_m'][0], [1e-6, 2e-6, -3e-6], rtol=1e-13)
    b = series(np.array([0., .001]), a[:, ::-1])
    np.testing.assert_allclose(m['tensor_m2'], b['tensor_m2'], rtol=1e-13)
    np.testing.assert_array_equal(m['components'], b['components'])

def test_original_32_statistics_agree():
    with h5py.File(BASE/'N32_REGRESSION.h5') as h:
        t = h['time_s'][:]; x = h['positions_m'][:]
    a, b = series(t, x), legacy(t, x)
    for k in ['centroid_m', 'tensor_m2', 'rg_m', 'aspect', 'directional_xy', 'principal_axis', 'components', 'largest_component', 'edges']:
        np.testing.assert_array_equal(a[k], b[k])

def test_wrong_ids_and_times_rejected():
    x = np.zeros((2, 64, 3))
    with pytest.raises(ValueError): validate([0, 1], x, np.arange(32))
    with pytest.raises(ValueError): validate([0, 0], x, np.arange(64))
    with pytest.raises(ValueError): validate([0, 1], x, np.arange(64)[::-1])

def test_events_keep_transient_edges():
    e = np.array([[False, False], [True, False], [False, True]])
    assert edge_events(np.array([0., .001, .002]), e) == {'0': [(.001, True), (.002, False)], '1': [(.002, True)]}

def test_preregistered_states_same_density_and_disconnected():
    from scripts.straight_tube_exp1.f10z_scaling.common import verify, read
    verify()
    audit = read(BASE/'INITIAL_STATE_AUDIT.json')
    densities = [a['number_density_m3'] for a in audit.values()]
    np.testing.assert_allclose(densities, densities[0], rtol=1e-14)
    with h5py.File(BASE/'INITIAL_STATES.h5') as h:
        for key in h:
            x = h[key][:]; assert x.shape == (int(key), 3)
            m = series(np.array([0.]), x[None])
            assert m['components'][0] == int(key)
            assert m['min_distance_m'][0] > 2*R

def test_convergence_refuses_mismatched_initial_and_detects_translation_error(tmp_path):
    from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
    from scripts.straight_tube_exp1.f10z_scaling.common import write, sha
    rng=np.random.default_rng(312)
    x=rng.normal(size=(64,3))*10e-6
    def fixture(name,states):
        out=tmp_path/name;out.mkdir()
        with h5py.File(out/'trajectory.h5','w') as h:
            h['positions_m']=states;h['time_s']=[0.,1e-6];h['particle_id']=np.arange(64);h.attrs['flow']='STRICT_ZERO'
        write(out/'RESULT.json',dict(status='PASS',N=64,accepted_steps=1,completed_s=1e-6,requested_s=1e-6,
            trajectory_sha256=sha(out/'trajectory.h5')))
        return out
    a=fixture('base',np.array([x,x]))
    identical=fixture('identical',np.array([x,x]))
    shifted_initial=fixture('badinitial',np.array([x+1e-10,x+1e-10]))
    displaced=fixture('different',np.array([x,x+np.array([1e-7,0,0])]))
    assert compare(a,identical)['admitted']
    assert compare(a,shifted_initial)['reason']=='INITIAL_STATE_OR_PARTICLE_COUNT_MISMATCH'
    result=compare(a,displaced)
    assert not result['admitted'] and not result['metric_gates']['COM_m']
    assert result['metric_gates']['Rg_m']
