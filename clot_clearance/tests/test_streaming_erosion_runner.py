"""Actual mean-flow integration and clock independence, no cyclic surrogate."""
import copy,json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.regularization.erosion_runner import run,validate

ROOT=Path(__file__).resolve().parents[1]


def small_config():
    c=json.loads((ROOT/'configs/streaming_erosion_demo_coarse.json').read_text())
    c['clot']['cells']=[4,3,3]
    c['streaming_erosion']['number_of_erosion_steps']=2
    return c


@pytest.fixture(scope='module')
def coupled(tmp_path_factory):
    root=tmp_path_factory.mktemp('erosion_coupled');outputs={}
    for name in ['baseline','legacy_steps','carrier','exposure']:
        c=small_config()
        if name=='legacy_steps':
            c['damage'].update(DeltaN=987654,cycles_per_update=123456)
            c['simulation']['represented_cycles']=999999
        if name=='carrier':
            c['streaming']['ultrasound_frequency_Hz']=7.3e6
            c['simulation']['representative_frequency_Hz']=1234.
        if name=='exposure':c['streaming_erosion']['erosion_exposure_step_s']=.02
        summary=run(c,root/name,quiet=True)
        with np.load(root/name/'states.npz') as z:data={k:z[k] for k in z.files}
        outputs[name]=(summary,data,root/name)
    return outputs


def test_small_mean_flow_safety_monotonicity_and_no_particle_deletion(coupled):
    summary,z,_=coupled['baseline'];fixed=z['fixed'];X=z['X']
    for k in ['x','v','damage','bond_D','erosion_driver_epsilon_plus']:assert np.isfinite(z[k]).all()
    assert summary['minimum_J_or_intrinsic_stretch']>.2
    assert np.min(np.diff(z['bond_D'],axis=0))>=0
    assert np.max(z['bond_D'])>0
    np.testing.assert_array_equal(z['x'][:,fixed],np.broadcast_to(X[fixed],z['x'][:,fixed].shape))
    np.testing.assert_array_equal(z['v'][:,fixed],0)
    assert z['x'].shape==(3,36,3)
    assert summary['all_particles_retained']
    assert summary['final']['broken_bond_count']==0
    assert summary['final']['detached_total_volume_m3']==0


def test_legacy_deltaN_and_represented_cycles_have_no_effect(coupled):
    for key in coupled['baseline'][1]:
        np.testing.assert_array_equal(coupled['baseline'][1][key],coupled['legacy_steps'][1][key])


def test_carrier_and_representative_frequency_metadata_have_no_effect(coupled):
    for key in coupled['baseline'][1]:
        np.testing.assert_array_equal(coupled['baseline'][1][key],coupled['carrier'][1][key])


def test_exposure_clock_never_changes_fixed_mechanical_window(coupled):
    base=coupled['baseline'][1];other=coupled['exposure'][1]
    np.testing.assert_allclose(base['mechanical_time_s'],[0,.002,.004],atol=1e-18)
    np.testing.assert_array_equal(base['mechanical_time_s'],other['mechanical_time_s'])
    np.testing.assert_allclose(base['erosion_exposure_time_s'],[0,.01,.02])
    np.testing.assert_allclose(other['erosion_exposure_time_s'],[0,.02,.04])
    np.testing.assert_array_equal(base['x'][1],other['x'][1])
    np.testing.assert_allclose(other['bond_D'][1],2*base['bond_D'][1],rtol=1e-15)


def test_every_update_uses_endpoint_strain_and_explicit_exposure(coupled):
    z=coupled['baseline'][1]
    np.testing.assert_array_equal(z['epsilon_plus'][1:],z['erosion_driver_epsilon_plus'][1:])
    expected=np.minimum(1.,z['bond_D'][:-1]+.01*.01*z['erosion_driver_epsilon_plus'][1:]/1e-4)
    np.testing.assert_allclose(z['bond_D'][1:],expected,atol=1e-18,rtol=1e-14)
    np.testing.assert_array_equal(z['bond_active'],z['bond_D']<1)
    np.testing.assert_array_equal(z['integrity'],1-z['bond_D'])


@pytest.mark.parametrize('change',[('damage','mode','COMBINED'),('streaming','oscillatory_amplitude',1.),
                                  ('streaming','mean_fraction',.5),('damage','D_break',.8)])
def test_runner_rejects_unsupported_mode_or_non_mean_loading(change):
    c=small_config();section,key,value=change;c[section][key]=value
    with pytest.raises(ValueError):validate(c)


def test_static_resolved_field_uses_strain_not_transient_amplitude(tmp_path):
    import pyvista as pv
    from scipy.spatial import Delaunay
    c=small_config();c['streaming_erosion']['number_of_erosion_steps']=1
    # Analytical steady incompressible shear, covering the entire moving clot.
    xyz=np.indices((2,2,2)).reshape(3,-1).T*.006-.003
    cells=Delaunay(xyz).simplices
    grid=pv.UnstructuredGrid(np.column_stack((np.full(len(cells),4),cells)).ravel(),np.full(len(cells),10,np.uint8),xyz)
    gradient=np.array([[0.,0.,10.],[0.,0.,0.],[0.,0.,0.]])
    grid['velocity']=xyz@gradient.T+[.02,0.,0.];grid['pressure']=np.zeros(len(xyz))
    path=tmp_path/'static_field.vtu';grid.save(path)
    c['streaming'].update(provider='resolved_field',field_path=str(path),metadata_path=dict(
        length_unit='m',velocity_unit='m/s',pressure_unit='Pa',time_unit='s',
        normal_convention='outward_solid',sample_domain='convex_hull'))
    summary=run(c,tmp_path/'resolved_run',quiet=True)
    with np.load(tmp_path/'resolved_run/states.npz') as z:
        endpoint=z['epsilon_plus'][-1]
        assert endpoint.max()>1e-10
        np.testing.assert_allclose(z['bond_D'][-1],.01*.01*endpoint/1e-4,atol=1e-18,rtol=1e-14)
    assert summary['minimum_J_or_intrinsic_stretch']>.2
    identity=json.loads((tmp_path/'resolved_run/IDENTITY.json').read_text())
    assert identity['flow_provider']=='ResolvedStreamingFieldProvider'
    assert identity['coupon_calibration_used'] is False
    settling=json.loads((tmp_path/'resolved_run/settling_windows.json').read_text())[0]
    assert len(settling['samples'])>=2
