"""The twelve requested integration gates; no old campaigns are rerun."""
import hashlib,json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[1]
def read(n):return json.loads((OUT/n).read_text())
def test_01_pd2_mechanics_cache_and_force_regression_unchanged():
    for p,h in read('PROTECTED_INPUT_HASHES.json').items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    a=read('PD2_TARGET_SUMMARIES.json')
    for k in ['shear_affine','nonaffine_axial','nonaffine_shear']:
        assert a[k]['previous_state_x_identical'] and a[k]['previous_force_identical'] and a[k]['previous_energy_identical']
    assert read('PD2_SOLVER_ADAPTER_CHECK.json')['PASS']
    # The historical normal-state expression differs by floating arithmetic
    # from the boundary surrogate. Evaluate its exact saved coordinates too.
    from verification.higher_order_baqp_boundary_reintegration_v1.run_audit import primary
    _,_,backend=primary()
    with np.load(ROOT/'verification/higher_order_baqp_energy_force_candidate_v1/AFFINE_AXIAL_STATE.npz') as old:
        new=backend.evaluate(old['x'])
        assert np.array_equal(new.force,old['force'])
        assert new.energy==old['energy_J']

def test_02_backend_api_and_tiny_historical_c1_conversion():
    from pd_clot.pfpd_candidate.kernel import Family
    from pd_clot.pfpd_candidate.energy import Material
    from pd_clot.pfpd_candidate.higher_order_rk.mechanics_candidate import Reference,GradientOperator,VariationalBAQP
    from pd_clot.pfpd_candidate.nonlocal_boundary.mechanics_backend import C1Backend,PD2Backend
    from pd_clot.pfpd_candidate.nonlocal_boundary.localized_conversion import convert_localized
    from pd_clot.pfpd_candidate.nonlocal_boundary.volume_force_conversion import convert_volume_force
    from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import BlockBoundary
    from pd_clot.pfpd_candidate.nonlocal_boundary.pd2_bulk_reference import field_for
    old=ROOT/'verification/higher_order_baqp_energy_force_candidate_v1'
    with np.load(old/'N7_REFERENCE_GEOMETRY.npz') as a:
        f=Family.from_arrays(a['X'],a['volume'],a['pairs'],.00025125,a['omega']);s,t=a['source'],a['target']
    with np.load(old/'N7_PD2_REFERENCE_CACHE.npz') as a:g=a['g']
    material=Material('svk',1000.,2000./3);field=field_for('shear_affine');x=field.position(f.X);zero=np.zeros_like(x)
    pd2=PD2Backend(VariationalBAQP(Reference.from_family(f),GradientOperator('PD2',s,t,g),material),f)
    out=convert_localized(pd2,x,zero,zero,f.volume)
    assert np.array_equal(out['f_BC_localized'],-pd2.force(x))
    assert np.isfinite(pd2.energy(x))
    assert pd2.node_F(x).shape==(len(x),3,3)
    assert pd2.bond_Ftilde(x).shape==(len(f.pairs),3,3)
    assert np.isfinite(pd2.psi_plus_diagnostic(x)).all()
    c1=C1Backend(f,material);geom=BlockBoundary(f.X.min(0),f.X.max(0))
    old_load=convert_volume_force(geom,f,field).nodal_force_N
    assert np.array_equal(old_load,convert_localized(c1,x,zero,zero,f.volume)['f_BC_localized'])

def test_03_affine_bulk_force():
    assert all(r['bulk_force_scaled']<=1e-10 for r in read('PD2_AFFINE_BULK_GATE.json'))
    assert read('PD2_FORCE_STENCIL_CLOSURE.json')['domains'][-1]['all_required_nodes_complete']

def test_04_nonaffine_nodewise_decomposition():
    summaries=read('PD2_TARGET_SUMMARIES.json')
    with np.load(OUT/'PD2_NODEWISE_BOUNDARY_DECOMPOSITION.npz') as a:
        for k in ['nonaffine_axial','nonaffine_shear']:
            get=lambda key:a[k+'__'+key]
            residual=get('f_FD')+get('f_BC_localized')+get('f_body_local')-get('f_model')
            assert np.linalg.norm(residual,axis=1).max()/summaries[k]['scales']['F_ref_N']<=1e-10
            assert np.linalg.norm(get('f_bulk_target'))>0
    assert all(r['PASS'] for r in read('PD2_ORACLE_VALIDATION.json'))

def test_05_fully_closed_center():
    a=read('PD2_CLOSED_CENTER_LOCALIZATION.json');closed=[r for r in a if r['fully_closed']]
    assert len(closed)==2 and all(r['boundary_force_scaled']<=1e-10 for r in closed)

def check_load(k):
    a=read('PD2_TARGET_SUMMARIES.json')[k]
    for name in ['load_checks']+(['normal_load_checks','tangential_load_checks'] if k=='shear_affine' else []):
        r=a[name]
        assert max(r['resultant_error'],r['moment_error'],r['virtual_work_error'])<=1e-6
def test_06_normal_force_moment_work():check_load('normal_affine')
def test_07_shear_force_moment_work():check_load('shear_affine')
def check_recovery(k):
    r=read('PD2_'+k.upper()+'.json')
    assert r['solver']['converged'] and r['force_residual']<=1e-8
    assert r['initial_perturbation_max_m']>0 and r['fixed_displacement_error_m']==0
    assert all(v['F_error_max_scaled']<=1e-5 and v['Psi_plus_error_max_scaled']<=1e-4 for v in r['fields'].values())
    if k in ('normal_affine','shear_affine'):
        a=read('PD2_AFFINE_LOAD_ROUNDOFF_ACCOUNTING.json')[k]
        assert a['external_minus_conversion_scaled']<=1e-10
        c=a['actual_external_load_checks']
        assert max(c['resultant_error'],c['moment_error'],c['virtual_work_error'])<=1e-6
def test_08_normal_equilibrium():check_recovery('normal_affine')
def test_09_shear_equilibrium():check_recovery('shear_affine')
def test_10_nonaffine_axial_recovery():check_recovery('nonaffine_axial')
def test_11_nonaffine_shear_recovery():check_recovery('nonaffine_shear')
def test_12_classical_shear_target_residual():
    assert read('PD2_CLASSICAL_TARGET_RESIDUAL.json')['error_scaled']<=1e-10
    s=read('PD2_TARGET_SUMMARIES.json')['nonaffine_shear']['scales']['F_ref_N']
    with np.load(OUT/'PD2_NODEWISE_BOUNDARY_DECOMPOSITION.npz') as a:
        residual=a['nonaffine_shear__f_FD']+a['nonaffine_shear__f_external_classical']
        assert np.linalg.norm(residual-a['nonaffine_shear__f_model'],axis=1).max()/s<=1e-10
    assert read('PD2_CLASSICAL_SHEAR.json')['force_residual']<=1e-8
