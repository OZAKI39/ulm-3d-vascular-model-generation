from common import *
from terms import metric,terms,RHO,MU,DT
import pyvista as pv

def test_production_binary_source_identity():
    a=read(HERE/'PRODUCTION_SOURCE_IDENTITY.json');assert a['source_identity']=='PRODUCTION_SOLVER_SOURCE_IDENTITY_ESTABLISHED'
    assert sha(a['actual_binary'])==a['actual_binary_sha256']
    manifest=read(HERE/'EXACT_DEPLOYED_SOURCE_MANIFEST.json')
    assert all(sha(Path(a['retained_patched_tree'])/p)==h for p,h in manifest.items())

def test_source_formula_parameters_and_flags():
    a=read(HERE/'STABILIZATION_SOURCE_AUDIT.json');assert a['JIA_ESMAILY_TRADITIONAL_FORM_MATCH']
    assert a['active_parameters']==dict(dt_s=.01,ctM=1.,ctC=36.,rho=1056.,mu=.00345312,nu=3.27e-6,rho_inf=.5,alpha_m=5/6,alpha_f=2/3,gamma=2/3)
    for r in a['references'].values():assert sha(SOURCE/'Code/Source/solver'/r['file'])==r['sha256']
    assert a['active_flags']['vmsFlag'] and not a['active_flags']['Taylor_Hood']

def test_metric_against_solver_cofactor_and_origin():
    x=np.array([[2.,0.,0.],[.3,3.,0.],[.2,.4,4.],[0.,0.,0.]])*1e-4
    J=(x[:3]-x[3]).T
    # Independent cofactor expression matching nn::gnn 3D, NOT another inverse.
    rows=np.stack([np.cross(J[:,1],J[:,2]),np.cross(J[:,2],J[:,0]),np.cross(J[:,0],J[:,1])]);rows/=np.dot(J[:,0],np.cross(J[:,1],J[:,2]))
    assert np.allclose(metric(x),rows.T@rows,rtol=5e-15,atol=1e-8)
    assert np.allclose(metric(x),metric(x[[1,0,2,3]]),rtol=1e-15,atol=1e-8)
    assert not np.allclose(metric(x),metric(x[[3,1,2,0]])) # Origin matters; no arbitrary permutation correction.
    assert np.allclose(metric(np.r_[np.eye(3),np.zeros((1,3))]),np.eye(3))

def test_metric_against_native_output_order():
    comp=read(PREVIOUS/'FIELD_COMPARABILITY.json');result={}
    for name in NAMES:
        mesh=pv.read(paths(name)[0]);native=pv.read(BASE/comp['fields'][name]['raw_payload']['native']);tet=mesh.cells.reshape(-1,5)[:,1:];nt=native.cells.reshape(-1,5)[:,1:]
        assert np.array_equal(np.sort(tet,axis=1),np.sort(nt,axis=1))
        ids=np.unique(np.load(PREVIOUS/f'{name}_accepted_1189.npz')['cell'])
        a=metric(np.asarray(mesh.points,float)[tet[ids]]);b=metric(np.asarray(native.points,float)[nt[ids]])
        assert np.allclose(a,b,rtol=1e-11,atol=1e-4)
        result[name]=float(np.max(np.abs(a-b)/np.max(np.abs(a),axis=(1,2))[:,None,None]))
    put('NATIVE_METRIC_VERIFICATION.json',result)

def test_formula_analytic_scale_and_density_units():
    K=np.array([np.diag([1e8,2e8,3e8])]);u=np.array([[.1,.2,.3]]);v=terms(K,u)[0]
    kt=40000.;ku=.01*1e8+.04*2e8+.09*3e8;ks=36*(MU/RHO)**2*14e16
    assert np.allclose(v[:3],[kt,ku,ks],rtol=2e-15)
    assert np.isclose(v[3]*RHO,1/np.sqrt(kt+ku+ks),rtol=1e-15)
    assert np.isclose(v[3]*v[4]*6e8,1)
    assert np.isclose(v[8]/v[3],1/np.sqrt(1-v[5]))

def test_same_cell_repeatability_and_all_wall_zero_u():
    z=np.load(HERE/'STABILIZATION_TERMS.npz');sel=(z['field']=='M0')&(z['pid']==1189)&(z['cell']==388769)
    K=z['Kxi_m_inv2'][sel];assert np.array_equal(K,np.broadcast_to(K[0],K.shape))
    v=z['values'][sel];u=z['u_m_s'][sel];assert np.array_equal(terms(K,u),v)
    sel=z['kind']=='cell_138584_centroid_gauss4';assert sel.sum()==5
    assert np.all(z['u_m_s'][sel]==0) and np.all(z['values'][sel,1]==0)

def test_frozen_event_positions_and_fraction_accounting():
    z=np.load(HERE/'STABILIZATION_TERMS.npz');assert z['valid'].all();assert len(z['cell'])==271473
    assert np.allclose(z['values'][:,5:8].sum(1),1,atol=5e-16,rtol=0)
    for name in NAMES:
        for pid in PIDS:
            x=np.load(CAUSE/f'original_{pid}.npz')['x_m'];sel=(z['field']==name)&(z['pid']==pid)&(z['kind']=='accepted');assert np.array_equal(x,z['x_m'][sel])

def test_stop_gate_no_new_CFD_or_patch():
    d=read(HERE/'DECISION_GATE.json');assert d['status']=='STOP_C' and d['new_pipe_CFD']==0 and d['new_BraVa_CFD']==0 and not d['new_solver_patch']
    for n in ['02_PIPE_BEFORE_AFTER.png','03_BRAVA_LOCAL_BEFORE_AFTER.png','04_CONTROL_VOLUME_BEFORE_AFTER.png']:assert not (HERE/n).exists()

def test_RT_analytic_constant_linear_divfree_nonconservative():
    from rt_reference import analytic_tests
    r=analytic_tests();assert r['incompatible_boundary_rejected'];assert r['constant']['centroid_max_component_difference']<1e-12
    assert r['linear_divergence_free']['centroid_max_component_difference']>.01 # RT0 cannot reproduce this general linear field.
    assert r['nonconservative']['max_cell_net_flux']<1e-12

def test_RT_ROI_mass_flux_and_boundary_semantics():
    d=read(HERE/'RT_REFERENCE.json');z=np.load(HERE/'RT_REFERENCE.npz')
    assert d['RT_REFERENCE_NOT_PRODUCTION_FIELD'] and d['cells']==468
    assert d['max_cell_Qnet_m3_s']<1e-21 and d['max_wall_flux_m3_s']<1e-21
    assert np.array_equal(np.load(PREVIOUS/'roi_1189.npz')['cells'],z['cells'])
    assert abs(d['old_boundary_Qnet_m3_s'])>1e-10 and d['cut_L1_change_m3_s']>1e-9
    assert d['cannot_preserve_all_original_cut_fluxes']
    for c in d['control_volumes']:assert abs(c['RT_Qnet_m3_s'])<1e-21

def test_all_recorded_inputs_unchanged():
    inputs=read(HERE/'INPUT_SHA256.json');assert all(sha(p)==h for p,h in inputs.items())
    for rel,root in [(CAUSE/'INPUT_SHA256.json',BASE),(PREVIOUS/'OUTPUT_SHA256.json',PREVIOUS)]:
        assert all(sha(root/p)==h for p,h in read(rel).items())

def test_figures_300dpi_and_no_nan_metadata():
    from PIL import Image
    for name in read(HERE/'FIGURE_METADATA.json')['figures']:
        im=Image.open(HERE/name);assert min(im.info['dpi'])>=299.9;assert im.width>3000
