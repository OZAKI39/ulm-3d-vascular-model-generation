"""Tests of THIS diagnostic only; these do not certify the CFD solution."""
from common import *
import numpy.testing as nt
from particle_3d.field import FrozenFEMField
from sample_fields import triangle_rule,strict_ids
from prepare_coordinates import cv_surface
from bounds_index import BoundsIndex
from analyze import quantities,pair_summary

def manufactured():
    points=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)*1e-3
    G=np.array([[1,2,-3],[-4,5,6],[7,-8,9]],float)
    u=points@G.T+np.array([.2,-.1,.3])
    return points,G,u,FrozenFEMField(points,np.array([[0,1,2,3]]),u,np.zeros(4))

def test_linear_vector_P1_gradient_curl_units():
    p,G,u,f=manufactured();x=p.mean(0);s=f.sample(x)
    nt.assert_allclose(s.velocity_m_s,u.mean(0),rtol=0,atol=1e-15)
    nt.assert_allclose(s.velocity_gradient_s_inv,G,rtol=0,atol=1e-12)
    nt.assert_allclose(s.vorticity_s_inv,[-14,-10,-6],rtol=0,atol=1e-12)
    assert np.linalg.norm(s.velocity_m_s)*1000>100 # m/s -> mm/s

def test_rigid_inverse_transform_vector_tensor():
    p,G,u,f=manufactured();R=np.array([[0.,-1,0],[1,0,0],[0,0,1]]);translation=np.array([.03,-.02,.01])
    rotated=FrozenFEMField(p@R.T+translation,f.tetra,u@R.T,np.zeros(4))
    sample=rotated.sample(p.mean(0)@R.T+translation)
    nt.assert_allclose(sample.velocity_m_s@R,u.mean(0),atol=1e-14)
    nt.assert_allclose(R.T@sample.velocity_gradient_s_inv@R,G,atol=1e-11)
    nt.assert_allclose((rotated.points_m-translation)@R,p,atol=1e-17)

def test_analytic_closed_flux_and_refinement():
    p,G,u,f=manufactured();cv=cv_surface(p,f.tetra,[0],set())
    expected=np.trace(G)*f.geometry.volume_m3[0]
    for level in [0,1,2,3,4,5]:
        bary=triangle_rule(level);flux=0
        for tri,n,area in zip(cv['triangles_m'],cv['normal_out'],cv['area_m2']):
            x=bary@tri;val=x@G.T+np.array([.2,-.1,.3]);flux+=float((val@n).mean()*area)
        nt.assert_allclose(flux,expected,rtol=0,atol=1e-21)

def test_piecewise_face_quadrature_refines_to_exact_integral():
    a=.231;exact=(1-a)**3/6;errors=[]
    for level in range(7):
        bary=triangle_rule(level);q=.5*np.maximum(bary[:,1]-a,0).mean();errors.append(abs(q-exact))
    assert errors[-1]<errors[0]/1000
    assert all(b<=a+1e-16 for a,b in zip(errors,errors[1:]))

def test_bounds_candidates_equal_bruteforce():
    rng=np.random.default_rng(173);p=rng.random((400,3))*.001+np.array([.03,-.04,.005]);tet=np.arange(400).reshape(-1,4)
    index=BoundsIndex(p,tet)
    for k in range(80):
        lo=rng.random(3)*.001+np.array([.03,-.04,.005]);hi=lo+(rng.random(3)*.0004 if k%2 else 1e-12)
        b=np.c_[lo,hi].ravel();expected=np.flatnonzero(np.all(index.lower<=hi,axis=1)&np.all(index.upper>=lo,axis=1))
        nt.assert_array_equal(index.candidates(b),expected)

def test_outside_is_not_extrapolated():
    p,G,u,f=manufactured();x=np.array([[.0001,.0001,.0001],[.002,0,0]])
    ids=strict_ids(f,x,np.array([0,-1]));nt.assert_array_equal(ids,[0,-1])
    assert not f.sample(x[1]).inside_lumen and np.isnan(f.sample(x[1]).velocity_m_s).all()

def test_missing_denominator_and_self_comparison():
    data=np.zeros((3,25));data[:,:3]=[1,2,3];n=np.tile([1.,0,0],(3,1));q=quantities(data,n)
    result=pair_summary(q,q,np.ones(3),np.ones(3,bool),1e-14)
    assert result['metrics']['u']['RMS_difference']==0
    assert all(v['RMS_difference']==0 for v in result['metrics'].values())
    other={k:v.copy() for k,v in q.items()};other['u'][1]=np.nan
    result=pair_summary(q,other,np.array([1.,2.,3.]),np.ones(3,bool),1e-14)
    assert result['total_points']==3 and result['valid_points']==2
    assert result['total_weight']==6 and result['valid_weight']==4 and result['missing_indices']==[1]

def test_real_all_segments_weight_and_822_identity():
    m=read(HERE/'FIELD_METRICS.json')
    for pid in PIDS:
        p=np.load(HERE/f'path_{pid}.npz')
        for order in [1,2,4]:
            ix=p['quadrature_order']==order
            nt.assert_allclose(p['arc_weight_m'][ix].sum(),p['original_arc_m'][-1],rtol=1e-13)
            nt.assert_allclose(p['time_weight_s'][ix].sum(),np.ptp(p['original_time_s']),rtol=1e-13)
        for name in NAMES:
            z=np.load(HERE/f'{name}_accepted_{pid}.npz')
            assert len(z['valid'])==len(p['original_time_s'])
    for name,r in m['822']['spatial_rise'].items():
        assert np.linalg.norm(r['endpoint_identity_error_m_s'])<1e-10
    assert len(m['822']['spatial_rise']['M0']['own_cells'])==4

def test_original_cell_flux_identity_and_residual_distinct():
    values=read(HERE/'M0_BALANCE.json')
    for row in values:assert abs(row['M0_identity_error_m3_s'])<1e-20
    negative=next(r for r in values if r['event']==1189 and r['layer']==0)
    assert negative['M0_exact_net_m3_s']<-5e-10

def test_figure_data_and_png_quality():
    from PIL import Image
    meta=read(HERE/'FIGURE_METADATA.json');assert len(meta)==5
    for name in meta:
        with Image.open(HERE/name) as im:
            assert min(im.size)>2000 and min(im.info.get('dpi',(0,0)))>=299.5
    peak=max(np.nanmax(norm(np.load(HERE/f'{name}_plane_1189.npz')['flow'][:,:3])) for name in NAMES)*1000
    assert meta['02_NEAR_WALL_FLOW.png']['speed_limits_mm_s'][1]==peak
    rows=list(__import__('csv').DictReader((HERE/'representative_values.csv').open()))
    metrics=read(HERE/'FIELD_METRICS.json')
    for row in rows:
        saved=metrics[row['event']]['representatives'][row['field']][row['endpoint']]['speed_mm_s']
        assert float(row['speed_mm_s'])==saved

def test_all_registered_inputs_unchanged():
    for relative,expected in read(HERE/'INPUT_SHA256.json').items():assert sha(BASE/relative)==expected,relative
    # Previous task's immutable original mesh, field, tracks and source ledger.
    for relative,expected in read(OLD/'INPUT_SHA256.json').items():assert sha(BASE/relative)==expected,relative

def test_spatial_gain_reconstruction():
    features=read(HERE/'SPATIAL_FEATURES.json');metrics=read(HERE/'FIELD_METRICS.json')
    for pid in map(str,PIDS):
        for name in NAMES:
            actual=features[pid]['fields'][name]['net_scalar_gain_mm_s']
            expected=metrics[pid]['representatives'][name]['scalar_speed_increase_mm_s']
            assert abs(actual-expected)<1e-8
