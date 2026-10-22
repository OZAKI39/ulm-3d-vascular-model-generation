import numpy as np
from scripts.straight_tube_exp1.f10ae.metrics import target_metrics,recovery
from scripts.straight_tube_exp1.f10ae.schedules import schedule
from scripts.straight_tube_exp1.f10ab.metrics import measure
def cloud():
    rng=np.random.default_rng(41)
    x=rng.normal(size=(64,3))*[9e-6,3e-6,2e-6]
    x-=x.mean(0)
    # Rotate to principal coordinates for an exact relative-aspect test.
    _,v=np.linalg.eigh(x.T@x);return x@v[:,::-1]
def test_T1_shape_COM_and_numeric_are_independent():
    x=cloud();y=x*[1.03,1,1]
    m=target_metrics(x,y+[1e-6,0,0],False)
    assert m['SHAPE_TARGET_PASS'] and not m['COM_TARGET_PASS'] and not m['NUMERICAL_PASS'] and not m['T1_DISCOVERY_PASS']
def test_aspect_gain_does_not_imply_absolute_extension():
    x=cloud();m=target_metrics(x,x*[.99,.98,.96])
    assert m['aspect_gain']>0 and m['axis_class']=='ANISOTROPIC_CONTRACTION' and not m['raw_extension_above_floor']
def test_recovery_fixed_ID_distinguishes_same_Rg_rearrangement():
    x=cloud();a=x*.95;b=x[::-1]
    r=recovery(x,a,b)
    assert r['raw_recovery_pass']['Rg'] and not r['raw_recovery_pass']['Kabsch_m']
def test_off_is_not_restoration():
    x=cloud();a=x*.95;r=recovery(x,a,a)
    assert not any(r['raw_recovery_pass'].values())
def test_two_segment_time_grid():
    a=dict(amplitudes=[0,0,700,700],phases=[0,0,0,np.pi],duration=.004);b=dict(a,segments=[dict(amplitudes=a['amplitudes'],phases=a['phases'],duration=.002)]*2)
    aa=list(schedule(a).steps(1e-6));bb=list(schedule(b).steps(1e-6))
    assert len(aa)==len(bb)==4000
    assert np.max(abs(np.array([t for t,h,s in aa])-np.array([t for t,h,s in bb])))<1e-15
    assert sum(h for t,h,s in bb)==sum(h for t,h,s in aa)
    assert sum(t==.002 for t,h,s in bb)==1
def test_rigid_translation_rotation_preserve_shape():
    x=cloud();q,_=np.linalg.qr(np.random.default_rng(9).normal(size=(3,3)))
    m=measure(x,x@q+[1e-6,0,0]);assert m['kabsch_rmsd_3d_m']<1e-18
    assert abs(m['similarity_scale_3d']-1)<1e-12
