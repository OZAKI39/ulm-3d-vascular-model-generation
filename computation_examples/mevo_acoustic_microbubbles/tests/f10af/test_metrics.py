import numpy as np
from scripts.straight_tube_exp1.f10af.metrics import spatial,rotation
def test_rigid_rotation_translation_and_factor_two():
    th=np.arange(64)*2*np.pi/64;x=np.c_[1e-4*np.cos(th),1e-4*np.sin(th),np.zeros(64)];v=np.c_[-3*x[:,1],3*x[:,0],np.zeros(64)]+[.1,-.2,0]
    r=rotation(x,v);assert abs(r['Omega_z_s']-3)<1e-10;assert abs(r['F10Y_omega_affine_s']-6)<1e-10;assert r['rigid_fraction']>1-1e-12
    assert abs(rotation(x,np.tile([.1,.2,0],(64,1)))['Omega_z_s'])<1e-10
def test_shear_not_rigid_rotation():
    rng=np.random.default_rng(42);x=rng.uniform(-.0004,.0004,(128,3));v=np.c_[3*x[:,1],np.zeros((128,2))];r=rotation(x,v);assert r['rigid_fraction']<.8
def test_random_cloud_not_pattern():
    rng=np.random.default_rng(71)
    for _ in range(20):
        x=rng.uniform([-.0005,-.0004,-.00002],[.0005,.0004,.00002],(64,3));r=spatial(x);assert not r['stripe_descriptor_pass'];assert not r['lattice_descriptor_pass'];assert not r['ring_descriptor_pass']
def test_synthetic_stripes_analysis_only():
    rng=np.random.default_rng(4);x=np.c_[np.repeat(np.arange(-2,3)*.0002,16)+rng.normal(0,2e-6,80),rng.uniform(-.0004,.0004,80),np.zeros(80)];r=spatial(x);assert r['stripe_descriptor_pass'];assert abs(r['stripe_spacing_m']-.0002)<1e-5
def test_synthetic_grid_and_ring_analysis_only():
    a=np.arange(-2,3)*.0002;xx,yy=np.meshgrid(a,a);x=np.c_[xx.ravel(),yy.ravel(),np.zeros(25)];assert spatial(x)['lattice_descriptor_pass']
    th=np.arange(64)*2*np.pi/64;x=np.c_[.0003*np.cos(th),.0003*np.sin(th),np.zeros(64)];assert spatial(x)['ring_descriptor_pass']
