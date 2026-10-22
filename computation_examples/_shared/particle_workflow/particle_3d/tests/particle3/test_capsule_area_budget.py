import numpy as np
from scipy.integrate import quad
from particle_3d.rbc_capillary_surrogate import oblate_area,capsule_area,capsule_length,area_feasible_interval
from particle_3d.particle_shapes import Capsule


def test_oblate_area_independent_surface_integral_and_sphere_limit(p3_geometries):
    assert oblate_area(2e-6,2e-6)==4*np.pi*(2e-6)**2
    for g in p3_geometries:
        a,c=float(g.a_m),float(g.c_m)
        ratio=c/a
        integral=quad(lambda t:np.sin(t)*np.sqrt(ratio**2*np.sin(t)**2+np.cos(t)**2),0,np.pi,epsabs=1e-12)[0]*2*np.pi*a*a
        assert abs(oblate_area(a,c)/integral-1)<1e-12


def test_capsule_volume_area_and_derived_radius_bounds(p3_geometries):
    for g in p3_geometries:
        lo,hi,budget=area_feasible_interval(g)
        assert 0<lo<hi and capsule_area(g.volume_m3,lo)<=budget
        assert abs(capsule_area(g.volume_m3,lo)/budget-1)<2e-12
        for r in np.linspace(lo,hi,9):
            cap=Capsule([0,0,0],[0,0,1],r,capsule_length(g.volume_m3,r))
            assert abs(cap.volume_m3/g.volume_m3-1)<1e-14
            assert abs(cap.area_m2/capsule_area(g.volume_m3,r)-1)<1e-14
            assert cap.area_m2<=budget*(1+1e-14)
