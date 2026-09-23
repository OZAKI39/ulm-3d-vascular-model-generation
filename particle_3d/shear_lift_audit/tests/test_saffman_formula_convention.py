import numpy as np
from audit_math import candidate

def test_radius_diameter_and_primary_coefficient_agree():
    a=1.3e-6;mu=.00345312;rho=1056.;s=2e-5;g=333.
    by_diameter=1.615*mu*(2*a)**2*s*np.sqrt(g/(mu/rho))
    assert np.isclose(candidate(a,s,g,mu,rho),by_diameter,rtol=2e-15,atol=0)
    # Independent conversion from Shi et al. Eq.(7), not the same code expression.
    d=2*a;Re=s*d/(mu/rho);Sr=g*d/s
    coefficient=18/np.pi**2*np.sqrt(Sr/Re)*2.254
    converted=coefficient*np.pi*d*d*rho*s*s/8
    assert np.isclose(candidate(a,s,g,mu,rho),converted,rtol=.001,atol=0)

def test_simple_shear_direction_sign():
    slip=np.array([1.,0,0]);curl=np.array([0.,0,-2.])
    assert np.array_equal(np.cross(slip,curl),[0,2,0])
