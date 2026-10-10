"""Physical boundary adapters consume the normalized provider payload."""
import numpy as np
from .api import RequestError,load_provider

def displacement_schedule(request,hemisphere=False):
 payload=load_provider(request)
 if payload['load_control']!='DISPLACEMENT_CONTROL':raise RequestError('CONTROL_MODE_NOT_AVAILABLE_FOR_THIS_DRIVER','The two demonstration drivers prescribe kinematics; use traction_external for a traction-controlled backend')
 d=np.asarray(payload['direction'],float);d/=np.linalg.norm(d);expected=np.array([1,0,0] if hemisphere else [0,1,0])
 if not np.allclose(d,expected,rtol=0,atol=1e-12):raise RequestError('UNSUPPORTED_BOUNDARY_DIRECTION','Mode I pulls along y; the hemisphere wall stretches along x')
 allowed={'HEMISPHERE_BOTTOM_DISK','Bonded bottom disk; ux=(lambda-1)*X, uy=uz=0'} if hemisphere else {'MODE_I_TOP_GRIP','Top grip y=5 mm; y=0 ligament and z=0 symmetry'}
 if payload['region'] not in allowed:raise RequestError('UNSUPPORTED_BOUNDARY_REGION','Select a supported physical boundary region')
 kin=payload['prescribed_kinematics'];key='stretch_ratio' if hemisphere else 'max_displacement_m'
 if key not in kin:raise RequestError('INVALID_KINEMATICS',f'This geometry requires {key}')
 maximum=kin[key]-1 if hemisphere else kin[key]/.01
 if maximum>(1. if hemisphere else 1.5):raise RequestError('DECLARED_LIMIT_EXCEEDED','Request exceeds the frozen demonstration envelope')
 return dict(maximum_parameter=maximum,max_steps=kin['max_steps'],nominal_loading_path=kin['loading_path'],adaptive_refinement=True,provider_payload=payload)

def traction_external(domain,V,facet_tags,tag,request,amplitude):
 """Reference dead traction on an explicitly selected physical facet tag.

 Direction is normalized in global coordinates, then projected into the surface
 tangent plane. Gaussian width is in reference metres. No dynamic inference.
 """
 import ufl
 from dolfinx import fem
 payload=load_provider(request)
 if payload['load_control']!='TRACTION_CONTROL':raise RequestError('INVALID_REQUEST','Traction payload required')
 z=payload['z_load']
 if not np.isfinite(amplitude) or amplitude<0 or amplitude>z['max_amplitude']:raise RequestError('DECLARED_LIMIT_EXCEEDED','Traction amplitude is outside the declared path')
 X=ufl.SpatialCoordinate(domain);N=ufl.FacetNormal(domain);dv=np.asarray(payload['direction'],float);dv/=np.linalg.norm(dv);direction=ufl.as_vector(dv)
 phi=ufl.exp(-sum((X[i]-z['center_m'][i])**2 for i in range(3))/(2*z['width_m']**2))
 pressure=fem.Constant(domain,np.float64(amplitude*z['p_Pa']));shear=fem.Constant(domain,np.float64(amplitude*z['tau_Pa']))
 T=phi*(-pressure*N+shear*(direction-ufl.dot(direction,N)*N))
 ds=ufl.Measure('ds',domain=domain,subdomain_data=facet_tags,metadata={'quadrature_degree':4})
 return ufl.dot(T,ufl.TestFunction(V))*ds(tag)
