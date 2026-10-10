"""One material card. UFL energy is the source of stress and consistent tangent."""
import numpy as np
CARD=dict(material_id='SVK_DEMO_1000PA_LAMBDA0',mu_Pa=1000.,lambda_Pa=0.,provenance='DEMONSTRATION_PARAMETERS_NOT_BIOCALIBRATED',energy='lambda/2 tr(E)^2 + mu E:E',Piola='d psi0 / dF',consistent_tangent='UFL derivative of the same weak energy gradient',Psi_plus='Green-strain positive eigenvalue diagnostic only; no phase/damage driving')
def ufl_energy(F):
    import ufl
    E=.5*(F.T*F-ufl.Identity(3));return 1000*ufl.inner(E,E)
def fields(F):
    E=.5*(np.swapaxes(F,-1,-2)@F-np.eye(3));psi=1000*np.sum(E*E,axis=(-2,-1));P=F@(2000*E);J=np.linalg.det(F)
    plus=1000*np.sum(np.maximum(np.linalg.eigvalsh(E),0)**2,axis=-1)
    sigma=(P@np.swapaxes(F,-1,-2))/J[...,None,None]
    return dict(F_minus_I=np.linalg.norm(F-np.eye(3),axis=(-2,-1)),J=J,J_minus_1=J-1,P=P,sigma=sigma,psi0=psi,Psi_plus=plus)
