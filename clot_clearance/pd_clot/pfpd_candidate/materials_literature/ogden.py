"""Isochoric one-term Ogden: analytic spectral matrix-function derivatives.

Eigendecomposition is used only to evaluate a matrix function. Its Frechet
 derivative uses divided differences with the analytic repeated-root limit;
 no differentiation through eigenvectors, jitter, clipping or abs(det F).
"""
import numpy as np
MU=1040.;ALPHA=7.21;BULK=100000.
CARD=dict(material_id='OGDEN1_BOVINE_LITERATURE_CANDIDATE_V1',mu_Pa=MU,alpha=ALPHA,K_Pa=BULK,Gc_J_per_m2=.96,source_DOI='10.1007/s13239-024-00733-2',elastic_parameter_origin='BOVINE_WHOLE_BLOOD_INVERSE_IDENTIFICATION_REPORTED',bulk_origin='LITERATURE_ASSUMPTION_NO_VOLUMETRIC_DATA',Gc_origin='LITERATURE_MODEL_DEPENDENT_INVERSE_FRACTURE_PARAMETER',calibration_scope='IMPORTED_NOT_RECALIBRATED_FOR_THIS_GEOMETRY_OR_HUMAN_CLOT',adaptation='EXPLICIT_STANDARD_ISOCHORIC_STRETCHES; printed paper does not explicitly show bars',Psi_plus='NOT_DEFINED_FOR_THIS_MATERIAL',units=dict(F='1',W0='Pa',W_iso='Pa',W_vol='Pa',P='Pa',tangent='Pa'))
def evaluate(F,tangent=False):
    F=np.asarray(F,dtype=np.float64)
    if F.shape[-2:]!=(3,3) or not np.isfinite(F).all():raise ValueError('Finite 3x3 F required')
    J=np.linalg.det(F)
    if np.any(J<=0):raise ValueError('INVALID_DEFORMATION: J must be positive, never replaced by abs(J)')
    C=F.swapaxes(-1,-2)@F;lam,Q=np.linalg.eigh(C)
    if np.any(lam<=0):raise ValueError('INVALID_DEFORMATION: C must be positive definite, no eigenvalue clipping')
    beta=ALPHA/2;r=beta-1;t=np.sum(lam**beta,axis=-1);b=J**(-ALPHA/3);S=(Q*(lam**r)[...,None,:])@Q.swapaxes(-1,-2);G=np.linalg.inv(F).swapaxes(-1,-2)
    dev=F@S-t[...,None,None]*G/3;Piso=(2*MU/ALPHA)*b[...,None,None]*dev;Pvol=BULK*(J*(J-1))[...,None,None]*G
    Wi=(2*MU/ALPHA**2)*(b*t-3);Wv=BULK/2*(J-1)**2
    out=dict(W0=Wi+Wv,W_iso=Wi,W_vol=Wv,P=Piso+Pvol,P_iso=Piso,P_vol=Pvol,J=J)
    if not tangent:return out
    aa=lam[..., :,None];bb=lam[...,None,:];delta=(aa-bb)/bb;close=np.abs(delta)<1e-7
    ratio=np.empty_like(delta);np.divide(np.expm1(r*np.log1p(delta)),delta,out=ratio,where=~close)
    L=np.where(close,r*((aa+bb)/2)**(r-1),bb**(r-1)*ratio)
    A=np.empty(F.shape[:-2]+(3,3,3,3));Qt=Q.swapaxes(-1,-2)
    for k in range(3):
        for l in range(3):
            H=np.zeros((3,3));H[k,l]=1.;dC=F.swapaxes(-1,-2)@H+H.T@F;dS=Q@(L*(Qt@dC@Q))@Qt
            dt=beta*np.sum(S*dC,axis=(-2,-1));q=G[...,k,l];dG=-G@H.T@G
            ddev=H@S+F@dS-dt[...,None,None]*G/3-t[...,None,None]*dG/3
            dPi=(2*MU/ALPHA)*b[...,None,None]*(ddev-(ALPHA/3)*q[...,None,None]*dev)
            dPv=BULK*((J*(2*J-1)*q)[...,None,None]*G+(J*(J-1))[...,None,None]*dG)
            A[...,k,l]=dPi+dPv
    out['A']=A;return out

def batch(F,tangent=False,block=4096):
    a=np.asarray(F);shape=a.shape[:-2];flat=a.reshape(-1,3,3);output={}
    for lo in range(0,len(flat),block):
        part=evaluate(flat[lo:lo+block],tangent)
        for k,v in part.items():
            if k not in output:output[k]=np.empty((len(flat),)+v.shape[1:])
            output[k][lo:lo+len(v)]=v
    return {k:v.reshape(shape+v.shape[1:]) for k,v in output.items()}
