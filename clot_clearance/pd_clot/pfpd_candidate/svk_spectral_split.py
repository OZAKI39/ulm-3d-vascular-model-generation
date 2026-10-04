"""Isolated Green--Lagrange/SVK Miehe-type spectral split.

Algebraic split: PhaseFieldX, commit a9714c122497f8829860efc68d776e4324d48eca,
Element/Phase_Field_Fracture/split_energy_stress_tangent_functions.py.
https://doi.org/10.1016/j.cma.2010.04.011 (split concept only).
The use of Green strain here does not imply equivalence to the small-strain FEM
solver. No AT2/history field, residual stiffness, eigenvalue smoothing or tangent
regularization is introduced. Batched NumPy eigh supplies the spectral routine.
"""
from dataclasses import dataclass
import numpy as np


def green_lagrange(F):
    F=np.asarray(F,dtype=float)
    if F.shape[-2:]!=(3,3) or not np.isfinite(F).all():
        raise ValueError('Finite 3x3 deformation gradient required')
    if np.any(np.linalg.det(F)<=0):raise ValueError('Positive determinant required')
    return .5*(F.swapaxes(-1,-2)@F-np.eye(3))


def spectral_green_split(E):
    E=np.asarray(E,dtype=float)
    if E.shape[-2:]!=(3,3) or not np.isfinite(E).all():raise ValueError('Finite symmetric 3x3 strain required')
    if not np.allclose(E,E.swapaxes(-1,-2),rtol=1e-12,atol=1e-14):raise ValueError('Symmetric Green strain required')
    eig,Q=np.linalg.eigh(E)
    plus=(Q*np.maximum(eig,0)[...,None,:])@Q.swapaxes(-1,-2)
    minus=(Q*np.minimum(eig,0)[...,None,:])@Q.swapaxes(-1,-2)
    return dict(eigenvalues=eig,E_plus=plus,E_minus=minus)


def _parts(E,lame_lambda,mu):
    if not np.isfinite([lame_lambda,mu]).all() or mu<=0 or lame_lambda+2*mu/3<=0:
        raise ValueError('Positive shear and bulk moduli required')
    E=np.asarray(E,float);d=spectral_green_split(E);tr=np.trace(E,axis1=-2,axis2=-1)
    tp=np.maximum(tr,0);tm=np.minimum(tr,0);ep=d['E_plus'];em=d['E_minus']
    d.update(Psi_plus=.5*lame_lambda*tp**2+mu*np.sum(ep*ep,axis=(-1,-2)),
        Psi_minus=.5*lame_lambda*tm**2+mu*np.sum(em*em,axis=(-1,-2)),
        S_plus=lame_lambda*tp[...,None,None]*np.eye(3)+2*mu*ep,
        S_minus=lame_lambda*tm[...,None,None]*np.eye(3)+2*mu*em)
    return d


def svk_split_energy(E,lame_lambda,mu):
    d=_parts(E,lame_lambda,mu);return d['Psi_plus'],d['Psi_minus']


def svk_split_second_pk(E,lame_lambda,mu):
    d=_parts(E,lame_lambda,mu);return d['S_plus'],d['S_minus']


@dataclass(frozen=True)
class SplitResponse:
    E: np.ndarray
    eigenvalues: np.ndarray
    E_plus: np.ndarray
    E_minus: np.ndarray
    Psi_plus: np.ndarray
    Psi_minus: np.ndarray
    Psi_total: np.ndarray
    S_plus: np.ndarray
    S_minus: np.ndarray
    S_total: np.ndarray
    P_total: np.ndarray


def svk_split_response(F,material):
    if material.model!='svk':raise ValueError('Spectral candidate is SVK only')
    F=np.asarray(F,float);E=green_lagrange(F);d=_parts(E,material.lame,material.mu)
    S=d['S_plus']+d['S_minus']
    return SplitResponse(E=E,**d,Psi_total=d['Psi_plus']+d['Psi_minus'],S_total=S,P_total=F@S)


@dataclass(frozen=True)
class DegradedSplitResponse:
    split: SplitResponse
    g: np.ndarray
    Psi: np.ndarray
    S: np.ndarray
    P: np.ndarray
    Y: np.ndarray


def degraded_svk_split_response(F,s,material):
    d=svk_split_response(F,material);s=np.asarray(s,float)
    if not np.isfinite(s).all() or np.any((s<0)|(s>=1)):raise ValueError('Early phase must satisfy 0 <= s < 1')
    s=np.broadcast_to(s,d.Psi_total.shape);g=(1-s)**2
    S=g[...,None,None]*d.S_plus+d.S_minus
    return DegradedSplitResponse(d,g,g*d.Psi_plus+d.Psi_minus,S,np.asarray(F)@S,d.Psi_plus)
