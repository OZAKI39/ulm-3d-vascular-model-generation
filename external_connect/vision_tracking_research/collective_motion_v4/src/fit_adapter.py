"""Reuse frozen F10U fit body; supply its exact linear Jacobian to SciPy.

No global monkeypatch: a private function object retains the original bytecode,
gates and logm conversion, with only the optimizer-call dependency substituted.
"""
import types,inspect
import numpy as np
from scipy.optimize import least_squares
from xcommon import fit_affine as historical_fit
def linear_jacobian_solver(fun,x0,**kwargs):
 values=inspect.getclosurevars(fun).nonlocals
 X=values['X'];assert X.ndim==2 and X.shape[1]==3 and np.asarray(x0).size==6
 # Preserve the frozen numerical solver whenever Huber's linear regime is
 # active. Sparse near-flat robust objectives can have non-unique optima;
 # an exact Jacobian then follows a different path despite identical gates.
 if kwargs.get('loss')!='linear' and np.max(np.abs(fun(x0)))>=kwargs.get('f_scale',1.):
  return least_squares(fun,x0,**kwargs)
 jac=np.kron(X,np.eye(2))
 # SciPy robust losses scale the returned Jacobian in place. A fresh array is
 # required at EVERY evaluation; returning the cached object changes the model.
 return least_squares(fun,x0,jac=lambda _:jac.copy(),**kwargs)
fit_affine=types.FunctionType(historical_fit.__code__,{**historical_fit.__globals__,'least_squares':linear_jacobian_solver},name='fit_affine_exact_jacobian',argdefs=historical_fit.__defaults__,closure=historical_fit.__closure__)
