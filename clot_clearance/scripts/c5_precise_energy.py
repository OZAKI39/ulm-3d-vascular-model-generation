"""Roundoff oracle: identical BAQP/spectral/energy code, extended arithmetic.

No constitutive or integration formula is copied. Function code objects are
rebound to a NumPy facade with long-double floating constructors. NumPy's
unsupported long-double linalg is evaluated in float64; its arguments are the
already computed Green strain (not F.T F minus identity). Thus the cancellation
in computing small Green strains and total energy differences is reduced.
This is an energy diagnostic only, never a stationary solver callback.
"""
import types
import numpy as np
from pd_clot.pfpd_candidate import baqp,svk_spectral_split as law
from pd_clot.pfpd_candidate.coupling_stability_audit import AuditProblem


class LinearAlgebra:
    def __getattr__(self,name):
        original=getattr(np.linalg,name)
        def run(a,*args,**kwargs):return original(np.asarray(a,dtype=float),*args,**kwargs)
        return run


class ExtendedNumPy:
    linalg=LinearAlgebra()
    def __getattr__(self,name):
        original=getattr(np,name)
        if name in ['zeros','ones','eye']:
            def run(*args,**kwargs):
                kwargs.setdefault('dtype',np.longdouble);return original(*args,**kwargs)
            return run
        return original


def rebind(fn,scope):return types.FunctionType(fn.__code__,scope,fn.__name__,fn.__defaults__,fn.__closure__)


def make_evaluator():
    facade=ExtendedNumPy();namespace=dict(vars(law));namespace.update(np=facade,float=np.longdouble)
    for name,value in vars(law).items():
        if isinstance(value,types.FunctionType) and value.__module__==law.__name__:namespace[name]=rebind(value,namespace)
    original_getter=baqp.State.energy.fget
    extended_getter=rebind(original_getter,dict(original_getter.__globals__,float=np.longdouble))
    class ExtendedState(baqp.State):energy=property(extended_getter)
    scope=dict(vars(baqp));scope.update(np=facade,float=np.longdouble,State=ExtendedState)
    evaluate=rebind(baqp.evaluate,scope)
    parts=rebind(AuditProblem.parts,dict(AuditProblem.parts.__globals__,float=np.longdouble))
    def energy(p,x,s):
        responses=[]
        def callback(F,mat):
            r=namespace['degraded_svk_split_response'](F,s,mat);responses.append(r);return r.Psi,r.P
        scope['constitutive']=callback
        state=evaluate(p.family,np.asarray(x,dtype=np.longdouble),p.material,g=np.ones(p.ns),h=np.ones(p.ns))
        state.W0=responses[0].Y;state.split_response=responses[0]
        return parts(p,state,np.asarray(s,dtype=np.longdouble))['total_J']/np.longdouble(p.energy_scale)
    return energy


energy=make_evaluator()


def directional_value(p,x,s,v,amplitude):
    xx=np.asarray(x,dtype=np.longdouble).copy();ss=np.asarray(s,dtype=np.longdouble).copy()
    xx[p.free]+=np.longdouble(amplitude)*np.longdouble(p.length_scale)*v[:p.nx].reshape(-1,3).astype(np.longdouble)
    ss+=np.longdouble(amplitude)*v[p.nx:].astype(np.longdouble)
    return energy(p,xx,ss)
