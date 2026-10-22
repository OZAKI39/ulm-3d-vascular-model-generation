"""Private compiled dependency in first_handoff_event only; unchanged source.

The original function's inner relative import receives a proxy ONLY in a
private function binding. No sys.modules/global monkeypatching is performed.
"""
from pathlib import Path
from types import MethodType, SimpleNamespace
import builtins, ctypes, hashlib, subprocess
import numpy as np

HERE=Path(__file__).resolve().parent
_LIB=None


def native_closest(point,triangles):
    global _LIB
    if _LIB is None:
        source=HERE/'closest_native.cpp'
        digest=hashlib.sha256(source.read_bytes()).hexdigest()[:16]
        path=HERE/f'closest_native_{digest}.so'
        if not path.exists():
            subprocess.run(['g++','-O2','-ffp-contract=off','-fno-fast-math','-shared','-fPIC',
                            str(source),'-o',str(path)],check=True)
        _LIB=ctypes.CDLL(str(path))
        ptr=ctypes.POINTER(ctypes.c_double)
        _LIB.closest.argtypes=[ptr,ptr,ctypes.c_int64,ptr,ptr]
        _LIB.closest.restype=None
    p=np.ascontiguousarray(point,dtype=np.float64)
    t=np.ascontiguousarray(triangles,dtype=np.float64)
    assert p.shape==(3,) and t.ndim==3 and t.shape[1:]==(3,3)
    out=np.empty((len(t),3));bary=np.empty_like(out)
    ptr=ctypes.POINTER(ctypes.c_double)
    _LIB.closest(p.ctypes.data_as(ptr),t.ctypes.data_as(ptr),len(t),out.ctypes.data_as(ptr),bary.ctypes.data_as(ptr))
    return out,bary


def import_binding(function,module_name,attributes):
    from particle_3d.particle6_stepper import bind_query_dependency
    def imports(name,globals=None,locals=None,fromlist=(),level=0):
        if name==module_name and level==1:
            actual=builtins.__import__(name,globals,locals,fromlist,level)
            return SimpleNamespace(**{item:attributes.get(item,getattr(actual,item)) for item in fromlist})
        return builtins.__import__(name,globals,locals,fromlist,level)
    return bind_query_dependency(function,{'__builtins__':dict(vars(builtins),__import__=imports)})


def accelerated_job(job):
    from particle_3d.particle9a_motion import Particle9AStepper
    from particle_3d.particle6_stepper import bind_query_dependency
    class CompiledGeometryStepper(Particle9AStepper):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs)
            method=self.advance_cached.__func__
            trial=method.__globals__['v1_trial']
            original=trial.__globals__['first_handoff_event']
            handoff=import_binding(original,'convex_triangle',{'triangle_closest_many':native_closest})
            trial=bind_query_dependency(trial,{'first_handoff_event':handoff})
            self.advance_cached=MethodType(bind_query_dependency(method,{'v1_trial':trial}),self)
    return import_binding(job,'particle9a_motion',{'Particle9AStepper':CompiledGeometryStepper})
