"""Optional compiled CPU sweep with the original sequential contact arithmetic.

No tree, BVH, force law or NNLS implementation is replaced. The original Python
reference stays intact. This candidate is opt-in and requires replay regression.
"""
import ctypes
import hashlib
from pathlib import Path
import shutil
import subprocess

import numpy as np

from mevo_mb.bubbles import contact


# Kept inside the Python source so existing source snapshots include the exact C.
SOURCE = r"""
#include <math.h>
#include <stdint.h>
typedef double (*ddot64)(int64_t, const double *, int64_t, const double *, int64_t);
int ordered_contact_sweep(double *x, const double *r,
                          const int64_t *pairs, int64_t count, ddot64 dot) {
    for (int64_t k=0; k<count; ++k) {
        int64_t i=pairs[2*k], j=pairs[2*k+1];
        double delta[3];
        for (int a=0; a<3; ++a) delta[a]=x[3*j+a]-x[3*i+a];
        double squared = dot ? dot(3,delta,1,delta,1) :
            (delta[0]*delta[0]+delta[1]*delta[1])+delta[2]*delta[2];
        double d=sqrt(squared);
        if (d==0.0) return 1;
        double amount=fmax(0.0,(r[i]+r[j])-d);
        double fraction=r[j]/(r[i]+r[j]);
        for (int a=0; a<3; ++a) x[3*i+a]-=((fraction*amount)*delta[a])/d;
        for (int a=0; a<3; ++a) x[3*j+a]+=(((1.0-fraction)*amount)*delta[a])/d;
    }
    return 0;
}
"""
FLAGS = ["-std=c99", "-O3", "-fno-fast-math", "-ffp-contract=off", "-fPIC", "-shared"]


def numpy_norm_implementation():
    """Reuse NumPy's exact ILP64 BLAS when its scalar norm calls that BLAS.

    NumPy releases differ here. Unknown implementations fail closed; this is
    not permission to substitute a merely close norm. Full replay remains a gate.
    """
    rng=np.random.default_rng(20261005)
    values=np.ascontiguousarray(rng.normal(size=(1000,3))*10.**rng.uniform(-10,1,(1000,1)))
    expected=np.array([np.linalg.norm(row) for row in values])
    scalar=np.sqrt((values[:,0]*values[:,0]+values[:,1]*values[:,1])+values[:,2]*values[:,2])
    if np.array_equal(expected,scalar):
        return None,None,dict(mode="unfused_ordered_three_term_sum",norm_probe_count=len(values))
    libraries=sorted((Path(np.__file__).parent.parent/"numpy.libs").glob("*openblas64*.so"))
    for path in libraries:
        library=ctypes.CDLL(str(path))
        function=getattr(library,"scipy_cblas_ddot64_",None)
        if function is None:
            continue
        function.argtypes=[ctypes.c_int64,np.ctypeslib.ndpointer(dtype=np.float64,flags="C_CONTIGUOUS"),
            ctypes.c_int64,np.ctypeslib.ndpointer(dtype=np.float64,flags="C_CONTIGUOUS"),ctypes.c_int64]
        function.restype=ctypes.c_double
        actual=np.array([np.sqrt(function(3,row,1,row,1)) for row in values])
        if np.array_equal(expected,actual):
            return library,ctypes.cast(function,ctypes.c_void_p),dict(mode="numpy_existing_scipy_cblas_ddot64_",
                library=str(path),library_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),norm_probe_count=len(values),
                numpy_version=np.__version__)
    raise RuntimeError("No bitwise match for NumPy scalar norm; compiled contact sweep unavailable")


class OrderedContactSweep:
    def __init__(self, cache_directory):
        compiler=shutil.which("cc")
        if compiler is None:
            raise RuntimeError("Optional ordered contact sweep needs a C compiler; Python remains the default")
        version=subprocess.check_output([compiler,"--version"],text=True)
        digest=hashlib.sha256((SOURCE+str(FLAGS)+version).encode()).hexdigest()
        folder=Path(cache_directory)/digest
        folder.mkdir(parents=True,exist_ok=True)
        source=folder/"ordered_contact_sweep.c"
        library=folder/"ordered_contact_sweep.so"
        source.write_text(SOURCE)
        command=[compiler,*FLAGS,str(source),"-o",str(library),"-lm"]
        if not library.exists():
            # Callers initialize candidates serially before any timed replay.
            subprocess.run(command,check=True,capture_output=True,text=True)
        self.library=ctypes.CDLL(str(library))
        self.function=self.library.ordered_contact_sweep
        self.norm_library,self.norm_pointer,norm_identity=numpy_norm_implementation()
        self.function.argtypes=[np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(dtype=np.int64,ndim=2,flags="C_CONTIGUOUS"),ctypes.c_int64,ctypes.c_void_p]
        self.function.restype=ctypes.c_int
        self.identity=dict(compiler=compiler,compiler_version=version,flags=FLAGS,command=command,
            source_sha256=hashlib.sha256(SOURCE.encode()).hexdigest(),
            library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
            numpy_norm=norm_identity,
            floating_point="float64, ordered pairs, no fast-math or fused multiply-add in contact updates; original NumPy norm implementation reused")

    def sweep(self,x,radii,pairs):
        if x.dtype!=np.float64 or not x.flags.c_contiguous or not x.flags.writeable or x.shape[1]!=3:
            raise ValueError("Ordered sweep requires writable contiguous float64 N by 3 coordinates")
        radii=np.ascontiguousarray(radii,dtype=np.float64)
        pairs=np.ascontiguousarray(pairs,dtype=np.int64)
        if len(radii)!=len(x) or pairs.shape!=(len(pairs),2):
            raise ValueError("Invalid radii or pairs shape")
        if len(pairs) and (pairs.min()<0 or pairs.max()>=len(x)):
            raise ValueError("Contact pair index out of bounds")
        if self.function(x,radii,pairs,len(pairs),self.norm_pointer):
            raise RuntimeError("Coincident centers are invalid; reduce timestep")

    def project(self,positions_m,radii_m,cylinder,tolerance_m,max_iterations,candidate_backend="brute",timings=None):
        # Identical control flow to contact.project_contacts; only its ordered
        # per-pair arithmetic moves across the Python/C boundary.
        x,wall_hit=cylinder.project(positions_m,radii_m)
        events=set()
        for iteration in range(max_iterations):
            pairs=contact.contact_candidates(x,radii_m,0.,candidate_backend,timings)
            if len(pairs):
                i,j=pairs.T
                pairs=pairs[radii_m[i]+radii_m[j]-np.linalg.norm(x[j]-x[i],axis=1)>tolerance_m]
            if not len(pairs):
                return x,wall_hit,sorted(events),iteration
            self.sweep(x,radii_m,pairs)
            events.update(map(tuple,pairs.tolist()))
            x,hit=cylinder.project(x,radii_m)
            wall_hit|=hit
        raise RuntimeError("Contact projection did not converge; no unconverged step accepted")
