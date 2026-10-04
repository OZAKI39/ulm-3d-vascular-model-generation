"""Isolated PETSc VI candidate. No production entry point imports this module.

The prior verified ResidualBlocks source supplies every derivative formula.
An inspected AST adapter replaces only its independent column loop, removes
projected-map row masking, and retains raw q/D derivatives. The generated source
is archived by the audit. PETSc, not this adapter, implements the VI algorithm.
"""
import ast
import inspect
import textwrap
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from scipy import sparse
from scripts import c5_stationarity_diagnostics as verified
from scripts.c5_activation_audit import q_fields, spectral, Problem


def bounds(p):
    return np.r_[np.full(p.nx,-np.inf),p.lower],np.r_[np.full(p.nx,np.inf),np.full(p.ns,p.sc)]


def residual(p,z,*,diagnostic_extension=False):
    x,s=p.unpack(z)
    if diagnostic_extension:
        # Exact algebraic continuation of the unchanged (1-s)^2 polynomial for
        # central differences at s=0. Never used for nonlinear solve callbacks.
        from . import svk_spectral_split as law
        tree=ast.parse(textwrap.dedent(inspect.getsource(law.degraded_svk_split_response)))
        fn=tree.body[0]
        guards=[node for node in fn.body if isinstance(node,ast.If)]
        assert len(guards)==1 and 'Early phase' in ast.unparse(guards[0])
        fn.body.remove(guards[0]);scope=dict(vars(law))
        exec(compile(ast.fix_missing_locations(tree),'<diagnostic polynomial continuation>','exec'),scope)
        with patch.object(spectral,'degraded_svk_split_response',scope[fn.name]):
            state=spectral.split_baqp(p.family,x,p.material,s)
        Y=state.W0;q=-2*(1-s)*Y+2*p.yc
        return np.r_[state.force[p.free].ravel()/p.config['mechanical_equilibrium']['force_floor_N'],q/np.maximum(Y,p.yc)]
    state=p.state(x,s);_,_,_,qs=q_fields(p,x,s)
    return np.r_[state.force[p.free].ravel()/p.config['mechanical_equilibrium']['force_floor_N'],qs]


def _source_nodes():
    fn=ast.parse(textwrap.dedent(inspect.getsource(verified.ResidualBlocks.__init__))).body[0]
    loop=[node for node in fn.body if isinstance(node,ast.For)]
    assert len(loop)==1 and ast.unparse(loop[0].iter)=='range(nx)'
    return fn,loop[0]


def _column_chunk(args):
    n,strain,lower,z,step,columns=args
    p=Problem(n,strain,lower);nx=p.nx;ns=p.ns
    # Exact old loop body, compact local output columns through a mapping proxy.
    class Columns:
        def __init__(self,rows):self.data=np.empty((rows,len(columns)));self.lookup={j:i for i,j in enumerate(columns)}
        def __setitem__(self,key,value):self.data[:,self.lookup[key[1]]]=value
    result=SimpleNamespace(M=Columns(nx),Yx=Columns(ns))
    _,loop=_source_nodes();loop.iter=ast.Name(id='columns',ctx=ast.Load())
    scope=dict(vars(verified));scope.update(self=result,p=p,z=z,step=step,nx=nx,ns=ns,
        floor=p.config['mechanical_equilibrium']['force_floor_N'],columns=columns)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[loop],type_ignores=[])),'<verified displacement columns>','exec'),scope)
    return columns,result.M.data,result.Yx.data


def _parallel_columns(obj,p,z,step,pool,workers):
    if pool is None:
        chunks=[_column_chunk((p.n,p.strain,p.lower,z,step,list(range(p.nx))))]
    else:
        tasks=[(p.n,p.strain,p.lower,z,step,c.tolist()) for c in np.array_split(np.arange(p.nx),workers) if len(c)]
        chunks=pool.map(_column_chunk,tasks)
    for ids,M,Yx in chunks:obj.M[:,ids]=M;obj.Yx[:,ids]=Yx


def adapted_constructor():
    fn,loop=_source_nodes();fn.name='build'
    i=fn.body.index(loop)
    fn.body[i]=ast.parse('_parallel_columns(self,p,z,step,pool,workers)').body[0]
    mask=[node for node in fn.body if ast.unparse(node).startswith('self.C[self.clipped] =')]
    assert len(mask)==1;fn.body.remove(mask[0])
    last=fn.body[-1];assert ast.unparse(last)=='self.diag = np.where(self.clipped, 1.0, self.Qs / self.D)'
    fn.body[-1]=ast.parse('self.diag = self.Qs / self.D').body[0]
    return ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))


def raw_blocks(p,z,pool=None,workers=1):
    tree=adapted_constructor();scope=dict(vars(verified));scope.update(_parallel_columns=_parallel_columns,pool=pool,workers=workers)
    exec(compile(tree,'<verified raw VI blocks>','exec'),scope)
    b=SimpleNamespace();scope['build'](b,p,z)
    return b


def sparse_jacobian(p,z,pool=None,workers=1):
    b=raw_blocks(p,z,pool,workers)
    # Explicit structural zeros keep AIJ preallocation invariant across Newton.
    nx,ns=p.nx,p.ns;n=nx+ns
    indptr=np.r_[np.arange(nx+1,dtype=np.int64)*n,nx*n+np.arange(1,ns+1,dtype=np.int64)*(nx+1)]
    phase_cols=np.column_stack([np.broadcast_to(np.arange(nx),(ns,nx)),nx+np.arange(ns)])
    indices=np.r_[np.tile(np.arange(n),nx),phase_cols.ravel()].astype(np.int32)
    data=np.r_[np.column_stack([b.M,b.B]).ravel(),np.column_stack([b.C,b.diag]).ravel()]
    return sparse.csr_matrix((data,indices,indptr),shape=(n,n))


def make_pool(workers):
    return ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn')) if workers>1 else None


class PETScAdapter:
    """Rank-zero physics with distributed vectors/AIJ and PETSc linear solves."""
    def __init__(self,p,comm,pool=None,workers=1):
        from petsc4py import PETSc
        self.PETSc=PETSc;self.p=p;self.comm=comm;self.mpi=comm.tompi4py();self.pool=pool;self.workers=workers
        self.v=PETSc.Vec().createMPI(p.nx+p.ns,comm=comm);self.f=self.v.duplicate()
        self.first,self.last=self.v.getOwnershipRange();self.jacobian_seconds=[]
        local=np.arange(self.first,self.last);n=p.nx+p.ns
        counts=np.where(local<p.nx,n,p.nx+1).astype(PETSc.IntType)
        if comm.size==1:
            self.J=PETSc.Mat().createAIJ(size=(n,n),nnz=counts,comm=comm)
        else:
            diag=np.where(local<p.nx,self.last-self.first,
                max(0,min(self.last,p.nx)-self.first)+(local>=p.nx)).astype(PETSc.IntType)
            self.J=PETSc.Mat().createAIJ(size=((self.last-self.first,n),(self.last-self.first,n)),nnz=(diag,counts-diag),comm=comm)
        self.J.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR,True)
        self.lo=self.v.duplicate();self.hi=self.v.duplicate();l,u=bounds(p);self.assign(self.lo,l);self.assign(self.hi,u)

    def assign(self,v,a):v.array[:]=a[self.first:self.last]
    def gather(self,v):return np.concatenate(self.mpi.allgather(v.array_r.copy()))
    def function(self,snes,v,f):
        z=self.gather(v);r=residual(self.p,z) if self.comm.rank==0 else None
        r=self.mpi.bcast(r,root=0);self.assign(f,r)
    def jacobian(self,snes,v,J,P):
        import time
        start=time.perf_counter();z=self.gather(v)
        matrix=sparse_jacobian(self.p,z,self.pool,self.workers) if self.comm.rank==0 else None
        matrix=self.mpi.bcast(matrix,root=0);local=matrix[self.first:self.last]
        J.setValuesCSR(local.indptr.astype(self.PETSc.IntType),local.indices.astype(self.PETSc.IntType),local.data)
        J.assemble();self.jacobian_seconds.append(time.perf_counter()-start)
        self.sparsity=dict(shape=list(matrix.shape),structural_nnz=int(matrix.nnz),numerical_nnz=int(np.count_nonzero(matrix.data)),
            CSR_bytes=int(matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes),
            block_nnz=dict(Jxx=self.p.nx**2,Jxs=self.p.nx*self.p.ns,Jsx=self.p.nx*self.p.ns,Jss=self.p.ns),
            PETSc_local_info=J.getInfo(),storage='AIJ arrowhead; no dense total or phase square matrix')
    def destroy(self):
        for obj in [self.v,self.f,self.J,self.lo,self.hi]:obj.destroy()
