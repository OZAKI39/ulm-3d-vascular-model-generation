"""C7 reduced phase VI, reusing the verified raw residual derivative source.

All excluded variables are expanded to the exact immutable lower bound before
any physical evaluation. The pool is computational, never an activation oracle.
"""
import ast, resource, time
from types import SimpleNamespace
import numpy as np
from scipy import sparse
from scripts import c5_activation_audit as old
from pd_clot.pfpd_candidate import petsc_vi_candidate as vi


class ReducedProblem:
    def __init__(self,full,ids):
        self.full=full;self.ids=np.asarray(sorted(set(ids)),int)
        for k in ['nx','n','strain','sc','yc','length_scale','energy_scale','config','free','family']:
            setattr(self,k,getattr(full,k))
        self.ns=len(self.ids);self.lower=full.lower[self.ids];self.weights=full.weights[self.ids]
    def expand(self,u):
        z=np.r_[u[:self.nx],self.full.lower.copy()];z[self.nx+self.ids]=u[self.nx:];return z
    def compress(self,z):return np.r_[z[:self.nx],z[self.nx+self.ids]]
    def unpack(self,u):
        x,s=self.full.unpack(self.expand(u));return x,s[self.ids]
    def state(self,x,s):
        phase=self.full.lower.copy();phase[self.ids]=s;state=self.full.state(x,phase)
        return SimpleNamespace(force=state.force,W0=state.W0[self.ids])


def neighbors(f,ids):
    ids=np.asarray(ids,int)
    if len(ids)==0:return ids
    nodes=np.zeros(len(f.X),bool);nodes[f.pairs[ids].ravel()]=True
    return np.flatnonzero(np.any(nodes[f.pairs],axis=1))


def initial_pool(p,x,s,previous_weak=(),threshold=1e-3):
    _,_,_,qs=old.q_fields(p,x,s)
    # Preserve all positive history, including a temporarily bound healed-free
    # variable; a superset of previous interior and weak variables is safe.
    seed=np.union1d(np.flatnonzero((s>1e-8)|(qs<=threshold)),np.asarray(previous_weak,dtype=int))
    return np.union1d(seed,neighbors(p.family,seed)).astype(int)


def excluded_violations(p,z,ids,tolerance=1e-8):
    x,s=p.unpack(z);qs=old.q_fields(p,x,s)[3];excluded=np.ones(p.ns,bool);excluded[ids]=False
    return np.flatnonzero(excluded&(qs < -tolerance))


def expand_pool(p,ids,violated):return np.union1d(ids,np.union1d(violated,neighbors(p.family,violated)))


def column_chunk(task):
    n,strain,lower,ids,u,step,columns=task
    from scripts.c7_early_branch import problem
    p=ReducedProblem(problem(n,strain,lower),ids);nx=p.nx;ns=p.ns;z=u
    class Columns:
        def __init__(self,rows):self.data=np.empty((rows,len(columns)));self.lookup={j:i for i,j in enumerate(columns)}
        def __setitem__(self,key,value):self.data[:,self.lookup[key[1]]]=value
    obj=SimpleNamespace(M=Columns(nx),Yx=Columns(ns));_,loop=vi._source_nodes();loop.iter=ast.Name(id='columns',ctx=ast.Load())
    scope=dict(vars(vi.verified));scope.update(self=obj,p=p,z=z,step=step,nx=nx,ns=ns,floor=p.config['mechanical_equilibrium']['force_floor_N'],columns=columns)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[loop],type_ignores=[])),'<verified reduced column loop>','exec'),scope)
    return columns,obj.M.data,obj.Yx.data,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024


def parallel_columns(obj,p,z,step,pool,workers):
    tasks=[(p.n,p.strain,p.full.lower,p.ids,z,step,c.tolist()) for c in np.array_split(np.arange(p.nx),workers) if len(c)]
    chunks=pool.map(column_chunk,tasks) if pool else map(column_chunk,tasks)
    obj.worker_peak_RSS=0
    for ids,M,Yx,rss in chunks:obj.M[:,ids]=M;obj.Yx[:,ids]=Yx;obj.worker_peak_RSS=max(obj.worker_peak_RSS,rss)


def jacobian(p,u,pool=None,workers=1):
    tree=vi.adapted_constructor();scope=dict(vars(vi.verified));scope.update(_parallel_columns=parallel_columns,pool=pool,workers=workers)
    exec(compile(tree,'<unchanged verified raw derivative constructor>','exec'),scope)
    b=SimpleNamespace();scope['build'](b,p,u)
    nx,ns=p.nx,p.ns;n=nx+ns
    indptr=np.r_[np.arange(nx+1,dtype=np.int64)*n,nx*n+np.arange(1,ns+1,dtype=np.int64)*(nx+1)]
    columns=np.column_stack([np.broadcast_to(np.arange(nx),(ns,nx)),nx+np.arange(ns)])
    indices=np.r_[np.tile(np.arange(n),nx),columns.ravel()].astype(np.int32)
    data=np.r_[np.column_stack([b.M,b.B]).ravel(),np.column_stack([b.C,b.diag]).ravel()]
    return sparse.csr_matrix((data,indices,indptr),shape=(n,n)),b.worker_peak_RSS


class Adapter(vi.PETScAdapter):
    def jacobian(self,snes,v,J,P):
        start=time.perf_counter();matrix,rss=jacobian(self.p,self.gather(v),self.pool,self.workers)
        J.setValuesCSR(matrix.indptr.astype(self.PETSc.IntType),matrix.indices.astype(self.PETSc.IntType),matrix.data);J.assemble()
        self.jacobian_seconds.append(time.perf_counter()-start);self.worker_peak_RSS=max(getattr(self,'worker_peak_RSS',0),rss)
        self.sparsity=dict(shape=list(matrix.shape),nnz=int(matrix.nnz),CSR_bytes=int(matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes),storage='Sparse AIJ reduced arrowhead; no full phase or full total matrix')


def solve_once(full,z,ids,folder,pool,workers):
    from petsc4py import PETSc
    from scripts import petsc_vi_audit as prior
    p=ReducedProblem(full,ids);a=Adapter(p,PETSc.COMM_SELF,pool,workers);a.assign(a.v,p.compress(z));sn=PETSc.SNES().create(PETSc.COMM_SELF)
    cfg=prior.config()['primary'];sn.setType(cfg['type']);sn.setFunction(a.function,a.f);sn.setJacobian(a.jacobian,a.J);sn.setVariableBounds(a.lo,a.hi)
    sn.setTolerances(atol=cfg['atol'],rtol=cfg['rtol'],stol=cfg['stol'],max_it=cfg['max_it']);sn.getLineSearch().setType(cfg['line_search'])
    ksp=sn.getKSP();ksp.setType(cfg['ksp']);pc=ksp.getPC();pc.setType(cfg['pc']);pc.setFactorSolverType(cfg['factor']);history=[]
    def monitor(sn,it,norm):
        zz=p.expand(a.gather(sn.getSolution()));x,s=full.unpack(zz);m=full.metrics(x,s)
        history.append(dict(iteration=int(it),norm=float(norm),force=m['force_residual'],KKT=m['projected_phase_KKT'],max_s=m['maximum_s']))
        old.save(folder/'ITERATIONS.json',history)
        print('RSLS',full.n,folder.name,it,m['force_residual'],m['projected_phase_KKT'],m['maximum_s'],flush=True)
    sn.setMonitor(monitor);start=time.perf_counter();error=None
    try:sn.solve(None,a.v)
    except PETSc.Error as exc:error=str(exc)
    z=p.expand(a.gather(a.v));r=dict(wall_seconds=time.perf_counter()-start,reason=int(sn.getConvergedReason()),iterations=int(sn.getIterationNumber()),linear_iterations=int(sn.getLinearSolveIterations()),error=error,jacobian_seconds=a.jacobian_seconds,worker_peak_RSS_bytes=getattr(a,'worker_peak_RSS',0),sparsity=getattr(a,'sparsity',{}),settings=cfg)
    old.save(folder/'PETSC_SUMMARY.json',r)
    xx,ss=full.unpack(z);np.savez_compressed(folder/'PETSC_STATE.npz',x=xx,s=ss,z=z,lower=full.lower,pool=p.ids)
    sn.destroy();a.destroy();return r,z
