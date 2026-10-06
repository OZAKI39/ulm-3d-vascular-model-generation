"""Parallel native-solver comparison on immutable exposed development QPs."""
import json
from pathlib import Path
from time import perf_counter
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import sparse
from mevo_mb.bubbles.constrained_motion import kkt_residual
from scripts.stage4l.dense_probe import write
from scripts.stage5a_r2.seal import BASE


def solve(task):
    folder, name, eps = task
    p=Path(folder)
    a,g=sparse.load_npz(p/'A.npz'),sparse.load_npz(p/'G.npz')
    data=dict(np.load(next(x for x in p.glob('*.npz') if x.name not in ['A.npz','G.npz'])))
    b,lo=data['b'],data['lower']; m,n=g.shape
    begin=perf_counter()
    if name == 'osqp':
        import osqp
        q=osqp.OSQP(); q.setup(P=a,q=-b,A=g,l=lo,u=np.full(m,np.inf),verbose=False,
                             eps_abs=eps,eps_rel=0,max_iter=300000,polishing=True,
                             polish_refine_iter=10, time_limit=10.)
        setup=perf_counter()-begin; result=q.solve(raise_error=False)
        u,lam=result.x,-result.y; status=result.info.status; it=result.info.iter
    elif name.startswith('prox'):
        import proxsuite as ps
        dense=name=='prox_dense'
        q=(ps.proxqp.dense if dense else ps.proxqp.sparse).QP(n,0,m)
        q.settings.eps_abs=eps;q.settings.eps_rel=0;q.settings.max_iter=1000
        q.settings.max_iter_in=1500;q.settings.nb_iterative_refinement=10;q.settings.compute_timings=True
        q.init(a.toarray() if dense else a,-b,None,None,g.toarray() if dense else g,lo,np.full(m,np.inf))
        setup=perf_counter()-begin;q.solve()
        u,lam=q.results.x,-q.results.z;status=str(q.results.info.status);it=q.results.info.iter
    else:
        import clarabel
        settings=clarabel.DefaultSettings();settings.verbose=False;settings.max_iter=500
        settings.tol_gap_abs=eps;settings.tol_gap_rel=eps;settings.tol_feas=eps
        settings.time_limit=10.
        q=clarabel.DefaultSolver(sparse.triu(a,format='csc'),-b,-g,-lo,[clarabel.NonnegativeConeT(m)],settings)
        setup=perf_counter()-begin;result=q.solve()
        u,lam=np.array(result.x),np.array(result.z);status=str(result.status);it=result.iterations
    seconds=perf_counter()-begin
    k,terms=kkt_residual(a,b,g,lo,u,lam)
    out=BASE/'diagnosis/native_v1';out.mkdir(exist_ok=True)
    np.savez_compressed(out/f'{p.name}_{name}_{eps:g}.npz',u=u,lam=lam)
    return dict(fixture=p.name,solver=name,eps=eps,status=status,kkt=k,terms=terms,
                objective=float(.5*u@(a@u)-b@u),seconds=seconds,setup_s=setup,iterations=int(it))


if __name__=='__main__':
    tasks=[(str(p),s,e) for p in sorted((BASE/'fixtures').iterdir()) if p.is_dir()
           for s in ['osqp','prox_dense','prox_sparse','clarabel'] for e in [1e-10,1e-12]]
    with ProcessPoolExecutor(max_workers=6) as pool:
        rows=list(pool.map(solve,tasks))
    write(BASE/'diagnosis/native_v1/SUMMARY.json',rows)
    for r in rows:
        print(r['fixture'],r['solver'],r['eps'],r['status'],r['kkt'],r['seconds'],flush=True)
