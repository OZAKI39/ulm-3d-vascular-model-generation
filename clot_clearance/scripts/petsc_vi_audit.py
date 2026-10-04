"""Reproducible, gated PETSc VI validation; no historical writers invoked."""
from pathlib import Path
import argparse,ast,ctypes,hashlib,json,os,sys,time,resource,subprocess
import numpy as np
from scripts import c5_activation_audit as prev
from scripts.c5_stationarity_diagnostics import full_acceptance
from pd_clot.pfpd_candidate import petsc_vi_candidate as vi
ROOT=prev.ROOT;OLD=prev.OUT;OUT=ROOT/'verification/pfpd_at1_petsc_vi_solver_v1'
read=prev.read;save=prev.save;Problem=prev.Problem
LABELS=['first_A1','continuation_01','continuation_04','offset_1e-05','same_load_as_original','offset_2e-05','offset_5e-05']


def config():return read(OUT/'provenance/PREDECLARED_CONFIG.json')
def state(label):return prev.ca.read_npz(OLD/f'reduced/{label}/STATE.npz')
def targets():return [dict(label=k,strain=read(OLD/f'reduced/{k}/SUMMARY.json')['strain']) for k in LABELS]
def initial():
    a=prev.ca.read_npz(OLD/'reproduction/C5_W_004.npz');r=read(OLD/'reproduction/C5_W_004.json')
    return a['x'],a['s'],r['strain']


def declare():
    from petsc4py import PETSc
    import petsc4py,mpi4py,platform
    from mpi4py import MPI
    info={k:int(v.split()[0])*1024 for k,v in (s.split(':',1) for s in Path('/proc/meminfo').read_text().splitlines())}
    physical=int(subprocess.check_output(['bash','-c','lscpu -p=CORE,SOCKET | sed "/^#/d" | sort -u | wc -l'],text=True))
    reserve=.2*info['MemTotal'];parent=.75*1024**3;per_worker=.6*1024**3
    memory_safe=max(1,int((info['MemAvailable']-reserve-parent)//per_worker));workers=min(physical-2,memory_safe)
    c=dict(task='PFPD_AT1_PETSC_VI_SOLVER_V1',force_tolerance=1e-9,KKT_tolerance=1e-8,
        agreement=dict(phase_Linf=1e-8,scaled_state_L2=1e-8,energy_relative=1e-9,reaction_relative=1e-8),
        jacobian_epsilons=[1e-5,3e-6,1e-6],jacobian_vectors=10,jacobian_seed=93121,jacobian_max_error=5e-5,
        primary=dict(type='vinewtonrsls',line_search='bt',atol=1e-12,rtol=1e-12,stol=1e-14,max_it=30,
            ksp='preonly',pc='lu',factor='mumps'),workers=workers,
        resources=dict(physical_cores=physical,memory_safe_worker_count=memory_safe,reserve_RAM_bytes=reserve,
            parent_budget_bytes=parent,worker_budget_bytes=per_worker,mem_available_bytes=info['MemAvailable'],
            BLAS_threads=1,FD_parallelism='ProcessPoolExecutor spawn; exact old displacement-column loop'),
        targets=targets(),oracle_use='Post-solve comparison only; diagnostic gate C reads a saved root, never passes it to any solver.',
        phase_F='raw q/max(Psi_plus,Yc), no projected clip',negative_phase_FD='Diagnostic-only exact polynomial continuation by removing original s-domain guard in a private cloned function; solver uses original guard.',
        derivative_reuse='AST clone of verified ResidualBlocks.__init__: parallelize its loop, remove C[clipped]=0, replace projected diag by Qs/D; all other derivative statements unchanged.',
        production_changes=False,stability_assessed=False)
    save(OUT/'provenance/PREDECLARED_CONFIG.json',c)
    package_names=['mumps','superlu_dist','suitesparse','hypre','metis','parmetis','scalapack','cuda','kokkos','kokkos_kernels','hip','sycl','opencl']
    env=dict(status='AVAILABLE',python=sys.executable,python_version=platform.python_version(),PETSc_version=PETSc.Sys.getVersion(),
        petsc4py_version=petsc4py.__version__,mpi4py_version=mpi4py.__version__,MPI_implementation=MPI.Get_library_version(),
        scalar_type=str(PETSc.ScalarType),integer_type=str(PETSc.IntType),packages={k:PETSc.Sys.hasExternalPackage(k) for k in package_names},
        CUDA_support=PETSc.Sys.hasExternalPackage('cuda'),VecCUDA_usable=PETSc.Sys.hasExternalPackage('cuda'),AIJCUSPARSE_usable=PETSc.Sys.hasExternalPackage('cuda'),
        config=petsc4py.get_config(),dependency_strategy='Existing PETSc runtime unchanged; numba and matching llvmlite wheels in task-local --target directory.')
    save(OUT/'provenance/PETSC_ENVIRONMENT.json',env)
    source=Path(petsc4py.get_config()['PETSC_DIR'])/'lib/petsc/conf/petscvariables'
    if source.exists():(OUT/'provenance/petscvariables.txt').write_text(source.read_text())
    (OUT/'provenance/ADAPTED_VERIFIED_CONSTRUCTOR.py.txt').write_text(ast.unparse(vi.adapted_constructor())+'\n')
    print('DECLARED',workers,'workers',env['PETSc_version'],flush=True)


def mapping_tests():
    from petsc4py import PETSc
    cases=[]
    for name,target,expected,sign in [('lower',-1.,0.,1),('upper',2.,1.,-1),('interior',.4,.4,0)]:
        x=PETSc.Vec().createSeq(1);x.set(.2);f=x.duplicate();lo=x.duplicate();lo.set(0);hi=x.duplicate();hi.set(1)
        J=PETSc.Mat().createAIJ([1,1],nnz=1);sn=PETSc.SNES().create();sn.setType('vinewtonrsls')
        def fun(sn,x,f):f.array[:]=x.array_r-target
        def jac(sn,x,J,P):J[0,0]=1;J.assemble()
        sn.setFunction(fun,f);sn.setJacobian(jac,J);sn.setVariableBounds(lo,hi)
        sn.getKSP().setType('preonly');sn.getKSP().getPC().setType('lu');sn.setTolerances(atol=1e-12)
        sn.solve(None,x);value=float(x.array[0]);q=value-target
        cases.append(dict(case=name,value=value,q=q,reason=int(sn.getConvergedReason()),passed=abs(value-expected)<1e-12 and (q*sign>0 if sign else abs(q)<1e-12)))
        sn.destroy()
        for o in [x,f,lo,hi,J]:o.destroy()
    p=Problem(5,targets()[0]['strain']);adapter=vi.PETScAdapter(p,PETSc.COMM_SELF)
    z=np.arange(p.nx+p.ns,dtype=float)/10000;adapter.assign(adapter.v,z);roundtrip=np.array_equal(adapter.gather(adapter.v),z)
    adapter.destroy();save(OUT/'mapping/SCALAR_AND_VECTOR.json',dict(cases=cases,roundtrip=roundtrip,passed=roundtrip and all(r['passed'] for r in cases)))
    assert roundtrip and all(r['passed'] for r in cases)


def gate():
    from petsc4py import PETSc
    cfg=config();pool=vi.make_pool(cfg['workers']);rows=[]
    try:
        ax,ss,eps=initial();states=[('A_last_s0',Problem(5,eps),ax,ss)]
        r=read(OLD/'activation/POINT.json');a=prev.ca.read_npz(OLD/'intact/activation_point/STATE.npz')
        states.append(('B_first_weak',Problem(5,r['strain']),a['x'],a['s']))
        r=read(OLD/'reduced/first_A1/SUMMARY.json');a=state('first_A1')
        states.append(('C_positive_root',Problem(5,r['strain'],a['lower']),a['x'],a['s']))
        for label,p,x,s in states:
            start=time.perf_counter();z=p.pack(x,s);a=vi.PETScAdapter(p,PETSc.COMM_SELF,pool,cfg['workers']);a.assign(a.v,z);a.jacobian(None,a.v,a.J,a.J)
            rng=np.random.default_rng(cfg['jacobian_seed']);checks=[];v=a.v.duplicate();w=v.duplicate()
            for i in range(cfg['jacobian_vectors']):
                direction=rng.normal(size=len(z));direction/=np.linalg.norm(direction);a.assign(v,direction);a.J.mult(v,w);actual=a.gather(w)
                for step in cfg['jacobian_epsilons']:
                    expected=(vi.residual(p,z+step*direction,diagnostic_extension=True)-vi.residual(p,z-step*direction,diagnostic_extension=True))/(2*step)
                    checks.append(dict(direction=i,epsilon=step,relative_error=float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected),1e-30))))
            result=dict(label=label,checks=checks,max_error=max(c['relative_error'] for c in checks),wall_seconds=time.perf_counter()-start,
                        sparsity=a.sparsity)
            result['passed']=result['max_error']<=cfg['jacobian_max_error'];save(OUT/f'jacobian/{label}.json',result);rows.append(result)
            if label=='A_last_s0':save(OUT/'JACOBIAN_SPARSITY.json',a.sparsity)
            v.destroy();w.destroy();a.destroy();print('JACOBIAN',label,result['max_error'],result['passed'],flush=True)
        result=dict(passed=all(r['passed'] for r in rows),max_error=max(r['max_error'] for r in rows),states=rows)
        save(OUT/'jacobian/GATE.json',result);assert result['passed'],'STOP: PETSC_JACOBIAN_GATE_FAILED'
    finally:
        if pool:pool.shutdown()


def acceptance(p,z,reason):
    # No PETSc reason code enters physical acceptance.
    x,s=p.unpack(z);interior=np.flatnonzero(s>p.lower+1e-10)
    result=full_acceptance(p,z,interior)
    result.update(PETSc_reason=int(reason),scientifically_accepted=result['valid'],interior_IDs=interior.tolist())
    return result


def compare(p,z,label):
    a=state(label);r=read(OLD/f'reduced/{label}/SUMMARY.json');x,s=p.unpack(z);metrics=p.metrics(x,s)
    reaction=-p.state(x,s).force[p.top,2].sum()/prev.v2.L**2
    errors=dict(phase_Linf=float(abs(s-a['s']).max()),scaled_state_L2=float(np.linalg.norm((x[p.free]-a['x'][p.free])/p.length_scale)),
        energy_relative=float(abs(metrics['total_J']/r['metrics']['total_J']-1)),reaction_relative=float(abs(reaction/r['reaction_stress_Pa']-1)))
    expected=set(a['active'].tolist());actual=set(np.flatnonzero(s>p.lower+1e-10).tolist())
    discovery=dict(expected_IDs=sorted(expected),actual_IDs=sorted(actual),precision=len(actual&expected)/max(len(actual),1),recall=len(actual&expected)/max(len(expected),1),oracle_access='post-solve only')
    return dict(errors=errors,agrees=all(errors[k]<=v for k,v in config()['agreement'].items()),discovery=discovery)


def solve(p,z0,folder,kind='vinewtonrsls',pool=None,workers=1,max_it=None):
    from petsc4py import PETSc
    from mpi4py import MPI
    comm=PETSc.COMM_WORLD;rank=comm.rank;folder.mkdir(parents=True,exist_ok=True)
    a=vi.PETScAdapter(p,comm,pool,workers);a.assign(a.v,z0);snes=PETSc.SNES().create(comm);snes.setType(kind)
    snes.setFunction(a.function,a.f);snes.setJacobian(a.jacobian,a.J);snes.setVariableBounds(a.lo,a.hi)
    cfg=config()['primary'];snes.setTolerances(atol=cfg['atol'],rtol=cfg['rtol'],stol=cfg['stol'],max_it=cfg['max_it'] if max_it is None else max_it)
    ksp=snes.getKSP();ksp.setType(cfg['ksp']);pc=ksp.getPC();pc.setType(cfg['pc']);pc.setFactorSolverType(cfg['factor'])
    opts=PETSc.Options();prefix='audit_';snes.setOptionsPrefix(prefix)
    settings=dict(snes_vi_monitor=True,snes_vi_monitor_residual=f'binary:{folder}/VI_RESIDUAL.dat',
        snes_vi_monitor_active=f'binary:{folder}/VI_ACTIVE.dat',snes_converged_reason=f'ascii:{folder}/CONVERGED_REASON.txt',snes_linesearch_type='bt')
    for k,v in settings.items():opts[prefix+k]=v
    snes.setFromOptions();history=[];last_set=set()
    # petsc4py lacks getLambda in 3.25; use the official C getter read-only.
    import petsc4py
    lib=ctypes.CDLL(str(Path(petsc4py.get_config()['PETSC_DIR'])/'lib/libpetsc.so'))
    get_lambda=lib.SNESLineSearchGetLambda;get_lambda.argtypes=[ctypes.c_void_p,ctypes.POINTER(ctypes.c_double)];get_lambda.restype=ctypes.c_int
    monitor_f=a.v.duplicate()
    get_active=lib.SNESVIGetActiveSetIS;get_active.argtypes=[ctypes.c_void_p]*3+[ctypes.POINTER(ctypes.c_void_p)]
    is_size=lib.ISGetLocalSize;is_size.argtypes=[ctypes.c_void_p,ctypes.POINTER(ctypes.c_int)]
    is_indices=lib.ISGetIndices;is_indices.argtypes=[ctypes.c_void_p,ctypes.POINTER(ctypes.POINTER(ctypes.c_int))]
    is_restore=lib.ISRestoreIndices;is_restore.argtypes=is_indices.argtypes
    is_destroy=lib.ISDestroy;is_destroy.argtypes=[ctypes.POINTER(ctypes.c_void_p)]
    def monitor(sn,it,norm):
        nonlocal last_set
        z=a.gather(sn.getSolution());lam=ctypes.c_double();code=get_lambda(ctypes.c_void_p(sn.getLineSearch().handle),ctypes.byref(lam));assert code==0
        # RSLS destroys IS_inact before SNESMonitor (virs.c:509,584).
        # Query a fresh PETSc active-set IS through its public C API instead.
        a.function(sn,sn.getSolution(),monitor_f)
        handle=ctypes.c_void_p();size=ctypes.c_int();indices=ctypes.POINTER(ctypes.c_int)()
        assert get_active(ctypes.c_void_p(sn.handle),ctypes.c_void_p(sn.getSolution().handle),ctypes.c_void_p(monitor_f.handle),ctypes.byref(handle))==0
        assert is_size(handle,ctypes.byref(size))==0
        assert is_indices(handle,ctypes.byref(indices))==0
        fixed=np.ctypeslib.as_array(indices,shape=(size.value,)).copy() if size.value else np.empty(0,np.int32)
        assert is_restore(handle,ctypes.byref(indices))==0
        assert is_destroy(ctypes.byref(handle))==0
        fixed=np.concatenate(a.mpi.allgather(fixed))
        free_phase=np.setdiff1d(np.arange(p.ns),fixed[fixed>=p.nx]-p.nx)
        if rank==0:
            x,s=p.unpack(z);r=vi.residual(p,z);metrics=p.metrics(x,s)
            atlo=s<=p.lower+1e-10;athi=s>=p.sc-1e-10;interior=~(atlo|athi);current=set(np.flatnonzero(interior).tolist())
            history.append(dict(iteration=int(it),residual_norm=float(norm),force=metrics['force_residual'],KKT=metrics['projected_phase_KKT'],
                lower_count=int(atlo.sum()),interior_count=int(interior.sum()),upper_count=int(athi.sum()),
                PETSc_free_phase_IDs=free_phase.tolist(),PETSc_free_phase_count=len(free_phase),active_set_changes=len(current^last_set),
                line_search_lambda=lam.value if it else None,linear_iterations=int(sn.getLinearSolveIterations()),reason=int(sn.getConvergedReason())))
            last_set=current;save(folder/'ITERATIONS.json',history)
            print(kind,'iteration',it,'force',metrics['force_residual'],'KKT',metrics['projected_phase_KKT'],'interior',interior.sum(),flush=True)
    snes.setMonitor(monitor);a.mpi.Barrier();start=time.perf_counter();error=None
    try:snes.solve(None,a.v)
    except PETSc.Error as exc:error=str(exc)
    elapsed=a.mpi.allreduce(time.perf_counter()-start,op=MPI.MAX);z=a.gather(a.v);reason=int(snes.getConvergedReason())
    result=None
    if rank==0:
        x,s=p.unpack(z)
        try:checked=acceptance(p,z,reason)
        except Exception as exc:checked=dict(scientifically_accepted=False,acceptance_error=str(exc))
        result=dict(strain=p.strain,solver=kind,ranks=comm.size,wall_seconds=elapsed,SNES_iterations=int(snes.getIterationNumber()),
            linear_iterations=int(snes.getLinearSolveIterations()),reason=reason,error=error,line_search=snes.getLineSearch().getType(),
            jacobian_seconds=a.jacobian_seconds,active_set_changes=sum(r['active_set_changes'] for r in history),
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,**checked)
        save(folder/'SUMMARY.json',result);np.savez_compressed(folder/'STATE.npz',x=x,s=s,z=z,z_initial=z0,lower=p.lower)
        if hasattr(a,'sparsity'):save(folder/'JACOBIAN_SPARSITY.json',a.sparsity)
    result=a.mpi.bcast(result,root=0);snes.destroy();monitor_f.destroy();a.destroy()
    for k in settings:del opts[prefix+k]
    return result,z


def run(mode,ranks=1):
    from petsc4py import PETSc
    assert read(OUT/'jacobian/GATE.json')['passed']
    rank=PETSc.COMM_WORLD.rank;workers=config()['workers'];pool=vi.make_pool(workers) if rank==0 else None
    try:
        if mode=='benchmark':
            x,s,eps=initial();t=targets()[0];p=Problem(5,t['strain'],s);z0=p.pack(prev.tracker.displacement_predictor(p,x,eps),s)
            r,z=solve(p,z0,OUT/f'mpi/ranks_{PETSc.COMM_WORLD.size}',pool=pool,workers=workers)
            if rank==0:r.update(compare(p,z,t['label']));save(OUT/f'mpi/ranks_{PETSc.COMM_WORLD.size}/SUMMARY.json',r)
        elif mode=='continuation':
            x,s,eps=initial();rows=[];start=time.perf_counter()
            for t in targets():
                p=Problem(5,t['strain'],s);z0=p.pack(prev.tracker.displacement_predictor(p,x,eps),s);folder=OUT/f'continuation/PETSC_RSLS/{t["label"]}'
                r,z=solve(p,z0,folder,pool=pool,workers=workers)
                if rank==0:r.update(compare(p,z,t['label']));save(folder/'SUMMARY.json',r)
                r=PETSc.COMM_WORLD.tompi4py().bcast(r,root=0);rows.append(dict(label=t['label'],**r))
                if rank==0:save(OUT/'continuation/PETSC_RSLS/ALL.json',dict(states=rows,wall_seconds=time.perf_counter()-start,passed=sum(a['scientifically_accepted'] and a['agrees'] for a in rows)))
                assert r['scientifically_accepted'],'STOP: C5 scientific acceptance failed'
                x,s=p.unpack(z);eps=t['strain']
        elif mode.startswith('ssls:'):
            label=mode.split(':')[1];a=prev.ca.read_npz(OUT/f'continuation/PETSC_RSLS/{label}/STATE.npz');t=next(t for t in targets() if t['label']==label)
            p=Problem(5,t['strain'],a['lower']);r,z=solve(p,a['z_initial'],OUT/f'ssls/{label}',kind='vinewtonssls',pool=pool,workers=workers)
            if rank==0:r.update(compare(p,z,label));save(OUT/f'ssls/{label}/SUMMARY.json',r)
    finally:
        if pool:pool.shutdown()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode');args=parser.parse_args()
    if args.mode=='declare':declare();mapping_tests()
    elif args.mode=='gate':gate()
    else:run(args.mode)
