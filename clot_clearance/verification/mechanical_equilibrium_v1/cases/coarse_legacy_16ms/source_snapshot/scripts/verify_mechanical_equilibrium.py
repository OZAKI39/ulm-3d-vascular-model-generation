"""Resource-limited parallel, no-damage window and SciPy benchmark campaign.

Run from clot_clearance: PYTHONPATH=. python -B scripts/verify_mechanical_equilibrium.py
All output directories must be new. ADR and erosion reruns remain gated by
the missing primary paper. The old 78 tests and all prior outputs are untouched.
"""
import os
for _key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[_key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import contextlib,copy,gzip,hashlib,json,csv,resource,shutil,time,traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/mechanical_equilibrium_v1'


def read(path):return json.loads(path.read_text())


def write(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)


def capacity(count):
    hardware=read(OUT/'HARDWARE.json')
    available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))*1024
    return max(1,min(count,hardware['available_physical_cores'],hardware['allowed_logical_cpu_count'],
                     int(.8*min(available,hardware['effective_available_memory_bytes']))//hardware['estimated_worker_memory_bytes']))


def work(job):
    out=OUT/'cases'/job['name'];out.mkdir(parents=True,exist_ok=False)
    os.environ['NUMBA_CACHE_DIR']=str(ROOT/f'build_mechanical_equilibrium_v1/numba_worker_{os.getpid()}')
    started=time.perf_counter();cpu=time.process_time()
    with (out/'run.log').open('x') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        try:
            import numpy as np
            from pd_clot.regularization.equilibrium import make_fixture,equilibrate,FrozenForceEvaluator
            c=read(ROOT/'configs/mechanical_equilibrium_verification_v1.json');c['name']=job['name']
            c['mechanical_equilibrium']['method']=job['method']
            if job['size']=='small':c['clot']['cells']=c['mechanical_equilibrium']['small_fixture_cells']
            if 'window_s' in job:c['mechanical_equilibrium']['legacy_window_s']=job['window_s']
            if job.get('provider')=='resolved':
                source=ROOT/'verification/streaming_erosion_v1/tests/test_static_resolved_field_use0/static_field.vtu'
                shutil.copy2(source,out/'static_field.vtu')
                c['streaming'].update(provider='resolved_field',field_path=str(out/'static_field.vtu'),metadata_path=dict(length_unit='m',velocity_unit='m/s',pressure_unit='Pa',time_unit='s',normal_convention='outward_solid',sample_domain='convex_hull'))
            write(out/'CONFIG.json',c)
            cloud=make_fixture(c);g=np.ones(len(cloud.pairs));x0=cloud.X.copy();v0=np.zeros_like(x0)
            source=[*ROOT.glob('pd_clot/*.py'),*ROOT.glob('pd_clot/regularization/*.py'),Path(__file__)]
            hashes={}
            for file in source:
                rel=file.relative_to(ROOT);dest=out/'source_snapshot'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(file,dest)
                hashes[str(rel)]=hashlib.sha256(file.read_bytes()).hexdigest()
            result=equilibrate(cloud,x0,v0,g,c,method=job['method'],verification_only=True)
            assert np.array_equal(g,np.ones(len(g)))
            force=FrozenForceEvaluator(cloud,g,c)(result.x)
            optional={} if result.v is None else dict(v=result.v)
            np.savez_compressed(out/'states.npz',X=cloud.X,x=result.x,x_initial=x0,pairs=cloud.pairs,
                volume=cloud.volume,fixed=cloud.fixed,integrity=g,bond_D=1-g,particle_damage=np.zeros(len(x0)),
                epsilon_plus=force.epsilon_plus,internal_force_N=force.internal_force_N,
                external_force_N=force.external_force_N,surface_traction_Pa=force.surface_traction_Pa,
                erosion_exposure_time_s=0.,**optional)
            columns={}
            for key in result.history[0]:
                values=[row.get(key) for row in result.history]
                if all(value is None for value in values):continue
                columns[key]=np.array(values) if isinstance(next(v for v in values if v is not None),str) else np.array([np.nan if v is None else v for v in values],float)
            np.savez_compressed(out/'history.npz',**columns)
            with gzip.open(out/'history.csv.gz','wt',newline='') as f:
                w=csv.DictWriter(f,list(result.history[0]));w.writeheader();w.writerows(result.history)
            write(out/'checks.json',result.checks)
            result.summary.update(job=job,wall_time_s=time.perf_counter()-started,cpu_time_s=time.process_time()-cpu,
                peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                particle_count=len(cloud.X),bond_count=len(cloud.pairs),free_particle_count=int((~cloud.fixed).sum()))
            write(out/'SUMMARY.json',result.summary)
            write(out/'IDENTITY.json',dict(config_sha256=hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest(),source_sha256=hashes,
                deterministic=True,randomness_used=False,history_sampling='Every mechanical step; SciPy accepted iterations plus final result',
                missing_numeric_diagnostics='NaN in history arrays means not applicable/not a check step; omitted if unavailable for every row',
                scipy_is_small_reference_only=True,damage_disabled=True,erosion_time_always_zero=True,
                force_units='N; existing NOSB kernel already includes particle volumes',
                RMS_definition='L2 of free particle force vector divided by sqrt(number of free particles)',
                runtime_threads={key:os.environ[key] for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']}))
            print(json.dumps(result.summary,indent=2),flush=True)
            return dict(name=job['name'],status=result.summary['status'],wall_time_s=result.summary['wall_time_s'])
        except Exception as exc:
            trace=traceback.format_exc();print(trace,flush=True)
            write(out/'FAILED.json',dict(error=repr(exc),traceback=trace,wall_time_s=time.perf_counter()-started))
            return dict(name=job['name'],status='FAILED',error=repr(exc))


def batch(jobs,label):
    n=capacity(len(jobs));start=time.perf_counter()
    print(json.dumps(dict(batch=label,parallel_workers=n,cases=[j['name'] for j in jobs])),flush=True)
    result=[]
    with ProcessPoolExecutor(max_workers=n,mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(work,j) for j in jobs]
        for future in as_completed(futures):
            row=future.result();result.append(row);print(json.dumps(row),flush=True)
    write(OUT/f'{label}_BATCH.json',dict(workers=n,wall_time_s=time.perf_counter()-start,jobs=result))
    return result


def compare():
    import numpy as np
    reference=read(OUT/'cases/small_scipy/SUMMARY.json')
    with np.load(OUT/'cases/small_scipy/states.npz') as z:ref={k:z[k] for k in z.files}
    ref_u=ref['x']-ref['X'];free=~ref['fixed'];umax=float(np.linalg.norm(ref_u,axis=1).max())
    acceptance=read(ROOT/'configs/mechanical_equilibrium_verification_v1.json')['mechanical_equilibrium']['reference_acceptance']
    comparisons=[]
    for ms in [2,4,8,16]:
        path=OUT/f'cases/small_legacy_{ms}ms';s=read(path/'SUMMARY.json')
        with np.load(path/'states.npz') as z:u=z['x']-z['X']
        diff=u-ref_u;rms=float(np.sqrt(np.mean(np.sum(diff[free]**2,axis=1))))
        rel_driver=abs(s['final']['maximum_epsilon_plus']-reference['final']['maximum_epsilon_plus'])/max(reference['final']['maximum_epsilon_plus'],1e-12)
        comparisons.append(dict(window_ms=ms,maximum_displacement_m=s['final']['maximum_displacement_m'],
            maximum_epsilon_plus=s['final']['maximum_epsilon_plus'],mean_epsilon_plus=s['final']['mean_epsilon_plus'],
            force_residual_L2_N=s['final']['force_residual_L2_N'],relative_force_residual=s['final']['relative_force_residual'],
            stored_energy_J=s['final']['stored_energy_J'],surface_traction_mean_Pa=s['final']['surface_traction_mean_Pa'],
            surface_traction_max_Pa=s['final']['surface_traction_max_Pa'],rms_displacement_difference_m=rms,
            maximum_displacement_difference_m=float(np.linalg.norm(diff[free],axis=1).max()),
            relative_L2_displacement_difference=float(np.linalg.norm(diff[free])/max(np.linalg.norm(ref_u[free]),1e-30)),
            rms_difference_over_reference_max=rms/max(umax,1e-30),relative_max_driver_difference=rel_driver,
            declared_comparison_passed=bool(reference['converged'] and s['converged'] and s['final']['relative_force_residual']<=acceptance['relative_force_residual'] and
                rel_driver<=acceptance['relative_max_driver_difference'] and rms/max(umax,1e-30)<=acceptance['rms_displacement_difference_over_reference_max'])))
    main=[read(OUT/f'cases/coarse_legacy_{ms}ms/SUMMARY.json') for ms in [2,4,8,16]]
    write(OUT/'COMPARISONS.json',dict(scipy_reference=reference,small_fixed_window_vs_scipy=comparisons,coarse_fixed_windows=main,
          adr_status='NEED_PAPER_INPUT',mechanical_equilibrium_status='NOT VERIFIED',
          pre_relaxation_status='NOT RUN: ADR/reference verification not complete',
          equilibrated_erosion_status='NOT RUN',slow_step_study_status='NOT RUN',
          note='Legacy-vs-SciPy agreement is limited to this small fixture; it cannot replace requested ADR verification.'))


def main():
    path=ROOT/'configs/mechanical_equilibrium_verification_v1.json'
    assert hashlib.sha256(path.read_bytes()).hexdigest()==read(OUT/'provenance/PREDECLARED_PLAN.json')['config_sha256']
    if (OUT/'cases').exists():raise FileExistsError('Refusing to overwrite existing benchmark cases')
    small=[dict(name='small_scipy',method='scipy_reference',size='small'),
           dict(name='small_scipy_resolved',method='scipy_reference',size='small',provider='resolved')]
    small += [dict(name=f'small_legacy_{ms}ms',method='legacy_fixed_damping',size='small',window_s=ms/1000) for ms in [2,4,8,16]]
    rows=batch(small,'SMALL')
    if any(row['status']=='FAILED' for row in rows):raise RuntimeError('Small-case execution failed; inspect individual logs before larger sweep')
    # No-damage window sensitivity remains allowed even if the algebraic oracle fails its independent criterion.
    rows=batch([dict(name=f'coarse_legacy_{ms}ms',method='legacy_fixed_damping',size='coarse',window_s=ms/1000) for ms in [2,4,8,16]],'COARSE')
    if any(row['status']=='FAILED' for row in rows):raise RuntimeError('Window case failed; outputs retained')
    compare()
    for path in (OUT/'cases').iterdir():
        write(path/'SHA256.json',{str(p.relative_to(path)):hashlib.sha256(p.read_bytes()).hexdigest() for p in path.rglob('*') if p.is_file()})


if __name__=='__main__':main()
