"""Parallel guarded ADR benchmarks, then explicitly gated erosion jobs."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']: os.environ[key]='1'
import argparse, contextlib, copy, hashlib, json, multiprocessing, resource, time, traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from pd_clot.regularization.equilibrium_adr import make_fixture, equilibrate
from pd_clot.regularization.erosion_equilibrium_runner import write, snapshot, save_equilibrium, run
from verify_mechanical_equilibrium import capacity

ROOT=Path(__file__).resolve().parents[1]
V=ROOT/'verification/mechanical_equilibrium_v1/adr_completion'
OLD=V.parent
CONFIG=ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json'


def read(path):return json.loads(path.read_text())


def work(job):
    os.environ['NUMBA_CACHE_DIR']=str(ROOT/f'build_mechanical_equilibrium_v1/adr_completion/worker_{os.getpid()}')
    c=read(CONFIG);start=time.perf_counter();cpu=time.process_time()
    out=(V/'cases'/job['name']) if job['kind']=='small' else ROOT/'results/streaming_erosion_equilibrium_v1'/job['name']
    out.parent.mkdir(parents=True,exist_ok=True)
    with (out.parent/(out.name+'.log')).open('x') as f,contextlib.redirect_stdout(f),contextlib.redirect_stderr(f):
        try:
            if job['kind']=='small':
                c['clot']['cells']=c['mechanical_equilibrium']['small_fixture_cells'];c['damage']['enabled']=False
                c['mechanical_equilibrium']['adr']['additional_iterations_after_convergence']=job.get('extra',0)
                if job.get('resolved'):
                    c['streaming'].update(provider='resolved_field',field_path=str(OLD/'cases/small_scipy_resolved/static_field.vtu'),
                        metadata_path=dict(length_unit='m',velocity_unit='m/s',pressure_unit='Pa',time_unit='s',normal_convention='outward_solid',sample_domain='convex_hull'))
                out.mkdir(exist_ok=False);snapshot(out,c)
                cloud=make_fixture(c);g=np.ones(len(cloud.pairs))
                result=equilibrate(cloud,cloud.X,np.zeros_like(cloud.X),g,c)
                save_equilibrium(out/'equilibrium',cloud,g,c,result)
                s=dict(result.summary,wall_time_s=time.perf_counter()-start,cpu_time_s=time.process_time()-cpu,
                       peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
                write(out/'SUMMARY.json',s)
            else:
                de=job['step'];c['streaming_erosion'].update(erosion_exposure_step_s=de,number_of_erosion_steps=round(.1/de))
                s=run(c,out)
            print(json.dumps(s,indent=2),flush=True)
            return dict(name=job['name'],status=s['status'],wall_time_s=time.perf_counter()-start)
        except Exception as exc:
            traceback.print_exc();out.mkdir(exist_ok=True)
            write(out/'FAILED.json',dict(error=repr(exc),traceback=traceback.format_exc()))
            return dict(name=job['name'],status='FAILED',error=repr(exc))


def batch(jobs,label):
    workers=capacity(len(jobs));start=time.perf_counter();rows=[]
    print('Workers',workers,flush=True)
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        for f in as_completed([pool.submit(work,j) for j in jobs]):
            row=f.result();rows.append(row);print(json.dumps(row),flush=True)
    write(V/(label+'_BATCH.json'),dict(workers=workers,jobs=rows,wall_time_s=time.perf_counter()-start))
    return rows


def compare_small():
    comparisons=[];threshold=read(CONFIG)['mechanical_equilibrium']['reference_acceptance']
    for name,refname in [('small_adr','small_scipy'),('small_adr_extra500','small_scipy'),('small_adr_resolved','small_scipy_resolved')]:
        p=V/'cases'/name;s=read(p/'SUMMARY.json');ref=read(OLD/'cases'/refname/'SUMMARY.json')
        z=np.load(p/'equilibrium/state.npz');rz=np.load(OLD/'cases'/refname/'states.npz');free=~z['fixed']
        assert np.array_equal(z['X'],rz['X']) and np.array_equal(z['pairs'],rz['pairs'])
        diff=(z['x']-rz['x'])[free];ref_u=(rz['x']-rz['X'])[free]
        rms=float(np.sqrt(np.mean(np.sum(diff**2,axis=1))));refmax=float(np.linalg.norm(ref_u,axis=1).max())
        driver=abs(s['final']['maximum_epsilon_plus']-ref['final']['maximum_epsilon_plus'])/ref['final']['maximum_epsilon_plus']
        comparisons.append(dict(name=name,reference=refname,rms_displacement_difference_m=rms,
            maximum_displacement_difference_m=float(np.linalg.norm(diff,axis=1).max()),
            relative_L2_displacement_difference=float(np.linalg.norm(diff)/np.linalg.norm(ref_u)),
            rms_difference_over_reference_max=rms/refmax,relative_max_driver_difference=driver,
            adr=s['final'],reference_final=ref['final'],passed=bool(s['converged'] and ref['converged'] and
                s['final']['relative_force_residual']<=threshold['relative_force_residual'] and
                rms/refmax<=threshold['rms_displacement_difference_over_reference_max'] and driver<=threshold['relative_max_driver_difference'])))
    write(V/'SMALL_REFERENCE_GATE.json',dict(passed=all(r['passed'] for r in comparisons),comparisons=comparisons))


def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['small','erosion']);a=p.parse_args()
    assert hashlib.sha256(CONFIG.read_bytes()).hexdigest()==read(V/'provenance/GUARDED_PLAN.json')['config_sha256']
    if a.phase=='small':
        rows=batch([dict(name='small_adr',kind='small'),dict(name='small_adr_extra500',kind='small',extra=500),
                    dict(name='small_adr_resolved',kind='small',resolved=True)],'SMALL')
        if any(r['status']=='FAILED' for r in rows):raise RuntimeError('Inspect retained failed cases')
        compare_small()
    else:
        if not read(V/'SMALL_REFERENCE_GATE.json')['passed']:raise RuntimeError('Independent equilibrium gate NOT passed')
        batch([dict(name=f'dt_erosion_{s}',kind='erosion',step=v) for s,v in [('0p02',.02),('0p01',.01),('0p005',.005)]],'EROSION')


if __name__=='__main__':main()
