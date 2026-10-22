"""Same frozen states, all registered backends, warm latency and independent gates."""
import time,traceback,os,resource
import numpy as np,h5py,warp as wp,psutil
from threadpoolctl import threadpool_limits,threadpool_info
from mevo_mb.fast_sim import f10g_contact as old
from mevo_mb.fast_sim.sequential_hard_sphere import sequential_step as original
from mevo_mb.fast_sim.accelerated_sequential import sequential_step as fast
from mevo_mb.fast_sim import accelerated_contact as opt
from mevo_mb.fast_sim.water_sequential import WaterSequential
from mevo_mb.fast_sim.accelerated_water import WaterAccelerated
from mevo_mb.fast_sim.accelerated_fields import WaterPacked
from mevo_mb.control.scheduled_pzt import Schedule
from scripts.straight_tube_exp1.f10z.run import Monitor
from .profile import Timers
from .common import *

def gpu_geometry(x,radius):
    import cupy as cp
    q=cp.asarray(x,dtype=cp.float64);i,j=cp.triu_indices(len(x),1);dx=q[j]-q[i]
    d=cp.linalg.norm(dx,axis=1);n=dx/d[:,None];gap=d-2*radius;e=gap/radius;mask=e<.1
    ii,jj,nn,ee=i[mask],j[mask],n[mask],e[mask]
    def phi(v):return .5/v-.45*cp.log(v)-(3/56)*v*cp.log(v)
    c=.5*(phi(cp.maximum(ee,.001))-phi(cp.asarray(.1)))
    blocks=c[:,None,None]*nn[:,:,None]*nn[:,None,:];A=cp.eye(3*len(x));axis=cp.arange(3);u=3*ii[:,None]+axis;v=3*jj[:,None]+axis
    cp.add.at(A,(u[:,:,None],u[:,None,:]),blocks);cp.add.at(A,(v[:,:,None],v[:,None,:]),blocks)
    A[u[:,:,None],v[:,None,:]]-=blocks;A[v[:,:,None],u[:,None,:]]-=blocks
    return A.get(),n.get(),gap.get()

def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='benchmark_replay');ap.add_argument('--reuse-fixed',default=None);a=ap.parse_args()
    verify();cfg=read(CONFIG);out=BASE/a.output;out.mkdir(exist_ok=False);start=time.perf_counter()
    monitor=Monitor();monitor.thread.start();records=[];geometry=[];full=[];linear=[]
    resume=BASE/a.reuse_fixed if a.reuse_fixed else None
    if resume:
      records=read(resume/'FIXED_BENCHMARK.json');geometry=read(resume/'GEOMETRY.json');linear=read(resume/'LINEAR_ALGEBRA.json')
      assert len(records)==96 and all(r['status']=='PASS' for r in records)
      write(out/'RESUME.json',dict(source=str(resume),sha256={p.name:sha(p) for p in resume.glob('*.json')}))
      write(out/'FIXED_BENCHMARK.json',records)
    c=context();p=policy();schedule=Schedule.from_record(read(SCALING/'RUN_SCHEDULES.json')['PAIR'])
    with h5py.File(BASE/'baseline/REFERENCE.h5') as h:
      for name,q in ([] if resume else h.items()):
        x,F,A=q['x'][:],q['force_N'][:],q['A'][:];z=q.attrs['zeta'];refx=q['reference/positions_m'][:]
        if name.endswith('_OFF'):settings=[('B0',1),('B1',1)]
        else:settings=[('B0',1),('B1',1),('B2_t2',2),('B3_CuPy',1)]
        for label,threads in settings:
          print('BEGIN',name,label,threads,flush=True)
          times=[];cpu=[];errors=[];status='PASS';failure=None;r=None
          with threadpool_limits(limits=threads):
            try:
              for k in range(cfg['warmup']+cfg['repeat_fixed']):
                t=time.perf_counter();tc=time.process_time()
                op=old.resistance(x,R) if label=='B0' else opt.resistance(x,R)
                r=original(x,F,R,z,1e-6,op['A'],p) if label=='B0' else fast(x,F,R,z,1e-6,op['A'],p,backend='cupy' if label=='B3_CuPy' else 'cpu')
                wall=time.perf_counter()-t;ct=time.process_time()-tc
                if k>=cfg['warmup']:times.append(wall);cpu.append(ct)
                assert r['converged']==q['reference'].attrs['converged'] and all(r['gates'].values())
                err=np.linalg.norm(r['positions_m']-refx,axis=1).max();errors.append(err);assert err<p['fixed_id_tolerance_m']
                By,_,gap=old.contact_geometry(r['positions_m'],R)
                balance=(z*A@r['velocity_m_s'].ravel()-F.ravel()-By.T@r['pair_force_N']).reshape(x.shape)
                assert np.linalg.norm(balance,axis=1).max()<=p['force_residual_atol_N']+p['force_residual_rtol']*r['force_scale_N']
            except Exception as ex:status='FAIL';failure=repr(ex);(out/(name+'_'+label+'_exception.txt')).write_text(traceback.format_exc())
            records.append(dict(case=name,N=len(x),backend=label,threads=threads,BLAS=threadpool_info(),status=status,failure=failure,
                wall=stats(times) if times else None,cpu=stats(cpu) if cpu else None,max_fixed_ID_m=float(max(errors,default=0)),
                active=r['active_count'] if r else None,iterations=r['qp_iterations'] if r else None,
                force_residual_N=r['force_equation_residual_N'] if r else None,complementarity_m=r['actual_complementarity_equivalent_gap_m'] if r else None))
        print(name,[(r['backend'],r['status'],round(r['wall']['p50']*1e3,3) if r['wall'] else None) for r in records if r['case']==name],flush=True)
        write(out/'FIXED_BENCHMARK.json',records)
        write(out/'GEOMETRY.json',geometry);write(out/'LINEAR_ALGEBRA.json',linear)
        if name.endswith('_ON'):
          row=dict(case=name,N=len(x),dense_bytes=q['B'].size*8,matrix_free_bytes=opt.contact_geometry(x,R)[0].nbytes,times={})
          for label,fn in [('dense_B',lambda:old.contact_geometry(x,R)),('matrix_free_B',lambda:opt.contact_geometry(x,R)),('old_resistance',lambda:old.resistance(x,R)),('new_resistance',lambda:opt.resistance(x,R)),('GPU_geometry_resistance_transfer',lambda:gpu_geometry(x,R))]:
            tt=[]
            try:
              for k in range(cfg['warmup']+cfg['repeat_fixed']):
                t=time.perf_counter();value=fn();elapsed=time.perf_counter()-t
                if k>=cfg['warmup']:tt.append(elapsed)
              row['times'][label]=stats(tt)
              if label.startswith('GPU'):row['GPU_A_error']=float(abs(value[0]-A).max());row['GPU_gap_error']=float(abs(value[2]-q['gap_m'][:]).max())
            except Exception as ex:row['times'][label]=dict(status='FAIL',reason=repr(ex))
          geometry.append(row)
        if name in ['N32_dispersed_ON','N64_historical_end_ON','N128_historical_end_ON','N128_fcc_ON']:
          rhs=q['B'][:8].T
          for backend in ['cpu','cupy']:
           for threads in ([1,2] if backend=='cpu' else [1]):
            fac=[];sol=[];error=[]
            with threadpool_limits(limits=threads):
             for k in range(cfg['warmup']+cfg['repeat_fixed']):
              t=time.perf_counter();f=opt.factorize(A,backend);tf=time.perf_counter()-t;t=time.perf_counter();v=f.solve(rhs);ts=time.perf_counter()-t
              if k>=cfg['warmup']:fac.append(tf);sol.append(ts);error.append(float(abs(A@v-rhs).max()))
            linear.append(dict(case=name,backend=backend,threads=threads,factor_s=stats(fac),solve_s=stats(sol),max_linear_residual=max(error),rhs_columns=8))
    write(out/'GEOMETRY.json',geometry);write(out/'LINEAR_ALGEBRA.json',linear)
    # Independent original Warp field/wall time including transfers; reset same initial state every time.
    with h5py.File(BASE/'STATES.h5') as h:
      for name in ['N32_dispersed','N64_historical_end','N128_dispersed','N128_historical_end']:
        x=h[name][:]
        for label,cls,kw in [('B0',WaterSequential,{}),('B1',WaterAccelerated,{}),('B3_CuPy',WaterAccelerated,{'backend':'cupy'}),('PACKED',WaterPacked,{})]:
          t=time.perf_counter();sim=cls(x,c['physics'],c['response'],c['basis'],c['cavity'],policy=p,**kw);sim.control(schedule.segments[0].event());cold=time.perf_counter()-t
          tt=[];ct=[];valid=True
          for k in range(cfg['warmup']+cfg['repeat_fixed']):
            sim.x=x.copy();t=time.perf_counter();tc=time.process_time();r=sim.step(1e-6);elapsed=time.perf_counter()-t;cpu=time.process_time()-tc;valid &= r['accepted']
            if k>=cfg['warmup']:tt.append(elapsed);ct.append(cpu)
          sim.x=x.copy()
          with Timers() as timer:sim.step(1e-6)
          kernel=[]
          def report(rows,**kwargs):
            for v in rows:kernel.append(dict(name=v.name,elapsed_ms=v.elapsed))
          sim.x=x.copy()
          with wp.ScopedTimer('GPU trace',print=True,cuda_filter=wp.TIMING_KERNEL,report_func=report):sim.step(1e-6)
          full.append(dict(case=name,N=len(x),backend=label,status='PASS' if valid else 'FAIL',cold_engine_s=cold,wall=stats(tt),cpu=stats(ct),transfer_bytes=timer.bytes,transfer_calls=timer.calls,kernels=kernel,RSS_MiB=psutil.Process().memory_info().rss/2**20))
          write(out/'FULL_STEP.json',full)
    # Explicit Warp/CuPy pointer and lifetime check, no mutation of physical state.
    import cupy as cp
    try:
      wp.synchronize();view=cp.asarray(sim.sim.x);interop=dict(pointer_equal=int(view.data.ptr)==sim.sim.x.ptr,exact_values=np.array_equal(view.get(),sim.sim.x.numpy()),owner_retained=True,stream_policy='wp.synchronize before borrowing; no external stream production pipeline')
    except Exception as ex:interop=dict(status='FAIL',reason=repr(ex))
    write(out/'INTEROP.json',interop);write(out/'GEOMETRY.json',geometry);write(out/'LINEAR_ALGEBRA.json',linear);write(out/'FULL_STEP.json',full)
    monitor.done.set();monitor.thread.join(4);write(out/'resources.json',monitor.rows)
    write(out/'COMPLETE.json',dict(status='COMPLETE',wall_s=time.perf_counter()-start,identity=identity(),reference_sha256=sha(BASE/'baseline/REFERENCE.h5')))
if __name__=='__main__':main()
