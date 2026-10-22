"""Reference-only profiling and immutable full-step fixtures, before any optimization."""
import time,os,io,contextlib,cProfile,pstats,resource,tracemalloc,subprocess,sys,json
import numpy as np,h5py,scipy,warp as wp
from scipy import linalg
from threadpoolctl import threadpool_info
from mevo_mb.fast_sim import f10g_contact as geom, sequential_hard_sphere as sqp, water_sequential as ws, water_b2 as wb
from mevo_mb.control.scheduled_pzt import Schedule
from line_profiler import LineProfiler
from .common import *

class Timers:
 def __init__(self):self.times={};self.original=[];self.bytes={'H2D':0,'D2H':0};self.calls={'H2D':0,'D2H':0}
 def wrap(self,obj,name,label,transfer=None):
  old=getattr(obj,name);self.original.append((obj,name,old))
  def fn(*a,**kw):
   begin=time.perf_counter();cpu=time.process_time();r=old(*a,**kw)
   self.times.setdefault(label,[]).append([time.perf_counter()-begin,time.process_time()-cpu])
   if transfer:
    v=r if transfer=='D2H' else a[1]
    if isinstance(v,np.ndarray):self.bytes[transfer]+=v.nbytes;self.calls[transfer]+=1
   return r
  setattr(obj,name,fn)
 def __enter__(self):
  for obj,name,label in [(geom,'normals','pair_normals'),(ws,'resistance','resistance'),(sqp,'contact_geometry','dense_contact_geometry'),(sqp,'active_set_qp','active_set_QP'),(ws,'sequential_step','SQP_total'),(linalg,'cho_factor','Cholesky_SPD'),(linalg,'cho_solve','linear_solve'),(wb.WaterB2,'fields','fields_total'),(wb.WaterB2,'wall_check','wall_total')]:self.wrap(obj,name,label)
  self.wrap(wp.array,'numpy','D2H_wait_and_copy','D2H');self.wrap(wp.array,'assign','H2D_upload','H2D')
  return self
 def __exit__(self,*args):
  for obj,name,old in reversed(self.original):setattr(obj,name,old)

def save_tree(h,obj):
 for k,v in obj.items():
  if isinstance(v,dict):save_tree(h.create_group(k),v)
  elif isinstance(v,list) and v and isinstance(v[0],dict):
   sub=h.create_group(k)
   for i,row in enumerate(v):save_tree(sub.create_group(str(i)),row)
  elif isinstance(v,(str,bool,int,float,np.number)) or v is None:h.attrs[k]='None' if v is None else v
  else:
   a=np.asarray(v)
   if a.dtype.kind not in 'OUS':h.create_dataset(k,data=a,compression='gzip' if a.ndim else None)
   else:h.attrs[k]=str(v)

def main():
 verify();out=BASE/'baseline';out.mkdir(exist_ok=False);start=time.perf_counter();c=context();cfg=read(CONFIG);meta=read(BASE/'STATE_RULES.json')
 import psutil
 state={}
 with h5py.File(BASE/'STATES.h5') as h:state={k:h[k][:] for k in h}
 schedules={k:Schedule.from_record(v) for k,v in read(SCALING/'RUN_SCHEDULES.json').items()}
 info=dict(identity=identity(),python=sys.version,packages={k:__import__('importlib.metadata').metadata.version(k) for k in ['numpy','scipy','warp-lang','cupy-cuda12x','threadpoolctl','line_profiler']},BLAS=threadpool_info(),cpu=subprocess.check_output(['lscpu'],text=True),gpu=subprocess.check_output(['nvidia-smi'],text=True),threads=psutil.Process().num_threads(),affinity=psutil.Process().cpu_affinity(),ps=subprocess.check_output(['ps','-eo','pid,pcpu,comm','--sort=-pcpu'],text=True).splitlines()[:12])
 stream=io.StringIO()
 with contextlib.redirect_stdout(stream):np.show_config();scipy.show_config()
 (out/'BLAS_CONFIG.txt').write_text(stream.getvalue());write(out/'ENVIRONMENT.json',info)
 records=[];profiles=[]
 with h5py.File(out/'REFERENCE.h5','w') as h:
  for name,x in state.items():
   begin=time.perf_counter();sim=ws.WaterSequential(x,c['physics'],c['response'],c['basis'],c['cavity'],policy=policy());cold=time.perf_counter()-begin
   sim.control(schedules['PAIR'].segments[0].event());acoustic=sim.fields(x)
   A=geom.resistance(x,R)['A'];zeta=sim.zeta
   forces={'ON':acoustic['total_force_N'],'OFF':np.zeros_like(x)}
   if meta[name]['force']=='compression':forces['ON']=(zeta*A@(-20*x).ravel()).reshape(x.shape)
   for action,F in forces.items():
    key=name+'_'+action;root=h.create_group(key)
    B,pairs,gap=geom.contact_geometry(x,R)
    begin=time.perf_counter();ref=sqp.sequential_step(x,F,R,zeta,1e-6,A,policy());elapsed=time.perf_counter()-begin
    save_tree(root,dict(x=x,radius_m=R,zeta=zeta,force_N=F,A=A,B=B,pairs=pairs,gap_m=gap,acoustic_fields=acoustic,reference=ref))
    rec=dict(case=key,N=len(x),force_kind=meta[name]['force'] if action=='ON' else 'OFF',cold_engine_s=cold,reference_step_s=elapsed,accepted=ref['converged'],active=ref['active_count'],iterations=ref['qp_iterations'],status=ref['solver_status'],gates=ref['gates'])
    records.append(rec);print(key,rec,flush=True)
   # Independent repeated physical step at frozen x; never advance a trajectory.
   if meta[name]['force']=='acoustic':
    sim.control(schedules['PAIR'].segments[0].event())
    for _ in range(cfg['warmup']):sim.x=x.copy();sim.step(1e-6)
    elapsed=[];cpu=[]
    for _ in range(cfg['repeat_fixed']):
     sim.x=x.copy();b=time.perf_counter();bc=time.process_time();sim.step(1e-6);elapsed.append(time.perf_counter()-b);cpu.append(time.process_time()-bc)
    sim.x=x.copy()
    with Timers() as timer:r=sim.step(1e-6)
    gpu_events=[]
    def report(rows,name):
     for v in rows:gpu_events.append({k:getattr(v,k) for k in ['name','elapsed'] if hasattr(v,k)})
    sim.x=x.copy()
    with wp.ScopedTimer('physical_step',print=False,cuda_filter=wp.TIMING_KERNEL,report_func=report):sim.step(1e-6)
    sim.x=x.copy();prof=cProfile.Profile();prof.enable();sim.step(1e-6);prof.disable();prof.dump_stats(str(out/(name+'.pstats')))
    with (out/(name+'_cprofile.txt')).open('w') as f:pstats.Stats(prof,stream=f).sort_stats('cumtime').print_stats(50)
    sim.x=x.copy();tracemalloc.start();sim.step(1e-6);mem=tracemalloc.get_traced_memory();tracemalloc.stop()
    profiles.append(dict(case=name,N=len(x),wall=stats(elapsed),CPU_s=stats(cpu),instrumented_inclusive={k:np.sum(v,axis=0) for k,v in timer.times.items()},transfer_bytes=timer.bytes,transfer_calls=timer.calls,kernels=gpu_events,tracemalloc_current_peak_bytes=mem,RSS_MiB=psutil.Process().memory_info().rss/2**20,threads=psutil.Process().num_threads(),BLAS=threadpool_info(),instrumentation_scope='Inclusive nested times; separate warmed uninstrumented total; fused acoustic/APECSS kernel cannot be split without changing implementation'))
   if name in ['N128_historical_end','N64_historical_end']:
    lp=LineProfiler(geom.normals,geom.resistance,geom.contact_geometry,sqp.sequential_step,sqp.active_set_qp)
    sim.x=x.copy();lp(sim.step)(1e-6)
    with (out/(name+'_lines.txt')).open('w') as f:lp.print_stats(stream=f)
 write(out/'FIXTURES.json',records);write(out/'PROFILING.json',profiles)
 write(out/'COMPLETE.json',dict(status='COMPLETE',wall_s=time.perf_counter()-start,reference_sha256=sha(out/'REFERENCE.h5'),identity=identity(),RSS_peak_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
if __name__=='__main__':main()
