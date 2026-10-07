"""Persistent Warp execution, CUDA graphs, event-aligned fixed macro steps."""
from time import perf_counter
import numpy as np
import warp as wp
from . import kernels as k


def complex_vectors(a):
    a=np.asarray(a); return np.ascontiguousarray(np.stack([a.real,a.imag],axis=-1))


class FastSimulator:
    def __init__(self, positions_m, physics, response, basis=None, device='cuda:0', graph_chunk=32):
        wp.init()
        self.device=wp.get_device(device)
        if not self.device.is_cuda:
            raise ValueError('FastSimulator requires the explicitly selected CUDA backend')
        self.p,self.response,self.basis=physics,response,basis
        x=np.asarray(positions_m,dtype=np.float64)
        if x.ndim==2: x=x[None,:,:]
        if x.ndim!=3 or x.shape[-1]!=3 or not np.isfinite(x).all(): raise ValueError('Expected finite [B,N,3]')
        self.initial=x.copy(); self.batch,self.n=x.shape[:2]; self.total=self.batch*self.n
        self.graph_chunk=int(graph_chunk); self.graphs={}; self.capture_s=0.0
        if not 1<=self.graph_chunk<=128: raise ValueError('Bounded graph chunk required')
        def arr(a,dtype): return wp.array(np.ascontiguousarray(a),dtype=dtype,device=self.device)
        self.x=arr(x.reshape(-1,3),wp.vec3d)
        self.tmp1=wp.empty_like(self.x); self.tmp2=wp.empty_like(self.x)
        self.force=wp.zeros_like(self.x); self.free=wp.zeros_like(self.x)
        self.vel1=wp.zeros_like(self.x); self.vel2=wp.zeros_like(self.x)
        self.harmonics=wp.zeros((self.total,8),dtype=wp.vec2d,device=self.device)
        self.pressure=wp.zeros(self.total,dtype=wp.vec2d,device=self.device)
        self.knots=arr(response.knots,wp.float64); self.values=arr(complex_vectors(response.values),wp.vec2d)
        self.valid=arr(response.valid,wp.int32); self.kv=arr(complex_vectors(response.knot_values),wp.vec2d)
        self.kval=arr(response.knot_valid,wp.int32)
        self.error=wp.zeros(1,dtype=wp.int32,device=self.device)
        self.diag=wp.zeros((self.total,10),dtype=wp.float64,device=self.device)
        self.uniform_p=wp.zeros(1,dtype=wp.vec2d,device=self.device)
        self.shape=(2,2,2) if basis is None else basis.shape
        self.origin=wp.vec3d(0.0) if basis is None else wp.vec3d(*basis.origin)
        self.spacing=wp.vec3d(1.0) if basis is None else wp.vec3d(*basis.spacing)
        self.grid=wp.zeros((int(np.prod(self.shape)),4),dtype=wp.vec2d,device=self.device)
        self.reset()
        begin=perf_counter(); wp.load_module(module=k,device=self.device); wp.synchronize_device(self.device)
        self.compile_s=perf_counter()-begin

    def reset(self):
        self.x.assign(self.initial.reshape(-1,3))
        diag=np.zeros((self.total,10)); diag[:,0]=np.inf
        self.diag.assign(diag); self.error.zero_()

    def control(self, event):
        if self.basis is None:
            p=float(event.get('pressure_Pa',0))*np.exp(1j*float(event.get('phase_rad',0)))
            self.uniform_p.assign(complex_vectors([p]))
            return {'time_s':event['time_s'],'uniform_pressure_real_Pa':float(p.real),'uniform_pressure_imag_Pa':float(p.imag)}
        field, grid, u=self.basis.combine(event['amplitude_Pa'],event['phase_rad'],event.get('enabled'))
        self.grid.assign(grid.reshape(-1,4,2))
        self.last_field=field
        return dict(event,control_real_Pa=u.real.tolist(),control_imag_Pa=u.imag.tolist())

    def launch(self, kernel, dim, args):
        wp.launch(kernel,dim=dim,inputs=args,device=self.device)

    def fields(self):
        self.launch(k.field_response,self.total,[self.x,self.grid,self.origin,self.spacing,*self.shape,
            int(self.basis is None),self.uniform_p,self.knots,self.values,self.valid,self.kv,self.kval,
            self.harmonics,self.force,self.pressure,self.diag,self.error])

    def pairs(self):
        p=self.p
        if p.secondary:
            self.launch(k.pair_force,(self.batch,self.n,self.n),[self.x,self.harmonics,self.force,self.n,
                p.radius_m,2*np.pi*p.frequency_Hz,p.density_kg_m3,int(p.near_field!='classical_far_field'),
                p.transition_over_sum,p.soft_core_over_sum,self.error])

    def step(self,dt):
        p=self.p; self.fields(); self.pairs()
        self.launch(k.free_velocity,self.total,[self.force,self.free,6*np.pi*p.viscosity_Pa_s*p.radius_m])
        velocity=self.free
        if p.local_lubrication:
            for it in range(p.mobility_iterations):
                out=self.vel1 if it%2==0 else self.vel2
                self.launch(k.local_mobility,self.total,[self.x,self.free,velocity,out,self.n,p.radius_m,
                    p.lubrication_activation,p.lubrication_floor])
                velocity=out
        self.launch(k.predict,self.total,[self.x,velocity,self.tmp1,dt,self.diag])
        x=self.tmp1
        for it in range(p.contact_iterations):
            out=self.tmp2 if it%2==0 else self.tmp1
            self.launch(k.contact,self.total,[x,out,self.n,p.radius_m,p.contact_relaxation,self.diag,self.error])
            x=out
        self.launch(k.accept,self.total,[x,self.x,self.n,p.radius_m,self.diag,self.error])

    def graph(self,dt,count):
        key=(float(dt),int(count))
        if key not in self.graphs:
            begin=perf_counter()
            with wp.ScopedCapture(device=self.device) as capture:
                for _ in range(count): self.step(dt)
            self.graphs[key]=capture.graph; self.capture_s+=perf_counter()-begin
        return self.graphs[key]

    def advance(self,dt,count,graphs=True):
        while count:
            size=min(count,self.graph_chunk)
            if graphs: wp.capture_launch(self.graph(dt,size))
            else:
                for _ in range(size): self.step(dt)
            count-=size

    def state(self):
        code=int(self.error.numpy()[0])
        if code: raise RuntimeError(f'GPU validity error bits={code}: field=1 table=2 coincident=4 nonfinite=8')
        return self.x.numpy().reshape(self.batch,self.n,3)

    def quality(self):
        d=self.diag.numpy()
        return dict(pressure_min_Pa=float(d[:,0].min()),pressure_max_Pa=float(d[:,1].max()),
                    max_predicted_penetration_m=float(d[:,2].max()),max_accepted_penetration_m=float(d[:,3].max()),
                    cumulative_correction_m=float(d[:,4].sum()),cumulative_physical_displacement_m=float(d[:,5].sum()),
                    max_contact_pass_correction_m=float(d[:,6].max()),max_speed_m_s=float(d[:,7].max()),
                    max_physical_substep_displacement_m=float(d[:,9].max()),error_bits=int(self.error.numpy()[0]),
                    correction_fraction=float(d[:,4].sum()/max(d[:,5].sum(),1e-300)))

    def run(self,duration_s,dt_s,events,sample_interval_s=1e-4,wall_budget_s=180,graphs=True,substeps=1):
        if duration_s<=0 or dt_s<=0 or substeps not in (1,2,4): raise ValueError('Invalid bounded time protocol')
        if not events or events[0]['time_s']!=0 or any(b['time_s']<=a['time_s'] for a,b in zip(events,events[1:])):
            raise ValueError('Controls must start at zero and be strictly ordered')
        if any(e['time_s']>=duration_s for e in events[1:]): raise ValueError('Control outside duration')
        # Compile and construct common CUDA graphs before timing; no trajectory is advanced during capture.
        self.control(events[0]); self.graph(dt_s/substeps,self.graph_chunk); self.graph(dt_s/substeps,1)
        wp.synchronize_device(self.device)
        if not np.array_equal(self.state(),self.initial): raise RuntimeError('Graph capture changed physical state')
        self.reset(); times=[0.0]; frames=[self.initial.copy()]; control_index=[0]; controls=[]
        macro=0; actual=0; started=perf_counter(); status='COMPLETE'; reason=None
        capture_before=self.capture_s; control_cost=0.0; sample_cost=0.0
        try:
            for seg,event in enumerate(events):
                start=float(event['time_s']); end=float(events[seg+1]['time_s']) if seg+1<len(events) else duration_s
                t0=perf_counter(); controls.append(self.control(event)); control_cost+=perf_counter()-t0
                if times[-1]==start: control_index[-1]=seg
                steps=int(np.floor((end-start)/dt_s+1e-10)); stride=max(1,int(round(sample_interval_s/dt_s)))
                done=0
                while done<steps:
                    count=min(stride,steps-done)
                    self.advance(dt_s/substeps,count*substeps,graphs)
                    done+=count; macro+=count; actual+=count*substeps
                    t0=perf_counter(); frame=self.state(); sample_cost+=perf_counter()-t0
                    frames.append(frame); times.append(start+done*dt_s); control_index.append(seg)
                    if perf_counter()-started>wall_budget_s: raise TimeoutError('Finite wall budget exhausted')
                    if self.quality()['max_accepted_penetration_m']>0.2*self.p.radius_m:
                        raise RuntimeError('Severe overlap: stop, retain failed attempt')
                remaining=end-(start+steps*dt_s)
                if remaining>max(1e-15,dt_s*1e-10):
                    self.advance(remaining/substeps,substeps,graphs); macro+=1; actual+=substeps
                    frames.append(self.state()); times.append(end); control_index.append(seg)
                else: times[-1]=end
        except (RuntimeError,TimeoutError) as exc:
            status='FAILED'; reason=str(exc)
        wp.synchronize_device(self.device); elapsed=perf_counter()-started
        return dict(time_s=np.array(times),positions_m=np.array(frames),control_index=np.array(control_index),
                    summary=dict(status=status,reason=reason,completed_physical_s=float(times[-1]),wall_s=elapsed,
                    macro_step_count=macro,actual_substep_count=actual,effective_steps_per_physical_second=actual/times[-1] if times[-1]>0 else None,
                    nominal_dt_s=dt_s,substeps_per_macro=substeps,compile_s=self.compile_s,graph_capture_s=self.capture_s,
                    in_run_graph_capture_s=self.capture_s-capture_before,control_update_s=control_cost,
                    sampled_download_s=sample_cost,batch_size=self.batch,particle_count_per_sample=self.n,
                    quality=self.quality(),controls=controls,backend='NVIDIA Warp',dtype='float64',model=self.p.metadata()))
