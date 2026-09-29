"""Shared-time CUDA population dynamics with explicit small-subset CPU authority.

FEM, state, finite triangle geometry, isolated P9-A resistance/contact, and
accepted state updates execute on GPU. Rare interacting components and hard
finite-edge path partitions use original CPU functions on those particles only.
No full population round-trip occurs in the dynamics path.
"""
from collections import Counter
from types import SimpleNamespace
import numpy as np
import torch
from ..particle_shapes import Sphere
from ..particle6_stepper import bind_query_dependency
from ..nearfield_regularization import NearFieldRegularizationV1
from ..nearfield_handoff import handoff_constraints,first_handoff_event,partitioned_wall_handoff_certificate,wall_handoff_certificate
from ..physical_time_refinement import TrialNeedsSubdivision,PhysicalTimeRefinementError
from .engine import PopulationEngine
from .gpu_state import GPUActivePool
from .gpu_field import GPUField
from .gpu_wall import GPUWall
from .gpu_neighbors import neighbors
from .gpu_resistance import assemble_isolated,solve_isolated
from .components import assemble_components,solve_components
from .pairs import candidates,sphere_gap
from .gpu_safety import GPUSafetyMonitor


def advance_quaternion(q,omega,dt):
    old=q/torch.linalg.vector_norm(q,dim=1)[:,None]
    magnitude=torch.linalg.vector_norm(omega,dim=1);half=magnitude*dt/2
    factor=(dt/2)*torch.sinc(half/torch.pi)
    scalar=torch.cos(half)[:,None];vector=omega*factor[:,None]
    out=torch.cat([scalar*old[:,:1]-(vector*old[:,1:]).sum(1)[:,None],
        scalar*old[:,1:]+old[:,:1]*vector+torch.linalg.cross(vector,old[:,1:])],dim=1)
    out/=torch.linalg.vector_norm(out,dim=1)[:,None]
    out=torch.where(((out*old).sum(1)<0)[:,None],-out,out)
    return torch.where((magnitude==0)[:,None],old,out)


class GPUPopulationEngine(PopulationEngine):
    def __init__(self,environment,*,dt=.0005,pairs=True,realization_id=0,recorder=None,
                 gpu_field=None,gpu_wall=None,cache_path=None):
        super().__init__(environment,dt=dt,pairs=pairs,realization_id=realization_id,recorder=recorder)
        self.pool=GPUActivePool();self.pairs=bool(pairs)
        self.safety=GPUSafetyMonitor(environment,self.dt,self.pool.device)
        self.field=gpu_field or GPUField(environment.field,cache_path=cache_path)
        self.wall=None if environment.wall is None else (gpu_wall or GPUWall(environment.wall))
        self.policy=NearFieldRegularizationV1();self.backend_counts=Counter();self.last_audit=None
        self.cached_field=None;self.cached_wall=None;self.device_output=[]
        self.prefetched=None
        if environment.classifier is not None:
            boxes=[]
            for role in sorted(set(environment.classifier.roles)):
                t=environment.classifier.xyz[environment.classifier.roles==role]
                boxes.append(np.r_[t.min((0,1)),t.max((0,1))])
            self.cap_boxes=torch.tensor(np.array(boxes),device=self.pool.device)
        else:self.cap_boxes=None

    def insert_accepted(self,particle,event):
        if event.get('realization_id',self.realization_id)!=self.realization_id:raise ValueError('Foreign realization')
        self.pool.insert(particle,event,self.time);self.births[particle.particle_id]=dict(event)
        self.safety.insert(self.pool.x[self.pool.n-1],self.time)
        self.counts['accepted_birth_count']+=1;self.maximum_active=max(self.maximum_active,self.pool.n)
        self.cached_field=None;self.cached_wall=None;self.prefetched=None;self.record()

    def admission_shapes(self,position,radius):
        p=self.pool;n=p.n
        if not n:return {}
        center=torch.tensor(position,device=p.device,dtype=torch.float64)
        given_radius=torch.tensor(radius,device=p.device,dtype=torch.float64)
        lower=torch.maximum(torch.full_like(p.radius[:n],2e-9),.001*torch.minimum(p.radius[:n],given_radius))
        scale=torch.maximum(p.x[:n].abs().max(1).values,center.abs().max())
        ro=512*torch.finfo(torch.float64).eps*torch.maximum(scale,torch.maximum(p.radius[:n],given_radius))
        potential=torch.nonzero(torch.linalg.vector_norm(p.x[:n]-center,dim=1)<=p.radius[:n]+radius+lower+ro).flatten()
        self.backend_counts['birth_pair_subset_particles_to_cpu']+=len(potential)
        return self._shapes(potential)

    def _shapes(self,indices):
        if not len(indices):return {}
        p=self.pool;ids=p.particle_id[indices].cpu().numpy();x=p.x[indices].cpu().numpy();r=p.radius[indices].cpu().numpy()
        return {int(i):Sphere(pos,float(a)) for i,pos,a in zip(ids,x,r)}

    def record(self):
        state=self.pool.device_snapshot(self.time)
        if self.recorder is None:
            # Explicit validation-only snapshots; production uses delayed output.
            self.snapshots.append({k:v.cpu().numpy() if isinstance(v,torch.Tensor) else v for k,v in state.items()})
        else:
            self.device_output.append(state)
            if len(self.device_output)>=64:self.flush_output()

    def flush_output(self):
        if not self.device_output:return
        rows=self.device_output;self.device_output=[];sizes=[len(r['particle_id']) for r in rows]
        keys=[k for k in rows[0] if k!='physical_time_s']
        host={k:torch.cat([r[k] for r in rows]).cpu().numpy() for k in keys}
        offset=0
        for row,size in zip(rows,sizes):
            self.recorder.snapshot(dict(physical_time_s=row['physical_time_s'],**{k:v[offset:offset+size] for k,v in host.items()}));offset+=size
        self.backend_counts['batched_output_transfers']+=1

    def _cpu_components(self,indices,free):
        p=self.pool;shapes=self._shapes(indices);ids=p.particle_id[indices].cpu().tolist()
        lf={int(i):v for i,v in zip(ids,free[indices].cpu().numpy())}
        system=assemble_components(shapes,lf,self.env.mu,self.env.wall,policy=self.policy,
            gradient_at=lambda x:self.env.field.sample(x).velocity_gradient_s_inv)
        constraints,_=bind_query_dependency(handoff_constraints,dict(all_pairs=candidates,pair_gap=sphere_gap))(shapes,self.env.wall,self.policy)
        solved=solve_components(system,particles=shapes,constraints=constraints)
        mapping={int(i):k for k,i in enumerate(system.ids)}
        result=np.array([solved.velocity[6*mapping[int(i)]:6*mapping[int(i)]+6] for i in ids])
        self.backend_counts['cpu_component_particle_solves']+=len(indices)
        return torch.tensor(result,device=p.device),solved.record

    def _outlets(self,velocity,dt):
        if self.cap_boxes is None:return dt,{}
        p=self.pool;x=p.x[:p.n];end=x+dt*velocity;pad=float(self.env.classifier.tolerance_m.max())
        lo=torch.minimum(x,end);hi=torch.maximum(x,end)
        possible=((lo[:,None,:]<=self.cap_boxes[None,:,3:]+pad)&(hi[:,None,:]>=self.cap_boxes[None,:,:3]-pad)).all(2).any(1)
        indices=torch.nonzero(possible).flatten();hits=[]
        if len(indices):
            starts=x[indices].cpu().numpy();ends=end[indices].cpu().numpy();ix=indices.cpu().tolist()
            for k,a,b in zip(ix,starts,ends):
                hit=self.env.classifier.first_event(a,b)
                if hit is not None:hits.append((hit.segment_fraction,k,hit.role))
        self.backend_counts['cpu_outlet_segment_queries']+=len(indices)
        if not hits:return dt,{}
        first=min(h[0] for h in hits)
        events={k:role for f,k,role in hits if abs(f-first)<=64*np.finfo(float).eps}
        return dt*first,events

    def _pair_root(self,velocity,dt):
        p=self.pool
        if not self.pairs or p.n<2:return None,None
        pairs,_,_=neighbors(p.x[:p.n],p.radius[:p.n],p.realization_id[:p.n],delta=dt*velocity)
        if not len(pairs):return None,pairs
        i,j=pairs[:,0],pairs[:,1];relative=p.x[i]-p.x[j];dv=velocity[i]-velocity[j]
        den=(dv*dv).sum(1);upper=torch.where(den>0,-(relative*dv).sum(1)/den,torch.zeros_like(den)).clamp(0.,dt)
        hmin=torch.linalg.vector_norm(relative+upper[:,None]*dv,dim=1)-p.radius[i]-p.radius[j]
        h0=torch.linalg.vector_norm(relative,dim=1)-p.radius[i]-p.radius[j]
        lower=torch.maximum(torch.full_like(hmin,2e-9),.001*torch.minimum(p.radius[i],p.radius[j]))
        scale=torch.stack([p.x[i].abs().max(1).values,p.x[j].abs().max(1).values,p.radius[i],p.radius[j]]).max(0).values
        pad=512*torch.finfo(torch.float64).eps*scale
        possible=pairs[(hmin-lower < -pad)&(h0-lower>pad)]
        if not len(possible):return None,pairs
        ix=torch.unique(possible);shapes=self._shapes(ix)
        ids=p.particle_id[ix].cpu().tolist();vs=velocity[ix].cpu().numpy()
        edge=[(int(p.particle_id[i]),int(p.particle_id[j])) for i,j in possible.cpu().tolist()]
        root=bind_query_dependency(first_handoff_event,dict(all_pairs=lambda s:edge))(
            shapes,{int(i):v for i,v in zip(ids,vs)},dt,None,self.policy)
        self.backend_counts['cpu_pair_sweep_particles']+=len(ix)
        return root,pairs

    def _pair_certificate(self,velocity,dt,pairs):
        if pairs is None or not len(pairs):return
        p=self.pool;i,j=pairs[:,0],pairs[:,1];relative=p.x[i]-p.x[j];delta=dt*(velocity[i]-velocity[j])
        denominator=(delta*delta).sum(1)
        t=torch.where(denominator>0,-(relative*delta).sum(1)/denominator,torch.zeros_like(denominator)).clamp(0.,1.)
        gap=torch.linalg.vector_norm(relative+t[:,None]*delta,dim=1)-p.radius[i]-p.radius[j]
        lower=torch.maximum(torch.full_like(gap,2e-9),.001*torch.minimum(p.radius[i],p.radius[j]))
        scale=torch.stack([p.x[i].abs().max(1).values,p.x[j].abs().max(1).values,
            (p.x[i]+dt*velocity[i]).abs().max(1).values,(p.x[j]+dt*velocity[j]).abs().max(1).values,
            p.radius[i],p.radius[j]]).max(0).values
        if bool((gap-lower < -512*torch.finfo(torch.float64).eps*scale).any()):
            raise TrialNeedsSubdivision('GPU_EXACT_SPHERE_SWEEP_BELOW_HANDOFF')

    def _remove_gpu_outlets(self,events):
        if not events:return
        p=self.pool;mask=torch.zeros(p.n,device=p.device,dtype=torch.bool)
        for index,role in sorted(events.items()):
            pid=int(p.particle_id[index]);birth=self.births[pid]
            row=dict(particle_id=pid,source_event_id=birth['source_event_id'],realization_id=self.realization_id,
                birth_time_s=birth['birth_time_s'],physical_time_s=self.time,residence_time_s=self.time-birth['birth_time_s'],
                status='EXITED',terminal_reason='FIRST_CENTER_CROSSING_OF_OFFICIAL_OPEN_OUTLET',outlet=role,exit_time_s=self.time)
            self.terminals.append(row);self.counts['exited_'+role]+=1;mask[index]=True
            if self.recorder is not None:self.recorder.terminal(row)
        self.safety.remove(mask);p.remove_mask(mask);self.cached_field=None;self.cached_wall=None;self.prefetched=None

    def _gpu_safety_terminal(self,reasons):
        p=self.pool;indices=torch.nonzero(reasons).flatten()
        if not len(indices):return False
        names={1:'ROUNDOFF_SCALE_STAGNATION_64_ACCEPTED_SUBSTEPS',2:'EXCESSIVE_SUBDIVISION_PROGRESS_GUARD',
            3:'REJECTED_TRIAL_PROGRESS_GUARD_512_EVALUATIONS_WITHOUT_ACCEPTED_STATE',4:'DETERMINISTIC_PROVIDER_CALL_GUARD',
            5:'EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS'}
        for index in indices.cpu().tolist():
            pid=int(p.particle_id[index]);birth=self.births[pid]
            supported,evidence=self.safety.support(index,float(p.radius[index]))
            row=dict(particle_id=pid,source_event_id=birth['source_event_id'],realization_id=self.realization_id,
                birth_time_s=birth['birth_time_s'],physical_time_s=self.time,residence_time_s=self.time-birth['birth_time_s'],
                status='SUPPORTED_STATIONARY' if supported else 'SOLVER_FAILURE',
                terminal_reason='MODEL_SUPPORTED_STATIONARY_TERMINATION' if supported else 'INTEGRATION_SAFETY_STOP_UNSUPPORTED',
                trigger=names[int(reasons[index])],outlet=None,exit_time_s=None,support_verified=supported,
                support_evidence=evidence,biological_trapping_claimed=False,dissolution_claimed=False)
            self.terminals.append(row);self.counts['supported_stationary_terminated_count' if supported else 'solver_failed_count']+=1
            if self.recorder is not None:self.recorder.terminal(row)
        mask=reasons!=0;self.safety.remove(mask);p.remove_mask(mask)
        self.cached_field=None;self.cached_wall=None;self.prefetched=None
        return True

    def _trial(self,target):
        p=self.pool;n=p.n;dt=float(target-self.time);start=self.time
        self.backend_counts['trial_count']+=1
        if self._gpu_safety_terminal(self.safety.provider()):return
        pre=self.prefetched;self.prefetched=None
        if pre is not None and (pre['time']!=self.time or pre['count']!=n):raise ValueError('STALE_BATCH_PREFETCH')
        sampled=pre['sampled'] if pre is not None else self.cached_field or self.field.sample(p.x[:n],p.tetra_id[:n])
        if bool((sampled[0]<0).any()):raise ValueError('GPU_CENTER_OUTSIDE_AUTHORITATIVE_FEM')
        data=sampled[1];free=torch.cat([data[:,:3],.5*data[:,13:16]],dim=1)
        if self.wall is not None:
            wall=pre['wall'] if pre is not None else self.cached_wall or self.wall.query(p.x[:n],p.radius[:n])
            if bool((wall[1][:,0]<torch.maximum(torch.full_like(p.radius[:n],2e-9),.001*p.radius[:n])-self.wall.roundoff).any()):
                raise ValueError('GPU_ACCEPTED_STATE_BELOW_WALL_HANDOFF')
        else:
            w=torch.zeros((n,14),device=p.device,dtype=torch.float64);w[:,0]=1.;w[:,7]=1.
            wall=(torch.full((n,),-1,device=p.device,dtype=torch.int64),w,None,
                torch.zeros((n,32,3),device=p.device,dtype=torch.float64),torch.zeros(n,device=p.device,dtype=torch.int64))
        if pre is None:
            matrix,rhs,diagonal,assembly=assemble_isolated(p.radius[:n],self.env.mu,wall[1],data[:,4:13].reshape(-1,3,3),free)
            velocity,hydro,audit=solve_isolated(matrix,rhs,diagonal,free,wall[3],wall[4])
        else:
            velocity,hydro,audit,assembly=pre['velocity'],pre['hydro'],pre['audit'],pre['assembly']
            self.backend_counts['batched_gpu_prefix_particle_updates']+=n
        fallback=audit['fallback']|(assembly['affine_solve_info']!=0)
        current_pairs=None;near=set();contacts=set();isolated=torch.ones(n,device=p.device,dtype=torch.bool)
        if self.pairs and n>1:
            current_pairs,geom,_=neighbors(p.x[:n],p.radius[:n],p.realization_id[:n])
            if len(current_pairs):
                ref=torch.minimum(p.radius[current_pairs[:,0]],p.radius[current_pairs[:,1]])
                lower=torch.maximum(torch.full_like(ref,2e-9),.001*ref)
                pair_scale=torch.stack([p.x[current_pairs[:,0]].abs().max(1).values,p.x[current_pairs[:,1]].abs().max(1).values,
                    p.radius[current_pairs[:,0]],p.radius[current_pairs[:,1]]]).max(0).values
                pair_pad=512*torch.finfo(torch.float64).eps*pair_scale
                active=(geom[:,0]<.05*ref)|(geom[:,0]<=lower+pair_pad)
                edges=current_pairs[active]
                if len(edges):
                    fallback[edges.flatten()]=True
                    isolated[edges.flatten()]=False
                    stable=p.particle_id[edges].cpu().tolist();near={tuple(sorted(e)) for e in stable}
                    touch=(geom[active,0]-lower[active]).abs()<=512*torch.finfo(torch.float64).eps*torch.maximum(p.x[edges[:,0]].abs().max(1).values,p.x[edges[:,1]].abs().max(1).values)
                    contacts={tuple(sorted(e)) for e in p.particle_id[edges[touch]].cpu().tolist()}
        ix=torch.nonzero(fallback).flatten();cpu_audit=None
        if len(ix):velocity[ix],cpu_audit=self._cpu_components(ix,free)
        self.backend_counts['gpu_particle_resistance_solves']+=n
        v=velocity[:,:3];dt,outlets=self._outlets(v,dt)
        if dt==0. and outlets:
            self._remove_gpu_outlets(outlets);return
        if dt<=0:raise ValueError('NONPOSITIVE_GPU_TRIAL_INTERVAL')
        pair_root,swept_pairs=self._pair_root(v,dt)
        handoff=None
        if self.wall is not None:
            path=self.wall.path(p.x[:n],v,p.radius[:n],dt)
            roots=torch.where(path[:,1]>=0,path[:,0],torch.full_like(path[:,0],float('inf')))
            root_time=float(roots.min())
            if root_time<dt:
                k=int(roots.argmin());handoff=dict(kind='WALL',particle_id=int(p.particle_id[k]),triangle_id=int(path[k,1]),time_from_start_s=root_time)
        if pair_root is not None and (handoff is None or pair_root['time_from_start_s']<handoff['time_from_start_s']):handoff=pair_root
        if handoff is not None:
            newtime=start+handoff['time_from_start_s']
            if newtime<=start:raise TrialNeedsSubdivision('GPU_HANDOFF_TIME_NOT_REPRESENTABLE')
            dt=newtime-start;outlets={}
        end=p.x[:n]+dt*v
        if self.wall is not None:
            path=self.wall.path(p.x[:n],v,p.radius[:n],dt)
            difficult=torch.nonzero(path[:,2]!=1).flatten()
            if len(difficult):
                shapes=self._shapes(difficult);ends=end[difficult].cpu().numpy();ids=p.particle_id[difficult].cpu().tolist()
                for pid,finish in zip(ids,ends):
                    s=shapes[int(pid)];lower=self.policy.lower_handoff_gap(s.radius_m)['h_lower_m']
                    safe,proof=wall_handoff_certificate(s,s.moved(finish),self.env.wall,lower)
                    if not safe and proof['proof']=='HANDOFF_SWEEP_NOT_CERTIFIED':safe,proof=partitioned_wall_handoff_certificate(s,s.moved(finish),self.env.wall,lower)
                    if not safe:raise TrialNeedsSubdivision(proof['proof'])
                self.backend_counts['cpu_partition_path_particles']+=len(difficult)
        self._pair_certificate(v,dt,swept_pairs)
        newfield=self.field.sample(end,sampled[0]);newwall=None if self.wall is None else self.wall.query(end,p.radius[:n])
        if bool((newfield[0]<0).any()):
            bad=set(torch.nonzero(newfield[0]<0).flatten().cpu().tolist())
            if not bad.issubset(outlets):raise TrialNeedsSubdivision('GPU_ENDPOINT_OUTSIDE_FEM_WITHOUT_OUTLET')
        if not bool(torch.isfinite(velocity).all()):raise ValueError('GPU_NONFINITE_VELOCITY')
        limiter=torch.zeros(n,device=p.device,dtype=torch.bool)
        if handoff is not None:
            responsible=[handoff[k] for k in ['particle_id','particle_j_id'] if k in handoff]
            for pid in responsible:limiter|=p.particle_id[:n]==pid
        self.safety.accepted(p.x[:n],end,v,start+dt,isolated,limiter)
        p.x[:n]=end;p.v[:n]=v;p.omega[:n]=velocity[:,3:];p.q[:n]=advance_quaternion(p.q[:n],velocity[:,3:],dt)
        self.time=start+dt if handoff is not None or outlets else float(target)
        p.age[:n]=self.time-p.birth_time[:n];p.tetra_id[:n]=newfield[0]
        if newwall is not None:p.wall_gap[:n]=newwall[1][:,0]
        self.cached_field=newfield;self.cached_wall=newwall
        self.counts['pair_nearfield_event_count']+=len(near-self._previous_pair_nearfield)
        self.counts['pair_collision_event_count']+=len(contacts-self._previous_pair_contacts)
        self._previous_pair_nearfield=near;self._previous_pair_contacts=contacts
        self.counts['maximum_simultaneous_pair_contacts']=max(self.counts['maximum_simultaneous_pair_contacts'],len(contacts))
        self.backend_counts['accepted_intervals']+=1;self.backend_counts['accepted_particle_updates']+=n
        self.backend_counts['accepted_gpu_authoritative_particle_updates']+=n-len(ix)
        self.backend_counts['accepted_cpu_component_authoritative_particle_updates']+=len(ix)
        self.last_audit=dict(physical_time_s=self.time,evaluated_at_time_s=start,accepted_dt_s=dt,
            handoff_event=handoff,pair_ids=sorted(near),pair_contact_ids=sorted(contacts),
            cpu_component_particle_count=len(ix),cpu_components=cpu_audit,
            maximum_gpu_backward_residual=float(audit['backward_residual'][~fallback].max()) if bool((~fallback).any()) else None,
            maximum_wall_contacts=int(wall[4].max()),gpu_particle_updates=n)
        self.record()
        if self.recorder is not None:self.recorder._row('interaction_audit',self.last_audit)
        self._remove_gpu_outlets(outlets)

    def advance(self,target):
        if not np.isfinite(target) or target<self.time:raise ValueError('Monotone finite global time required')
        def interval(end,depth=0):
            events=0
            while self.time<end:
                if not self.pool.n:self.time=float(end);return
                start=self.time
                try:self._trial(end)
                except TrialNeedsSubdivision as error:
                    mid=start+(end-start)/2
                    if depth>=48 or mid<=start or mid>=end:
                        raise PhysicalTimeRefinementError('GPU_PHYSICAL_TIME_REFINEMENT_LIMIT',(start,end),state=str(error)) from error
                    self.backend_counts['rejected_trials']+=1
                    interval(mid,depth+1);interval(end,depth+1);return
                if self.time==start:
                    if not self.pool.n:continue
                    # Zero-time outlet removal is the only admissible no-time event.
                    events+=1
                    if events>48:raise RuntimeError('GPU_ZERO_TIME_EVENT_LIMIT')
                else:events=0
        interval(float(target))
        if self.pool.n:self._gpu_safety_terminal(self.safety.nominal(self.time,self.pool.x[:self.pool.n]))
        self.accounting()
