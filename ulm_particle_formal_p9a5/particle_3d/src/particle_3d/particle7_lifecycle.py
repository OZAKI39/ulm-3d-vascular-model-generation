"""Event-split inlet lifecycle. Plug transport is an explicit validation fixture.

The infrastructure accepts an external certified dynamics mover. It does not
introduce RBC lubrication or replace the P6.5 resistance solver.
"""
from copy import deepcopy
from dataclasses import replace
from collections import deque
import math
import numpy as np
from .injection_population import InjectionScheduler
from .lammps_neighbors import standalone_candidates


class PopulationLifecycle:
    def __init__(self,scheduler,admission,classifier,*,bridge=None,mover=None):
        self.scheduler=scheduler; self.admission=admission; self.classifier=classifier; self.bridge=bridge; self.mover=mover
        self.time_s=0.; self.pending={s:deque() for s in ['MB','RBC']}; self.active={}; self.births={}; self.events=[]; self.exits=[]; self.timeline=[]
        self.neighbor_comparisons=0; self.neighbor_mismatches=0
    def sync(self):
        if self.bridge is not None:
            self.bridge.write(self.active.values()); actual=self.bridge.rebuild()
            expected=standalone_candidates({i:p.shape() for i,p in self.active.items()},self.bridge.policy)
            self.neighbor_comparisons+=1; self.neighbor_mismatches+=len(set(actual)^set(expected))
            if actual!=expected: raise AssertionError('Dynamic LAMMPS neighbor mismatch')
    def retry(self):
        # FIFO independently in each species, as persisted pending_MB/pending_RBC.
        # A blocked RBC does not suppress actual MB admission attempts.
        for species in ('MB','RBC'):
            queue=self.pending[species]
            while queue:
                event=queue[0]
                if event['scheduled_time_s']>self.time_s: break
                p=self.admission.attempt(event,self.active)
                if p is None: break
                queue.popleft(); event['admitted_time_s']=self.time_s
                event['birth_position_m']=p.position.tolist(); event['birth_shape_mode']=p.shape().mode
                self.active[p.particle_id]=p; self.births[p.particle_id]=deepcopy(event)
                if self.bridge is not None: self.bridge.insert(p)
                self.sync()
    def advance(self,target):
        if target<self.time_s: raise ValueError('Physical time decreased')
        if target==self.time_s: return
        while self.time_s<target:
            old=self.time_s; dt=target-old
            # A mover must return a certified endpoint with center first-outlet
            # events along its path. Default is exact constant-velocity validation.
            proposed={i:replace(p,position=p.position+p.velocity*dt) for i,p in self.active.items()} if self.mover is None else self.mover(self.active,old,target)
            if isinstance(proposed,tuple):
                proposed,transport_end=proposed
                if not old<transport_end<=target: raise ValueError('Mover endpoint time outside requested interval')
                dt=transport_end-old
            hits={}
            for i,p in self.active.items():
                hit=self.classifier.first_event(p.position,proposed[i].position)
                if hit is not None:
                    if not hit.role.startswith('OUTLET_'): raise ValueError('Non-outlet boundary in lifecycle outlet classifier')
                    hits[i]=hit
            fraction=min((h.segment_fraction for h in hits.values()),default=1.)
            stop=old+dt*fraction
            self.active={i:(proposed[i] if fraction==1. else replace(p,position=p.position+fraction*(proposed[i].position-p.position))) for i,p in self.active.items()}
            self.time_s=stop
            for i,h in sorted(hits.items()):
                if h.segment_fraction!=fraction: continue
                p=self.active[i]; birth=self.births[i]
                self.exits.append(dict(particle_id=i,particle_type=birth['species'],birth_scheduled_time=birth['scheduled_time_s'],
                    birth_admitted_time=birth['admitted_time_s'],exit_time=stop,lifetime=stop-birth['admitted_time_s'],outlet_role=h.role,
                    geometry_provenance=deepcopy(birth['provenance']),volume_m3=birth['volume_m3'],
                    trajectory_summary=dict(birth_position_m=birth['birth_position_m'],exit_position_m=h.position_m.tolist(),
                    displacement_m=float(np.linalg.norm(h.position_m-np.asarray(birth['birth_position_m']))))))
                del self.active[i]
                if self.bridge is not None: self.bridge.remove(i)
            self.sync(); self.retry()
            if self.time_s>=target: break
            if not hits: raise ValueError('Mover stopped early without a boundary event')
            if stop==old and not any(h.segment_fraction==fraction for h in hits.values()): raise ValueError('No lifecycle progress')
        self.time_s=target
    def step_to(self,target):
        if target<self.time_s: raise ValueError('Time cannot decrease')
        while self.scheduler.next_time()<=target:
            t=self.scheduler.next_time(); self.advance(t)
            event=self.scheduler.pop(); self.events.append(deepcopy(event)); self.pending[event['species']].append(event); self.retry()
        self.advance(target); self.scheduler.time_s=target
        row=self.accounting(); self.timeline.append(row); return row
    def accounting(self):
        result=dict(time_s=self.time_s,mb_expected=self.scheduler.mb_clock.cumulative(self.time_s),rbc_target_m3=self.scheduler.rbc_clock.cumulative(self.time_s))
        for species in ('MB','RBC'):
            key=species.lower(); scheduled=[e for e in self.events if e['species']==species]
            admitted=[e for e in self.births.values() if e['species']==species]
            pending=list(self.pending[species]); active=[self.births[i] for i in self.active if self.births[i]['species']==species]
            exited=[e for e in self.exits if e['particle_type']==species]
            for label,items in [('scheduled',scheduled),('admitted',admitted),('pending',pending),('active',active),('exited',exited)]:
                result[f'{label}_{key}_count']=len(items)
                if species=='RBC': result[f'{label}_rbc_volume_m3']=math.fsum(e['volume_m3'] for e in items)
            if len(scheduled)!=len(admitted)+len(pending) or len(admitted)!=len(active)+len(exited): raise AssertionError('Particle conservation failed')
        return result
    def state(self):
        return dict(time_s=self.time_s,scheduler=self.scheduler.state(),pending={k:list(v) for k,v in self.pending.items()},
            births=self.births,events=self.events,exits=self.exits,timeline=self.timeline,
            active=[p.to_dict() for p in self.active.values()],admission_counts=dict(self.admission.counts),
            admission_candidates=self.admission.candidates,neighbor_comparisons=self.neighbor_comparisons,neighbor_mismatches=self.neighbor_mismatches)
