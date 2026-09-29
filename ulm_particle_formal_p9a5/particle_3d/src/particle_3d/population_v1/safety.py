"""Explicit inherited progress guards; support requires the original evidence."""
from collections import deque
from types import SimpleNamespace
import numpy as np
from ..formal_dynamics_p9a5 import AuditSink
from ..particle_shapes import Sphere
from ..particle65_motion import assemble_v1
from ..particle9a_motion import augment_planar_system
from ..nearfield_handoff import handoff_constraints
from ..nearfield_regularization import NearFieldRegularizationV1
from ..resistance_solver import solve_resistance
from ..wall_gap import wall_gap


POLICY=dict(name='INHERITED_PROGRESS_GUARDS_AND_SUPPORTED_STATIONARY_REMOVAL_V1',
    unchanged_nominal_intervals=32,roundoff_position_window=64,roundoff_multiplier=16,
    progress_check_provider_interval=128,rejected_provider_checks_without_acceptance=4,
    minimum_64_interval_progress_fraction=.01,maximum_provider_calls=128000,
    supported_window=32,supported_rank=3,stationary_policy='TERMINATE_AND_REMOVE',
    supported_terminal_reason='MODEL_SUPPORTED_STATIONARY_TERMINATION',
    unsupported_stop='SOLVER_FAILURE_NOT_STATIONARY_NOT_OUTLET',
    nominal_clock='GLOBAL_REALIZATION_GRID',support_scope='ISOLATED_PARTICLE_ORIGINAL_SINGLE_MB_CONTRACT',
    physiological_trapping_claimed=False,dissolution_claimed=False)


class PopulationSafetyEvent(RuntimeError):
    def __init__(self,pid,reason):self.pid=int(pid);self.reason=reason;super().__init__(reason)


class CPUSafetyMonitor:
    def __init__(self,env,dt):self.env=env;self.dt=dt;self.rows={};self.last_tick=None
    def insert(self,p,time):
        self.rows[p.particle_id]=dict(calls=0,count=1,last_count=0,no_acceptance=0,
            history=deque([(time,p.position.copy())],maxlen=64),support=deque(maxlen=32),
            previous_nominal=p.position.copy(),unchanged=0)
    def provider(self,pid):
        r=self.rows[pid];r['calls']+=1
        if r['calls']>128000:raise PopulationSafetyEvent(pid,'DETERMINISTIC_PROVIDER_CALL_GUARD')
        if r['calls']%128:return
        r['no_acceptance']=r['no_acceptance']+1 if r['count']==r['last_count'] else 0;r['last_count']=r['count']
        if r['no_acceptance']>=4:raise PopulationSafetyEvent(pid,'REJECTED_TRIAL_PROGRESS_GUARD_512_EVALUATIONS_WITHOUT_ACCEPTED_STATE')
        if len(r['history'])==64:
            positions=np.array([v[1] for v in r['history']]);pad=0. if self.env.wall is None else self.env.wall.roundoff_m
            isolated=len(r['support'])==32 and all(bool(v['solver']) for v in r['support'])
            if isolated and np.linalg.norm(np.ptp(positions,axis=0))<=16*pad:raise PopulationSafetyEvent(pid,'ROUNDOFF_SCALE_STAGNATION_64_ACCEPTED_SUBSTEPS')
            if r['history'][-1][0]-r['history'][0][0]<self.dt*.01:raise PopulationSafetyEvent(pid,'EXCESSIVE_SUBDIVISION_PROGRESS_GUARD')
    def accepted(self,world):
        a=world.projection;components=a.get('components',[dict(particle_ids=tuple(sorted(world.shapes)),**a)])
        byid={i:c for c in components for i in c['particle_ids']}
        for pid,shape in world.shapes.items():
            r=self.rows[pid];r['count']+=1;r['history'].append((world.time_s,shape.center_m.copy()))
            component=byid[pid]
            r['support'].append(dict(solver=component if len(component['particle_ids'])==1 else {},velocity=world.velocities[pid].tolist()))
    def nominal(self,time,shapes):
        ratio=time/self.dt;tick=round(ratio)
        if abs(ratio-tick)>8*abs(np.spacing(ratio)) or tick==self.last_tick:return []
        self.last_tick=tick;events=[]
        for pid,s in shapes.items():
            r=self.rows[pid];r['unchanged']=r['unchanged']+1 if np.array_equal(s.center_m,r['previous_nominal']) else 0
            r['previous_nominal']=s.center_m.copy()
            isolated=len(r['support'])==32 and all(bool(v['solver']) for v in r['support'])
            if isolated and r['unchanged']>=32:events.append((pid,'EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS'))
        return events
    def support(self,pid):
        rows=list(self.rows[pid]['support']);ok=AuditSink.supported(SimpleNamespace(last=rows))
        return bool(ok),rows


def verify_original_support(env,positions,velocities,radius):
    """Rare terminal audit: original isolated CPU solve at all 32 actual old states."""
    if len(positions)!=32:return False,[]
    policy=NearFieldRegularizationV1();rows=[]
    for x,v in zip(positions,velocities):
        shapes={1:Sphere(x,radius)};sample=env.field.sample(x)
        free={1:np.r_[sample.velocity_m_s,.5*sample.vorticity_s_inv]}
        system=assemble_v1(shapes,free,env.mu,env.wall,policy=policy)
        system=augment_planar_system(system,shapes,env.mu,env.wall,
            lambda p:env.field.sample(p).velocity_gradient_s_inv,wall_gap)
        constraints,_=handoff_constraints(shapes,env.wall,policy)
        solved=solve_resistance(system,particles=shapes,constraints=constraints)
        rows.append(dict(solver=solved.record,velocity=v.tolist(),position_m=x.tolist(),
            independently_recomputed_velocity_m_s=solved.velocity[:3].tolist()))
    ok=AuditSink.supported(SimpleNamespace(last=rows))
    ok=ok and all(np.max(np.abs(r['independently_recomputed_velocity_m_s']))<=r['solver'].get('contact_kkt',{}).get('velocity_budget_m_s',0.) for r in rows)
    return bool(ok),rows
