"""Event-addressed Poisson source and single-shot active-population admission."""
from dataclasses import replace
import numpy as np
from ..continuous_infusion import poisson_interval
from ..particle6_stepper import bind_query_dependency
from ..injection_admission import FiniteSizeAdmission
from .pairs import sphere_gap
from .provenance import content_sha

ARRIVAL,DIAMETER,POSITION,ORIENTATION=940,941,942,943


def rng(seed,realization,event,role):
    for value in [seed,realization,event,role]:
        if isinstance(value,bool) or not isinstance(value,(int,np.integer)) or value<0:
            raise ValueError('Nonnegative integer RNG address required')
    if event<1 or role not in (ARRIVAL,DIAMETER,POSITION,ORIENTATION):raise ValueError('Invalid event/role')
    return np.random.default_rng([int(seed),int(realization),int(event),int(role)])


class PerfusionSource:
    def __init__(self,env,*,concentration_m3,seed,realization_id=0):
        self.env=env;self.seed=int(seed);self.realization_id=int(realization_id)
        self.rate=float(concentration_m3)*float(env.sampler.Q_m3_s)
        if not np.isfinite(self.rate) or self.rate<0:raise ValueError('Finite nonnegative source rate required')
        rng(seed,realization_id,1,ARRIVAL)
        self.next_sid=1;self.previous_time=0.;self.pending=None;self.stopped_at=None
        self.identity=dict(model='GLOBAL_TIME_SINGLE_SHOT_CONTINUOUS_PERFUSION_V1',
            flow_sha256=env.flow_sha256,distribution_contract=env.distribution.contract,
            distribution_contract_sha256=env.distribution.contract_sha256,
            concentration_m3=float(concentration_m3),Q_in_m3_s=float(env.sampler.Q_m3_s),
            lambda_source_s_inv=self.rate,master_seed=seed,realization_id=realization_id,
            rng='PCG64_SEED_SEQUENCE_MASTER_REALIZATION_EVENT_ROLE',
            arrival_role=ARRIVAL,diameter_role=DIAMETER,position_role=POSITION,orientation_role=ORIENTATION,
            admission='ONE_SIZE_ONE_POSITION_ONE_CHECK_ACTIVE_WALL_AND_PAIRS_NO_RETRY')
        self.identity_sha256=content_sha(self.identity)
        # Private binding changes only sphere geometry evaluation, never the
        # original checker or any global module. No retry method is called.
        self.check=bind_query_dependency(FiniteSizeAdmission.check,dict(pair_gap=sphere_gap))

    def marks(self,sid):
        get=lambda role:rng(self.seed,self.realization_id,sid,role)
        arrival_u=float(get(ARRIVAL).random());dt=poisson_interval(arrival_u,self.rate)
        d,draw=self.env.distribution.sample(get(DIAMETER))
        x,tri=self.env.sampler.sample(get(POSITION))
        q=get(ORIENTATION).standard_normal(4);q/=np.linalg.norm(q)
        return dict(source_event_id=int(sid),realization_id=self.realization_id,particle_id=None,
            arrival_uniform=arrival_u,arrival_delta_t_s=dt,diameter_m=float(d),radius_m=float(d/2),
            diameter_um=float(d*1e6),position_m=x[0].tolist(),inlet_triangle_id=int(tri[0]),
            q=q.tolist(),diameter_source_draw=draw,diameter_draw_count=1,position_draw_count=1,
            source_identity_sha256=self.identity_sha256)

    def peek(self):
        if self.stopped_at is not None or self.rate==0:return None
        if self.pending is None:
            self.pending=self.marks(self.next_sid)
            birth=self.previous_time+self.pending['arrival_delta_t_s']
            if not np.isfinite(birth) or birth<=self.previous_time:
                raise FloatingPointError('POISSON_CLOCK_FLOAT64_RESOLUTION_EXHAUSTED')
            self.pending['birth_time_s']=birth
        return self.pending

    def admit(self,engine):
        row=dict(self.peek())
        if engine.time!=row['birth_time_s']:raise ValueError('ADMISSION_REQUIRES_EXACT_BIRTH_TIME')
        if engine.realization_id!=self.realization_id:raise ValueError('Foreign realization')
        provisional=engine.counts['accepted_birth_count']+1
        event=dict(species='MB',particle_id=provisional,radius_m=row['radius_m'],q=row['q'])
        active=engine.admission_shapes(row['position_m'],row['radius_m']) if hasattr(engine,'admission_shapes') else engine.pool.shapes()
        p,status,detail=self.check(self.env.checker,event,np.array(row['position_m']),active)
        row.update(admission_check_count=1,admission_status=status,admission_detail=detail,
            status='ACCEPTED_BIRTH' if p is not None else 'REJECTED_SOURCE_EVENT',
            rejection_reason=None if p is not None else status)
        engine.counts['source_event_count']+=1
        if p is not None:
            if status!='ACCEPTED':raise ValueError('INCONSISTENT_ADMISSION')
            row['particle_id']=provisional
            row['clearance_m']=detail.get('wall_gap_m',detail.get('wall_gap_lower_bound_m',0.))
            engine.insert_accepted(p,row)
        else:
            key='pair_rejected_inlet_count' if status.startswith('PAIR_') else 'wall_rejected_source_count' if status.startswith('WALL_') else 'other_rejected_source_count'
            engine.counts[key]+=1
        engine.source_ledger.append(row)
        if engine.recorder is not None and hasattr(engine.recorder,'source'):engine.recorder.source(row)
        self.previous_time=row['birth_time_s'];self.next_sid+=1;self.pending=None
        engine.accounting()
        return row

    def stop(self,time):
        self.stopped_at=float(time)
        # A future precomputed mark is NOT a realized source event or rejection.
        self.pending=None


def perfuse_until(engine,source,end_time):
    """Finite-window helper for tests/CLI; production coverage uses the same events."""
    while (row:=source.peek()) is not None and row['birth_time_s']<=end_time:
        engine.step_to(row['birth_time_s']);source.admit(engine)
    engine.step_to(end_time)
    return engine.accounting()
