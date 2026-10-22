"""Fixed input window orchestration only; inherited particle physics is unchanged."""
from collections import Counter
from decimal import Decimal
import gzip,hashlib
import numpy as np
from .output import PopulationOutput
from .pairs import candidates,sphere_gap
from ..particle_shapes import Sphere
from ..nearfield_regularization import NearFieldRegularizationV1

SCHEMA='FIXED_WINDOW_POPULATION_EXPERIMENT_V1'
STOPPING='FIXED_INFUSION_AND_FIXED_POST_INFUSION_OBSERVATION'
CENSORED='RIGHT_CENSORED_OBSERVATION_END'


def validate_contract(c,source=None):
    required=['schema','flow_sha256','T_infusion_s','T_observation_s','T_end_s','concentration_m3',
        'Q_in_m3_s','lambda_source_s_inv','source_distribution_sha256','master_seed','realization_id',
        'nominal_dt_s','stopping_policy','outcome_dependent_stop','right_censor_policy']
    if any(k not in c for k in required):raise ValueError('Incomplete fixed-window contract')
    if c['schema']!=SCHEMA or c['stopping_policy']!=STOPPING or c['outcome_dependent_stop'] is not False or c['right_censor_policy']!='ACTIVE_AT_T_END':
        raise ValueError('Outcome-dependent or unknown protocol is forbidden')
    for k in ['T_infusion_s','T_observation_s','T_end_s','concentration_m3','Q_in_m3_s','lambda_source_s_inv','nominal_dt_s']:
        if isinstance(c[k],bool) or not np.isfinite(c[k]) or c[k]<0:raise ValueError('Invalid '+k)
    if c['T_observation_s']!=10. or c['T_end_s']!=float(Decimal(str(c['T_infusion_s']))+Decimal('10.0')) or c['nominal_dt_s']<=0:
        raise ValueError('Observation must be exactly 10 s; endpoint must equal infusion + 10 s')
    if not np.isclose(c['lambda_source_s_inv'],c['concentration_m3']*c['Q_in_m3_s'],rtol=1e-14,atol=0):raise ValueError('Source rate mismatch')
    for key in ['master_seed','realization_id']:
        if isinstance(c[key],bool) or not isinstance(c[key],int) or c[key]<0:raise ValueError('Invalid '+key)
    if source is not None:
        pairs=[('flow_sha256','flow_sha256'),('source_distribution_sha256','distribution_contract_sha256'),
            ('master_seed','master_seed'),('realization_id','realization_id'),('concentration_m3','concentration_m3')]
        for a,b in pairs:
            if c[a]!=source.identity[b]:raise ValueError('Runtime source contract mismatch: '+a)
        for a,b in [('Q_in_m3_s','Q_in_m3_s'),('lambda_source_s_inv','lambda_source_s_inv')]:
            if not np.isclose(c[a],source.identity[b],rtol=1e-13,atol=0):raise ValueError('Runtime source mismatch: '+a)
    return c


def terminal_category(row):
    status=row['status']
    if status=='EXITED':return 'EXITED_O'+str(int(row['outlet'].split('_')[-1]))
    return {'SUPPORTED_STATIONARY':'MODEL_SUPPORTED_STATIONARY_TERMINATION',
        'SOLVER_FAILURE':'NUMERICAL_FAILURE',CENSORED:CENSORED}[status]


class FixedWindowOutput(PopulationOutput):
    """Separate exogenous proposals from downstream-dependent admission results."""
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.streams['source_proposals']=gzip.open(self.folder/'source_proposals.jsonl.gz','wb',compresslevel=3)
        self.hashes['source_proposals']=hashlib.sha256()
    def proposal(self,row):self._row('source_proposals',row)
    def terminal(self,row):self._row('terminals',dict(row,terminal_category=terminal_category(row)))


def _host(a):
    return a.detach().cpu().numpy() if hasattr(a,'detach') else np.asarray(a)


def censor_remaining(engine,end):
    """Snapshot final accepted state and final-geometry pair status; never extrapolate."""
    if engine.time!=end:raise ValueError('Censor requires exact fixed observation endpoint')
    p=engine.pool;n=p.n
    if not n:return []
    arrays={k:_host(v[:n]).copy() for k,v in p.arrays().items()}
    shapes={int(arrays['particle_id'][k]):Sphere(arrays['x'][k],arrays['radius'][k]) for k in range(n)}
    policy=NearFieldRegularizationV1();byid={pid:[] for pid in shapes}
    for i,j in candidates(shapes):
        gap=sphere_gap(shapes[i],shapes[j],i,j);ref=min(shapes[i].radius_m,shapes[j].radius_m)
        lower=policy.lower_handoff_gap(ref)['h_lower_m'];weight=policy.activation_weight(gap.gap_m/ref)
        state=policy.interaction_state(gap.gap_m,ref,lower,gap.roundoff_budget_m)
        if weight>0 or gap.gap_m<=lower+gap.roundoff_budget_m:
            for a,b in [(i,j),(j,i)]:byid[a].append(dict(partner_particle_id=b,gap_m=gap.gap_m,
                interaction_state=state,nearfield_weight=weight,geometry_time_s=end))
    rows=[]
    for k,pid0 in enumerate(arrays['particle_id']):
        pid=int(pid0);birth=engine.births[pid];row=dict(particle_id=pid,source_event_id=birth['source_event_id'],
            realization_id=engine.realization_id,birth_time_s=birth['birth_time_s'],physical_time_s=end,terminal_time_s=end,
            age_s=end-birth['birth_time_s'],residence_time_s=end-birth['birth_time_s'],
            status=CENSORED,terminal_category=CENSORED,terminal_reason='FIXED_OBSERVATION_WINDOW_ENDED',outlet=None,exit_time_s=None,
            last_position_m=arrays['x'][k].tolist(),last_velocity_m_s=arrays['v'][k].tolist(),
            last_tetra_id=int(arrays['tetra_id'][k]),last_wall_gap_m=float(arrays['wall_gap'][k]),
            active_pair_state=byid[pid],pair_state_definition='FRESH_EXACT_SPHERE_GEOMETRY_AT_T_END_USING_UNCHANGED_V1_THRESHOLDS',
            last_velocity_definition='HELD_ON_PRECEDING_ACCEPTED_INTERVAL')
        rows.append(row);engine.terminals.append(row);engine.counts['censored_count']+=1
        if engine.recorder is not None:engine.recorder.terminal(row)
    if hasattr(p,'remove_mask'):
        import torch
        remove=torch.ones(n,device=p.device,dtype=torch.bool);engine.safety.remove(remove);p.remove_mask(remove)
    else:
        p.remove(shapes)
        for pid in shapes:engine.safety.rows.pop(pid,None)
    return rows


def run_fixed_window(engine,source,contract,progress=None):
    """One exact global clock; outcomes are diagnostics and never clock inputs."""
    validate_contract(contract,source)
    if engine.time!=0 or engine.realization_id!=contract['realization_id'] or engine.dt!=contract['nominal_dt_s']:
        raise ValueError('Fresh realization, matching dt and realization ID required')
    tin=float(contract['T_infusion_s']);end=float(contract['T_end_s']);timeline=[];proposals=[]
    def mark(phase):
        row=dict(engine.accounting(),phase=phase,source_on=engine.time<tin,
            T_infusion_s=tin,T_end_s=end);timeline.append(row)
        if progress:progress(row)
    mark('INFUSION_ON' if tin>0 else 'INFUSION_OFF')
    if tin==0:source.stop(0.)
    while engine.time<end:
        source_on=engine.time<tin;row=source.peek() if source_on else None
        birth=row['birth_time_s'] if row is not None and row['birth_time_s']<tin else float('inf')
        r=engine.time/engine.dt;tick=(int(np.floor(r+8*abs(np.spacing(r))))+1)*engine.dt
        if tick<=engine.time:tick=engine.time+engine.dt
        # Empty populations still reach fixed times; no zero-state dynamics are needed.
        target=min(tick if engine.pool.n else float('inf'),birth,tin if source_on else end,end)
        engine.advance(float(target))
        if engine.time!=target:raise RuntimeError('Integrator did not reach exact requested time')
        if source_on and target==tin:
            source.stop(tin);mark('INFUSION_OFF')
        elif row is not None and target==birth:
            proposal=dict(row);proposals.append(proposal)
            if engine.recorder is not None and hasattr(engine.recorder,'proposal'):engine.recorder.proposal(proposal)
            source.admit(engine)
        if progress:progress(dict(engine.accounting(),phase='INFUSING' if engine.time<tin else 'OBSERVING',source_on=engine.time<tin))
    mark('OBSERVATION_END_BEFORE_CENSOR')
    censor_remaining(engine,end)
    if hasattr(engine,'flush_output'):engine.flush_output()
    mark('EXPERIMENT_END')
    categories=Counter(terminal_category(r) for r in engine.terminals)
    if sum(categories.values())!=engine.counts['accepted_birth_count']:raise AssertionError('Exclusive terminal accounting failed')
    if len({r['particle_id'] for r in engine.terminals})!=len(engine.terminals):raise AssertionError('Duplicate terminal ID')
    if len(proposals)!=engine.counts['source_event_count']:raise AssertionError('Proposal accounting failed')
    if any(r['birth_time_s']>=tin for r in proposals):raise AssertionError('Post-infusion birth')
    return dict(status='COMPLETED_FIXED_WINDOW',**engine.accounting(),terminal_categories=dict(categories),
        previous_protocol='OUTCOME_DEPENDENT_COVERAGE_STOP',current_protocol=SCHEMA,
        T_infusion_s=tin,T_observation_s=10.,T_end_s=end,infusion_stop_time_s=tin,
        source_stop_rule_is_data_dependent=False,expected_source_count=source.rate*tin,
        actual_source_count=len(proposals),
        source_rejection_count=len(proposals)-engine.counts['accepted_birth_count'],
        source_proposals=proposals,population_timeline=timeline,
        no_numerical_failure_reported=engine.counts['solver_failed_count']==0,
        cfd_physical_accuracy_validated_by_protocol=False)
