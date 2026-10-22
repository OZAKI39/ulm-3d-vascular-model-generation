"""Protocol regression only: small synthetic trajectories and source-only ensembles."""
import os,json,copy
from pathlib import Path
from decimal import Decimal
from types import SimpleNamespace
from collections import Counter
import numpy as np
import pytest
from particle_3d.population_v1.fixed_window import *
from particle_3d.population_v1.source import PerfusionSource,rng,ARRIVAL
from particle_3d.population_v1.provenance import canonical,write_json
from particle_3d.population_v1.engine import PopulationEngine
from particle_3d.continuous_infusion import poisson_interval
from particle_3d.field import FrozenFEMField
from particle_3d.injection_admission import FiniteSizeAdmission
from test_phase4 import source_env
from test_phase2 import particle,event


def save(name,row):write_json(Path(os.environ.get('FIXED_VALIDATION_DIR','/tmp/fixed_window_validation'))/(name+'.json'),row)

def contract(source,dt=.37,tin=1.240815660):
    return dict(schema=SCHEMA,flow_sha256=source.env.flow_sha256,T_infusion_s=tin,T_observation_s=10.,
        T_end_s=float(Decimal(str(tin))+Decimal('10')),concentration_m3=source.identity['concentration_m3'],
        Q_in_m3_s=source.env.sampler.Q_m3_s,lambda_source_s_inv=source.rate,
        source_distribution_sha256=source.env.distribution.contract_sha256,master_seed=source.seed,
        realization_id=source.realization_id,nominal_dt_s=dt,stopping_policy=STOPPING,
        outcome_dependent_stop=False,right_censor_policy='ACTIVE_AT_T_END')

class EmptyPool:
    n=0

class DiagnosticEngine:
    """Scheduling-only fixture; synthetic immediate terminal outcomes, no dynamics."""
    def __init__(self,mode,dt=.37):
        self.time=0.;self.dt=dt;self.realization_id=0;self.pool=EmptyPool();self.counts=Counter();self.recorder=None
        self.terminals=[];self.source_ledger=[];self.mode=mode;self.calls=[]
    def advance(self,t):self.time=t;self.calls.append(t)
    def accounting(self):return dict(self.counts,physical_time_s=self.time,active_count=0,conservation_pass=True)

class DiagnosticSource:
    def __init__(self,env,mode,rate=1200,seed=17):
        self.raw=PerfusionSource(env,concentration_m3=rate/env.sampler.Q_m3_s,seed=seed);self.mode=mode
    def __getattr__(self,k):return getattr(self.raw,k)
    def stop(self,t):self.raw.stop(t)
    def admit(self,e):
        row=dict(self.peek());pid=row['source_event_id'];e.counts['source_event_count']+=1;e.counts['accepted_birth_count']+=1
        outlet=f'OUTLET_{(pid-1)%3+1:02d}' if self.mode=='all_fast' else 'OUTLET_01'
        e.counts['exited_'+outlet]+=1;e.terminals.append(dict(particle_id=pid,status='EXITED',outlet=outlet,physical_time_s=e.time))
        e.source_ledger.append(row);self.raw.previous_time=row['birth_time_s'];self.raw.next_sid+=1;self.raw.pending=None
        return row


def test_A_B_early_coverage_and_downstream_independent_proposals():
    env=source_env();results=[]
    for mode in ['all_fast','one_outlet_only']:
        s=DiagnosticSource(env,mode);e=DiagnosticEngine(mode);c=contract(s)
        result=run_fixed_window(e,s,c);assert e.time==c['T_end_s'];assert s.stopped_at==c['T_infusion_s']
        results.append((result,e))
    a,e=results[0];b,_=results[1];coverage_time=e.terminals[59]['physical_time_s']
    assert coverage_time<c['T_infusion_s'];assert len(a['source_proposals'])>60
    assert canonical(a['source_proposals'])==canonical(b['source_proposals'])
    assert any(r['birth_time_s']>coverage_time for r in a['source_proposals'])
    save('A_B',dict(PASS=True,scope='SCHEDULER_FIXTURE_NO_PHYSICS',old_20_each_outlet_coverage_time_s=coverage_time,
        T_infusion_s=c['T_infusion_s'],T_end_s=e.time,proposal_count=len(a['source_proposals']),source_proposals_bit_identical=True,
        last_birth_time_s=a['source_proposals'][-1]['birth_time_s']))


def slow_environment():
    env=source_env();old=env.field
    env.field=FrozenFEMField(old.geometry.points,old.geometry.tetra,np.tile([1e-6,0.,0.],(4,1)),np.zeros(4))
    env.checker=FiniteSizeAdmission(None,None,wall=env.wall,field=env.field)
    return env


def test_C_D_E_H_actual_cpu_censor_exact_time_and_accounting():
    env=slow_environment();s=PerfusionSource(env,concentration_m3=2/env.sampler.Q_m3_s,seed=4);c=contract(s)
    e=PopulationEngine(env,dt=c['nominal_dt_s'],pairs=True)
    result=run_fixed_window(e,s,c);assert result['actual_source_count']>0
    assert all(t['status']==CENSORED for t in e.terminals)
    assert all(t['terminal_time_s']==c['T_end_s'] for t in e.terminals)
    assert all(t['age_s']==c['T_end_s']-t['birth_time_s'] for t in e.terminals)
    assert all(t['last_tetra_id']==0 for t in e.terminals)
    assert max(r['physical_time_s'] for r in e.snapshots)==c['T_end_s']
    assert any(r['physical_time_s']==c['T_infusion_s'] for r in e.snapshots)
    assert all(r['birth_time_s']<c['T_infusion_s'] for r in result['source_proposals'])
    save('C_D_E_H',dict(PASS=True,scope='ACTUAL_CPU_INTEGRATOR_SYNTHETIC_CONSTANT_FLOW',result=result,terminals=e.terminals))


def test_E_birth_exactly_at_stop_excluded_and_empty_reaches_end():
    env=source_env();s=DiagnosticSource(env,'all_fast',rate=1);first=s.peek()['birth_time_s'];e=DiagnosticEngine('all_fast')
    c=contract(s,tin=first);result=run_fixed_window(e,s,c)
    assert result['actual_source_count']==0 and e.time==c['T_end_s'] and s.stopped_at==first
    save('E_boundary',dict(PASS=True,event_exactly_T_infusion_excluded=True,empty_population_reaches_end=True,T_infusion_s=first,T_end_s=e.time))


def test_contract_guards():
    env=source_env();s=PerfusionSource(env,concentration_m3=1,seed=2);c=contract(s)
    for key,value in [('outcome_dependent_stop',True),('T_observation_s',12.),('T_end_s',14.),('nominal_dt_s',0.)]:
        bad=dict(c);bad[key]=value
        with pytest.raises(ValueError):validate_contract(bad,s)


def test_F_actual_cuda_cpu_and_output_batch_reproducibility(tmp_path):
    torch=pytest.importorskip('torch')
    if not torch.cuda.is_available():pytest.skip('Actual CUDA required')
    from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
    env=slow_environment();ledgers=[];terminals=[];hashes=[];backend=[]
    for cls,chunk in [(PopulationEngine,1),(GPUPopulationEngine,1),(GPUPopulationEngine,1000)]:
        s=PerfusionSource(env,concentration_m3=2/env.sampler.Q_m3_s,seed=4);c=contract(s)
        e=cls(env,dt=c['nominal_dt_s'],pairs=True);folder=tmp_path/f'run{len(ledgers)}'
        writer=FixedWindowOutput(folder,dict(experiment_contract=c,scope='SYNTHETIC_SMALL_VALIDATION'),chunk_rows=chunk);e.recorder=writer
        result=run_fixed_window(e,s,c);m=writer.close({k:v for k,v in result.items() if k not in ['source_proposals','population_timeline']})
        ledgers.append(canonical(result['source_proposals']));terminals.append(e.terminals);hashes.append(m['ledger_scientific_hashes']);backend.append(dict(getattr(e,'backend_counts',{})))
    assert ledgers[0]==ledgers[1]==ledgers[2]
    assert hashes[0]['source']==hashes[1]['source']==hashes[2]['source']
    error=max(np.max(abs(np.array(a['last_position_m'])-b['last_position_m'])) for a,b in zip(terminals[0],terminals[1]))
    assert error<1e-12
    assert all(r['terminal_time_s']==c['T_end_s'] for group in terminals for r in group)
    save('F',dict(PASS=True,source_proposals_bit_identical=True,admission_ledgers_bit_identical_for_this_fixture=True,
        chunks=[1,1,1000],position_error_m=error,backend_counts=backend,hashes=hashes,terminals=terminals))


def test_G_source_only_poisson_fixed_window():
    # Exact same event-addressed arrival RNG and interval function; no particle trajectories.
    rate=1200.;tin=1.240815660;expected=rate*tin;n=1024;counts=[]
    for rid in range(n):
        t=0.;sid=1
        while True:
            t+=poisson_interval(float(rng(2026092909,rid,sid,ARRIVAL).random()),rate)
            if t>=tin:break
            sid+=1
        counts.append(sid-1)
    values=np.array(counts);mean=values.mean();var=values.var(ddof=1)
    mean_se=np.sqrt(expected/n);variance_se=np.sqrt((expected+2*expected**2)/(n-1))
    assert abs(mean-expected)<6*mean_se and abs(var-expected)<6*variance_se
    folder=Path(os.environ.get('FIXED_VALIDATION_DIR','/tmp/fixed_window_validation'));folder.mkdir(parents=True,exist_ok=True)
    np.savetxt(folder/'G_source_only_counts.csv',np.c_[np.arange(n),values],delimiter=',',header='realization_id,source_count',comments='',fmt='%d')
    save('G',dict(PASS=True,realizations=n,no_dynamics_executed=True,lambda_s_inv=rate,T_infusion_s=tin,expected=expected,
        mean=mean,sample_variance=var,mean_se=mean_se,variance_se=variance_se,mean_z=(mean-expected)/mean_se,variance_z=(var-expected)/variance_se,
        criterion='Predeclared 6 standard errors; arrival-only events use same RNG addresses and poisson_interval as actual source; not 1024 production runs'))


def test_H_terminal_categories_and_legacy_cli_rejected():
    rows=[dict(status='EXITED',outlet=f'OUTLET_{k:02d}') for k in [1,2,3]]
    rows += [dict(status=s) for s in ['SUPPORTED_STATIONARY','SOLVER_FAILURE',CENSORED]]
    categories=[terminal_category(r) for r in rows]
    assert len(set(categories))==6
    import subprocess,sys
    result=subprocess.run([sys.executable,'-m','particle_3d.population_v1.run','--environment','unused.json','--output','unused-output','--contract','unused-contract.json','--git-commit','test','--minimum-exits','20'],capture_output=True,text=True)
    assert result.returncode!=0 and 'unrecognized arguments: --minimum-exits 20' in result.stderr
    save('H_categories_CLI_guard',dict(PASS=True,mutually_distinct_categories=categories,legacy_cli_exit_code=result.returncode))


def test_F_event_mark_batch_sizes():
    env=source_env();source=PerfusionSource(env,concentration_m3=1200/env.sampler.Q_m3_s,seed=2026092909)
    outputs=[]
    for batch in [1,17,256]:
        # Reverse dispatch within each batch, then restore immutable event-ID order.
        marks={}
        for start in range(1,1601,batch):
            for sid in reversed(range(start,min(1601,start+batch))):marks[sid]=source.marks(sid)
        t=0.;rows=[]
        for sid in range(1,1601):
            row=marks[sid];t+=row['arrival_delta_t_s']
            if t>=1.240815660:break
            row['birth_time_s']=t;rows.append(row)
        assert t>=1.240815660
        outputs.append(canonical(rows))
    assert outputs[0]==outputs[1]==outputs[2]
    save('F_mark_batches',dict(PASS=True,batch_sizes=[1,17,256],dispatch_order='REVERSED_WITHIN_BATCH_THEN_EVENT_ID_SORTED',source_proposal_count=len(rows),bit_identical=True))
