from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pyvista as pv
from scipy.stats import kstest, chisquare
from particle_3d.inlet_flux import InletFluxSampler
from particle_3d.injection_method_c import TruncatedSonoVue
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.validation_boundary import ValidationBoundaryClassifier
from particle_3d.population_v1.source import PerfusionSource,rng,ARRIVAL,perfuse_until
from particle_3d.population_v1.engine import PopulationEngine
from particle_3d.population_v1.provenance import canonical
from particle_3d.continuous_infusion import poisson_interval
from test_phase2 import particle,event,save,shear_environment


def source_env():
    env=shear_environment()
    tri=np.array([[[-1e-4,-2e-5,1e-5],[-1e-4,0,1e-5],[-1e-4,-2e-5,3e-5]],
                  [[-1e-4,0,1e-5],[-1e-4,2e-5,1e-5],[-1e-4,0,3e-5]]])
    env.sampler=InletFluxSampler(tri,np.array([[.01,.01,.03],[.01,.01,.03]]))
    from human_vessel_flow.paths import ROOT
    found=ROOT/'inputs/distributions/sonovue'
    env.distribution=TruncatedSonoVue(found);env.flow_sha256='SYNTHETIC_TEST_EXPLICITLY_NOT_PRODUCTION'
    env.checker=FiniteSizeAdmission(None,None,wall=env.wall,field=env.field)
    return env


def test_source_statistics_and_parallel_reproducibility():
    env=source_env();s=PerfusionSource(env,concentration_m3=4e9,seed=2026092901)
    n=30000;rows=[s.marks(i) for i in range(1,n+1)]
    times=np.array([r['arrival_delta_t_s'] for r in rows]);d=np.array([r['diameter_m'] for r in rows])
    poisson=kstest(times,'expon',args=(0,1/s.rate));diameter=kstest(d,env.distribution.cdf)
    counts=np.bincount([r['inlet_triangle_id'] for r in rows],minlength=2)
    theoretical=np.asarray(env.sampler.weights)/env.sampler.Q_m3_s
    flux=chisquare(counts,n*theoretical)
    x=np.array([r['position_m'] for r in rows]);mean_error=x[0]+(x-x[0]).mean(0)-env.sampler.expectation()
    # Predeclared broad 6-sigma statistical acceptance, not exact empirical matching.
    assert abs(1/times.mean()/s.rate-1)<6/np.sqrt(n)
    assert poisson.pvalue>1e-5 and diameter.pvalue>1e-5 and flux.pvalue>1e-5
    assert np.all(abs(mean_error)<6*x.std(0)/np.sqrt(n)+1e-17)
    with ThreadPoolExecutor(4) as executor:parallel=list(executor.map(s.marks,range(1,201)))
    assert canonical(parallel)==canonical(rows[:200])
    assert canonical(PerfusionSource(env,concentration_m3=4e9,seed=2026092901,realization_id=1).marks(1))!=canonical(rows[0])
    save('phase4_source_statistics',dict(PASS=True,n=n,lambda_s_inv=s.rate,empirical_rate_s_inv=1/times.mean(),
        poisson_ks_statistic=poisson.statistic,poisson_pvalue=poisson.pvalue,diameter_ks_statistic=diameter.statistic,
        diameter_pvalue=diameter.pvalue,flux_counts=counts,flux_expected=n*theoretical,flux_pvalue=flux.pvalue,
        flux_centroid_error_m=mean_error,workers_1_vs_4_identical=True,diameter_contract=env.distribution.contract))


def test_single_shot_occupied_inlet_rejection():
    env=source_env();s=PerfusionSource(env,concentration_m3=4e9,seed=99)
    row=s.peek();engine=PopulationEngine(env,pairs=True);engine.advance(row['birth_time_s'])
    p=particle(env,1,row['position_m']);engine.insert_accepted(p,event(p,row['birth_time_s']))
    # Account explicitly for the preexisting fixture birth; it is not this source's draw.
    engine.counts['source_event_count']=1
    result=s.admit(engine)
    assert result['particle_id'] is None and result['rejection_reason'].startswith('PAIR_')
    assert all(result[k]==1 for k in ['diameter_draw_count','position_draw_count','admission_check_count'])
    assert engine.pool.n==1 and s.next_sid==2 and len(engine.source_ledger)==1
    assert engine.accounting()['conservation_pass']
    save('phase4_occupied_inlet',dict(PASS=True,event=result,accounting=engine.accounting()))


def test_dynamic_birth_split_and_accounting():
    env=source_env();s=PerfusionSource(env,concentration_m3=2e14,seed=4)
    engine=PopulationEngine(env,pairs=True)
    result=perfuse_until(engine,s,.003)
    assert result['source_event_count']>0
    for birth in engine.births.values():
        first=next(r for r in engine.snapshots if birth['particle_id'] in r['particle_id'])
        assert first['physical_time_s']==birth['birth_time_s']
        k=list(first['particle_id']).index(birth['particle_id']);assert first['age'][k]==0.
    save('phase4_event_split',dict(PASS=True,accounting=result,source_ledger=engine.source_ledger,
         concentration_role='SYNTHETIC_FIXTURE_ONLY_NOT_BRAVA_PRODUCTION'))


def test_multiple_outlets_ties_and_no_id_reuse():
    env=shear_environment()
    faces=np.array([[3,0,1,2],[3,0,2,3]]).ravel();caps={}
    for k,yrange in enumerate([[-2e-4,-.5e-4],[.5e-4,2e-4]],1):
        pts=np.array([[2e-6,yrange[0],0],[2e-6,yrange[1],0],[2e-6,yrange[1],1e-4],[2e-6,yrange[0],1e-4]])
        caps[f'OUTLET_{k:02d}']=pv.PolyData(pts,faces)
    env.classifier=ValidationBoundaryClassifier(caps);engine=PopulationEngine(env,pairs=True)
    for pid,center in [(2,[0,-1e-4,1e-5]),(7,[0,1e-4,1e-5]),(11,[-1e-6,1e-4,1e-5])]:
        p=particle(env,pid,center,.375e-6);engine.insert_accepted(p,event(p))
    engine.step_to(.001)
    assert engine.pool.n==0
    times={r['particle_id']:r['exit_time_s'] for r in engine.terminals}
    assert abs(times[2]-.0002)<1e-14 and times[2]==times[7] and abs(times[11]-.0003)<1e-14
    try:engine.insert_accepted(particle(env,2,[0,0,1e-5]),event(particle(env,2,[0,0,1e-5]),engine.time))
    except ValueError as error:assert str(error)=='STABLE_ID_NEVER_REUSED'
    else:raise AssertionError('Stable ID was reused')
    save('phase4_outlets',dict(PASS=True,terminals=engine.terminals,accounting=engine.accounting()))
