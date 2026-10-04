"""Run the frozen three-step study and audit the saved states independently.

Run from clot_clearance with PYTHONPATH=.:
  python scripts/verify_streaming_erosion.py --run
  python scripts/verify_streaming_erosion.py --analyze
Existing runs are never overwritten. The 0.01 s directory aliases coarse.
"""
import argparse,copy,csv,hashlib,json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from pd_clot.regularization.erosion_runner import run

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results/streaming_erosion_demo_v1'
VERIFY=ROOT/'verification/streaming_erosion_v1'
CASES=[('dt_erosion_0p02',.02),('coarse',.01),('dt_erosion_0p005',.005)]
METRICS={
    'mean_damage':1e-12,'maximum_damage':1e-12,
    'damage_stored_energy_drop_J':1e-30,'stored_energy_J':1e-30,
    'maximum_displacement_m':1e-15,'maximum_epsilon_plus':1e-12,
    'fraction_of_damage_inside_influence_region':1e-12,
    'damage_weighted_distance_from_source_m':1e-15}


def read(path):return json.loads(path.read_text())


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)


def execute():
    path=ROOT/'configs/streaming_erosion_demo_coarse.json';c=read(path)
    plan=read(VERIFY/'provenance/PREDECLARED_PLAN.json')
    assert hashlib.sha256(path.read_bytes()).hexdigest()==plan['new_erosion_config_sha256']
    for name,_ in CASES:
        if (RESULTS/name).exists():raise FileExistsError(RESULTS/name)
    for name,step in CASES:
        config=copy.deepcopy(c)
        config['name']='streaming_erosion_demo_'+name
        config['streaming_erosion'].update(erosion_exposure_step_s=step,number_of_erosion_steps=round(.1/step))
        run(config,RESULTS/name)
    (RESULTS/'dt_erosion_0p01').symlink_to('coarse',target_is_directory=True)
    write(RESULTS/'CASE_INDEX.json',dict(cases=CASES,dt_erosion_0p01='Relative symlink to coarse; one computation, not an independent replicate'))


def audit(path):
    c=read(path/'CONFIG.json');identity=read(path/'IDENTITY.json');summary=read(path/'SUMMARY.json');h=read(path/'history.json')
    with np.load(path/'states.npz') as data:z={k:data[k] for k in data.files}
    n=c['streaming_erosion']['number_of_erosion_steps'];step=c['streaming_erosion']['erosion_exposure_step_s']
    assert n*step==.1
    np.testing.assert_allclose(z['erosion_exposure_time_s'],np.arange(n+1)*step,rtol=0,atol=1e-15)
    np.testing.assert_allclose(z['mechanical_time_s'],np.arange(n+1)*.002,rtol=0,atol=1e-15)
    for key,value in z.items():assert np.isfinite(value).all(),key
    D=z['bond_D'];assert np.all(np.diff(D,axis=0)>=0) and np.all((D>=0)&(D<=1))
    np.testing.assert_array_equal(z['integrity'],1.-D);np.testing.assert_array_equal(z['bond_active'],D<1.)
    np.testing.assert_array_equal(z['D_before_update'][1:],D[:-1])
    x=z['x'];pairs=z['pairs'];X=z['X'];fixed=z['fixed'];volume=z['volume']
    length=np.linalg.norm(X[pairs[:,1]]-X[pairs[:,0]],axis=1)
    eps=np.maximum(np.linalg.norm(x[:,pairs[:,1]]-x[:,pairs[:,0]],axis=2)/length-1.,0.)
    np.testing.assert_array_equal(z['erosion_driver_epsilon_plus'][1:],eps[1:])
    expected=np.minimum(1.,D[:-1]+step*.01*eps[1:]/1e-4*(D[:-1]<1))
    np.testing.assert_allclose(D[1:],expected,rtol=1e-14,atol=1e-18)
    np.testing.assert_array_equal(x[:,fixed],np.broadcast_to(X[fixed],x[:,fixed].shape))
    np.testing.assert_array_equal(z['v'][:,fixed],0)
    for i,row in enumerate(h):
        degree=np.bincount(pairs.ravel(),minlength=len(X))
        pd=np.bincount(pairs.ravel(),weights=np.repeat(D[i],2),minlength=len(X))/degree
        np.testing.assert_allclose(z['damage'][i],pd,rtol=1e-11,atol=1e-16)
        p=pairs[D[i]<1]
        graph=coo_matrix((np.ones(2*len(p)),(p.ravel(),p[:,::-1].ravel())),shape=(len(X),len(X)))
        _,labels=connected_components(graph,directed=False)
        attached=np.isin(labels,np.unique(labels[fixed]));np.testing.assert_array_equal(attached,z['attached'][i])
        assert row['broken_bond_count']==int((D[i]==1).sum())
        assert row['detached_total_volume_m3']==float(volume[~attached].sum())
        assert row['particle_count']==len(X)
        assert abs(row['relative_numerical_energy_residual'])<=.05 or abs(row['numerical_energy_residual_J'])<=1e-15
        assert row['mechanical_time_s']==z['mechanical_time_s'][i]
        assert row['erosion_exposure_time_s']==z['erosion_exposure_time_s'][i]
        np.testing.assert_allclose(z['bond_stored_energy_drop_J'][i].sum(),row['damage_stored_energy_drop_J'],rtol=1e-12,atol=1e-30)
    for rel,sha in identity['source_sha256'].items():
        assert hashlib.sha256((path/'source_snapshot'/rel).read_bytes()).hexdigest()==sha
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==sha
    assert c['streaming']['oscillatory_amplitude']==0 and c['streaming']['mean_fraction']==1
    for key,value in dict(C_t_s_inv=.01,m_t=1.,epsilon_ref=1e-4,threshold_ratio=0.,D_break=1.).items():assert c['streaming_erosion'][key]==value
    plan=read(VERIFY/'provenance/PREDECLARED_PLAN.json')
    assert hashlib.sha256((ROOT/'configs/streaming_erosion_demo_coarse.json').read_bytes()).hexdigest()==plan['new_erosion_config_sha256']
    return dict(status='PASS',case=path.name,checked_snapshots=len(h),particles=len(X),bonds=len(pairs),
                checks=['finite','analytic_rate_from_saved_endpoint','two_clocks','irreversible','D_break_1','connectivity',
                        'all_particles_retained','fixed_exact','energy_budget','source_hashes','frozen_parameters'],
                minimum_J=summary['minimum_J_or_intrinsic_stretch'],max_energy_residual=max(abs(row['relative_numerical_energy_residual']) for row in h),
                source_hash=hashlib.sha256(json.dumps(identity['source_sha256'],sort_keys=True).encode()).hexdigest())


def analyze():
    audits=[audit(RESULTS/name) for name,_ in CASES]
    write(VERIFY/'comparisons/STATE_AUDIT.json',audits)
    histories={name:read(RESULTS/name/'history.json') for name,_ in CASES}
    for rows in histories.values():
        for row in rows:row['stored_energy_J']=row['elastic_energy_J']+row['stabilization_energy_J']
    reference=histories['dt_erosion_0p005'];comparisons=[]
    for name,step in CASES[:-1]:
        for t in np.arange(1,6)*.02:
            row=min(histories[name],key=lambda h:abs(h['erosion_exposure_time_s']-t))
            ref=min(reference,key=lambda h:abs(h['erosion_exposure_time_s']-t))
            assert abs(row['erosion_exposure_time_s']-t)<1e-14 and abs(ref['erosion_exposure_time_s']-t)<1e-14
            for key,floor in METRICS.items():
                a,b=row[key],ref[key]
                comparisons.append(dict(case=name,step_s=step,erosion_exposure_time_s=float(t),metric=key,value=a,
                    reference=b,absolute_difference=abs(a-b),relative_difference=abs(a-b)/max(abs(b),floor),denominator_floor=floor))
    worst=max(comparisons,key=lambda r:r['relative_difference'])
    final=[r for r in comparisons if abs(r['erosion_exposure_time_s']-.1)<1e-14]
    summaries={name:read(RESULTS/name/'SUMMARY.json') for name,_ in CASES}
    settled=all(s['mechanical_settling_status']=='DRIVER_SETTLING_CHECK_PASSED' for s in summaries.values())
    status='INCONCLUSIVE' if not settled else 'CONSISTENT' if worst['relative_difference']<=.05 else 'NOT YET CONVERGED'
    result=dict(status=status,reason='Any failed settling makes slow-step accuracy inconclusive; no window or parameters retuned' if not settled else '5% descriptive comparison tolerance; this is not mesh convergence',
        common_slow_exposure_times_s=[.02,.04,.06,.08,.1],reference_step_s=.005,
        relative_difference_definition='abs(case-reference)/max(abs(reference), metric_floor); no initial zero state',metric_floors=METRICS,
        maximum_relative_difference=worst,maximum_final_relative_difference=max(final,key=lambda r:r['relative_difference']),
        by_metric={k:dict(maximum_common_time_relative_difference=max(r['relative_difference'] for r in comparisons if r['metric']==k),
                        maximum_final_relative_difference=max(r['relative_difference'] for r in final if r['metric']==k)) for k in METRICS},
        cases=[dict(name=name,step_s=step,mechanical_settling_status=summaries[name]['mechanical_settling_status'],
                    maximum_driver_relative_drift=summaries[name]['maximum_driver_relative_drift'],final=summaries[name]['final']) for name,step in CASES])
    write(VERIFY/'comparisons/TIME_STEP_STUDY.json',result)
    with (VERIFY/'comparisons/time_step_differences.csv').open('x',newline='') as f:
        writer=csv.DictWriter(f,list(comparisons[0]));writer.writeheader();writer.writerows(comparisons)
    for name,_ in CASES:
        path=RESULTS/name
        write(path/'SHA256.json',{str(p.relative_to(path)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(path.rglob('*')) if p.is_file()})
    print(json.dumps({k:result[k] for k in ['status','maximum_relative_difference','maximum_final_relative_difference']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',action='store_true');p.add_argument('--analyze',action='store_true');args=p.parse_args()
    if args.run:execute()
    if args.analyze:analyze()
