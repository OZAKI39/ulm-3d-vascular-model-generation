"""Independent observable judgments and clearly labelled historical supplements."""
import math
from pathlib import Path
import numpy as np
from py_scripts.fluid_physics.analysis import correlation_time_samples
from .equilibration_analysis import metric_rules, sufficient_statistics


OBSERVABLES = ['TEMPERATURE_STATIONARITY','TEMPERATURE_TARGET_MATCH','PRESSURE_STATIONARITY',
               'STRUCTURE_STATIONARITY','SAMPLING_SUFFICIENCY','RESTART_VALIDITY','OVERALL_EQUILIBRIUM']


def not_run():
    return {key:{'status':'NOT_TESTED','execution_status':'NOT_RUN','window_star':None,
                 'description':None,'steady_statistics':None} for key in OBSERVABLES}


def describe(t, x, window):
    lo,hi=window;t=np.asarray(t);x=np.asarray(x)
    mask=(t>lo+1e-12)&(t<=hi+1e-12)&np.isfinite(x)
    t,x=t[mask],x[mask]
    if not len(x):return {'sample_count':0,'window_star':window,'mean':None}
    halves=[x[t<=(lo+hi)/2+1e-12],x[t>(lo+hi)/2+1e-12]]
    quarters=[x[(t>a+1e-12)&(t<=b+1e-12)] for a,b in zip(np.linspace(lo,hi,5)[:-1],np.linspace(lo,hi,5)[1:])]
    return {'sample_count':len(x),'window_star':window,'mean':float(x.mean()),
            'std':float(x.std(ddof=1)) if len(x)>1 else None,'min':float(x.min()),'max':float(x.max()),
            'half_means':[float(a.mean()) if len(a) else None for a in halves],
            'quarter_means':[float(a.mean()) if len(a) else None for a in quarters],
            'scope':'All fixed-window samples, descriptive only; not a steady-state estimate or CI.'}


def hac_trend(t,x,lags=(20,40,80)):
    """Newey-West slope diagnostic, not an equilibrium autocorrelation or mean CI.

    A fixed bandwidth sensitivity analysis on one specified interval. Finite-box
    linear trend + weakly dependent residual assumptions are explicitly retained.
    """
    t,x=np.asarray(t,float),np.asarray(x,float)
    X=np.column_stack([np.ones(len(t)),t-t.mean()]);inv=np.linalg.inv(X.T@X)
    beta=inv@X.T@x;e=x-X@beta;z=X*e[:,None];results=[]
    for L in lags:
        if L>=len(x)//2:continue
        meat=z.T@z
        for j in range(1,L+1):
            cross=z[j:].T@z[:-j];meat+=(1-j/(L+1))*(cross+cross.T)
        se=math.sqrt(max(0.,float((inv@meat@inv)[1,1])))*math.sqrt(len(x)/(len(x)-2))
        results.append({'lag_samples':L,'slope_per_star':float(beta[1]),'slope_SE':se,
                        'approximate_95pct_slope_interval':[float(beta[1]-1.96*se),float(beta[1]+1.96*se)]})
    direction='DECREASING' if results and all(r['approximate_95pct_slope_interval'][1]<0 for r in results) else \
              'INCREASING' if results and all(r['approximate_95pct_slope_interval'][0]>0 for r in results) else 'DIRECTION_UNRESOLVED'
    return {'direction':direction,'sensitivity':results,
            'scope':'SUPPLEMENTAL only: linear trend with correlated residuals, asymptotic HAC bands. Bandwidth sensitivity is not a stationarity test, independent replication, or replacement of frozen acceptance rules.',
            'measured_steady_ACF':None}


def observable_analysis(data,plan,structure,restart=None,complete=True):
    out=not_run();t=np.asarray(data.get('time_star',[]));lo,hi=plan['formal_window_star']
    if not len(t):return out
    c,s=plan['criteria'],plan['parameters'];mask=(t>lo+1e-12)&(t<=hi+1e-12);ts=t[mask]
    full=(complete and len(t)==plan['steps']//s['snapshot_every']+1 and t[0]==0
          and np.allclose(np.diff(t),plan['sample_interval_star'],atol=1e-12,rtol=0)
          and math.isclose(float(t[-1]),plan['duration_star'],abs_tol=1e-12)
          and len(ts)==round((hi-lo)/plan['sample_interval_star']))
    conservation=True
    for key,expected in [('N',s['expected_N']),('total_mass_star',s['expected_N']*s['m_star'])]:
        if key not in data or not np.all(np.asarray(data[key])==expected):conservation=False
    full=full and conservation
    metrics={};block_counts=[]
    for key,(scale,limit,qlimit) in metric_rules(plan).items():
        desc=describe(t,data[key],[lo,hi]);item={'status':'INCONCLUSIVE','description':desc,
                'window_star':[lo,hi],'steady_statistics':None,'point_screen':'NOT_EVALUATED'}
        metrics[key]=item
        if not full or len(ts)<16 or not np.isfinite(np.asarray(data[key])[mask]).all() or any(q is None for q in desc.get('quarter_means',[])):
            item['reason']='Fixed final interval incomplete; no window relocation or qualified CI.';continue
        x=np.asarray(data[key])[mask];halves=desc['half_means'];q=desc['quarter_means']
        drift=abs(halves[1]-halves[0])/scale;qr=np.ptp(q)/scale
        item.update(two_half_change_over_scale=float(drift),quarter_range_over_scale=float(qr),
                    two_half_limit=limit,quarter_limit=qlimit)
        # Preserve the existing coherent-drift criterion separately from HAC.
        chunks=[x[(ts>a+1e-12)&(ts<=b+1e-12)] for a,b in zip(np.linspace(lo,hi,5)[:-1],np.linspace(lo,hi,5)[1:])]
        noise=np.sqrt(np.mean(np.concatenate([a-a.mean() for a in chunks])**2))/scale
        coherent=(np.all(np.diff(q)>0) or np.all(np.diff(q)<0)) and drift>c['coherent_drift_over_within_quarter_rms_min']*noise and (drift>limit or qr>qlimit)
        item['coherent_drift']=bool(coherent)
        if coherent:item.update(status='FAIL',point_screen='FAIL',reason='TRANSIENT_PERSISTS under the existing coherent material trend rule.');continue
        if drift>limit or qr>qlimit:
            item.update(point_screen='FAIL',reason='Frozen point screen exceeded; this alone does not establish persistent drift.');continue
        item['point_screen']='PASS'
        tau=[correlation_time_samples(a) for a in [x,x[ts<=(lo+hi)/2+1e-12],x[ts>(lo+hi)/2+1e-12]]]
        block=max(plan['minimum_block_samples'],math.ceil(c['ACF_multiplier']*max(tau)))
        stat=sufficient_statistics(ts,x,plan['minimum_block_duration_star'],block,c['min_blocks'])
        hs=[sufficient_statistics(ts[m],x[m],plan['minimum_block_duration_star'],block,c['min_half_blocks']) for m in
            [ts<=(lo+hi)/2+1e-12,ts>(lo+hi)/2+1e-12]]
        item.update(conditional_ACF_samples=tau,block_samples=block,block_count=stat['block_count'],
                    tail_samples=stat['discarded_tail_samples'],conditional_block_statistics=stat)
        block_counts.append(stat['block_count'])
        if stat['status']!='SUFFICIENT' or any(h['status']!='SUFFICIENT' for h in hs):
            item['reason']='Insufficient complete blocks; no qualified CI.';continue
        upper=(abs(hs[1]['mean']-hs[0]['mean'])+sum(h['ci95_halfwidth'] for h in hs))/scale
        item.update(drift_upper_bound=upper,reason='Conservative sum of fixed half-mean t interval halfwidths, with complete-block mean/CI samples identical.')
        if upper<=limit:item.update(status='PASS',steady_statistics=stat)
    mapping={'TEMPERATURE_STATIONARITY':'kBT_COM_star','PRESSURE_STATIONARITY':'pressure_star'}
    for name,key in mapping.items():out[name]=metrics[key]
    temp=out['TEMPERATURE_STATIONARITY'];target=s['kBT_star']
    match={'status':'INCONCLUSIVE','window_star':[lo,hi],'description':temp['description'],
           'steady_statistics':None,'reason':'Steady temperature and its uncertainty have not been established.'}
    if temp['status']=='PASS':
        stats=temp['steady_statistics'];err=abs(stats['mean']/target-1);hw=stats['ci95_halfwidth']/target
        match.update(status='PASS' if err+hw<=.02 else 'FAIL' if err-hw>.02 else 'INCONCLUSIVE',
                     steady_statistics=stats,error_including_uncertainty=err+hw,
                     reason='Frozen 2% band including the complete-block mean uncertainty.')
    out['TEMPERATURE_TARGET_MATCH']=match
    structure_metrics={k:v for k,v in metrics.items() if k.startswith('kernel_')}
    structure_pass=structure.get('status')=='PASS' and all(v['status']=='PASS' for v in structure_metrics.values())
    out['STRUCTURE_STATIONARITY']={'status':'PASS' if structure_pass else 'INCONCLUSIVE','window_star':[lo,hi],
        'description':structure,'kernel_metrics':structure_metrics,'steady_statistics':None,
        'reason':'Fixed distributions and snapshot quantiles, not frozen particle positions. Three late frames do not provide independent ensemble replication.'}
    sampling_pass=all(v['status']=='PASS' for v in metrics.values())
    out['SAMPLING_SUFFICIENCY']={'status':'PASS' if sampling_pass else 'INCONCLUSIVE','window_star':[lo,hi],
        'description':{'formal_samples':len(ts),'max_complete_blocks_at_physical_minimum':len(ts)//plan['minimum_block_samples'],
                       'validated_independent_information':None if not sampling_pass else min(block_counts),
                       'physical_minimum_block_samples':plan['minimum_block_samples'],
                       'min_blocks':c['min_blocks'],'min_half_blocks':c['min_half_blocks']},
        'reason':'Block capacity is an upper bound, not established independence. Transient ACF is not a steady-state correlation time.'}
    out['RESTART_VALIDITY']={'status':(restart or {}).get('status','NOT_TESTED'),'description':restart or {'execution_status':'NOT_RUN'}}
    transient=any(v.get('coherent_drift') for v in metrics.values())
    branch='TRANSIENT_PERSISTS' if transient else 'STATIONARY_THERMAL_BIAS' if temp['status']=='PASS' and match['status']=='FAIL' else \
           'EQUILIBRIUM_SCREEN_PASS' if sampling_pass and structure_pass and match['status']=='PASS' else 'SAMPLING_INCONCLUSIVE'
    out['OVERALL_EQUILIBRIUM']={'status':'PASS' if branch=='EQUILIBRIUM_SCREEN_PASS' else 'FAIL' if branch in ['TRANSIENT_PERSISTS','STATIONARY_THERMAL_BIAS'] else 'INCONCLUSIVE',
        'branch':branch,'window_star':[lo,hi],'selection':None,'conservation_check':'PASS' if conservation else 'FAIL_OR_MISSING',
        'equilibrated_state_ready':branch=='EQUILIBRIUM_SCREEN_PASS' and out['RESTART_VALIDITY']['status']=='PASS',
        'reason':'Per-observable progress retained; no claim about viscosity, EOS validation, walls or biological components.'}
    return out


def deterministic_virial(directory, frames, spec):
    """CPU fC/fD reconstruction on saved pre-force frames; total stress is untouched."""
    results=[];root=Path(directory);box=np.asarray(spec['domain_star']);rc=spec['rc_star'];V=float(np.prod(box))
    for frame in frames:
        prefix=root/f"snapshot_{frame['index']:02d}"
        x=np.load(str(prefix)+'_positions.npy').astype(float);v=np.load(str(prefix)+'_pre_velocities.npy').astype(float)
        den=np.load(str(prefix)+'_kernel_number.npy').reshape(-1).astype(float)
        stress=np.load(str(prefix)+'_stresses.npy').astype(float)
        # Same bounded NumPy all-pairs traversal as the existing native snapshot
        # auditor; do not add/install SciPy into the frozen environment.
        chunks=[]
        for begin in range(0,len(x),128):
            delta=x[begin:begin+128,None,:]-x[None,:,:];delta-=box*np.rint(delta/box)
            dist2=np.sum(delta*delta,axis=2)
            ii,jj=np.nonzero((dist2<=rc*rc)&(dist2>=1e-6))
            ii+=begin;keep=ii<jj;chunks.append(np.column_stack([ii[keep],jj[keep]]))
        pairs=np.concatenate(chunks)
        i,j=pairs.T;dr=x[i]-x[j];dr-=box*np.rint(dr/box);r2=np.sum(dr*dr,axis=1);valid=r2>=1e-6
        i,j,dr,r2=i[valid],j[valid],dr[valid],r2[valid];r=np.sqrt(r2);e=dr/r[:,None]
        derivative=20*(21/(2*math.pi))/rc**4*(r/rc)*(r/rc-1)**3
        inv_i,inv_j=1/den[i]**2,1/den[j]**2;cand=spec['candidate']
        pi=cand['sound_speed_star']**2*(spec['m_star']*den[i]-cand['rho_0_star'])
        pj=cand['sound_speed_star']**2*(spec['m_star']*den[j]-cand['rho_0_star'])
        fC=-(inv_i*pi+inv_j*pj)*derivative
        fD=5*cand['viscosity_mu_star']*np.minimum(0,(inv_i+inv_j)*derivative/r)*np.sum(e*(v[i]-v[j]),axis=1)
        pc=float(np.sum(r*fC)/(3*V));pd=float(np.sum(r*fD)/(3*V));total=float(stress[:,[0,3,5]].sum()/(3*V))
        results.append({'step':frame['step'],'conservative_virial_star':pc,'dissipative_virial_star':pd,
                        'total_native_virial_star':total,'stochastic_plus_reconstruction_residual_star':total-pc-pd})
    return {'frames':results,'scope':'Float64 CPU reconstruction of frozen native fC and fD formulas on identical saved pre-integration frames. Residual includes stochastic virial plus reconstruction roundoff, not an independently regenerated RNG. Six frames cannot quantify full-window variance components.'}


def historical_supplement(data, plan, structure, directory, frames, elapsed):
    t=data['time_star'];lo,hi=plan['formal_window_star'];mask=(t>lo+1e-12)&(t<=hi+1e-12)
    details={}
    for key in ['kBT_COM_star','pressure_star','kernel_number_mean','kernel_number_variance']:
        details[key]={'fixed_formal_description':describe(t,data[key],[lo,hi]),
                      'fixed_tenths':[describe(t,data[key],[a,a+.1]) for a in [0.,.1,.2,.3]],
                      'supplemental_HAC_trend':hac_trend(t[mask],data[key][mask])}
    s=plan['parameters'];sigma=float(np.std(data['pressure_star'][mask],ddof=1));thermal=s['task']['n_star']*s['kBT_star']
    nmean=math.ceil((1.96*sigma/thermal)**2);nhalves=math.ceil(8*(1.96*sigma/thermal)**2)
    condition=[]
    for name,n in [('mean_95pct_halfwidth_at_thermal_pressure',nmean),('sum_of_two_half_mean_halfwidths_at_thermal_pressure',nhalves)]:
        condition.append({'criterion':name,'optimistic_iid_samples':n,'duration_star':n*plan['sample_interval_star'],
                          'physical_duration_s':n*plan['sample_interval_star']*s['locked_units']['t0'],
                          'proportional_wall_s_at_old_cadence':n*s['snapshot_every']*elapsed/plan['steps']})
    components=deterministic_virial(directory,frames,s)
    differences=[abs(a['conservative_virial_star']-b['conservative_virial_star']) for a,b in zip(components['frames'],structure['snapshots'])]
    if len(differences)!=len(frames) or max(differences)>.05:
        raise ValueError('CONSERVATIVE_PRESSURE_RECONSTRUCTION_DISAGREES_WITH_OLD_CPU_AUDIT')
    components['crosscheck_against_previous_CPU_density_reconstruction']={'status':'PASS','max_absolute_difference_star':max(differences),
        'tolerance_star':.05,'reason':'Current reconstruction uses saved native density; earlier audit recomputed density on CPU. Difference is retained, never subtracted from the native total pressure.'}
    late=[f for f in components['frames'] if f['step']*s['dt_star']>=lo-1e-12]
    components['late_conservative_virial_change_pa']=(late[-1]['conservative_virial_star']-late[0]['conservative_virial_star'])*s['locked_units']['si_per_star']['pressure']
    components['late_interpretation']='Only three frames: conservative mechanical background rises as structure evolves, while stochastic residual masks it in instantaneous total pressure. This supports a relaxation hypothesis but does not establish a full-window total-pressure drift or its variance decomposition.'
    return {'evidence_kind':'HISTORICAL_CPU_SUPPLEMENT_NOT_NEW_GPU_MEASUREMENT','source_formal_window_star':[lo,hi],
        'old_result_unchanged':'STATIONARITY_OR_SAMPLING_INCONCLUSIVE','details':details,
        'temperature_interpretation':'All four fixed quarter means decline; negative HAC slope persists across lags 20/40/80. This supports continuing cooling under the diagnostic assumptions, although the old 3-RMS coherent-drift rule did not pass. No steady temperature bias is identified.',
        'pressure_interpretation':'Nonmonotone quarter means and HAC slope intervals crossing zero: persistent drift is not established. Instantaneous stochastic/dissipative virial is large. Existing same-frame reductions and phase checks passed; no evidence here requires deleting random stress or substituting input EOS.',
        'structure_interpretation':'Kernel mean/variance fall systematically and nearest q10 rises across the three late frames. Distribution evolution remains plausible; three spatially correlated snapshots cannot prove steady distribution or crystallization.',
        'pressure_sampling':{'measurement':'ONE_INSTANTANEOUS_SAMPLE_EVERY_200_STEPS','is_200_step_average':False,
            'stress_contains':['conservative','dissipative','stochastic'],'pressure_scale_pa':s['locked_units']['si_per_star']['pressure'],
            'thermal_pressure_tolerance_pa':thermal*s['locked_units']['si_per_star']['pressure'],
            'formal_sigma_pa':sigma*s['locked_units']['si_per_star']['pressure'],
            'formal_mean_CI95':None,'optimistic_iid_mean_halfwidth_pa':1.96*sigma/math.sqrt(int(mask.sum()))*s['locked_units']['si_per_star']['pressure'],
            'conditional_costs':condition,'conditional_assumptions':'Independent samples, unchanged noise variance, normal approximation and zero actual drift. These are optimistic planning scenarios, NOT a measured steady ACF or a convergence forecast. Correlation or nonzero drift increases cost.',
            'PROPOSED':{'status':'PROPOSED_NOT_APPLIED','sampling':'Benchmark non-invasive every-step mechanical-pressure accumulation into the unchanged 200-step output bins; keep all conservative/dissipative/random stresses and report interval coverage. Measure cost and correlation before approval.',
                        'criteria':'Keep original threshold and result. Separately discuss whether absolute thermal-pressure drift precision serves the intended physical question. No relaxed acceptance is implemented.'}},
        'pressure_components':components,
        'sampling':{'actual_formal_samples':int(mask.sum()),'physical_minimum_block_samples':plan['minimum_block_samples'],
                    'upper_bound_blocks_before_measured_correlation':int(mask.sum())//plan['minimum_block_samples'],
                    'validated_independent_blocks':None,'measured_steady_ACF':None,
                    'extended_window_capacity_at_minimum_block':1000//plan['minimum_block_samples'],
                    'condition':'The proposed 1000-sample window doubles capacity, not guaranteed independent information or equilibrium.'}}
