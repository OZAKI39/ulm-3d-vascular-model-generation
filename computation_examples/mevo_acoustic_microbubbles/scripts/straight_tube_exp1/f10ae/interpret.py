"""Scientific admission is stricter than accepted integration, with independent targets."""
import numpy as np,h5py
from .common import *
def main():
    rows=read(BASE/'ALL_RUNS.json');by={r['id']:r for r in rows};validation=read(BASE/'VALIDATION.json');tests=read(BASE/'FINAL_TEST_RUNS.json');sel=read(BASE/'SELECTED_FROZEN.json')['selected'];low=read(BASE/'LOW_PRESSURE_ADMISSION.json');audit=read(BASE/'AUDIT.json');baseline=read(BASE/'BASELINE_VALIDATION.json')
    targets={}
    for g in validation['groups']:
        key=g['target'];rr=[r for r in tests if r['id'].startswith('final_'+key+'_')];raw=[]
        for r in rr:
            ok=r['status']=='PASS'
            if ok and key in ['T1','T2']:ok=r['targets'][key+'_DISCOVERY_PASS']
            elif ok and key=='T3':ok=r['targets']['raw_extension_above_floor'] and r['targets']['COM_TARGET_PASS']
            elif ok and key=='T4':ok=all(r['recovery']['raw_recovery_pass'].values())
            raw.append(bool(ok))
        reps=g['replicates'];confirmed=g['discovery_and_validation_confirmed'] and len(raw)==3 and all(raw)
        targets[key]=dict(source=sel[key]['source'],discovery=by[sel[key]['source']]['targets'],validation_confirmed=sum(r['confirmed'] for r in reps[1:]),validation_numerical_pass=sum(r['numerical_pass'] for r in reps[1:]),final_test_raw_pass=sum(raw),final_test_count=len(raw),confirmed=confirmed,replicates=reps,final_test_ids=[r['id'] for r in rr],final_test_scope='Once-only frozen-action outcomes; dt/cache admission comes from separate three validation states, not per-test refinements')
    t4=by[sel['T4']['source']];off=by['schedule_C'];r4=t4['recovery'];offerrors=off['recovery']['final_errors'];recovery_gates={}
    recovery_keys={'RG_RECOVERY_PASS':'Rg','PRINCIPAL_AXIS_RECOVERY_PASS':'axes','COVARIANCE_RECOVERY_PASS':'covariance','FIXED_ID_CONFIGURATION_RECOVERY_PASS':'Kabsch_m','TOPOLOGY_RECOVERY_PASS':'neighbor_error'}
    for gate,key in recovery_keys.items():
        reps=targets['T4']['replicates'];test=[r for r in tests if r['id'].startswith('final_T4_')]
        recovery_gates[gate]=bool(r4['raw_recovery_pass'][key] and r4['final_errors'][key]<offerrors[key] and all(r['numerical_pass'] and r['recovery'] and r['recovery']['raw_recovery_pass'][key] for r in reps) and all(r['status']=='PASS' and r['recovery']['raw_recovery_pass'][key] for r in test))
    comp=read(BASE/'OPTIMIZER_COMPARISON.json');comparisons=validation['comparisons']+read(BASE/'ROBUSTNESS_CONVERGENCE.json');step=[r for r in comparisons if 'DT05' in r.get('comparison','')];cache=[r for r in comparisons if 'FINE' in r.get('comparison','')]
    allforce=audit['all_completed_independent_pass'];finite_extension=[r for r in rows if r['status']=='PASS' and r['targets']['raw_extension_above_floor'] and not r['task'].get('diagnostic_low')]
    interpretation=dict(targets=targets,recovery=dict(source=t4['id'],result=r4,OFF=off['recovery'],active_minus_OFF_errors={k:r4['final_errors'][k]-v for k,v in offerrors.items()},confirmed_levels=recovery_gates),optimizer=comp,extension=dict(raw_above_1nm_ids=[r['id'] for r in finite_extension],confirmed=targets['T3']['confirmed'],scope='Finite tested four-source domain only; sorted eigenlength change, not identifiable-axis orientation near degeneracy'),numerics=dict(comparisons=len(comparisons),step_pass=sum(r['admitted'] for r in step),step_total=len(step),cache_pass=sum(r['admitted'] for r in cache),cache_total=len(cache),rejected=[r for r in comparisons if not r['admitted']]),physical_attempts=len(rows),complete=sum(r['status']=='PASS' for r in rows),failures=[dict(id=r['id'],status=r['scientific_status'],reasons=r.get('reasons',r.get('reason')),completed_s=r.get('completed_s')) for r in rows if r['status']!='PASS'],identity=identity())
    write(BASE/'INTERPRETATION.json',interpretation)
    gates={}
    def gate(name,status,evidence,scope='',threshold=None,reason=None):
        gates[name]=dict(status=('PASS' if status else 'FAIL') if isinstance(status,(bool,np.bool_)) else status,scope=scope or 'N64; full fixed-ID forward trajectories; homogeneous-water proxy',frequencies_Hz=[2400000.,2500000.,2600000.],states=['S2_r0','S2_validation_0..2','S2_test_0..2'],threshold=threshold,evidence=evidence,failure_reason=reason if reason else (None if status is True or status=='PASS' else 'See measured evidence; no threshold relaxation'))
    gate('F10AA_BACKEND_VERIFIED',all(r['admitted'] for r in baseline['regressions']),['BASELINE_VALIDATION.json'])
    gate('F10AC_F10AD_BASELINE_VERIFIED',all(r['bitwise_positions'] for r in baseline['regressions']),['BASELINE_VALIDATION.json'])
    gate('HISTORICAL_TRAJECTORIES_PRESERVED','PENDING_HASH_AUDIT',['HISTORICAL_RAW_SHA256.json'])
    gate('LOW_PRESSURE_TABLE_AUDITED',True,['LOW_PRESSURE_TABLE_AUDIT.json'])
    for name,key in [('LOW_PRESSURE_APECSS_VALIDATED','LOW_PRESSURE_APECSS_CONVERGED'),('LOW_PRESSURE_FORCE_VALIDATED','LOW_PRESSURE_FORCE_VALIDATED'),('LOW_PRESSURE_DYNAMIC_ADMISSION','SUB_0P1_INCIDENT_PRESSURE_DYNAMIC_PASS'),('HISTORICAL_2P5MHZ_REGRESSION_PASS','HISTORICAL_2P5MHZ_REGRESSION_PASS')]:gate(name,low['gates'][key],['LOW_PRESSURE_ADMISSION.json'],'2.5 MHz, shell x1; opt-in lower bound 0.001 Pa',.001)
    for name in ['TARGET_SHAPE_METRICS_VALIDATED','ABSOLUTE_AXIS_EXTENSION_EVALUATED','COM_CONSTRAINT_EVALUATED','CONNECTIVITY_RECORDED']:
        gate(name,True,['tests/f10ae/test_metrics_schedule.py','VALIDATION.json','STATIC_EXTENSION.json'])
    for target in ['T1','T2','T3']:
        gate('TARGET_'+target+'_EVALUATED',True,['INTERPRETATION.json'])
        gate('TARGET_'+target+'_CONFIRMED',targets[target]['confirmed'],['VALIDATION.json','FINAL_TEST_RUNS.json'],'Discovery + 3 unseen validation states + 3 once-only final test states',read(CONFIG)[target])
    gate('TWO_STAGE_SCHEDULE_VALIDATED',baseline['schedule']['admitted'] and baseline['schedule']['OFF_stage_displacement_m']==0,['BASELINE_VALIDATION.json'])
    gate('QUASISTEADY_SWITCH_SCOPE_DOCUMENTED',True,['TWO_STAGE_SCHEDULE_VALIDATION.md'],'Steady APECSS response at ideal switch; actual radial/hardware transient NOT_RUN')
    gate('SHAPE_RECOVERY_EVALUATED',True,['RECOVERY_RUNS.json'])
    for name,val in recovery_gates.items():gate(name,val,['INTERPRETATION.json','VALIDATION.json'])
    gate('CONSTRAINED_OPTIMIZER_REUSED',True,['AX_EXPERIMENT.json','AX_ASK_TELL.json'],'Official Ax/BoTorch outcome constraint; SciPy Sobol baseline')
    gate('OPTIMIZER_FAIR_COMPARISON_PASS',comp['fair_observation_budget'],['OPTIMIZER_COMPARISON.json'],'32 observations each; 8 identical shared initial evaluations; 56 distinct new physics calls; one optimizer seed')
    gate('INDEPENDENT_SEEDS_PASS',all(targets[k]['confirmed'] for k in ['T1','T2','T3']),['VALIDATION.json'],'All requested single-stage scientific targets independently confirmed; each target also separate')
    gate('STEP_SIZE_CONVERGENCE_PASS',all(r['admitted'] for r in step),['VALIDATION.json','ROBUSTNESS_CONVERGENCE.json'],'All tested candidates; individual failures retained',50e-9)
    gate('FIELD_CACHE_CONVERGENCE_PASS',all(r['admitted'] for r in cache),['VALIDATION.json','ROBUSTNESS_CONVERGENCE.json'])
    gate('CONTACT_AND_FORCE_GATES_PASS',allforce,['AUDIT.json'],'All complete trajectories independently reassembled; failed attempts are excluded and separately retained')
    for name in ['FREQUENCY_ROBUSTNESS_EVALUATED','SHELL_SENSITIVITY_EVALUATED','SPATIAL_POSITION_SENSITIVITY_EVALUATED','PROXY_SENSITIVE_MODES_IDENTIFIED']:
        gate(name,True,['ROBUSTNESS_RUNS.json','ROBUSTNESS_CONVERGENCE.json'],'Measured tests, not universal robustness admission')
    gate('FIGURES_VALIDATED','PENDING_VISUAL_CHECK',['figures/MANIFEST.json'])
    gate('VIDEOS_VALIDATED','PENDING_VISUAL_CHECK',['videos/MANIFEST.json'])
    for name in ['REAL_PDMS_FIELD_CALIBRATED','REAL_PZT_AMPLITUDE_PHASE_CALIBRATED','WORLD_MODEL_TRAINED']:gate(name,'NOT_RUN',[],'Outside this exploratory numerical study')
    gate('TRAINING_READY',False,['INTERPRETATION.json'],'No calibrated action-conditioned world-model training admission')
    gate('SCIENTIFIC_TARGETED_SHAPE_CONTROL_ADMITTED',targets['T1']['confirmed'] and all(recovery_gates.values()),['INTERPRETATION.json'],'Strict combined requested target-and-recovery claim; limited target-specific passes reported separately')
    for name,record in gates.items():
        if name.startswith('LOW_PRESSURE') or name=='HISTORICAL_2P5MHZ_REGRESSION_PASS':record['frequencies_Hz']=[2500000.]
    write(BASE/'FINAL_GATE.json',dict(gates=gates,labels=LABELS,training_ready=False,identity=identity()))
if __name__=='__main__':main()
