"""Independent read-only checks of archived states, clocks, events and identities."""
import time
import subprocess
import numpy as np
import h5py
from mevo_mb.control.scheduled_pzt import Schedule, Segment
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from .common import BASE, ROOT, OLD, Q, CONFIG, read, write, sha, identity, policy, context


def complex_array(values):
    return np.array([complex(v['real'],v['imag']) for v in values])


def main():
    begin=time.perf_counter();checks={};evidence={};p=policy();cfg=read(CONFIG)
    def check(name,value):
        checks[name]=bool(value)
        if not value:raise AssertionError(name)
    baseline=read(BASE/'cleanup/BASELINE.json')
    changed=[s for s,h in baseline['protected_sha256'].items() if not (ROOT/s).exists() or sha(ROOT/s)!=h]
    check('950_frozen_source_config_test_dependency_files',not changed and len(baseline['protected_sha256'])==950)
    # Detect modification/removal of any pre-existing tracked history, not just code.
    diff=subprocess.check_output(['git','diff','--name-status',cfg['public_baseline'],'--','.'],cwd=ROOT,text=True)
    check('no_modified_or_removed_baseline_tracked_files',all(s.startswith('A\t') for s in diff.splitlines()))
    remote=read(BASE/'REMOTE_MANIFEST.json')
    check('every_returned_file_SHA256',all(sha(BASE/k)==v['sha256'] for k,v in remote['files'].items()))
    check('preregistered_input_hashes',all(sha(ROOT/k)==v for k,v in read(BASE/'PREREGISTRATION.json')['input_sha256'].items()))
    test=read(BASE/'local/PROTOCOL_TESTS.json')
    check('protocol_tests_match_final_module',test['returncode']==0 and test['protocol_source_sha256']==sha(ROOT/'src/mevo_mb/control/scheduled_pzt.py'))
    with h5py.File(BASE/'pilot_v1/runs/X_COMPAT_DT1/trajectory.h5') as a,h5py.File(Q/'b_navigation_v1/runs/NAV_X_DT1/trajectory.h5') as b:
        check('full_precision_X_replay_all_steps_bitwise_identical',np.array_equal(a['positions_m'][:],b['positions_m'][:]))
    check('inherited_CPU_GPU_APECSS_wall_checks',read(BASE/'pilot_v1/CHECKS.json')['status']=='PASS')
    ctx=context();bases={'coarse':ctx['basis'],'fine':CachedWaterBasis(OLD/'cache/fine.h5')}
    initial=dict(np.load(BASE/'INITIAL_STATES.npz'));rows=[]
    ii,jj=np.triu_indices(32,1)
    for path in sorted(BASE.glob('*/runs/*/RESULT.json')):
        r=read(path);folder=path.parent;inp=read(folder/'INPUT.json');tag=r['name']
        schedule=Schedule.from_record(inp['schedule']);events=read(folder/'SWITCHES.json')
        with h5py.File(folder/'trajectory.h5') as h:
            t=h['time_s'][:];x=h['positions_m'][:];sm=h['step_metrics'][:];ex=h['extra_metrics'][:];clock=h['step_time_duration_segment'][:]
            ids=h['particle_id'][:];force_time=h['force_samples/time_s'][:];reactions=h['pair_force_N'][:]
            check(tag+'_strict_zero_flow',h.attrs['flow']=='STRICT_ZERO')
            check(tag+'_source_commit',h.attrs['execution_commit']==inp['identity']['execution_commit'])
            check(tag+'_fixed_32_IDs',np.array_equal(ids,np.arange(32)) and x.shape[1:]==(32,3))
        check(tag+'_trajectory_SHA256',sha(folder/'trajectory.h5')==r['trajectory_sha256'])
        check(tag+'_registered_initial_exact',np.array_equal(x[0],initial[inp['initial']]))
        check(tag+'_finite_states',np.isfinite(x).all())
        check(tag+'_no_missing_or_repeated_time',np.all(np.diff(t)>0) and np.allclose(np.diff(t),clock[sm[:,0].astype(bool),1],rtol=0,atol=1e-15))
        check(tag+'_ideal_clock_start_and_end',t[0]==0 and abs(t[-1]-r['completed_s'])<1e-15)
        check(tag+'_scheduled_steps',np.allclose(clock,[(t,h,s.segment_id) for t,h,s in schedule.steps(inp['dt_s'])][:len(clock)],rtol=0,atol=1e-15))
        overlap=float(np.maximum(2*p['radius_m']-np.linalg.norm(x[:,ii]-x[:,jj],axis=-1),0).max())
        check(tag+'_independent_true_sphere_overlap',overlap<=min(p['contact_overlap_m'],p['hard_geometry_m']))
        radial=float(np.linalg.norm(x[:,:,1:],axis=-1).max()+p['radius_m'])
        axial=float(abs(x[:,:,0]).max()+p['radius_m'])
        check(tag+'_inside_main_tube_and_sealed_ends',radial<.0008 and axial<.01625)
        accepted=sm[sm[:,0].astype(bool)]
        for label,col,limit in [('mobility',1,p['mobility_residual_m']),('iteration_correction',3,p['contact_correction_m']),
                                ('total_correction',4,p['contact_correction_m']),('wall',5,0),('primal',12,p['primal_violation_m']),
                                ('dual',13,p['dual_negative_N']),('linear_complementarity',15,p['complementarity_equivalent_gap_m']),
                                ('true_complementarity',17,p['actual_geometry_complementarity_m'])]:
            check(tag+'_'+label,float(accepted[:,col].max(initial=0))<=limit)
        force_scale=np.maximum.reduce([ex[:,5],abs(reactions).max(axis=1),np.full(len(sm),1e-18)])
        check(tag+'_true_force_equation',np.all(accepted[:,14]<=p['force_residual_atol_N']+p['force_residual_rtol']*force_scale[sm[:,0].astype(bool)]))
        # Original table queried strictly, including exact zero; no clipped pressure accepted.
        ctx['response'].query(ex[:,3:5])
        max_event_error=0.
        for j,event in enumerate(events):
            s=schedule.segments[j];k=np.argmin(abs(t-s.start_s));xx=np.array(event['positions_m'])
            check(tag+f'_event_{j}_time_and_state',event['ideal_command_time_s']==s.start_s and event['effective_time_s']==s.start_s and np.array_equal(xx,x[k]))
            pb,_=bases[inp['cache']].query_basis(xx)
            for when,control in [('before',event['before_control']),('after',event['after_control'])]:
                expected=np.zeros(32,complex)
                if control is not None:
                    e=Segment(**control).event();expected=bases[inp['cache']].vector(e['amplitude_Pa'],e['phase_rad'],e['enabled'])@pb
                measured=complex_array(event[when+'_pressure_Pa'])
                err=float(abs(expected-measured).max());max_event_error=max(max_event_error,err)
                check(tag+f'_event_{j}_{when}_actual_pressure',err<1e-8)
        rows.append(dict(name=tag,frames=len(t),steps=len(sm),max_independent_overlap_m=overlap,
                         min_wall_margin_m=min(.0008-radial,.01625-axial),max_event_pressure_error_Pa=max_event_error,
                         last_force_sample_time_s=float(force_time[-1]),status=r['status']))
    schedules=read(BASE/'RUN_SCHEDULES.json');original=cfg['original_X_control']
    for shift in cfg['phase_full_trajectory_deg']:
        seg=Schedule.from_record(schedules[f'X_PHASE_{shift:+d}']).segments[0]
        check(f'phase_{shift:+d}_only_Xplus_changes',list(seg.amplitude_Pa_equivalent)==original['amplitude_Pa'] and list(seg.enabled)==original['enabled'] and list(seg.phase_rad[:3])==original['phase_rad'][:3] and abs(seg.phase_rad[3]-original['phase_rad'][3]-np.deg2rad(shift))<1e-15)
    evidence.update(trajectories=rows,remote_files_verified=len(remote['files']),preserved_baseline_files=len(baseline['protected_sha256']),
                    response_intervals_Pa=ctx['response'].metadata['intervals_Pa'],valid_response_intervals=ctx['response'].metadata['valid'],
                    caveat='Table endpoints and strict simulator admission verified; no new experiment or contact-physics calibration')
    elapsed=time.perf_counter()-begin
    record=dict(identity=identity(),status='PASS',checks=checks,evidence=evidence,wall_s=elapsed)
    write(BASE/'local/INDEPENDENT_VERIFICATION.json',record)
    gate=read(BASE/'analysis/FINAL_GATE.json');gate['ENGINEERING_PASS']='PASS'
    gate['engineering_evidence']='local/INDEPENDENT_VERIFICATION.json';write(BASE/'analysis/FINAL_GATE.json',gate)
    cost=read(BASE/'analysis/COSTS.json')
    cost['local_protocol_tests_s']=test['wall_s'];cost['independent_verification_s']=elapsed
    cost['total_scientific_s']=sum(cost[k] for k in ['remote_campaign_wall_s','local_field_diagnostics_s','analysis_s','local_protocol_tests_s','independent_verification_s'])
    check('budget_600s',cost['total_scientific_s']<=cfg['total_scientific_budget_s'])
    cost['budget_scope']='All four remote processes (including checks, output and loading), field-only diagnostics/refinement/branch review, actual analysis, protocol tests and independent numerical verification. File transfer, figure rendering and engineering/development time separately excluded.'
    write(BASE/'analysis/COSTS.json',cost)
    print('PASS:',len(checks),'checks;',len(rows),'trajectories;',len(remote['files']),'remote files; science seconds',cost['total_scientific_s'])


if __name__=='__main__':main()
