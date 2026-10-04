"""Conservative synthesis of map, variational, domain and continuation evidence."""
import numpy as np
from scipy.linalg import solve
from pd_clot.pfpd_candidate.coupling_stability_audit import (
    OLD,OUT,audit_config,read_json,save_json,read_npz,load_audit_state,deterministic_starts)


def spectrum_summary():
    a=audit_config();summary={};alignments=[];relaxed=[];probes=[]
    paths=read_json(OUT/'energy_path/SAVED_PATHS.json')
    failed=[r for r in paths['loading_unloading'] if r['load_step']==10 and r['kind']=='phase_after']
    for key in ['A','B','C']:
        rows=[read_json(OUT/f'spectrum/{key}/eps_{eps:g}/SUMMARY.json') for eps in a['jacobian_relative_epsilon']]
        matrices=[read_npz(OUT/f'spectrum/{key}/eps_{eps:g}/MATRICES.npz') for eps in a['jacobian_relative_epsilon']]
        rhos=np.array([r['selected_branch_rho'] for r in rows]);change=float(np.ptp(rhos)/max(abs(rhos).max(),1e-30))
        stable=change<=a['spectral_epsilon_relative_tolerance']
        ordinary=all(r['ordinary_Phi_spectrum_valid'] for r in rows) and stable
        summary[key]=dict(ordinary_status='RESOLVED' if ordinary else 'UNRESOLVED',
            ordinary_rho=float(np.median(rhos)) if ordinary else None,
            conditional_branch_rho=float(np.median(rhos)),conditional_epsilon_relative_spread=change,
            conditional_epsilon_stable=stable,rows=rows)
        center=matrices[1];lam=center['eigenvalues'][0];v=center['eigenvectors'][:,0]
        summary[key]['conditional_dominant']=dict(real=float(lam.real),imag=float(lam.imag),magnitude=float(abs(lam)))
        for omega in a['relaxation_omega_diagnostic']:
            relaxed.append(dict(state=key,omega=omega,conditional=not ordinary,
                lambda_real=float(lam.real),lambda_imag=float(lam.imag),
                relaxed_real=float(((1-omega)+omega*lam).real),relaxed_imag=float(((1-omega)+omega*lam).imag),
                relaxed_magnitude=float(abs((1-omega)+omega*lam))))
        if key=='C' and abs(lam.imag)<1e-10:
            v=v.real;chosen=read_npz(OLD/'loading_unloading/load_010/iteration_006_phase.npz')
            if v@(chosen['s_new']-chosen['s_old'])<0:v=-v
            for r in failed:
                state=read_npz(r['source']);delta=state['s_new']-state['s_old']
                cosine=float(v@delta/(np.linalg.norm(v)*np.linalg.norm(delta))) if np.linalg.norm(delta) else None
                alignments.append(dict(iteration=r['iteration'],delta_L2=r['delta_s_L2'],L2_ratio=r['L2_ratio'],
                    Linf_ratio=r['Linf_ratio'],successive_cosine=r['increment_cosine'],mode_cosine=cosine,
                    conditional_lambda_magnitude=float(abs(lam))))
            p,data=load_audit_state(key);base=p.phi(data['s'],data['x_start']);v/=np.max(abs(v))
            for eps in a['jacobian_relative_epsilon']:
                h=eps*max(.1,float(data['s'].max()))
                for sign in [-1,1]:
                    direction=sign*v;r=p.phi(data['s']+h*direction,data['x_start'])
                    row=dict(state=key,epsilon=eps,sign=sign,status=r['status'])
                    if r['status']=='OK':
                        measured=(r['s_output']-base['s_output'])/h;predicted=lam.real*direction
                        row.update(relative_JVP_error=float(np.linalg.norm(measured-predicted)/np.linalg.norm(predicted)),
                            observed_directional_gain=float(np.linalg.norm(measured)/np.linalg.norm(direction)),
                            selector_switches=int(np.sum(r['active_set']!=base['active_set'])),
                            minimum_Green_strain=r['before']['minimum_Green_strain'])
                    probes.append(row)
    save_json(OUT/'spectrum/EPSILON_STABILITY.json',summary)
    save_json(OUT/'spectrum/MODE_ALIGNMENT.json',dict(scope='Selected frozen branch at C; actual Phi is nonsmooth',rows=alignments,directional_probes=probes))
    middle=[r for r in alignments if 4<=r['iteration']<=20]
    save_json(OUT/'spectrum/OBSERVED_MODE_INTERPRETATION.json',dict(
        status='OBSERVED_DIVERGENCE_IS_CONSISTENT_WITH_DOMINANT_FIXED_POINT_MODE',
        scope='CONDITIONAL frozen branch only; not proof of an ordinary Phi Jacobian',
        iteration_range=[4,20],minimum_mode_cosine=min(r['mode_cosine'] for r in middle),
        maximum_relative_amplification_difference=max(abs(r['L2_ratio']/r['conditional_lambda_magnitude']-1) for r in middle),
        ordinary_Phi_spectrum='UNRESOLVED',directional_derivative_is_not_two_sided=True))
    save_json(OUT/'spectrum/RELAXATION_THEORY_ONLY.json',relaxed)
    # Independent cross-check of the smooth trial-state derivative, using the
    # variational blocks and the mechanical implicit-function derivative.
    cross=[]
    for key in ['B','C']:
        p,data=load_audit_state(key);b=read_npz(OUT/f'hessian/{key}/eps_3e-06/BLOCKS.npz')
        fd=read_npz(OUT/f'spectrum/{key}/eps_3e-06/MATRICES.npz')
        state=p.state(data['x'],data['s']);Y=state.W0
        multiplier=p.yc/Y**2*p.energy_scale/(2*p.weights*(1-data['s']))
        trial_H=multiplier[:,None]*(b['H_sx']@solve(b['H_xx'],b['H_xs'],assume_a='gen'))
        relative=float(np.linalg.norm(trial_H-fd['trial_derivative'])/np.linalg.norm(fd['trial_derivative']))
        cross.append(dict(state=key,trial_derivative_Hessian_vs_FD_relative_error=relative,
            interpretation='Smooth stationary trial branch only; not a replacement for switched Phi'))
    save_json(OUT/'comparisons/HESSIAN_MAP_CROSSCHECK.json',cross)
    return summary,alignments


def curvature_and_reference():
    a=audit_config();rows=read_json(OUT/'hessian/BASE_STATES.json')+read_json(OUT/'hessian/OPTIMIZER_STATES.json')
    grouped={}
    for r in rows:grouped.setdefault(r['state'],[]).append(r)
    curvature={}
    for key,values in grouped.items():
        eigs=np.array([r['minimum_reduced_eigenvalue'] for r in values])
        stable=bool(np.ptp(eigs)/max(abs(eigs).max(),1e-30)<=a['spectral_epsilon_relative_tolerance'])
        negative=stable and all(r['status']=='NEGATIVE_CURVATURE_FOUND' for r in values)
        positive=stable and all(r['status']=='LOCALLY_POSITIVE' and r['weak_lower_DOFs']==0 for r in values)
        feasible=all(any(all(p.get('bounds',False) and p.get('in_tensile_domain',False) for p in probe['points'])
                         for probe in r['minimum_mode_probes']) for r in values)
        curvature[key]=dict(minimum_eigenvalues=eigs.tolist(),epsilon_stable=stable,
            feasible_negative_mode=negative and feasible,
            status='NEGATIVE_CURVATURE_FOUND' if negative and feasible else 'LOCALLY_POSITIVE' if positive else 'INCONCLUSIVE',
            min_Hxx_eigenvalue=min(r['minimum_mechanical_eigenvalue'] for r in values),
            Hxx_condition_max=max(r['Hxx_condition_number'] for r in values),
            min_Schur_eigenvalue=min([r['minimum_Schur_eigenvalue'] for r in values if r['minimum_Schur_eigenvalue'] is not None],default=None))
    save_json(OUT/'hessian/EPSILON_STABILITY.json',curvature)
    optimizers=read_json(OUT/'monolithic/ALL_OPTIMIZERS.json');p,starts=deterministic_starts()
    first_x,first_s=starts['failed_iteration_000_phase_after'];late_x,late_s=starts['failed_iteration_023_phase_after']
    first_energy=p.metrics(first_x,first_s)['total_J'];late_energy=p.metrics(late_x,late_s)['total_J']
    previous_x,previous_s=starts['previous_accepted_with_affine_boundary_lift']
    prev_energy=p.metrics(previous_x,previous_s)['total_J'];comparison=[]
    for row in optimizers:
        key=row['method']+'_'+row['start'];curv=curvature.get(key)
        accepted=bool(row['first_order_candidate_passed'] and curv and curv['status']=='LOCALLY_POSITIVE')
        state=read_npz(OUT/'monolithic'/row['method']/row['start']/'STATE.npz')
        comparison.append(dict(method=row['method'],start=row['start'],strain=row['strain'],
            optimizer_success=row['optimizer_success'],first_order_passed=row['first_order_candidate_passed'],
            curvature=None if not curv else curv['status'],valid_local_minimum=accepted,metrics=row['metrics'],
            domain_exit=row['exit'] is not None,exit_details=row['exit'],
            recorded_point_kind='optimizer_endpoint' if row['optimizer_success'] else 'last_domain_valid_evaluation_not_a_solution',
            relative_to_previous_affine_lift_J=row['metrics']['total_J']-prev_energy if row['strain']==p.strain else None,
            relative_to_failed_first_phase_J=row['metrics']['total_J']-first_energy if row['strain']==p.strain else None,
            relative_to_failed_late_phase_J=row['metrics']['total_J']-late_energy if row['strain']==p.strain else None,
            max_x_distance_to_previous_lift_m=float(abs(state['x']-previous_x).max()),
            max_s_distance_to_previous_load=float(abs(state['s']-previous_s).max())))
    save_json(OUT/'monolithic/ACCEPTANCE_AND_COMPARISON.json',dict(previous_affine_lift_energy_J=prev_energy,
        failed_first_phase_energy_J=first_energy,failed_late_phase_energy_J=late_energy,rows=comparison))
    previous_problem,previous_data=load_audit_state('B')
    references=[]
    for label,x,s in [('previous_accepted_affine_lift_at_target',previous_x,previous_s),
                      ('first_failed_phase_after',first_x,first_s),('late_failed_phase_after',late_x,late_s)]:
        references.append(dict(state=label,strain=p.strain,metrics=p.metrics(x,s),
            max_x_distance_to_previous_lift_m=float(abs(x-previous_x).max()),
            max_s_distance_to_previous_load=float(abs(s-previous_s).max())))
    save_json(OUT/'energy_path/REFERENCE_STATES.json',dict(
        previous_accepted_at_own_load=previous_problem.metrics(previous_data['x'],previous_data['s']),
        previous_accepted_strain=previous_problem.strain,target_strain=p.strain,target_load_states=references,
        caution='Subtracting the old-load energy is a plotting baseline, not a same-load stability comparison. The affine lift is an initial guess, not an accepted solution.'))
    return curvature,comparison


def ordering_and_continuation():
    a=audit_config();paths=read_json(OUT/'energy_path/SAVED_PATHS.json');ordering={}
    for name,step in [('loading_unloading',10),('loading_unloading_extra_iterations',8)]:
        rows=[r for r in paths[name] if r['load_step']==step];phase=[r for r in rows if r['kind']=='phase_after']
        gate=[r for r in rows if r['kind']=='mechanical_gate_stop']
        amplified=[r for r in phase if r.get('L2_ratio') is not None and r['L2_ratio']>1+a['rho_margin'] and r['phase_increment']>1e-9]
        start=min([r['iteration'] for r in amplified],default=None);end=gate[0]['iteration'] if gate else None
        # Midrange excludes roundoff-scale increments and late nonlinear growth.
        middle=[r for r in phase if 1e-9<r['phase_increment']<1e-4 and r.get('L2_ratio') is not None]
        def limits(key):return [min(r[key] for r in middle),max(r[key] for r in middle)] if middle else None
        ordering[name]=dict(status='FOLLOWS_DIVERGENCE' if start is not None and end is not None and start<end else 'INCONCLUSIVE',
            first_clearly_amplifying_iteration=start,first_domain_exit_iteration=end,
            L2_amplification_range=limits('L2_ratio'),Linf_amplification_range=limits('Linf_ratio'),successive_cosine_range=limits('increment_cosine'),
            scope='Temporal ordering supports a later domain violation; it does not establish causality',
            max_phase_energy_increase_J=max(r['phase_energy_change_J'] for r in phase),
            max_mechanical_energy_increase_J=max(r['mechanical_energy_change_J'] for r in phase))
    save_json(OUT/'energy_path/TENSION_GATE_ORDERING.json',ordering)
    reference=read_npz(OUT/'monolithic/L-BFGS-B/failed_iteration_000_phase_after/STATE.npz')
    rows=[]
    for r in read_json(OUT/'continuation/ALL_SUBSTEPS.json'):
        state=read_npz(OUT/f"continuation/substeps_{r['substeps']}/FINAL_STATE.npz")
        rows.append(dict(substeps=r['substeps'],endpoint_reached=r['endpoint_reached'],status=r['status'],reason=r['reason'],
            accepted_substeps=len(r['accepted_substeps']),last_accepted_strain=r['accepted_substeps'][-1]['strain'] if r['accepted_substeps'] else a['subdivision_start'],
            attempted_strain=r['stopped_or_last_strain'],
            endpoint_max_s=float(state['s'].max()) if r['endpoint_reached'] else None,
            endpoint_energy_J=r['last_computed_metrics']['total_J'] if r['endpoint_reached'] else None,
            endpoint_max_delta_x_to_stationary_reference_m=float(abs(state['x']-reference['x']).max()) if r['endpoint_reached'] else None,
            endpoint_max_delta_s_to_stationary_reference=float(abs(state['s']-reference['s']).max()) if r['endpoint_reached'] else None,
            last_audited_inner_state=r['last_computed_metrics'],
            reference_is_a_saddle_not_a_valid_local_minimum=True))
    save_json(OUT/'continuation/ENDPOINT_COMPARISON.json',dict(label='STEP_SUBDIVISION_DIAGNOSTIC_ONLY',rows=rows))
    return ordering


def main():
    spectra,alignments=spectrum_summary();curvature,optimizers=curvature_and_reference();ordering=ordering_and_continuation()
    failure_cases=[r for r in optimizers if r['strain']==audit_config()['subdivision_endpoint']]
    negative_stationary=[r for r in failure_cases if r['first_order_passed'] and r['curvature']=='NEGATIVE_CURVATURE_FOUND']
    valid_minima=[r for r in failure_cases if r['valid_local_minimum']]
    root='COUPLED_BRANCH_LOSS_OR_LOCALIZATION_INSTABILITY' if negative_stationary else 'INCONCLUSIVE'
    status={'COUPLING_STABILITY_AUDIT':'COMPLETE','FIXED_POINT_MAP':read_json(OUT/'fixed_point/REPRODUCTION.json')['status'],
        'FIXED_POINT_SPECTRAL_RADIUS':'RESOLVED' if all(r['ordinary_status']=='RESOLVED' for r in spectra.values()) else 'UNRESOLVED',
        'MONOLITHIC_REFERENCE':'LOCAL_MINIMUM_FOUND' if valid_minima else 'NO_VALID_LOCAL_MINIMUM_FOUND' if negative_stationary else 'INCONCLUSIVE',
        'REDUCED_HESSIAN':'NEGATIVE_CURVATURE_FOUND' if negative_stationary else 'INCONCLUSIVE',
        'ROOT_CAUSE':root,'PURE_TENSION_GATE':'FOLLOWS_DIVERGENCE' if all(r['status']=='FOLLOWS_DIVERGENCE' for r in ordering.values()) else 'INCONCLUSIVE',
        'STEP_SUBDIVISION':'DIAGNOSTIC_ONLY','PFPD_AT1_EVOLUTION_V1':'REMAINS_NOT_VERIFIED','PRODUCTION MODEL':'UNCHANGED'}
    save_json(OUT/'STATUS.json',status)
    print(status,flush=True)


if __name__=='__main__':main()
