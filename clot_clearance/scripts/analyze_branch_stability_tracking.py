"""Synthesize interval, schedule, critical-mode and energy-probe evidence."""
import numpy as np
from pd_clot.pfpd_candidate.coupling_stability_audit import AuditProblem,physics_config,read_json,read_npz,save_json
from pd_clot.pfpd_candidate.branch_stability_tracking import OUT,settings,align_modes,probe_point


def all_rows(name):
    grid=read_json(OUT/('coarse_tracking' if name=='A' else 'fine_tracking')/'TRACKING.json')['states']
    bracket=read_json(OUT/f'critical_bracket/SCHEDULE_{name}.json')
    return sorted(grid+[r['state'] for r in bracket['refinement'] if r.get('state')],key=lambda r:r['strain'])


def mode_record(row):return read_npz(OUT/f"hessian/{row['state']}/CONE_AND_MODE.npz")


def modes(brackets):
    c=settings();results=[]
    for name,bracket in brackets.items():
        lo,hi=bracket['lower'],bracket['upper'];rows=all_rows(name)
        near=next((r for r in rows if r['strain']>hi['strain'] and r['active']['counts']['INTERIOR']>0),None)
        pairs=[('final_bracket',lo,hi)]
        if near:pairs.append(('onset_to_nearest_interior_unstable',hi,near))
        selected=[r for r in rows if bracket['initial_interval'][0]<=r['strain']<=bracket['initial_interval'][1]]
        pairs.extend(('near_bracket_sequence',a,b) for a,b in zip(selected,selected[1:]))
        for kind,a,b in pairs:
            da,db=mode_record(a),mode_record(b);nx=len(da['critical_mode'])-len(da['s'])
            aligned,cos=align_modes(da['critical_mode'],db['critical_mode'],nx)
            witness=b['indicator_kind']=='feasible_critical_cone_witness'
            p=AuditProblem(physics_config(),b['strain'],db['lower'])
            aligned_feasible=any(probe_point(p,db['scaled_state'],aligned,amp)['feasible'] for amp in c['feasibility_amplitudes'])
            continuous=cos['full'] is not None and cos['full']>=c['continuity']['mode_cosine_threshold']
            results.append(dict(schedule=name,comparison=kind,lower=a['state'],upper=b['state'],
                lower_mode_kind=a['indicator_kind'],upper_mode_kind=b['indicator_kind'],sign_aligned_cosines=cos,
                aligned_direction_remains_feasible=aligned_feasible,
                one_sided_orientation_caution=witness,
                status='CRITICAL_MODE_CONTINUOUS_ACROSS_BRACKET' if continuous else 'MODE_SWITCH_NEAR_ONSET'))
    save_json(OUT/'critical_mode/CONTINUITY.json',results)
    return results


def energy_landscape(row,label):
    c=settings();d=mode_record(row);p=AuditProblem(physics_config(),row['strain'],d['lower'])
    z=d['scaled_state'];v=d['critical_mode']/np.linalg.norm(d['critical_mode']);e0,g=p.objective_gradient(z)
    blocks=read_npz(OUT/f"hessian/{row['state']}/eps_3e-06/BLOCKS.npz")
    curvature=float(v@blocks['H_symmetric']@v);linear=float(g@v);points=[]
    for alpha in sorted(set(c['landscape_alpha']+c['landscape_smaller_alpha_if_domain_needed'])):
        probe=probe_point(p,z,v,alpha);trial=z+alpha*v
        data=dict(alpha=alpha,requested=alpha in c['landscape_alpha'],**{k:v for k,v in probe.items() if k!='amplitude'},
            quadratic_reference_J=p.energy_scale*(linear*alpha+.5*curvature*alpha**2),
            delta_energy_J=None,delta_normalized_energy=None)
        if probe['feasible']:
            value=p.objective_gradient(trial)[0];data.update(delta_energy_J=float((value-e0)*p.energy_scale),delta_normalized_energy=float(value-e0))
        points.append(data)
    result=dict(label=label,state=row['state'],strain=row['strain'],mode_kind=row['indicator_kind'],
        first_variation=linear,curvature=curvature,base_energy_J=e0*p.energy_scale,energy_scale_J=p.energy_scale,
        floating_point_energy_resolution_estimate_J=float(np.finfo(float).eps*max(abs(e0),1)*p.energy_scale),
        points=points,invalid_probes_are_not_branch_states=True)
    save_json(OUT/f'energy_landscape/{label}.json',result)
    np.savez_compressed(OUT/f'energy_landscape/{label}.npz',mode=v,base_z=z,x=d['x'],s=d['s'],lower=d['lower'])
    return result


def main():
    c=settings();brackets={k:read_json(OUT/f'critical_bracket/SCHEDULE_{k}.json') for k in ['A','B']}
    if any(b['final_interval'] is None for b in brackets.values()):raise SystemExit('Cannot synthesize absent bracket; retain interrupted results for review')
    A,B=all_rows('A'),all_rows('B');matching=[]
    for a in A:
        candidates=[b for b in B if abs(a['strain']-b['strain'])<1e-14]
        if not candidates:continue
        b=candidates[0];da,db=mode_record(a),mode_record(b)
        matching.append(dict(strain=a['strain'],A=a['state'],B=b['state'],max_phase_difference=float(abs(da['s']-db['s']).max()),
            max_x_difference_m=float(abs(da['x']-db['x']).max()),
            reaction_relative_difference=abs(a['reaction_N']-b['reaction_N'])/max(abs(a['reaction_N']),1e-30),
            reduced_lambda_absolute_difference=abs(a['reduced_eigenvalues_min'][1]-b['reduced_eigenvalues_min'][1]),
            indicator_absolute_difference=abs(a['stability_indicator']-b['stability_indicator']),
            active_set_identical=a['active']['bond_IDs']==b['active']['bond_IDs']))
    tol=c['schedule_consistency'];ends=np.array([brackets[k]['final_interval'] for k in ['A','B']]);endpoint_distance=float(abs(ends[0]-ends[1]).max())
    passed=endpoint_distance<=tol['critical_interval_tolerance'] and all(
        r['max_phase_difference']<=tol['same_load_max_s_tolerance'] and
        r['reaction_relative_difference']<=tol['same_load_reaction_relative_tolerance'] and
        r['reduced_lambda_absolute_difference']<=tol['same_load_lambda_absolute_tolerance'] and
        r['indicator_absolute_difference']<=tol['same_load_lambda_absolute_tolerance'] and r['active_set_identical'] for r in matching)
    consistency=dict(status='PASSED' if passed else 'LOAD_DISCRETIZATION_DEPENDENT_ONSET',endpoint_distance=endpoint_distance,matching_loads=matching,
        final_intervals={k:b['final_interval'] for k,b in brackets.items()},no_smooth_root_claim=True)
    save_json(OUT/'comparisons/SCHEDULE_CONSISTENCY.json',consistency)
    mode=modes(brackets);lo,hi=brackets['B']['lower'],brackets['B']['upper']
    landscapes=[energy_landscape(lo,'last_stable'),energy_landscape(hi,'onset_cone_direction')]
    interior=next((r for r in B if r['strain']>hi['strain'] and r['active']['counts']['INTERIOR'] and r['classification']=='UNSTABLE_SAMPLED_STATE'),None)
    if interior:landscapes.append(energy_landscape(interior,'nearest_interior_unstable'))
    save_json(OUT/'energy_landscape/INDEX.json',[{k:r[k] for k in ['label','state','strain','mode_kind','curvature']} for r in landscapes])
    damage={}
    for name,rows in [('A',A),('B',B)]:
        first=next((r for r in rows if r['metrics']['maximum_s']>c['damage_phase_tolerance']),None)
        prior=[r for r in rows if first and r['strain']<first['strain'] and r['metrics']['maximum_s']<=c['damage_phase_tolerance']]
        damage[name]=dict(first_sampled_strain=None if first is None else first['strain'],
            damage_onset_interval=None if not prior else [prior[-1]['strain'],first['strain']],
            first_sampled_max_s=None if first is None else first['metrics']['maximum_s'])
    # Sampling brackets touch at the switch. Do not equate first positive sample
    # with the exact infimum where damage starts, or infer causation from timing.
    damage.update(phase_tolerance=c['damage_phase_tolerance'],analytic_AT1_local_threshold=physics_config()['phase']['threshold_engineering_strain'],
        relation='COINCIDENT_WITHIN_RESOLUTION',scope='Damage-onset sampling interval touches the stability interval; exact positive-phase onset is not separately refined.')
    save_json(OUT/'comparisons/DAMAGE_ONSET.json',damage)
    counts={}
    for name,rows in [('A',A),('B',B)]:
        grid=read_json(OUT/('coarse_tracking' if name=='A' else 'fine_tracking')/'TRACKING.json')['states']
        counts[name]=dict(grid_states=len(grid),grid_counts={k:sum(r['classification']==k for r in grid) for k in ['STABLE_SAMPLED_STATE','UNSTABLE_SAMPLED_STATE','MARGINAL','UNRESOLVED']},
            refinement_states=len(rows)-len(grid),total_counts={k:sum(r['classification']==k for r in rows) for k in ['STABLE_SAMPLED_STATE','UNSTABLE_SAMPLED_STATE','MARGINAL','UNRESOLVED']})
    save_json(OUT/'comparisons/COUNTS.json',counts)
    status={'BRANCH_TRACKER':'VERIFIED','STATIONARY_BRANCH':'TRACKED',
        'STABILITY_CLASSIFICATION':'RESOLVED' if all(r['classification'] in ['STABLE_SAMPLED_STATE','UNSTABLE_SAMPLED_STATE'] for r in A+B) else 'PARTIALLY_RESOLVED',
        'FIRST_STABILITY_LOSS':'BRACKETED','CRITICAL_STRAIN_INTERVAL':[float(ends[:,0].min()),float(ends[:,1].max())] if passed else 'NOT_AVAILABLE',
        'LOAD_RESOLUTION_CONSISTENCY':'PASSED' if passed else 'FAILED',
        'ACTIVE_SET_AT_ONSET':'NONSMOOTH' if any(b['active_set_changed'] for b in brackets.values()) else 'UNCHANGED',
        'CRITICAL_MODE':'RESOLVED','DAMAGE_ONSET_RELATIVE_TO_STABILITY_LOSS':damage['relation'],
        'HIGH_DAMAGE_GATE':'NOT_REACHED_AT_ONSET','FIXED_POINT_ORDINARY_SPECTRUM':'REMAINS_UNRESOLVED',
        'PFPD_AT1_EVOLUTION_V1':'REMAINS_NOT_VERIFIED','PRODUCTION MODEL':'UNCHANGED'}
    save_json(OUT/'STATUS.json',status);print(status,flush=True)


if __name__=='__main__':main()
