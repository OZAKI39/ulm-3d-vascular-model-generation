"""Observe saved trajectories and existing solver diagnostics; never advance dynamics."""
from pathlib import Path
import csv,gzip,json,hashlib,sys,time
from collections import Counter
import numpy as np
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R.parents[1]/'src'))
from runner import make_environment,safe
from particle_3d.validation_boundary import ValidationBoundaryClassifier

def read(p):return json.loads(Path(p).read_text())
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2,default=safe,ensure_ascii=False,allow_nan=False)+'\n')
def csvwrite(p,rows):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
        w.writerows({k:json.dumps(v,default=safe,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()} for row in rows)
def stats(a):
    a=np.asarray(a,float)
    return dict(n=len(a),minimum=float(a.min()),median=float(np.median(a)),mean=float(a.mean()),maximum=float(a.max())) if len(a) else dict(n=0,minimum=None,median=None,mean=None,maximum=None)
def count_nonfinite(x):
    if isinstance(x,dict):return sum(count_nonfinite(v) for v in x.values())
    if isinstance(x,list):return sum(count_nonfinite(v) for v in x)
    return int(not np.isfinite(x)) if isinstance(x,float) else 0
def outcome(r):return r['outlet'] if r['status']=='COMPLETED' else ('STATIONARY' if r['status']=='STATIONARY' else r['status'])

def certificate_leaves(c):
    if c['proof']=='CONTINUOUS_ORIGINAL_SUPPORT_PLANE_MINUS_H_LOWER':
        return [c]
    assert c['proof']=='UNION_OF_ORIGINAL_CONTINUOUS_HANDOFF_CERTIFICATES','Unknown accepted certificate schema'
    assert c['proof_partition_only'] and c['held_velocity_path_unchanged']
    children=c['subcertificates'];assert children and children[0]['t0_fraction']==0 and children[-1]['t1_fraction']==1
    assert all(s['t0_fraction']<s['t1_fraction'] for s in children)
    assert all(a['t1_fraction']==b['t0_fraction'] for a,b in zip(children[:-1],children[1:]))
    return [leaf for s in children for leaf in certificate_leaves(s['certificate'])]

def one(label,e,env):
    pid=e['particle_id'];stem=R/'outputs'/label/'trajectories'/f'mb_{pid:06d}'
    m=read(stem.with_suffix('.json'));a=np.load(stem.with_suffix('.npz'))['samples']
    with gzip.open(stem.with_suffix('.audit.jsonl.gz'),'rt') as f:obs=[json.loads(line) for line in f]
    accepted=[r for r in obs if r['accepted']];solves=[r for r in accepted if 'solver' in r]
    stationary='STATIONARY' in str(m['failure_detail']) or 'ROUNDOFF_SCALE_STAGNATION' in str(m['failure_detail'])
    last=solves[-32:]
    supported=stationary and len(last)==32 and all(r['solver'].get('contact_redundancy',{}).get('rank_after')==3
        and max(abs(np.array(r['velocity'][:3])))<=r['solver'].get('contact_kkt',{}).get('velocity_budget_m_s',0)
        and np.all(np.asarray(r['solver'].get('multipliers',[-1]))>=0) for r in last)
    outlet=None;classification_verified=False;inlet_escape=None
    if m['completed']:
        status='COMPLETED';hit=env.classifier.first_event(a[-2,1:4],a[-1,1:4]);classification_verified=hit is not None and hit.role==m['exit_outlet']
        outlet='O'+str(int(m['exit_outlet'].split('_')[-1]))
    elif stationary:status='STATIONARY'
    elif m['end_reason']=='PHYSICAL_RESIDENCE_HORIZON_REACHED':status='TIME_LIMIT'
    else:
        status='SOLVER_FAIL'
        if 'CENTER_OUTSIDE' in str(m['failure_detail']):
            inlet=ValidationBoundaryClassifier({'INLET':env.boundaries['INLET']})
            hits=[(i,inlet.first_event(x,y)) for i,(x,y) in enumerate(zip(a[:-1,1:4],a[1:,1:4]))]
            hits=[(i,h) for i,h in hits if h is not None]
            if hits:
                status='OTHER';outlet='INLET_ESCAPE';i,h=hits[-1]
                inlet_escape=dict(segment_index=i,segment_start=a[i,1:4],segment_end=a[i+1,1:4],hit_position=h.position_m,
                    initial_position=a[0,1:4],initial_velocity=a[0,4:7],first_diagnostic=obs[0],last_diagnostics=obs[-5:])
    tol=env.wall.roundoff_m
    penetration=int(np.count_nonzero(a[:,14]<-tol));handoff=int(np.count_nonzero(a[:,15]<-tol))
    topcerts=[c for r in accepted for c in r['continuous_certificates']]
    certs=[leaf for c in topcerts for leaf in certificate_leaves(c)]
    certificate_violations=sum(c['minimum_g_nf_bound_m'] < -c['roundoff_m'] for c in certs)
    nonfinite=int((~np.isfinite(a)).sum());diagnostic_nonfinite=sum(count_nonfinite(r) for r in accepted)
    steps=a[1:,18];its=[r['solver']['contact_kkt']['iterations'] for r in solves if 'contact_kkt' in r['solver']]
    contacts=[r['solver'].get('contact_count',0) for r in solves]
    contact_onsets=sum(v>0 and (i==0 or contacts[i-1]==0) for i,v in enumerate(contacts))
    corrections=[]
    for r in solves:
        for d in r.get('planar',[]):
            n=np.array(d['wall_normal_xyz']);bulk=np.array(d['bulk_velocity_xyz']);target=np.array(d['target_tangential_velocity_xyz'])
            corrections.append(float(d['wall_weight']*np.linalg.norm(target-(bulk-n*np.dot(n,bulk)))))
    final_sample=env.field.sample(a[-1,1:4]);finalspeed=float(np.linalg.norm(final_sample.velocity_m_s)) if final_sample.inside_lumen else None
    row=dict(bubble_id=pid,flow=label,diameter_um=e['diameter_um'],radius_um=e['radius_m']*1e6,
        initial_x_m=float(a[0,1]),initial_y_m=float(a[0,2]),initial_z_m=float(a[0,3]),status=status,outlet=outlet,
        transit_time_s=m['residence_time_s'] if m['completed'] else None,physical_simulated_time_s=float(a[-1,0]),
        absolute_acquisition_final_time_s=m['last_physical_time_s'],path_length_m=float(np.linalg.norm(np.diff(a[:,1:4],axis=0),axis=1).sum()),
        minimum_wall_gap_m=float(a[:,14].min()),minimum_h_over_a=float(a[:,14].min()/e['radius_m']),minimum_g_nf_m=float(a[:,15].min()),
        maximum_lubrication_activation=max([b['w'] for r in solves for b in r.get('interactions',[])],default=0),
        max_planar_wall_tangential_target_correction_m_s=max(corrections,default=0),
        max_unconstrained_hydrodynamic_velocity_correction_m_s=max([float(np.linalg.norm(np.array(r['hydro'][:3])-np.array(r['free'][:3]))) for r in solves],default=0),
        handoff_events=sum(r.get('handoff_event') is not None for r in accepted),contact_events=contact_onsets,
        contact_accepted_steps=sum(c>0 for c in contacts),max_simultaneous_contacts=max(contacts,default=0),
        contact_active_set_iterations=stats(its),minimum_timestep_s=float(steps.min()),maximum_timestep_s=float(steps.max()),
        nearwall_exposure_h_over_a_le_0p1_s=float(np.sum(steps[a[:-1,14]/e['radius_m']<=.1])),
        final_x_m=float(a[-1,1]),final_y_m=float(a[-1,2]),final_z_m=float(a[-1,3]),
        penetration_count=penetration,handoff_violation_count=handoff,continuous_certificate_violations=certificate_violations,
        nan_inf_count=nonfinite,accepted_diagnostic_nan_inf_count=diagnostic_nonfinite,
        unclassified_solver_corruption=int(status=='OTHER' and outlet!='INLET_ESCAPE')+int(m['completed'] and not classification_verified),
        outlet_classification_verified=classification_verified,stationary_three_contact_supported=bool(supported),
        local_final_fluid_speed_m_s=finalspeed,accepted_steps=m['accepted_steps'],rejected_trials=m['rejected_trials'],
        maximum_refinement_depth=m['maximum_refinement_depth'],end_reason=m['end_reason'],failure_detail=m['failure_detail'],wall_seconds=m['wall_seconds'])
    audit=dict(bubble_id=pid,flow=label,status=status,stationary=stationary,stationary_supported=bool(supported),
        position_m=a[-1,1:4],diameter_um=e['diameter_um'],local_fluid_speed_m_s=finalspeed,wall_gap_m=float(a[-1,14]),
        last_accepted_solve=last[-1] if last else None,last_32_contact_certificate_pass=bool(supported),
        failure_detail=m['failure_detail'],inlet_escape=inlet_escape,
        rejection_reasons=dict(Counter(str(r.get('error')) for r in obs if not r['accepted'])),
        continuous_certificate_count=len(certs),continuous_top_level_certificate_count=len(topcerts),continuous_certificate_violations=certificate_violations,
        scope='Saved states and original diagnostics only; no time advancement or radius experiment')
    return row,audit

def main():
    start=time.time();events=read(R/'data/paired_events.json');rows={};audits=[]
    for label in ['OLD','NEW']:
        env=make_environment(R/'server_bundle',label);rr=[]
        for e in events:
            row,ad=one(label,e,env);rr.append(row);audits.append(ad)
        rows[label]={r['bubble_id']:r for r in rr};csvwrite(R/'data'/f'{label.lower()}_paired_metrics.csv',rr)
        dump(R/'data'/f'{label.lower()}_paired_metrics.json',rr)
    transitions=[];pairedtransit=[]
    for e in events:
        pid=e['particle_id'];o,n=rows['OLD'][pid],rows['NEW'][pid]
        change=o['status']==n['status']=='COMPLETED' and o['outlet']!=n['outlet']
        transitions.append(dict(bubble_id=pid,diameter_um=e['diameter_um'],old_status=o['status'],new_status=n['status'],
            old_outlet=o['outlet'],new_outlet=n['outlet'],old_outcome=outcome(o),new_outcome=outcome(n),outlet_changed=change,
            old_transit_time_s=o['transit_time_s'],new_transit_time_s=n['transit_time_s'],old_min_gap_m=o['minimum_wall_gap_m'],new_min_gap_m=n['minimum_wall_gap_m']))
        if o['status']==n['status']=='COMPLETED':pairedtransit.append(dict(bubble_id=pid,old_s=o['transit_time_s'],new_s=n['transit_time_s'],
            delta_s=n['transit_time_s']-o['transit_time_s'],ratio=n['transit_time_s']/o['transit_time_s']))
    csvwrite(R/'data/paired_outcome_transition.csv',transitions);dump(R/'data/paired_outcome_transition.json',transitions)
    dump(R/'data/paired_transit_time.json',pairedtransit);dump(R/'data/saved_state_safety_audit.json',audits)
    pointrows=[]
    for e in events:
        p=read(R/'outputs/NEW_POINT'/f'point_{e["particle_id"]:06d}.json');n=rows['NEW'][e['particle_id']]
        assert p['initial_center_m']==e['birth_center_m']
        po='O'+str(int(p['outlet'].split('_')[-1])) if p['outlet'] else p['end_reason']
        pointrows.append(dict(bubble_id=e['particle_id'],point_outcome=po,mb_outcome=outcome(n),
            outcome_differs=po!=outcome(n),both_exited_route_differs=n['status']=='COMPLETED' and po in ['O1','O2','O3'] and po!=n['outlet'],point_end_reason=p['end_reason']))
    csvwrite(R/'data/point_vs_mb.csv',pointrows);dump(R/'data/point_vs_mb.json',pointrows)
    flow=read(R/'data/old_new_flow_contract.json');flux=read(R/'data/reintegrated_flow_flux.json');cohort=read(R/'data/cohort_provenance.json');limit=read(R/'data/known_limitation_p1_network_bc_local_conservation_baseline.json')
    old=list(rows['OLD'].values());new=list(rows['NEW'].values());so=[r['bubble_id'] for r in old if r['status']=='STATIONARY'];sn=[r['bubble_id'] for r in new if r['status']=='STATIONARY']
    counts={label:dict(Counter(r['status'] for r in rr.values())) for label,rr in rows.items()}
    matrix=dict(Counter(t['old_outcome']+' -> '+t['new_outcome'] for t in transitions));pm=dict(Counter(p['point_outcome']+' -> '+p['mb_outcome'] for p in pointrows))
    numerical_fail=any(r['status']=='SOLVER_FAIL' or any(r[k] for k in ['penetration_count','handoff_violation_count','continuous_certificate_violations','nan_inf_count','accepted_diagnostic_nan_inf_count','unclassified_solver_corruption']) for r in new)
    review=any(r['status'] in ['OTHER','TIME_LIMIT'] or (r['status']=='STATIONARY' and not r['stationary_three_contact_supported']) for r in new)
    status='NETWORK_FLOW_MB_VALIDATION_NUMERICAL_FAIL' if numerical_fail else ('NETWORK_FLOW_MB_VALIDATION_PHYSICS_REVIEW_REQUIRED' if review else ('NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY' if sn else 'NETWORK_FLOW_MB_VALIDATION_PASS'))
    summary=dict(old_flow_path=flow['paths']['OLD'],old_flow_sha=flow['sha256']['OLD'],new_flow_path=flow['paths']['NEW'],new_flow_sha=flow['sha256']['NEW'],
        geometry_identical=flow['geometry_identical'],geometry_bitwise_identical=flow['geometry_bitwise_identical'],old_fluid_split=flux['OLD']['split'],new_fluid_split=flux['NEW']['split'],
        cohort_source=cohort['source'],cohort_manifest_sha=cohort['manifest_sha256'],N=30,old_reused=True,
        old_completed=counts['OLD'].get('COMPLETED',0),old_stationary=len(so),old_fail=counts['OLD'].get('SOLVER_FAIL',0),
        new_completed=counts['NEW'].get('COMPLETED',0),new_stationary=len(sn),new_fail=counts['NEW'].get('SOLVER_FAIL',0),old_status_counts=counts['OLD'],new_status_counts=counts['NEW'],
        old_outlet_counts={o:sum(r['status']=='COMPLETED' and r['outlet']==o for r in old) for o in ['O1','O2','O3']},
        new_outlet_counts={o:sum(r['status']=='COMPLETED' and r['outlet']==o for r in new) for o in ['O1','O2','O3']},
        outlet_changed_count=sum(t['outlet_changed'] for t in transitions),outlet_changed_ids=[t['bubble_id'] for t in transitions if t['outlet_changed']],
        outcome_changed_count=sum(t['old_outcome']!=t['new_outcome'] for t in transitions),transition_matrix=matrix,
        stationary_old_ids=so,stationary_new_ids=sn,stationary_resolved_ids=[i for i in so if rows['NEW'][i]['status']=='COMPLETED'],new_stationary_ids=sorted(set(sn)-set(so)),
        mean_transit_time_old=stats([r['transit_time_s'] for r in old if r['status']=='COMPLETED'])['mean'],
        mean_transit_time_new=stats([r['transit_time_s'] for r in new if r['status']=='COMPLETED'])['mean'],
        paired_transit_time_statistics=dict(delta_s=stats([r['delta_s'] for r in pairedtransit]),new_over_old=stats([r['ratio'] for r in pairedtransit]),
            shorter_count=sum(r['delta_s']<0 for r in pairedtransit),longer_count=sum(r['delta_s']>0 for r in pairedtransit)),
        min_gap_old_statistics=stats([r['minimum_wall_gap_m'] for r in old]),min_gap_new_statistics=stats([r['minimum_wall_gap_m'] for r in new]),
        min_h_over_a_old_statistics=stats([r['minimum_h_over_a'] for r in old]),min_h_over_a_new_statistics=stats([r['minimum_h_over_a'] for r in new]),
        nearwall_exposure_old_s=stats([r['nearwall_exposure_h_over_a_le_0p1_s'] for r in old]),nearwall_exposure_new_s=stats([r['nearwall_exposure_h_over_a_le_0p1_s'] for r in new]),
        penetration_count_new=sum(r['penetration_count'] for r in new),handoff_violation_count_new=sum(r['handoff_violation_count'] for r in new),
        continuous_certificate_violations_new=sum(r['continuous_certificate_violations'] for r in new),solver_fail_count_new=counts['NEW'].get('SOLVER_FAIL',0),
        nan_inf_count_new=sum(r['nan_inf_count'] for r in new),accepted_diagnostic_nan_inf_count_new=sum(r['accepted_diagnostic_nan_inf_count'] for r in new),
        unclassified_solver_corruption_new=sum(r['unclassified_solver_corruption'] for r in new),inlet_escape_ids_new=[r['bubble_id'] for r in new if r['outlet']=='INLET_ESCAPE'],
        point_vs_MB_transition=pm,point_vs_MB_outcome_difference_count=sum(r['outcome_differs'] for r in pointrows),point_vs_MB_both_exited_route_difference_count=sum(r['both_exited_route_differs'] for r in pointrows),
        new_flow_local_flux_max_error=limit['root_section_max_error'],new_flow_local_flux_rms_error=limit['root_section_RMS_error'],
        science_files_changed=None,taylor_hood_termination_status='TAYLOR_HOOD_VALIDATION_STOPPED_RESOURCE_COST',
        large_sample_recommendation='NOT_READY_FOR_LARGE_SAMPLE' if numerical_fail or review else 'READY_FOR_LARGE_SAMPLE_REVIEW',final_status=status,
        scientific_scope='N=30 paired validation only; not population flow split or physiological ground truth',
        analysis_elapsed_seconds=time.time()-start)
    dump(R/'data/analysis_summary.json',summary)
    print(json.dumps(summary,indent=2,default=safe),flush=True)
if __name__=='__main__':main()
