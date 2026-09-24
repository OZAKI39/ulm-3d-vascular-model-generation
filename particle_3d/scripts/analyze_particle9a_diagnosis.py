#!/usr/bin/env python3
"""Read-only numerical analysis of bounded diagnostic replays; no integration."""
from collections import Counter
import csv, gzip, json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/'particle_3d/reports/particle9a_2mmps_diagnosis'
DATA=REPORT/'data'
_INLET_CLASSIFIER=None

def dump(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def readrows(path):
    with gzip.open(path,'rt') as f:return [json.loads(line) for line in f]

def csvwrite(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for row in rows:
            w.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(dict,list)) else v for k,v in row.items()})

def outcome(meta,arr):
    if meta['completed']:return 'COMPLETED'
    reason=str(meta.get('failure_detail'))
    if reason=='CENTER_OUTSIDE_FROZEN_LUMEN':
        # A global fitted cap plane is not exact at sub-picometre displacement.
        # Use finite original triangles, never the fitted-plane sign alone.
        hit=_INLET_CLASSIFIER.first_event(arr[0,1:4],arr[-1,1:4])
        if hit is not None and hit.role=='INLET':return 'INLET_ESCAPE'
        return 'OUTSIDE_OTHER_UNRESOLVED'
    if 'SUBDIVISION_PROGRESS' in reason:return 'SUBDIVISION_STOP'
    if meta['end_reason']=='PHYSICAL_RESIDENCE_HORIZON_REACHED':return 'TIME_LIMIT'
    return 'OTHER_UNRESOLVED'

def load(model,pid,base='replays'):
    stem=DATA/base/model/'trajectories'/f'mb_{pid:06d}'
    return (json.loads(stem.with_suffix('.json').read_text()),np.load(stem.with_suffix('.npz'))['samples'],
            readrows(stem.with_suffix('.trials.jsonl.gz')),readrows(stem.with_suffix('.states.jsonl.gz')))

def trial_stats(rows):
    solved=[r for r in rows if 'resistance_condition_estimate' in r]
    accepted=[r for r in rows if r['trial_accepted']]
    conditions=[r['resistance_condition_estimate'] for r in solved]
    active=[r['handoff_active'] for r in solved]
    kkt=[r['solver_record']['contact_kkt'] for r in solved if 'contact_kkt' in r['solver_record']]
    angles=[r['normal_change_angle_from_previous_deg'] for r in rows]
    # 10/45/90 deg are reported descriptive bins, never model/acceptance gates.
    activity={}
    for flag in [False,True]:
        selected=[r for r in solved if r['handoff_active']==flag]
        if selected:
            activity[str(flag)]=dict(count=len(selected),rejected=sum(not r['trial_accepted'] for r in selected),
                min_start_g_nf_m=min(r['g_nf_m'] for r in selected),max_start_g_nf_m=max(r['g_nf_m'] for r in selected),
                normal_speed_after_min_m_s=min(r['normal_velocity_after_constraint'] for r in selected),
                normal_speed_after_max_m_s=max(r['normal_velocity_after_constraint'] for r in selected),
                feature_counts=dict(Counter(r['nearest_wall_feature'] for r in selected)))
    return dict(trials=len(rows),accepted=len(accepted),rejections=dict(Counter(r['trial_rejection_reason'] for r in rows if not r['trial_accepted'])),
        handoff_activity_groups=activity,
        maximum_contact_kkt_condition=max((r['condition'] for r in kkt),default=None),
        maximum_contact_velocity_budget_m_s=max((r['velocity_budget_m_s'] for r in kkt),default=None),
        maximum_scaled_resistance_backward_residual=max((r['solver_record']['relative_backward_residual'] for r in solved),default=None),
        minimum_scaled_resistance_eigenvalue=min((r['solver_record']['minimum_scaled_eigenvalue'] for r in solved),default=None),
        triangle_switch_count=sum(r['wall_triangle_changed'] for r in rows),maximum_normal_change_deg=max(angles,default=0.),
        normal_angle_counts_descriptive_only={str(a):sum(v>=a for v in angles) for a in (10,45,90)},
        handoff_active_trials=sum(active),handoff_activation_deactivation_count=sum(a!=b for a,b in zip(active,active[1:])),
        minimum_start_g_nf_m=min((r['g_nf_m'] for r in rows),default=None),
        minimum_rejected_certificate_g_nf_m=min((r.get('handoff_certificate',{}).get('minimum_g_nf_m',np.inf) for r in rows),default=np.inf),
        condition_min=min(conditions,default=None),condition_max=max(conditions,default=None),
        max_depth=max((r['subdivision_depth'] for r in rows),default=0),
        accepted_physical_time_s=sum(r['accepted_dt_s'] for r in accepted),
        nominal_requested_dt_s=rows[-1]['requested_dt_s'] if rows else None,
        accepted_time_over_one_nominal_dt=sum(r['accepted_dt_s'] for r in accepted)/rows[-1]['requested_dt_s'] if rows else None,
        min_normal_velocity_before_constraint=min((r['normal_velocity_before_constraint'] for r in solved),default=None),
        max_abs_normal_velocity_after_constraint=max((abs(r['normal_velocity_after_constraint']) for r in solved),default=None))

def main():
    global _INLET_CLASSIFIER
    from particle_3d.particle81_simulation import environment
    from particle_3d.validation_boundary import ValidationBoundaryClassifier
    env=environment();_INLET_CLASSIFIER=ValidationBoundaryClassifier({'INLET':env.boundaries['INLET']})
    audit=json.loads((DATA/'original_smoke_audit.json').read_text()); originals={r['particle_id']:r for r in audit['rows']}
    ab=[];cache={};parity=[];inlets=[];stalls=[];closure=[];timecheck=[];inlet_summary=[];stall_summary=[];divergence=[];contact_anomalies=[]
    checkpoints={}
    for pid in originals:
        for model in ['P65','P9A']:
            m,a,t,s=load(model,pid); cache[(pid,model)]=(m,a,t,s)
            ab.append(dict(particle_id=pid,model=model,birth_time_s=m['birth_time_s'],radius_m=m['radius_m'],
                birth_metadata_sha256=m['birth_metadata_sha256'],outcome=outcome(m,a),outlet=m['exit_outlet'],
                travel_time_s=m['residence_time_s'],minimum_gap_m=m['minimum_original_wall_gap_m'],
                minimum_accepted_g_nf_m=float(a[:,15].min()),accepted_substeps=m['accepted_steps'],
                rejected_trials=m['rejected_trials'],max_subdivision_depth=max(r['subdivision_depth'] for r in t),
                final_time_s=m['last_elapsed_time_s'],final_position_m=a[-1,1:4].tolist(),
                end_reason=m['end_reason'],failure_detail=m['failure_detail']))
            if model=='P9A':
                old=np.load(ROOT/f'particle_3d/outputs/particle9a_2mmps/trajectories/mb_{pid:06d}.npz')['samples']
                parity.append(dict(particle_id=pid,same_shape=a.shape==old.shape,
                    bitwise_all_samples=bool(np.array_equal(a,old)),max_abs_difference=float(abs(a-old).max()) if a.shape==old.shape else None))
                accepted=[r for r in t if r['trial_accepted']]
                for i in np.unique(np.linspace(0,len(accepted)-1,min(25,len(accepted))).astype(int)):
                    r=dict(accepted[i],outcome=originals[pid]['outcome'],selection='EVEN_ACCEPTED_INDEX_INCLUDING_FIRST_LAST',
                        state_role='ACCEPTED_TRIAL_START; EXACT_LOCATION_OF_SAVED_SOLVE_AND_CLOSURE')
                    closure.append(r)
        assert cache[(pid,'P65')][0]['birth_metadata_sha256']==cache[(pid,'P9A')][0]['birth_metadata_sha256']
        assert cache[(pid,'P65')][0]['integration_config']==cache[(pid,'P9A')][0]['integration_config']
        assert np.array_equal(cache[(pid,'P65')][1][0],cache[(pid,'P9A')][1][0])
        orig=originals[pid]['outcome']
        if orig=='INLET_ESCAPE':
            p65=cache[(pid,'P65')];p9=cache[(pid,'P9A')]
            cap_id=p9[0]['birth_metadata']['anchor_triangle'];tri=_INLET_CLASSIFIER.xyz[cap_id]
            n=_INLET_CLASSIFIER.normal[cap_id].copy()
            if n@np.array(p9[0]['inlet_plane']['inward_normal'])<0:n*=-1
            def exact_plane(row):
                r=dict(row,best_fit_inlet_plane_distance_m=row['distance_to_inlet_plane_m'],
                    distance_to_inlet_plane_m=float((np.array(row['position_xyz'])-tri[0])@n),
                    inlet_distance_role='ORIGINAL_ADMISSION_ANCHOR_CAP_TRIANGLE_PLANE',inlet_cap_triangle_id=cap_id)
                if 'velocity_xyz' in r:r['velocity_dot_inlet_inward_normal']=float(np.array(r['velocity_xyz'])@n)
                if 'FEM_velocity_xyz' in r:r['FEM_inward_velocity_m_s']=float(np.array(r['FEM_velocity_xyz'])@n)
                r['wall_inlet_normal_angle_deg']=float(np.degrees(np.arccos(np.clip(np.array(r['wall_normal_xyz'])@n,-1,1))))
                return r
            for model in ['P65','P9A']:
                m,a,t,s=cache[(pid,model)]
                inlets.extend(dict(exact_plane(r),record_kind='ACCEPTED_STATE') for r in s)
                inlets.extend(dict(exact_plane(r),record_kind='TRIAL_VELOCITY_EVALUATION') for r in t)
            a,b=exact_plane(p65[2][0]),exact_plane(p9[2][0])
            comparison_time=float(p9[1][1,0]);a_at_b=np.array([np.interp(comparison_time,p65[1][:,0],p65[1][:,j]) for j in [1,2,3]])
            hit=_INLET_CLASSIFIER.first_event(p9[1][0,1:4],p9[1][-1,1:4])
            record=dict(particle_id=pid,P65_outcome=outcome(p65[0],p65[1]),P9A_outcome=outcome(p9[0],p9[1]),
                first_velocity_divergence_evaluation_time_s=0.,first_position_divergence_time_s=float(p9[1][1,0]),
                first_velocity_difference_m_s=float(np.linalg.norm(np.array(a['velocity_xyz'])-b['velocity_xyz'])),
                first_position_difference_m=float(np.linalg.norm(a_at_b-p9[1][1,1:4])),
                position_comparison_role='P65_PIECEWISE_EULER_PATH_AT_FIRST_P9A_SAVED_ENDPOINT',
                cap_crossing_triangle_id=hit.role_triangle_id,cap_crossing_segment_fraction=hit.segment_fraction,
                inlet_cap_origin_m=tri[0].tolist(),inlet_cap_inward_normal=n.tolist(),
                P65_inward_velocity_m_s=a['velocity_dot_inlet_inward_normal'],
                P9A_inward_velocity_m_s=b['velocity_dot_inlet_inward_normal'],
                P9A_last_signed_distance_m=exact_plane(p9[3][-1])['distance_to_inlet_plane_m'],
                P65_first_signed_distance_m=exact_plane(p65[3][1])['distance_to_inlet_plane_m'])
            for k in ['FEM_inward_velocity_m_s','bulk_t_speed','linear_shear_speed','wall_weight','nearest_wall_feature',
                'gap_ratio','distance_to_inlet_plane_m','wall_triangle_touches_inlet_rim','wall_inlet_normal_angle_deg',
                'wall_point_distance_to_inlet_rim_m','target_tangential_speed','tangential_speed','closure_angle_deg',
                'closure_relative_mismatch','nearest_wall_triangle_id','wall_target_velocity_xyz','blended_target_velocity_xyz']:
                record[k]=b[k]
            record['blended_target_inward_m_s']=float(np.array(b['blended_target_velocity_xyz'])@n)
            inlet_summary.append(record)
            knots=np.unique(np.r_[p65[1][:,0],p9[1][:,0]])
            for time in knots[knots<=min(p65[1][-1,0],p9[1][-1,0])]:
                positions=[];velocities=[]
                for samples in [p65[1],p9[1]]:
                    positions.append(np.array([np.interp(time,samples[:,0],samples[:,j]) for j in [1,2,3]]))
                    idx=max(1,int(np.searchsorted(samples[:,0],time,side='left')))
                    velocities.append(samples[idx,4:7])
                divergence.append(dict(particle_id=pid,time_s=float(time),physical_time_s=float(time+p9[0]['birth_time_s']),
                    P65_position_m=positions[0].tolist(),P9A_position_m=positions[1].tolist(),
                    position_difference_m=float(np.linalg.norm(positions[0]-positions[1])),
                    velocity_difference_m_s=float(np.linalg.norm(velocities[0]-velocities[1])),
                    P65_inward_m_s=float(velocities[0]@n),P9A_inward_m_s=float(velocities[1]@n),
                    inward_sign_difference=bool((velocities[0]@n)*(velocities[1]@n)<0),
                    sampling='UNION_OF_ACCEPTED_TIMES; EXACT_LINEAR_EULER_PATH; LEFT_HELD_VELOCITY_EXCEPT_AT_BIRTH'))
        if orig=='SUBDIVISION_STOP':
            m,a,t,s=cache[(pid,'P9A')]
            window=t[-100:]
            valid_boundaries=[float(v) for v in a[:,0] if v>0 and v/.00025==round(v/.00025) and v<window[0]['time_s']]
            # Select an already saved nominal boundary preceding the window.
            checkpoint=valid_boundaries[-1] if valid_boundaries else 0.
            checkpoints[str(pid)]=checkpoint
            summary=dict(particle_id=pid,P65_outcome=outcome(*cache[(pid,'P65')][:2]),P9A_outcome=outcome(m,a),
                         checkpoint_time_s=checkpoint,P9A_last100=trial_stats(window),P9A_all=trial_stats(t),
                         minimum_accepted_g_nf_m=float(a[:,15].min()))
            # Exact first bisection and first depth >=7: latter means <1% nominal
            # per trial, the existing guard scale, NOT a new physics criterion.
            summary['first_rejected_trial_time_s']=next((r['time_s'] for r in t if not r['trial_accepted']),None)
            summary['first_depth_ge_7_time_s']=next((r['time_s'] for r in t if r['subdivision_depth']>=7),None)
            for row in window:
                sr=row.get('solver_record',{})
                if sr.get('contact_count',0)>1:
                    ids=[c[3] for c in sr['contact_ids'] if c[1]=='WALL']
                    shared=set(env.wall.global_node_ids[ids[0]])
                    for tri_id in ids[1:]:shared.intersection_update(env.wall.global_node_ids[tri_id])
                    contact_anomalies.append(dict(particle_id=pid,trial_index=row['trial_index'],time_s=row['time_s'],
                        contact_triangle_ids=ids,shared_global_node_ids_zero_based=sorted(map(int,shared)),
                        shared_finite_edge=len(shared)==2,solver_record=sr,trial_rejection_reason=row['trial_rejection_reason'],
                        handoff_certificate=row.get('handoff_certificate')))
            for model in ['P65','P9A']:
                stalls.extend(dict(r,window='LAST_100_TRIALS_FROM_BIRTH_'+model) for r in cache[(pid,model)][2][-100:])
            path=DATA/'checkpoint_replays/checkpoint/trajectories'/f'mb_{pid:06d}.json'
            if path.exists():
                cm,ca,ct,cs=load('checkpoint',pid,'checkpoint_replays')
                count=cm['checkpoint']['accepted_prefix_samples']
                equal=bool(np.array_equal(ca[:count],a[:count]))
                assert equal and ca[count-1,0]==checkpoint
                selected=[r for r in ct if r['time_s']>=checkpoint]
                summary.update(checkpoint_P65_outcome=outcome(cm,ca),checkpoint_prefix_bitwise_equal=equal,
                    checkpoint_prefix_samples=count,checkpoint_P65_final_time_s=float(ca[-1,0]),
                    checkpoint_P65_last100=trial_stats(selected[-100:]),
                    checkpoint_snapshot=cm['checkpoint'])
                stalls.extend(dict(r,window='LAST_100_TRIALS_P65_AFTER_SAME_P9A_CHECKPOINT') for r in selected[-100:])
            stall_summary.append(summary)
        if orig=='TIME_LIMIT':
            for model in ['P65','P9A']:
                m,a,t,s=cache[(pid,model)];dts=np.diff(a[:,0]);distance=np.linalg.norm(np.diff(a[:,1:4],axis=0),axis=1)
                duration=a[-1,0];accepted=[r for r in t if r['trial_accepted']]
                handoff_time=sum(r['accepted_dt_s'] for r in accepted if r.get('handoff_active',False))
                last=a[1:,0]>.9*duration
                timecheck.append(dict(particle_id=pid,model=model,outcome=outcome(m,a),elapsed_time_s=float(duration),
                    path_length_m=float(distance.sum()),mean_speed_m_s=float(distance.sum()/duration),
                    last_10percent_time_mean_speed_m_s=float(distance[last].sum()/dts[last].sum()),
                    minimum_gap_m=float(a[:,14].min()),maximum_gap_m=float(a[:,14].max()),
                    handoff_active_time_fraction=float(handoff_time/duration),
                    fully_weighted_wall_time_fraction=sum(r['accepted_dt_s'] for r in accepted if r.get('wall_weight')==1)/duration,
                    final_speed_m_s=float(np.linalg.norm(a[-1,4:7])),final_position_m=a[-1,1:4].tolist()))
    csvwrite(DATA/'ab_outcomes.csv',ab);csvwrite(DATA/'inlet_escape_states.csv',inlets)
    csvwrite(DATA/'handoff_stall_trials.csv',stalls);csvwrite(DATA/'closure_consistency_states.csv',closure)
    csvwrite(DATA/'time_limit_check.csv',timecheck);csvwrite(DATA/'inlet_summary.csv',inlet_summary)
    csvwrite(DATA/'inlet_divergence.csv',divergence)
    dump(DATA/'contact_solver_edge_audit.json',contact_anomalies)
    dump(DATA/'original_p9a_replay_parity.json',parity);dump(DATA/'checkpoint_times.json',checkpoints)
    # JSON supports no Infinity; absence is explicit and includes the reason.
    def finite(x):
        if isinstance(x,dict):return {k:finite(v) for k,v in x.items()}
        if isinstance(x,list):return [finite(v) for v in x]
        if isinstance(x,float) and not np.isfinite(x):return None
        return x
    closure_groups={}
    for group in dict.fromkeys(originals[pid]['outcome'] for pid in originals):
        rows=[r for r in closure if r['outcome']==group]
        closure_groups[group]=dict(n=len(rows),median_speed_ratio=float(np.median([r['closure_speed_ratio'] for r in rows])),
            min_speed_ratio=min(r['closure_speed_ratio'] for r in rows),max_speed_ratio=max(r['closure_speed_ratio'] for r in rows),
            median_mismatch=float(np.median([r['closure_relative_mismatch'] for r in rows])),
            max_angle_deg=max(r['closure_angle_deg'] or 0 for r in rows))
    summary=finite(dict(counts={m:dict(Counter(r['outcome'] for r in ab if r['model']==m)) for m in ['P65','P9A']},
        inlet=inlet_summary,stall=stall_summary,time_limit=timecheck,closure_groups=closure_groups,
        all_original_p9a_bitwise_equal=all(r['bitwise_all_samples'] for r in parity)))
    dump(DATA/'analysis_summary.json',summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
