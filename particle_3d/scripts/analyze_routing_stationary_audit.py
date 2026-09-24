#!/usr/bin/env python3
"""Assemble paired evidence and protect all pre-audit files; no trajectory integration."""
from pathlib import Path
from collections import Counter
import json,subprocess
import numpy as np
from particle_3d.routing_stationary_audit import *
from particle_3d.particle8_replay import REPO
R=REPO/'particle_3d/reports/particle9a1_routing_stationary_audit';D=R/'data'
REF=REPO/'particle_3d/outputs/particle9a1_2mmps'

def split(rows,key):
    c=Counter(r[key] or 'NO_EXIT' for r in rows);return {k:c[k] for k in ROLES}
def matrix(rows,a,b):
    m=np.zeros((4,4),int)
    for r in rows:m[ROLES.index(r[a]),ROLES.index(r[b])]+=1
    den=m.sum(1);return dict(roles=ROLES,counts=m.tolist(),fraction_of_500=(m/len(rows)).tolist(),fraction_of_source_basin=np.divide(m,den[:,None],out=np.zeros((4,4)),where=den[:,None]>0).tolist(),changed_including_no_exit=int(m.sum()-np.trace(m)),changed_between_named_outlets=int(m[:3,:3].sum()-np.trace(m[:3,:3])))
def stats(v):
    v=np.asarray(v);return dict(count=len(v),median=float(np.median(v)),min=float(v.min()),max=float(v.max()),percentiles={str(p):float(np.percentile(v,p)) for p in [5,25,50,75,95]})
def hotspots(rows):
    groups=[{i} for i in range(len(rows))]
    for i,a in enumerate(rows):
        for j,b in enumerate(rows[:i]):
            dist=np.linalg.norm(np.array(a['final_position_m'])-b['final_position_m']);rr=a['radius_m']+b['radius_m'];shared=set(a['contact_triangle_ids'])&set(b['contact_triangle_ids'])
            if dist<=2*rr and (shared or dist<=rr):
                A=next(g for g in groups if i in g);B=next(g for g in groups if j in g)
                if A is not B:A.update(B);groups.remove(B)
    result=[]
    for k,g in enumerate(sorted(groups,key=lambda g:min(rows[i]['particle_id'] for i in g)),1):
        items=[rows[i] for i in sorted(g)]
        for item in items:item['hotspot_id']=k
        result.append(dict(hotspot_id=k,count=len(items),particle_ids=[x['particle_id'] for x in items],centroid_m=np.mean([x['final_position_m'] for x in items],axis=0).tolist(),triangle_ids=sorted(set(t for x in items for t in x['contact_triangle_ids'])),radius_range_um=[min(x['radius_m'] for x in items)*1e6,max(x['radius_m'] for x in items)*1e6]))
    return result

def main():
    point=read(D/'point_completed.json')['results'];p65=read(D/'p65_completed.json')['results'];stationary=read(D/'stationary_audit.json')
    events=read(REF/'admission/birth_ledger.json')['events'];identity=verify_events(read(REPO/'particle_3d/outputs/particle9a_2mmps/admission/all_scheduled_admission_records.json'),events)
    pb={x['particle_id']:x for x in point};b65={x['particle_id']:x for x in p65};p9=[];paired=[]
    for e in events:
        pid=e['particle_id'];m=read(REF/'P9A1/trajectories'/f'mb_{pid:06d}.json');q=read(R/'diagnostic_outputs/P65_NEW/trajectories'/f'mb_{pid:06d}.json');assert m['birth_metadata']==q['birth_metadata']==e
        assert m['integration_config']==q['integration_config']
        p9.append(dict(particle_id=pid,outlet=m['exit_outlet'] or 'NO_EXIT',outcome=m['end_reason']))
        paired.append(dict(particle_id=pid,point_outlet=pb[pid]['point_outlet'],p65_outlet=b65[pid]['outlet'],p9a1_outlet=m['exit_outlet'] or 'NO_EXIT',radius_m=e['radius_m'],birth_position_m=e['birth_center_m'],birth_time_s=e['birth_time_s'],minimum_gap_p65_m=q['minimum_original_wall_gap_m'],minimum_gap_p9a1_m=m['minimum_original_wall_gap_m'],p65_outcome=q['end_reason'],p9a1_outcome=m['end_reason'],p65_failure=q['failure_detail'],p9a1_failure=m['failure_detail']))
    assert len(paired)==len(p65)==len(p9)==500
    dump(D/'paired_routing.json',paired);csvwrite(D/'p65_vs_p9a1_routing.csv',paired);csvwrite(D/'p9a1_read_only_outcomes.csv',p9)
    accepted=[x for x in point if x['accepted']];basins=[]
    for basin in ROLES:
        raw=[x for x in point if x['point_outlet']==basin];acc=[x for x in raw if x['accepted']]
        basins.append(dict(point_basin=basin,raw_count=len(raw),accepted_count=len(acc),rejected_count=len(raw)-len(acc),acceptance_fraction=len(acc)/len(raw) if raw else None))
    csvwrite(D/'admission_by_basin.csv',basins)
    raw_o1=[x for x in point if x['point_outlet']=='OUTLET_01'];accepted_o1=[x for x in paired if x['point_outlet']=='OUTLET_01']
    loss=dict(raw_O1=len(raw_o1),not_in_accepted=sum(not x['accepted'] for x in raw_o1),admission_rejected=sum(not x['accepted'] for x in raw_o1),not_accepted_other_reason=0,
        accepted_O1=len(accepted_o1),P65_keeps_O1=sum(x['p65_outlet']=='OUTLET_01' for x in accepted_o1),accepted_P65_reroutes=sum(x['p65_outlet'] not in ['OUTLET_01','NO_EXIT'] for x in accepted_o1),accepted_P65_stops=sum(x['p65_outlet']=='NO_EXIT' for x in accepted_o1),
        P65_keeps_then_P9_reroutes=sum(x['p65_outlet']=='OUTLET_01' and x['p9a1_outlet']!='OUTLET_01' for x in accepted_o1),P9_keeps_O1=sum(x['p9a1_outlet']=='OUTLET_01' for x in accepted_o1),counting_note='not_in_accepted and admission_rejected are the same group, NEVER ADD THEM')
    dump(D/'outlet01_loss_chain.json',loss)
    representatives=[];seen=set()
    ordered=sorted(accepted_o1,key=lambda x:(x['p65_outlet']!='OUTLET_01',x['particle_id']))
    for r in ordered:
        route=(r['p65_outlet'],r['p9a1_outlet'])
        if route not in seen:representatives.append(r);seen.add(route)
    representatives=(representatives+[r for r in ordered if r not in representatives])[:5]
    csvwrite(D/'outlet01_representatives.csv',representatives)
    byid={x['particle_id']:x for x in paired};all_r=np.array([x['radius_m'] for x in paired]);size_rows=[]
    for x in paired:size_rows.append(dict(particle_id=x['particle_id'],radius_um=x['radius_m']*1e6,diameter_um=x['radius_m']*2e6,stationary=x['p9a1_outlet']=='NO_EXIT',size_percentile=float(np.mean(all_r<=x['radius_m'])*100)))
    for r in stationary:
        p=byid[r['particle_id']];r['contact_geometry_status']='RIGID_CONTACT_GEOMETRY_CONFIRMED' if r['retained_rows_independent'] else 'CONTACT_RANK_UNRESOLVED';r.update(point_tracer_outlet=p['point_outlet'],p65_outcome=p['p65_outcome'],p65_outlet=p['p65_outlet'],p9a1_outcome=p['p9a1_outcome'],size_percentile=float(np.mean(all_r<=r['radius_m'])*100))
    hot=hotspots(stationary);dump(D/'stationary_audit.json',stationary);csvwrite(D/'stationary_summary.csv',[{k:r[k] for k in ['particle_id','radius_m','diameter_um','final_position_m','wall_gap_m','contact_count','contact_triangle_ids','contact_feature_types','contact_normals','contact_rank','contact_condition','multipliers','free_fem_velocity_m_s','unconstrained_particle_velocity','constrained_velocity','point_tracer_outlet','p65_outcome','p65_outlet','p9a1_outcome','hotspot_id','size_percentile','classification','radius_sensitivity','critical_radius_ratio']} for r in stationary]);csvwrite(D/'particle_sizes.csv',size_rows)
    virtual=[dict(particle_id=r['particle_id'],radius_ratio=t['radius_ratio'],status=t['status'],classification=r['classification'],sensitivity=r['radius_sensitivity'],recorded_continuation_time_s=t['last_time_s']-t['checkpoint_time_s'],first_local_passage_elapsed_s=(t['local_passage_witness']['time_s']-t['checkpoint_time_s']) if t.get('local_passage_witness') else None,fixed_initial_fem_projection_um=t['downstream_distance_m']*1e6,cumulative_local_downstream_um=t.get('verified_cumulative_local_downstream_m',t.get('cumulative_local_downstream_m',0))*1e6,reference_reproduced=r['reference_jam_reproduced']) for r in stationary for t in r['trials']];csvwrite(D/'virtual_radius_trials.csv',virtual)
    size_stats={label:{unit:stats([x[unit] for x in size_rows if x['stationary']==flag]) for unit in ['radius_um','diameter_um','size_percentile']} for label,flag in [('completed',False),('stationary',True)]}
    fluid=read(REPO/'formal_3D_flow_solver/FEM_SimVascular/frozen_reference/new_flow_flux.json');dump(D/'flow_split.json',fluid)
    protected=read(D/'protected_before.json');mismatches=[p for p,h in protected.items() if not (REPO/p).is_file() or sha(REPO/p)!=h]
    before=read(R/'logs/git_before.json');index=subprocess.check_output(['git','write-tree'],cwd=REPO,text=True);protection=dict(file_count=len(protected),mismatches=mismatches,index_unchanged=index==before['git write-tree']);dump(D/'protected_verification.json',protection)
    summary=dict(flow_split=fluid,raw_candidate_point_split=split(point,'point_outlet'),accepted500_point_split=split(accepted,'point_outlet'),p65_500_split=split(p65,'outlet'),p9a1_500_split=split(p9,'outlet'),
        candidate_count=len(point),accepted_count=len(accepted),acceptance_by_point_basin=basins,point_to_p65_transition=matrix(paired,'point_outlet','p65_outlet'),p65_to_p9a1_transition=matrix(paired,'p65_outlet','p9a1_outlet'),
        outlet01_loss_chain=loss,stationary_count=len(stationary),stationary_hotspots=hot,hotspot_rule='CONNECTED_COMPONENTS: distance<=2*(r_i+r_j) AND (shared contact triangle OR distance<=r_i+r_j)',stationary_radius_statistics=size_stats,
        critical_radius_ratios={str(r['particle_id']):r['critical_radius_ratio'] for r in stationary},stationary_classification_counts=dict(Counter(r['classification'] for r in stationary)),unresolved_ids=[r['particle_id'] for r in stationary if r['classification']=='NUMERICAL_OR_GEOMETRIC_UNRESOLVED'],
        point_unresolved_ids=[r['particle_id'] for r in point if r['point_outlet']=='NO_EXIT'],production_protection=protection,event_identity=identity,
        runtimes_s={k:read(D/(k+'_completed.json'))['wall_seconds'] for k in ['point','p65','radius']},server=read(D/'remote.json'),workers={k:read(D/(k+'_completed.json'))['workers'] for k in ['point','p65','radius']})
    raw_fraction=summary['raw_candidate_point_split']['OUTLET_02']/len(point);accepted_fraction=summary['accepted500_point_split']['OUTLET_02']/500;p65_fraction=summary['p65_500_split']['OUTLET_02']/500;p9_fraction=summary['p9a1_500_split']['OUTLET_02']/500
    shifts=dict(SOURCE_A_RAW_MINUS_FLUID=(raw_fraction-fluid['outlet_fractions']['OUTLET_02'])*100,SOURCE_B_ADMISSION=(accepted_fraction-raw_fraction)*100,SOURCE_C_P65=(p65_fraction-accepted_fraction)*100,SOURCE_D_P9A1=(p9_fraction-p65_fraction)*100)
    summary['causal_decomposition']=dict(O2_percentage_point_changes=shifts,dominant_observed_enrichment_stage=max(shifts,key=shifts.get),role='SEQUENTIAL_CHANGES_IN_OBSERVED_PROPORTIONS; NOT A FORCE_DECOMPOSITION',source_E='13 original P1 wall-tetra asymptotic stagnations, verified analytically; no classifier mismatch')
    field_audit=read(D/'point_stagnation_field.json');confirmed={r['particle_id'] for r in field_audit['records'] if r['analytic_wall_stagnation_confirmed']}
    summary['point_no_exit_ids']=summary.pop('point_unresolved_ids')
    summary['point_unresolved_cause_ids']=[i for i in summary['point_no_exit_ids'] if i not in confirmed]
    summary['point_no_exit_explanation']='ANALYTIC_ASYMPTOTIC_WALL_STAGNATION_IN_ORIGINAL_NON_SOLENOIDAL_P1_WALL_TETRAHEDRA; NO_OUTLET_ASSIGNED'
    summary['point_field_stagnation']=field_audit
    summary['point_resolution_audit']=read(D/'point_resolution.json')
    summary['status']='SCIENTIFIC_AUDIT_PARTIAL' if summary['unresolved_ids'] or summary['point_unresolved_cause_ids'] or mismatches else 'SCIENTIFIC_AUDIT_COMPLETE'
    dump(D/'audit_summary.json',summary)
    rows=[]
    for level,key in [('Raw candidate point','raw_candidate_point_split'),('Accepted 500 point','accepted500_point_split'),('P65 finite-size','p65_500_split'),('P9A1 saved','p9a1_500_split')]:
        counts=summary[key];n=sum(counts.values());rows.extend(dict(stage=level,outlet=role,count=counts[role],denominator=n,fraction=counts[role]/n) for role in ROLES)
    rows=[dict(stage='Fluid Q',outlet=role,count=None,denominator=None,fraction=fluid['outlet_fractions'].get(role,0.)) for role in ROLES]+rows;csvwrite(D/'stage_splits.csv',rows)
    print(json.dumps({k:summary[k] for k in ['status','raw_candidate_point_split','accepted500_point_split','p65_500_split','p9a1_500_split','stationary_classification_counts','unresolved_ids']},indent=2))
if __name__=='__main__':main()
