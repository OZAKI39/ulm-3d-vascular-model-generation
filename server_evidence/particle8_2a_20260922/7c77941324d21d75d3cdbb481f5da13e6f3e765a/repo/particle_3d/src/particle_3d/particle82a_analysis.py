"""Quantitative paired-cohort analysis; unresolved events stay visible."""
from pathlib import Path
from collections import Counter
import argparse, csv, gzip, json
import numpy as np
from scipy.stats import kstest, ks_2samp, wasserstein_distance
from .particle82a_trajectories import iter_rows
from .particle82a_admission import context
from .particle82a_pipeline import BASINS, CAUSES, atomic_npz
from .particle82_provenance import atomic_json, sha256


def descriptive(values):
    a=np.asarray(values,float)
    return dict(n=len(a),mean=float(a.mean()),min=float(a.min()),p10=float(np.quantile(a,.1)),
        p25=float(np.quantile(a,.25)),median=float(np.median(a)),p75=float(np.quantile(a,.75)),
        p90=float(np.quantile(a,.9)),p95=float(np.quantile(a,.95)),max=float(a.max())) if len(a) else dict(n=0)


def distribution_distance(values, distribution):
    a=np.sort(np.asarray(values,float))
    if not len(a):return dict(n=0)
    nodes=np.unique(np.r_[a,distribution.low,distribution.high]);left=nodes[:-1];right=nodes[1:]
    empirical=np.searchsorted(a,(left+right)/2,side='right')/len(a)
    y0=empirical-distribution.cdf(left);y1=empirical-distribution.cdf(right)
    area=(np.abs(y0)+np.abs(y1))*.5*(right-left)
    crossing=y0*y1<0
    area[crossing]=(right[crossing]-left[crossing])*(y0[crossing]**2+y1[crossing]**2)/(2*np.abs(y1[crossing]-y0[crossing]))
    ks=kstest(a,distribution.cdf)
    return dict(n=len(a),KS_statistic=float(ks.statistic),KS_pvalue=float(ks.pvalue),
        Wasserstein_um=float(area.sum()),reference='EXACT_FROZEN_CONTINUOUS_SONOVUE_CDF',
        pvalue_role='DESCRIPTIVE_ONLY_NOT_A_MULTIPLE_TESTING_DISCOVERY_CLAIM')


def fractions(labels):
    c=Counter(labels);n=len(labels)
    return {b:c[b]/n if n else 0. for b in BASINS}


def tv(left,right):return float(sum(abs(left[b]-right[b]) for b in BASINS)/2)


def csv_write(path, rows):
    rows=list(rows)
    if not rows:return
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def admission_analysis(admission, report, *, development_partial=False):
    admission=Path(admission);report=Path(report);report.mkdir(parents=True,exist_ok=True)
    rows=list(iter_rows(admission));events=[r['event'] for r in rows];n=len(events)
    if n<100000 and not development_partial:raise ValueError('Formal common ledger requires at least 100000 events')
    c=context();distribution=c.distribution
    anchors=np.array([e['anchor_m'] for e in events]);uv=c.geometry.coordinates(anchors)*1e6
    sizes=np.array([e['first_diameter_um'] for e in events]);labels=[e['point_tracer_basin'] for e in events]
    before=fractions(labels);methods={};flat=[];invariants={}
    csv_write(report/'common_inlet_audit_ledger.csv',[
        dict(event_id=e['event_id'],x_m=e['anchor_m'][0],y_m=e['anchor_m'][1],z_m=e['anchor_m'][2],
             basin=e['point_tracer_basin'],first_diameter_um=e['first_diameter_um'],
             Dmax_um=e['Dmax_um'],geometry_pass_probability=e['geometry_pass_probability'],
             aperture_passable=e['aperture_passable'],center_on_plane_accepted=e['center_on_plane_accepted'],
             first_cause=e['current_first_cause']) for e in events])
    # Full RNG states remain in immutable common-event records, exported verbatim.
    with gzip.open(report/'COMMON_INLET_AUDIT_LEDGER.jsonl.gz','wt') as f:
        for e in events:f.write(json.dumps(e,ensure_ascii=False)+'\n')
    baseline_geometry={}
    for basin in BASINS:
        mask=np.array([b==basin for b in labels]);es=[e for e,m in zip(events,mask) if m]
        baseline_geometry[basin]=dict(scheduled=len(es),
            Dmax_um=descriptive([e['Dmax_um'] for e in es]),
            amax_um=descriptive([e['Dmax_um']/2 for e in es]),
            geometry_pass_probability=float(np.mean([e['geometry_pass_probability'] for e in es])) if es else None,
            geometry_probability_MC_standard_error=float(np.std([e['geometry_pass_probability'] for e in es],ddof=1)/np.sqrt(len(es))) if len(es)>1 else None,
            open_aperture_passable_fraction=float(np.mean([e['aperture_passable'] for e in es])) if es else None,
            center_on_plane_acceptance=float(np.mean([e['center_on_plane_accepted'] for e in es])) if es else None)
    arrays=dict(anchors_m=anchors,uv_um=uv,first_diameter_um=sizes,
        anchor_basin=np.array([BASINS.index(b) for b in labels]),Dmax_um=np.array([e['Dmax_um'] for e in events]))
    edges=[np.linspace(uv[:,k].min()-1e-5,uv[:,k].max()+1e-5,33) for k in range(2)]
    before_hist=np.histogram2d(uv[:,0],uv[:,1],bins=edges)[0]/n
    for method in 'ABC':
        records=[r['methods'][method] for r in rows];accepted=[a for a in records if a['accepted']]
        mask=np.array([a['accepted'] for a in records]);ds=np.array([a['diameter_um'] for a in accepted])
        after_labels=[a['birth_point_basin'] for a in accepted];anchor_after=[a['anchor_basin'] for a in accepted]
        failure=Counter()
        for a in records:failure.update(a.get('failure_counts',{}))
        final=np.array([a['birth_center_m'] for a in accepted]);finaluv=c.geometry.coordinates(final)*1e6
        hist=np.histogram2d(finaluv[:,0],finaluv[:,1],bins=edges)[0]/max(1,len(finaluv))
        group={}
        for basin in BASINS:
            selected=[a for a in accepted if a['anchor_basin']==basin]
            original_at_accepted=[e['first_diameter_um'] for e,a in zip(events,records) if a['accepted'] and a['anchor_basin']==basin]
            scheduled=[a for a in records if a['anchor_basin']==basin]
            group[basin]=dict(scheduled=len(scheduled),accepted=len(selected),
                acceptance_fraction=len(selected)/len(scheduled) if scheduled else None,
                original_diameters_um=descriptive(original_at_accepted),
                accepted_diameters_um=descriptive([a['diameter_um'] for a in selected]),
                size_distance=distribution_distance([a['diameter_um'] for a in selected],distribution),
                inward_distance_um=descriptive([a.get('s_birth_m',0.)*1e6 for a in selected]),
                inward_distance_over_radius=descriptive([a.get('s_birth_m',0.)/a['radius_m'] for a in selected]))
        methods[method]=dict(scheduled=n,accepted=len(accepted),acceptance_fraction=len(accepted)/n,
            status_counts=dict(Counter(a['status'] for a in records)),
            attempt_count=descriptive([a['attempt_count'] for a in records]),
            attempts_per_accepted_event=descriptive([a['attempt_count'] for a in accepted]),
            rejection_counts={key:failure[key] for key in CAUSES if key!='ACCEPTED'},
            total_candidate_attempts=sum(a['attempt_count'] for a in records),
            original_anchor_basin_fraction=before,accepted_anchor_basin_fraction=fractions(anchor_after),
            accepted_birth_point_basin_fraction=fractions(after_labels),
            accepted_birth_point_basin_counts=dict(Counter(after_labels)),
            anchor_to_birth_basin_counts={b:dict(Counter(a['birth_point_basin'] for a in accepted if a['anchor_basin']==b)) for b in BASINS},
            basin_total_variation_vs_original_anchor=tv(before,fractions(after_labels)),
            accepted_anchor_total_variation_vs_original_anchor=tv(before,fractions(anchor_after)),
            inlet_density_total_variation_32x32=float(np.abs(hist-before_hist).sum()/2) if method!='C' else None,
            density_comparison_role='C_ANCHOR_FIDELITY_SEPARATE_FROM_CURVED_INWARD_BIRTH_CENTER',
            accepted_size_um=descriptive(ds),size_distance=distribution_distance(ds,distribution),
            paired_original_size_um=descriptive(sizes[mask]),
            paired_mean_diameter_change_um=float(np.mean(ds-sizes[mask])),
            paired_size_KS=float(ks_2samp(ds,sizes[mask]).statistic),
            paired_size_Wasserstein_um=float(wasserstein_distance(ds,sizes[mask])),
            by_anchor_basin=group)
        if method=='C':
            methods[method].update(inward_distance_um=descriptive([a['s_birth_m']*1e6 for a in accepted]),
                inward_distance_over_radius=descriptive([a['s_birth_m']/a['radius_m'] for a in accepted]),
                accepted_but_anchor_aperture_fails=sum(not a['anchor_aperture_passable'] for a in accepted),
                accepted_and_anchor_aperture_passes=sum(a['anchor_aperture_passable'] for a in accepted),
                normal_fallback_count=sum(a['normal_fallback'] for a in records),
                geometry_evaluations=descriptive([a['evaluations'] for a in records]),
                unresolved_detail_counts=dict(Counter(a.get('detail','UNSPECIFIED') for a in records if a['status']=='NUMERICAL_GEOMETRY_UNRESOLVED')),
                physical_crossing_certified=False,
                important_limitation='A valid inner birth center alone cannot justify bypassing a solid-wall/perimeter conflict at the anchor.')
            birth_geometry=[]
            for a in accepted:
                wall,perimeter,cap=c.geometry.distances(a['birth_center_m'],full=True)
                birth_geometry.append(dict(event_id=a['event_id'],anchor_m=a['anchor_m'],birth_center_m=a['birth_center_m'],
                    radius_m=a['radius_m'],s_birth_m=a['s_birth_m'],wall_gap_m=wall-a['radius_m'],
                    open_cap_clearance_m=cap-a['radius_m'],full_domain_margin_m=a['full_margin_m'],
                    first_legal_bracket_m=a['bracket_m'],minimality_tolerance_m=a['minimality_tolerance_m'],
                    anchor_basin=a['anchor_basin'],birth_basin=a['birth_point_basin'],
                    anchor_aperture_passable=a['anchor_aperture_passable']))
            atomic_json(report/'METHOD_C_BIRTH_GEOMETRY.json',dict(rows=birth_geometry,
                entry_path_is_representation_not_physical_swept_sphere=True))
            methods[method]['birth_wall_gap_nm']=descriptive([a['wall_gap_m']*1e9 for a in birth_geometry])
            methods[method]['birth_open_cap_clearance_nm']=descriptive([a['open_cap_clearance_m']*1e9 for a in birth_geometry])
        invariants[method]=dict(common_input_hashes_match=all(a['common_event_sha256']==rows[i]['methods']['A']['common_event_sha256'] for i,a in enumerate(records)),
            anchor_records_unchanged=all(a['anchor_m']==e['anchor_m'] for a,e in zip(records,events)),
            size_draw_unchanged=all(a['diameter_um']==e['first_diameter_um'] for a,e in zip(records,events)) if method!='B' else None,
            B_accepted_positions_equal_anchors=all(a['birth_center_m']==a['anchor_m'] for a in accepted) if method=='B' else None)
        arrays[method+'_accepted']=mask;arrays[method+'_diameter_um']=np.array([a['diameter_um'] if a['accepted'] else np.nan for a in records])
        arrays[method+'_birth_m']=np.array([a['birth_center_m'] if a['accepted'] else [np.nan]*3 for a in records])
        arrays[method+'_s_birth_m']=np.array([a.get('s_birth_m',np.nan) for a in records])
        arrays[method+'_status']=np.array([a['status'] for a in records])
        arrays[method+'_birth_basin']=np.array([BASINS.index(a['birth_point_basin']) if a['accepted'] else -1 for a in records])
        for e,a in zip(events,records):
            pos=a['birth_center_m'] or [None]*3
            flat.append(dict(event_id=e['event_id'],method=method,anchor_basin=e['point_tracer_basin'],
                accepted=a['accepted'],status=a['status'],attempts=a['attempt_count'],
                first_diameter_um=e['first_diameter_um'],accepted_diameter_um=a['diameter_um'] if a['accepted'] else None,
                birth_x_m=pos[0],birth_y_m=pos[1],birth_z_m=pos[2],birth_point_basin=a['birth_point_basin'],
                s_birth_m=a.get('s_birth_m'),entry_transition_time_s=a.get('entry_transition_time_s'),
                common_event_sha256=a['common_event_sha256']))
    csv_write(report/'A_B_C_admission_events.csv',flat)
    atomic_npz(report/'admission_plot_data.npz',**arrays)
    official=dict(mean_um=distribution.target_mean_um(),quantiles_um=dict(zip(['D10','D25','D50','D75','D90','D95'],distribution.inverse_cdf([.1,.25,.5,.75,.9,.95]).tolist())),
        histogram_sha256=c.contract['histogram_sha256'],sampler_source_sha256=c.contract['sampler_source_sha256'])
    summary=dict(common_event_count=n,development_partial=development_partial,
        result_role='DEVELOPMENT_PARTIAL_NOT_FORMAL' if development_partial else 'FORMAL_100000_COMMON_EVENTS',
        point_tracer_basin_counts=dict(Counter(labels)),original_size=descriptive(sizes),
        original_draw_distance=distribution_distance(sizes,distribution),official_sonovue=official,
        geometry_by_basin=baseline_geometry,methods=methods,invariants=invariants,
        current_first_failure_causes=dict(Counter(e['current_first_cause'] for e in events)),
        full_domain_at_anchor_count=sum(e['full_domain_at_anchor'] for e in events),
        ledger_sha256=sha256(report/'COMMON_INLET_AUDIT_LEDGER.jsonl.gz'),
        shared_basin_warning='Unresolved point paths remain a separate group; labels never drive admission or size sampling.')
    atomic_json(report/'ADMISSION_COMPARISON.json',summary)
    return summary


def trajectory_analysis(root, report):
    root=Path(root);report=Path(report);methods={};flat=[]
    for method in 'ABC':
        records=[json.loads(p.read_text()) for p in sorted((root/'formal'/method/'trajectories').glob('mb_*.json')) if not p.name.endswith('.receipt.json')]
        completed=[r for r in records if r['completed']]
        states=Counter()
        for r in records:states.update(r.get('nearfield_states',{}))
        methods[method]=dict(admitted_integrated=len(records),completed=len(completed),
            outlet_counts={b:sum(r['exit_outlet']==b for r in records) for b in BASINS[:3]},
            end_reasons=dict(Counter(r['end_reason'] for r in records)),
            safety_stops=sum(r['end_reason']=='INTEGRATION_SAFETY_STOP' for r in records),
            stop_rate=sum(r['end_reason']=='INTEGRATION_SAFETY_STOP' for r in records)/max(1,len(records)),
            residence_time_s=descriptive([r['residence_time_s'] for r in completed]),
            completed_path_length_um=descriptive([r['path_length_m']*1e6 for r in completed]),
            all_saved_path_length_um=descriptive([r['path_length_m']*1e6 for r in records]),
            diameter_um=descriptive([r['diameter_um'] for r in records]),
            completed_diameter_um=descriptive([r['diameter_um'] for r in completed]),
            minimum_wall_gap_nm=descriptive([r['minimum_original_wall_gap_m']*1e9 for r in records if 'minimum_original_wall_gap_m' in r]),
            nearfield_saved_sample_counts=dict(states),
            represented_birth_residence_definition='FROM_FIRST_REPRESENTED_CENTER_TO_FIRST_CENTER_OUTLET_CROSSING',
            computational_stops_not_physical_capture=True)
        if method=='C':
            methods[method]['by_anchor_aperture_passability']={str(flag):dict(
                count=sum(r['birth_metadata']['entry_aperture_passable']==flag for r in records),
                outlet_counts={b:sum(r['birth_metadata']['entry_aperture_passable']==flag and r['exit_outlet']==b for r in records) for b in ['OUTLET_01','OUTLET_02','OUTLET_03']},
                safety_stops=sum(r['birth_metadata']['entry_aperture_passable']==flag and r['end_reason']=='INTEGRATION_SAFETY_STOP' for r in records)) for flag in [True,False]}
        for r in records:
            flat.append(dict(method=method,event_id=r['particle_id'],diameter_um=r['diameter_um'],
                anchor_basin=r['birth_metadata']['point_tracer_basin'],birth_basin=r['birth_metadata']['birth_point_basin'],
                completed=r['completed'],outlet=r['exit_outlet'],end_reason=r['end_reason'],
                residence_time_s=r['residence_time_s'],path_length_m=r['path_length_m'],
                minimum_wall_gap_m=r.get('minimum_original_wall_gap_m'),failure_detail=r['failure_detail'],
                sample_count=r['sample_count']))
    csv_write(report/'A_B_C_full_trajectory_catalog.csv',flat)
    atomic_json(report/'FULL_TRAJECTORY_COMPARISON.json',methods);return methods


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['admission','trajectories'])
    p.add_argument('--input',required=True);p.add_argument('--report',required=True)
    p.add_argument('--development-partial',action='store_true');a=p.parse_args()
    if a.action=='admission':admission_analysis(a.input,a.report,development_partial=a.development_partial)
    else:trajectory_analysis(a.input,a.report)


if __name__=='__main__':main()
