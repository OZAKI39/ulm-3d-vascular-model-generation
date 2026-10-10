"""Time-resolved causal audit on saved native states; no added mechanics cases."""
import json,math,time
import numpy as np
from scipy.spatial import cKDTree
from .config import dump
from .failure_map import csv_out
from .fragment_analysis import read_dump
from .contact_counterfactual import O,BASE,geometry

def first_time(t,mask):
    ids=np.flatnonzero(mask)
    return float(t[ids[0]]) if len(ids) else None

def category_intervals(t,modes,mode):
    # Left-hold duration on the saved grid. Final sample has zero extrapolated time.
    durations=np.diff(t,append=t[-1]);flag=modes==mode
    return dict(first_saved_time_s=first_time(t,flag),saved_grid_duration_s=float(durations[flag].sum()),saved_frames=int(flag.sum()),right_censored=bool(flag[-1]))

def summarize_case(name,path,factor):
    out=O/'diagnostics'/name;out.mkdir(parents=True,exist_ok=True)
    r=json.loads((path/'request.json').read_text());result=json.loads((path/'result.json').read_text());g=geometry();n=g['clot_count'];N=len(g['X']);b=g['bonds'];typ=g['bond_types'];keys=b[:,0]*N+b[:,1]
    z=np.load(path/'accepted_trajectory.npz');X=z['position_m'];V=z['velocity_m_s'];step=z['step'];t=step*r['dt_s'];labels=z['component_id'];attached=z['attached'];m=np.genfromtxt(path/'fragment_metrics.csv',delimiter=',',names=True,dtype=None,encoding='utf8')
    e=np.load(path/'broken_bonds.npz');events=e['events'];indices=e['initial_bond_indices'];alive=np.ones(len(b),bool);cursor=0;mass=r['rho_kg_m3']*r['dx_m']**3;dx=r['dx_m'];delta=r['delta_m'];kpair=18*r['K_Pa']/(np.pi*delta**4);kphysical=15*kpair*dx**6/delta
    rows=[];components=[];native=iter(read_dump(path/'trajectory.lammpstrj'));force_prefix=[]
    first_break=min(v for v in [result['first_bulk_break_time'],result['first_interface_break_time']] if v is not None)
    for k,s in enumerate(step):
        while cursor<len(events) and events[cursor,1]<=s:alive[indices[cursor]]=False;cursor+=1
        ns,cols,a=next(native);assert ns==s;a=a[np.argsort(a[:,0])]
        Fload=a[:,[cols.index('v_load_f'+d) for d in 'xyz']]
        if t[k]<first_break:force_prefix.append(a[:,[cols.index('f'+d) for d in 'xyz']])
        # Exclude every currently surviving bond, matching special_bonds 0 1 1.
        pairs=cKDTree(X[k]).query_pairs(1.35*dx,output_type='ndarray');pk=pairs[:,0]*N+pairs[:,1];loc=np.searchsorted(keys,pk);valid=loc<len(keys);bonded=np.zeros(len(pk),bool);ii=np.flatnonzero(valid);bonded[ii]=(keys[loc[ii]]==pk[ii])&alive[loc[ii]]
        pairs=pairs[(~bonded)&(pairs[:,0]<n)];dist=np.linalg.norm(X[k,pairs[:,1]]-X[k,pairs[:,0]],axis=1)
        deep=pairs[dist<.5*dx];affected=np.unique(deep.ravel());affected=affected[affected<n]
        active=dist<factor*dx if factor is not None else np.zeros(len(dist),bool)
        energy_contact=float(.5*kphysical*np.square(factor*dx-dist[active]).sum()) if factor is not None else 0.
        ids,inv,count=np.unique(labels[k],return_inverse=True,return_counts=True);nc=len(ids);vcom=np.stack([np.bincount(inv,weights=V[k,:n,d],minlength=nc)/count for d in range(3)],axis=1);com=np.stack([np.bincount(inv,weights=X[k,:n,d],minlength=nc)/count for d in range(3)],axis=1)
        att=np.bincount(inv,weights=attached[k].astype(float),minlength=nc)>0;free=np.flatnonzero(~att);comKE=.5*mass*count*np.sum(vcom*vcom,axis=1)
        for j in sorted(free,key=lambda j:count[j],reverse=True)[:5]:
            components.append(dict(step=int(s),time_s=float(t[k]),component_min_node_id=int(ids[j]),node_count=int(count[j]),mass_kg=float(mass*count[j]),COM_x_m=com[j,0],COM_y_m=com[j,1],COM_z_m=com[j,2],COM_vx_m_s=vcom[j,0],COM_vy_m_s=vcom[j,1],COM_vz_m_s=vcom[j,2],COM_KE_J=float(comKE[j]),identity_rule='minimum original node id; a daughter may inherit its parent label'))
        row={key:m[k][key].item() if hasattr(m[k][key],'item') else m[k][key] for key in m.dtype.names}
        row.update(contact_case=name,bulk_survival=1-row['bulk_broken_fraction'],interface_survival=1-row['interface_broken_fraction'],largest_clot_component_fraction=row['largest_free_fraction'],free_component_count=len(free),free_resolved_component_count=int((count[free]>=3).sum()),clot_KE_J=float(.5*mass*np.sum(V[k,:n]**2)),free_component_COM_KE_J=float(comKE[free].sum()),external_power_W=float(np.sum(Fload*V[k])),numerical_drag_power_W=float(-r['numerical_drag_rate_s_inv']*mass*np.sum(V[k,:n]**2)),nonbonded_pairs_under_half_dx=len(deep),clot_fraction_in_deep_overlap=len(affected)/n,minimum_nonbonded_distance_dx=float(dist.min()/dx) if len(dist) else None,contact_active_nonbonded_pairs=int(active.sum()),reconstructed_contact_overlap_potential_J=energy_contact,max_clot_speed_m_s=float(np.linalg.norm(V[k,:n],axis=1).max()))
        rows.append(row)
    assert next(native,None) is None
    for power,work in [('external_power_W','sampled_external_work_J'),('numerical_drag_power_W','sampled_drag_work_J')]:
        p=np.array([a[power] for a in rows]);w=np.r_[0,np.cumsum(.5*(p[:-1]+p[1:])*np.diff(t))]
        for a,val in zip(rows,w):a[work]=float(val)
    csv_out(out/'timeseries.csv',rows);csv_out(out/'key_independent_components.csv',components)
    np.savez_compressed(out/'prebreak_total_force.npz',force_N=force_prefix)
    summary=dict(case=name,source_path=str(path),contact_distance_factor=factor,first_bulk_break_time_s=result['first_bulk_break_time'],first_interface_break_time_s=result['first_interface_break_time'],first_any_break_time_s=first_break,
        first_independent_component_saved_time_s=first_time(t,np.array([a['free_component_count']>0 for a in rows])),first_independent_resolved_component_saved_time_s=first_time(t,np.array([a['free_resolved_component_count']>0 for a in rows])),
        categories={mode:category_intervals(t,m['failure_mode'],mode) for mode in np.unique(m['failure_mode'])},final=rows[-1],peak_clot_KE_J=max(a['clot_KE_J'] for a in rows),max_deep_overlap_fraction=max(a['clot_fraction_in_deep_overlap'] for a in rows),runtime=result['runtime'],postprocess_wall_s=result['postprocess_wall_s'],scope='DIAGNOSTIC_ONLY' if factor is None else 'EXPLORATORY_ONLY')
    dump(out/'summary.json',summary);return summary,rows

def compare():
    start=time.perf_counter();criteria=json.loads((O/'PREDECLARED_CAUSAL_CRITERIA.json').read_text());cfg=json.loads((O/'FROZEN_BASELINE.json').read_text())['request']
    cases=[('CONTACT_ON',BASE,1.35),('CONTACT_OFF',O/'CONTACT_OFF',None)]
    if (O/'CONTACT_MODIFIED/result.json').exists():cases.append(('CONTACT_MODIFIED',O/'CONTACT_MODIFIED',1.0))
    summaries={};series={};allrows=[]
    for name,path,f in cases:
        summaries[name],series[name]=summarize_case(name,path,f);allrows.extend(series[name]);print('DIAGNOSED',name,summaries[name]['final']['failure_mode'],flush=True)
    on=np.load(BASE/'accepted_trajectory.npz');off=np.load(O/'CONTACT_OFF/accepted_trajectory.npz');assert np.array_equal(on['step'],off['step']);t=on['step']*cfg['dt_s'];dx=cfg['dx_m'];cs=np.sqrt(cfg['G_Pa']/cfg['rho_kg_m3'])
    earliest=min(summaries[x]['first_any_break_time_s'] for x in ['CONTACT_ON','CONTACT_OFF']);pre=t<earliest
    dp=np.linalg.norm(on['position_m']-off['position_m'],axis=2).max(axis=1)/dx;dv=np.linalg.norm(on['velocity_m_s']-off['velocity_m_s'],axis=2).max(axis=1)/cs
    aa=np.load(O/'diagnostics/CONTACT_ON/prebreak_total_force.npz')['force_N'][:pre.sum()];bb=np.load(O/'diagnostics/CONTACT_OFF/prebreak_total_force.npz')['force_N'][:pre.sum()]
    first_delta=abs(summaries['CONTACT_ON']['first_any_break_time_s']-summaries['CONTACT_OFF']['first_any_break_time_s'])/cfg['dt_s']
    parity=pre.sum()>=criteria['minimum_prebreak_saved_frames'] and dp[pre].max()<=criteria['pre_break_position_tolerance_dx'] and dv[pre].max()<=criteria['pre_break_velocity_tolerance_shear_wave_speed'] and first_delta<=criteria['first_any_break_tolerance_steps']+1e-8
    parity_data=dict(status='PASS' if parity else 'FAIL',saved_frames=int(pre.sum()),last_prebreak_saved_time_s=float(t[pre][-1]),max_position_error_dx=float(dp[pre].max()),max_velocity_error_over_shear_wave_speed=float(dv[pre].max()),max_total_force_error_N=float(abs(aa-bb).max()),first_any_break_delta_steps=first_delta)
    bon=np.array([x['bulk_broken_fraction'] for x in series['CONTACT_ON']]);boff=np.array([x['bulk_broken_fraction'] for x in series['CONTACT_OFF']]);don=np.array([x['debris_fraction'] for x in series['CONTACT_ON']]);doff=np.array([x['debris_fraction'] for x in series['CONTACT_OFF']])
    onset_on=summaries['CONTACT_ON']['categories'].get('PULVERIZATION',{}).get('first_saved_time_s');onset_off=summaries['CONTACT_OFF']['categories'].get('PULVERIZATION',{}).get('first_saved_time_s');ratio=onset_off/onset_on if onset_off is not None and onset_on else None;drop=float(don[-1]-doff[-1])
    contrast=first_time(t,don-doff>=criteria['primary_terminal_debris_absolute_reduction'])
    if onset_on is not None and (onset_off is None or onset_off>onset_on):contrast=min(contrast,onset_on) if contrast is not None else onset_on
    deep=np.array([x['clot_fraction_in_deep_overlap']>=criteria['severe_penetration_clot_fraction'] for x in series['CONTACT_OFF']]);sustained=np.convolve(deep.astype(int),np.ones(criteria['severe_penetration_consecutive_saved_frames'],int),mode='valid')>=criteria['severe_penetration_consecutive_saved_frames'];deep_start=first_time(t[:len(sustained)],sustained);confound=deep_start is not None and (contrast is None or deep_start<=contrast)
    primary=drop>=criteria['primary_terminal_debris_absolute_reduction'] or (onset_on is not None and (onset_off is None or ratio>=criteria['primary_pulverization_onset_ratio']))
    if not parity or confound:role='INCONCLUSIVE'
    elif primary:role='PRIMARY_CONTRIBUTOR_SUPPORTED'
    elif drop>=criteria['moderate_debris_reduction'] or (ratio is not None and ratio>=criteria['moderate_onset_ratio']):role='PARTIAL_CONTRIBUTOR'
    else:role='NOT_PRIMARY'
    comparison=dict(scope='EXPLORATORY_ONLY',pre_break_parity=parity_data,contact_role=role,stage_B_allowed=role=='PRIMARY_CONTRIBUTOR_SUPPORTED',terminal_debris_absolute_reduction=drop,pulverization_time_ratio=ratio,
        first_clear_causal_contrast_time_s=contrast,first_sustained_deep_overlap_time_s=deep_start,penetration_confounds_comparison=confound,
        first_position_divergence_above_0p01dx_s=first_time(t,dp>.01),first_bulk_broken_fraction_difference_gt_0p001_s=first_time(t,abs(bon-boff)>.001),
        cases=summaries,predeclared_criteria=criteria,analysis_wall_s=time.perf_counter()-start,energy_semantics='Physical nodal KE and COM KE in J; external/drag work from saved-grid quadrature; contact overlap potential reconstructed by integrating repulsive force to onset. No closed energy balance, no inferred fracture energy, no native pe used as dissipation.')
    csv_out(O/'CAUSAL_COMPARISON.csv',allrows);dump(O/'CAUSAL_COMPARISON.json',comparison);dump(O/'checks/PRE_BREAK_PARITY.json',parity_data)
    (O/'CONTACT_ROLE_DECISION.md').write_text(f"# Frozen-threshold Stage A decision\n\nCONTACT_ROLE={role}; PRE_BREAK_PARITY={parity_data['status']}. Terminal debris reduction={drop:.8f}; first pulverization on/off={onset_on}/{onset_off} s; onset ratio={ratio}. First clear contrast={contrast} s; sustained deep-overlap onset={deep_start} s; confound={confound}.\n\nStage B allowed: {comparison['stage_B_allowed']}. Criteria were saved and hashed before CONTACT_OFF was run. CONTACT_OFF remains DIAGNOSTIC_ONLY and cannot enter production. Deep-overlap counts are saved-time geometric diagnostics at the frozen 0.5dx/1%/3-frame thresholds, not a proof of continuum nonpenetration.\n")
    print('CAUSAL_DECISION',role,'debris_reduction',drop,'parity',parity_data,'overlap',deep_start,flush=True);return comparison
