"""Raw field comparisons, scale profiles and predeclared causal-control metrics."""
import csv
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from scripts.activation_root_cause import ROOT,OUT,RESULTS,STRAINS,read,save,arrays
from pd_clot.pfpd_candidate.activation_boundary_control import BoundaryProblem,FULL,AXIAL


def stats(values):
    a=np.asarray(values,float)
    if len(a)==0:return dict(count=0,mean=None,median=None,P95=None,max=None,CV=None)
    mean=float(a.mean())
    return dict(count=len(a),mean=mean,median=float(np.median(a)),P95=float(np.percentile(a,95)),
                max=float(a.max()),CV=float(a.std()/abs(mean)) if mean else 0.)


def csv_save(path, rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def discrepancy(a,b):
    return float(abs(a-b)/(.5*(abs(a)+abs(b))))


def historical_ids(n):
    if n==5:
        d=read(ROOT/'verification/pfpd_at1_c5_activation_stationarity_audit_v1/activation/POINT.json')
        return np.asarray(d['IDs']['WEAK_LOWER'],int)
    d=read(ROOT/'verification/pfpd_at1_c7_early_branch_v1/activation/WEAK.json')
    return np.asarray(d['q_IDs']['WEAK_LOWER'],int)


def analyze(n,e,mode):
    p=BoundaryProblem(n,e,mode);f=p.family
    folder=OUT/'equilibria'/f'C{n}/{mode}/{e:.8f}'
    a=arrays(folder/'STATE.npz');x=a['x'];s=a['s']
    eq=p.state(x,s);sp=eq.split_response.split
    aff=p.state(p.affine(),p.lower);asp=aff.split_response.split
    cfg=read(OUT/'CONFIG.json');Y=sp.Psi_plus;Ya=asp.Psi_plus
    safe=Ya>cfg['A_psi_denominator_floor_fraction_Yc']*p.yc
    amp=np.divide(Y,Ya,out=np.full_like(Y,np.nan),where=safe)
    epos=np.maximum(sp.eigenvalues.max(axis=1),0)
    epos_a=np.maximum(asp.eigenvalues.max(axis=1),0)
    ae=np.divide(epos,epos_a,out=np.full_like(Y,np.nan),where=epos_a>1e-15)
    mid=f.X[f.pairs].mean(axis=1);orientation=f.xi/f.length[:,None]
    d=.000125-abs(mid[:,2]);dh=d/p.cloud.spacing;dd=d/f.horizon
    axial=abs(orientation[:,2])>=.95
    first=axial&np.isclose(dh,.5,rtol=0,atol=1e-8)
    i,j=f.pairs.T;neighbor=np.bincount(f.pairs.ravel(),minlength=len(f.X))
    condition=np.linalg.cond(eq.K);endpoint_support=f.omega0[f.pairs].mean(axis=1)
    endpoint_condition=condition[f.pairs].max(axis=1)
    ids=historical_ids(n);historical_mask=np.zeros(p.ns,bool);historical_mask[ids]=True
    rows=[]
    for b in np.flatnonzero(axial):
        rows.append(dict(bond_ID=int(b),historical_activation=bool(historical_mask[b]),i=int(i[b]),j=int(j[b]),
            midpoint_x_m=float(mid[b,0]),midpoint_y_m=float(mid[b,1]),midpoint_z_m=float(mid[b,2]),
            d_clamp_m=float(d[b]),d_over_h=float(dh[b]),d_over_delta=float(dd[b]),abs_nz=float(abs(orientation[b,2])),
            Psi_plus_affine_Pa=float(Ya[b]),Psi_plus_eq_Pa=float(Y[b]),Psi_minus_eq_Pa=float(sp.Psi_minus[b]),
            A_psi=float(amp[b]) if safe[b] else None,delta_Psi_Pa=float(Y[b]-Ya[b]),
            Psi_normalized=float(Y[b]/p.yc),q_over_Yc=float(2*(p.yc-Y[b])/p.yc),
            max_positive_E_eq=float(epos[b]),max_positive_E_affine=float(epos_a[b]),A_E=float(ae[b]),
            actual_bond_stretch=float(np.linalg.norm(x[j[b]]-x[i[b]])/f.length[b]-1),
            endpoint_support_mean=float(endpoint_support[b]),endpoint_shape_condition_max=float(endpoint_condition[b])))
    csv_save(RESULTS/'tables'/f'C{n}_{mode}_{e:.4f}_AXIAL_BONDS.csv',rows)
    profile=[]
    # Physical midpoint layers shared across orientation classes; axial-only profile.
    for layer in np.unique(np.round(dh[axial],10)):
        mask=axial&np.isclose(dh,layer,rtol=0,atol=1e-8)
        profile.append(dict(d_clamp_m=float(d[mask].mean()),d_over_h=float(layer),d_over_delta=float(dd[mask].mean()),
            count=int(mask.sum()),mean_A_psi=float(amp[mask].mean()),P95_A_psi=float(np.percentile(amp[mask],95)),
            max_A_psi=float(amp[mask].max()),mean_Psi_normalized=float(Y[mask].mean()/p.yc),
            mean_q_over_Yc=float((2*(p.yc-Y[mask])/p.yc).mean()),affine_Psi_stats=stats(Ya[mask]),equilibrium_Psi_stats=stats(Y[mask])))
    save(RESULTS/'tables'/f'C{n}_{mode}_{e:.4f}_PROFILE.json',profile)
    un=x-p.affine();u=x-f.X;layers=[]
    for z in np.unique(f.X[:,2]):
        mask=f.X[:,2]==z;distance=.000125-abs(z)
        layers.append(dict(z_m=float(z),distance_over_h=float(distance/p.cloud.spacing),
            nonaffine_RMS_m=float(np.sqrt(np.mean(np.sum(un[mask]**2,axis=1)))),
            mean_uz_deviation_m=float(un[mask,2].mean()),mean_lateral_displacement_m=float(np.linalg.norm(u[mask,:2],axis=1).mean()),
            max_lateral_displacement_m=float(np.linalg.norm(u[mask,:2],axis=1).max()),
            mean_ux_m=float(u[mask,0].mean()),mean_uy_m=float(u[mask,1].mean())))
    csv_save(RESULTS/'tables'/f'C{n}_{mode}_{e:.4f}_NONAFFINE_LAYERS.csv',layers)
    support_rows=[]
    activation_nodes=np.unique(f.pairs[ids]);zvalues=f.X[:,2];most_interior=np.min(.000125-abs(f.X),axis=1)
    central=np.isclose(most_interior,most_interior.max(),atol=1e-15)
    for k in range(len(f.X)):
        support_rows.append(dict(particle_ID=k,x_m=float(f.X[k,0]),y_m=float(f.X[k,1]),z_m=float(f.X[k,2]),
            neighbor_count=int(neighbor[k]),weighted_family_measure=float(f.omega0[k]),
            shape_condition=float(condition[k]),shape_trace=float(np.trace(eq.K[k])),
            moment_eigen_min=float(np.linalg.eigvalsh(eq.K[k])[0]),
            normalization_reciprocal_omega0=float(1/f.omega0[k]),
            activation_endpoint=bool(k in activation_nodes),most_interior_available=bool(central[k])))
    if e==.003 and mode==FULL:csv_save(RESULTS/'tables'/f'C{n}_PARTICLE_SUPPORT.csv',support_rows)
    peers=[]
    for k in activation_nodes:
        same=zvalues==zvalues[k]
        peers.append(dict(particle_ID=int(k),position_m=f.X[k],neighbor_count=int(neighbor[k]),
            weighted_family_measure=float(f.omega0[k]),shape_condition=float(condition[k]),
            same_layer_support=stats(f.omega0[same]),same_layer_condition=stats(condition[same]),
            same_layer_neighbor_count=stats(neighbor[same]),
            support_percentile_in_layer=float(np.mean(f.omega0[same]<=f.omega0[k])*100),
            condition_percentile_in_layer=float(np.mean(condition[same]<=condition[k])*100)))
    correlations={}
    for name,v in [('weighted_support',endpoint_support),('shape_condition',endpoint_condition),('d_over_h',dh)]:
        rho,prob=spearmanr(amp[axial&safe],v[axial&safe]);correlations[name]=dict(rho=float(rho),p_value=float(prob))
    orientation_groups=[]
    for lo,hi in zip([0,.25,.5,.75,.95],[.25,.5,.75,.95,1.0000001]):
        mask=(abs(orientation[:,2])>=lo)&(abs(orientation[:,2])<hi)
        orientation_groups.append(dict(abs_nz_interval=[lo,min(hi,1.)],affine_Psi_plus=stats(Ya[mask]),
                                       equilibrium_Psi_plus=stats(Y[mask]),A_psi=stats(amp[mask&safe])))
    expected=np.diag([1,1,1+e]);Eexpected=.5*(expected.T@expected-np.eye(3))
    affine_check=dict(Fbar_max_error=float(abs(aff.Fbar-expected).max()),Ftilde_max_error=float(abs(aff.Ftilde-expected).max()),
        E_max_error=float(abs(asp.E-Eexpected).max()),expected_F=expected,expected_E=Eexpected,
        expected_principal_Green_strains=np.linalg.eigvalsh(Eexpected),
        principal_Green_min_over_bonds=asp.eigenvalues.min(axis=0),principal_Green_max_over_bonds=asp.eigenvalues.max(axis=0),
        Psi_plus_all=stats(Ya),Psi_plus_axial=stats(Ya[axial]),Psi_minus_max=float(asp.Psi_minus.max()),
        q_min=float((2*(p.yc-Ya)).min()),q_max=float((2*(p.yc-Ya)).max()),
        relative_driver_span=float(np.ptp(Ya)/Ya.mean()),
        free_force_residual=float(np.linalg.norm(aff.force[p.free])/p.config['mechanical_equilibrium']['force_floor_N']))
    # Existing-field diagnostic, not another solve: resolve the centerline gaps
    # responsible for the strongest axial bond strain.
    central_ids=np.flatnonzero(np.max(abs(f.X[:,:2]),axis=1)<1e-15)
    central_ids=central_ids[np.argsort(f.X[central_ids,2])]
    centerline=[]
    for left,right in zip(central_ids[:-1],central_ids[1:]):
        gap=float((x[right,2]-f.X[right,2])-(x[left,2]-f.X[left,2]))
        centerline.append(dict(left_ID=int(left),right_ID=int(right),z_left_m=float(f.X[left,2]),
            z_right_m=float(f.X[right,2]),reference_gap_m=float(f.X[right,2]-f.X[left,2]),
            displacement_gap_m=gap,engineering_strain=gap/float(f.X[right,2]-f.X[left,2]),
            fraction_of_global_displacement=gap/(e*.00025)))
    csv_save(RESULTS/'tables'/f'C{n}_{mode}_{e:.4f}_CENTERLINE_GAPS.csv',centerline)
    affine_force_layers=[]
    for z in np.unique(f.X[:,2]):
        mask=f.X[:,2]==z
        affine_force_layers.append(dict(z_m=float(z),affine_axial_force_sum_N=float(aff.force[mask,2].sum()),
            affine_force_L2_N=float(np.linalg.norm(aff.force[mask])),equilibrium_axial_force_sum_N=float(eq.force[mask,2].sum())))
    top=np.flatnonzero(Y>=Y.max()*(1-1e-8));peak=[]
    for k in top:
        peak.append(dict(bond_ID=int(k),endpoints=f.pairs[k],midpoint_m=mid[k],Ftilde=eq.Ftilde[k],
            principal_Green_strains=sp.eigenvalues[k],Psi_plus_Pa=float(Y[k]),Psi_minus_Pa=float(sp.Psi_minus[k]),
            A_psi=float(amp[k]),A_E=float(ae[k]),d_over_h=float(dh[k]),d_over_delta=float(dd[k]),
            actual_bond_stretch=float(np.linalg.norm(x[j[k]]-x[i[k]])/f.length[k]-1)))
    summary=dict(level=n,strain=e,BC_MODE=mode,affine=affine_check,
        equilibrium=read(folder/'SUMMARY.json'),A_psi_axial=stats(amp[axial&safe]),
        A_psi_all=stats(amp[safe]),A_E_axial=stats(ae[axial]),delta_Psi_axial=stats((Y-Ya)[axial]),
        first_layer_Psi_normalized=stats(Y[first]/p.yc),first_layer_A_psi=stats(amp[first]),
        first_layer_A_E=stats(ae[first]),correlations=correlations,orientation_groups=orientation_groups,
        nonaffine_RMS_m=float(np.sqrt(np.mean(np.sum(un**2,axis=1)))),layers=layers,profile=profile,
        centerline_gaps=centerline,affine_force_layers=affine_force_layers,
        activation_endpoint_comparison=peers,most_interior_support=stats(f.omega0[central]),
        most_interior_condition=stats(condition[central]),peak_bonds=peak,
        Ftilde_bond_kinematic_identity_max_error=float(abs(np.einsum('bij,bj->bi',eq.Ftilde,f.xi)-(x[j]-x[i])).max()),
        near_zero_affine_denominator_excluded=int((~safe).sum()))
    save(OUT/'analysis'/f'C{n}_{mode}_{e:.4f}.json',summary)
    return summary


def collapse(profiles):
    cfg=read(OUT/'CONFIG.json');rows=[]
    for e in STRAINS:
        levels=[profiles[(n,e,FULL)]['profile'] for n in (5,7)]
        result=dict(strain=e,metrics={})
        for metric in ['mean_A_psi','P95_A_psi','max_A_psi','mean_Psi_normalized','mean_q_over_Yc']:
            entry={}
            for coord,label in [('d_over_h','h'),('d_over_delta','delta')]:
                xs=[np.array([r[coord] for r in level]) for level in levels]
                ys=[np.array([r[metric] for r in level]) for level in levels]
                grid=np.linspace(max(x[0] for x in xs),min(x[-1] for x in xs),cfg['profile']['common_grid_points'])
                values=[np.interp(grid,x,y) for x,y in zip(xs,ys)]
                entry['RMSE_'+label]=float(np.sqrt(np.mean((values[0]-values[1])**2)))
                entry['domain_'+label]=[float(grid[0]),float(grid[-1])]
                csv_save(RESULTS/'tables'/f'COLLAPSE_{e:.4f}_{metric}_{label}.csv',
                         [dict(coordinate=float(g),C5=float(v),C7=float(w)) for g,v,w in zip(grid,*values)])
            result['metrics'][metric]=entry
        rows.append(result)
    primary=[r['metrics']['mean_A_psi'] for r in rows]
    label=('GRID_LAYER_SUPPORTED' if all(r['RMSE_h']<.5*r['RMSE_delta'] for r in primary) else
           'PHYSICAL_HORIZON_SUPPORTED' if all(r['RMSE_delta']<.5*r['RMSE_h'] for r in primary) else 'MIXED')
    out=dict(classification=label,primary_metric='mean_A_psi',rows=rows,fit='none',interpolation=cfg['profile']['interpolation'])
    save(OUT/'PROFILE_SCALING.json',out);return out


def decision(profiles,scaling):
    cfg=read(OUT/'CONFIG.json');comparison=[]
    for n in (5,7):
        full=profiles[(n,.003,FULL)];axial=profiles[(n,.003,AXIAL)]
        hotspot_reduction=1-axial['first_layer_Psi_normalized']['max']/full['first_layer_Psi_normalized']['max']
        comparison.append(dict(level=n,first_layer_hotspot_reduction=float(hotspot_reduction),
            baseline_max_driver=full['equilibrium']['max_driver_normalized'],control_max_driver=axial['equilibrium']['max_driver_normalized'],
            baseline_first_layer=full['first_layer_Psi_normalized'],control_first_layer=axial['first_layer_Psi_normalized'],
            baseline_max_Apsi=full['A_psi_all']['max'],control_max_Apsi=axial['A_psi_all']['max'],
            baseline_reaction=full['equilibrium']['reaction_stress_Pa'],control_reaction=axial['equilibrium']['reaction_stress_Pa']))
    db=discrepancy(*(r['baseline_max_driver'] for r in comparison));dc=discrepancy(*(r['control_max_driver'] for r in comparison))
    reduction=1-dc/db
    clamp=('STRONG' if all(r['first_layer_hotspot_reduction']>=.2 for r in comparison) and reduction>=.3 else
           'WEAK' if all(abs(r['first_layer_hotspot_reduction'])<.05 for r in comparison) and abs(reduction)<.05 else 'MIXED')
    rho={name:[profiles[(n,.003,FULL)]['correlations'][name]['rho'] for n in (5,7)] for name in ['weighted_support','shape_condition']}
    support=('STRONG' if any(all(abs(v)>=.7 for v in values) for values in rho.values()) else
             'WEAK' if all(abs(v)<.3 for values in rho.values() for v in values) else 'MIXED')
    affine='ABSENT' if all(r['affine']['Psi_plus_all']['CV']<=1e-8 and r['affine']['relative_driver_span']<=1e-8 for r in profiles.values()) else 'PRESENT'
    hotspot='PRESENT' if all(profiles[(n,.003,FULL)]['A_psi_all']['max']>=1.2 for n in (5,7)) else 'ABSENT'
    if affine=='PRESENT':root='BAQP_RECONSTRUCTION_BOUNDARY_EFFECT'
    elif clamp=='STRONG':root='FULL_CLAMP_KINEMATIC_BOUNDARY_LAYER'
    elif clamp=='WEAK' and support=='STRONG':root='NONLOCAL_SUPPORT_TRUNCATION'
    elif clamp=='MIXED' and support=='STRONG' and all(r['first_layer_hotspot_reduction']>=.2 for r in comparison):root='COMBINED_CLAMP_AND_SUPPORT'
    else:root='INCONCLUSIVE'
    predictions=[]
    for n in (5,7):
        for mode in (FULL,AXIAL):
            y=[profiles[(n,e,mode)]['equilibrium']['max_Psi_plus_Pa'] for e in STRAINS]
            slope=(y[1]-y[0])/(STRAINS[1]-STRAINS[0]);prediction=STRAINS[1]+(cfg['Yc_Pa']-y[1])/slope
            historical=read(ROOT/'verification'/('pfpd_at1_c5_activation_stationarity_audit_v1/activation/POINT.json' if n==5 else 'pfpd_at1_c7_early_branch_v1/activation/WEAK.json'))['strain']
            predictions.append(dict(level=n,mode=mode,epsilon_prediction=float(prediction),local_linear_slope_Pa=slope,
                true_baseline_activation=historical if mode==FULL else None,
                relative_prediction_error=float(prediction/historical-1) if mode==FULL else None,
                interpretation='two-point pre-activation mechanical extrapolation, not an activation result'))
    baseline_eps=[r['true_baseline_activation'] for r in predictions if r['mode']==FULL]
    result=dict(MODEL=cfg['model'],BASELINE_ACTIVATION_OBJECTIVITY='NOT_ESTABLISHED',AFFINE_HOTSPOT=affine,
        EQUILIBRIUM_HOTSPOT=hotspot,PROFILE_SCALING=scaling['classification'],SUPPORT_ASSOCIATION=support,
        CLAMP_ASSOCIATION=clamp,ROOT_CAUSE_PRIMARY=root,CANDIDATE_FIX='NONE',CANDIDATE_IMPLEMENTED='NO',
        C5_C7_ACTIVATION_DISCREPANCY_BASELINE=discrepancy(*baseline_eps),C5_C7_ACTIVATION_DISCREPANCY_CANDIDATE=None,
        BASELINE_SHIFT_RELATIVE_TO_C5=float(baseline_eps[1]/baseline_eps[0]-1),CANDIDATE_EFFECT='NOT_RUN',
        PHASE_SUPPORT_LOCALIZATION='NOT_ASSESSED',C7_STABILITY='NOT_RUN',C9='NOT_RUN',
        clamp_comparison=comparison,baseline_driver_discrepancy=db,control_driver_discrepancy=dc,
        driver_discrepancy_reduction=reduction,support_rho=rho,predictions=predictions,
        phase_A_complete=True,root_cause_rule='predeclared control and association thresholds; correlation is not causal proof',
        stop_reason='DIAGNOSIS_INCONCLUSIVE' if root=='INCONCLUSIVE' else 'PHASE_A_COMPLETE_CANDIDATE_REVIEW')
    save(OUT/'PHASE_A_DECISION.json',result);return result


def main():
    profiles={}
    for n in (5,7):
        for e in STRAINS:
            for mode in (FULL,AXIAL):profiles[(n,e,mode)]=analyze(n,e,mode)
    scaling=collapse(profiles);result=decision(profiles,scaling)
    print({k:v for k,v in result.items() if k.isupper()},flush=True)


if __name__=='__main__':main()
