"""Complete-window fixed-ID, shape and event comparisons, with honest cost limits."""
import numpy as np
import h5py
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from .common import *
from .metrics import series, validate, edge_events

def load(folder):
    r = read(folder/'RESULT.json')
    assert sha(folder/'trajectory.h5') == r['trajectory_sha256']
    with h5py.File(folder/'trajectory.h5') as h:
        t, x = validate(h['time_s'][:], h['positions_m'][:], h['particle_id'][:])
        assert h.attrs['flow']=='STRICT_ZERO'
    assert x.shape[1] == r.get('N', x.shape[1])
    assert len(t) == r['accepted_steps'] + 1 and abs(t[-1]-r['completed_s'])<1e-13
    if r['status']=='PASS': assert abs(t[-1]-r['requested_s'])<1e-13
    return r, t, x

def summarize(folder):
    r, t, x = load(folder); m = series(t, x)
    with h5py.File(folder/'trajectory.h5') as h:
        steps = h['step_metrics'][:]
        cols = h['step_metrics'].attrs['columns'].split(',')
        if 'max_pressure_Pa' in cols:
            p = steps[:, cols.index('max_pressure_Pa')]
            pmin = steps[:, cols.index('min_pressure_Pa')]
        else:
            extra = h['extra_metrics'][:]
            pmin, p = extra[:, 3], extra[:, 4]
        force = h['force_samples/total_force_N'][:]
    # Pressure domain is checked per particle by the engine and independently by each step extrema.
    legal_pressure = bool(np.all(((pmin==0)&(p==0)) | ((pmin>=.1)&(p<=40000))))
    ppolicy=policy()
    independently_checked={}
    if 'force_scale_N' in cols:
        limits={'residual_m':ppolicy['mobility_residual_m'],
            'overlap_m':min(ppolicy['hard_geometry_m'],ppolicy['contact_overlap_m']),
            'max_iteration_m':ppolicy['contact_correction_m'], 'total_contact_m':ppolicy['contact_correction_m'],
            'wall_m':0., 'primal_m':ppolicy['primal_violation_m'], 'dual_N':ppolicy['dual_negative_N'],
            'complementarity_m':ppolicy['complementarity_equivalent_gap_m'],
            'actual_complementarity_m':ppolicy['actual_geometry_complementarity_m']}
        independently_checked={k:bool(np.all(steps[:,cols.index(k)]<=v)) for k,v in limits.items()}
        independently_checked['force_balance']=bool(np.all(steps[:,cols.index('force_residual_N')]<=
            ppolicy['force_residual_atol_N']+ppolicy['force_residual_rtol']*steps[:,cols.index('force_scale_N')]))
        independently_checked['all_attempts_accepted']=bool(np.all(steps[:,0]==1))
        if r['status']=='PASS': assert all(independently_checked.values()), independently_checked
    zero_motion = float(np.linalg.norm(x-x[0], axis=-1).max())
    events = edge_events(t, m['edges'])
    persistence = {}
    for threshold in [.0002, .002]:
        count = 0
        for k in range(m['edges'].shape[1]):
            e = m['edges'][:, k]
            starts = np.flatnonzero(e & ~np.r_[False, e[:-1]])
            ends = np.flatnonzero(~e & np.r_[False, e[:-1]])
            for start in starts:
                if start==0: continue
                later = ends[ends>start]; end = t[later[0]] if len(later) else t[-1]
                count += int(end-t[start]>=threshold-1e-13)
        persistence[str(threshold)] = count
    r.update(N=x.shape[1], input_path=str(folder.relative_to(ROOT)), Rg_initial_final_m=m['rg_m'][[0,-1]],
        Rg_change_m=float(m['rg_m'][-1]-m['rg_m'][0]), Rg_change_percent=float(100*(m['rg_m'][-1]/m['rg_m'][0]-1)),
        components_initial_final=m['components'][[0,-1]], largest_initial_final=m['largest_component'][[0,-1]],
        aspect_initial_final=m['aspect'][[0,-1]], max_fixed_ID_motion_m=zero_motion,
        min_distance_m=float(m['min_distance_m'].min()),
        pressure_legal=legal_pressure, pressure_range_Pa=[float(pmin.min()),float(p.max())],
        independent_step_gate_check=independently_checked,
        sampled_mean_force_N=float(np.linalg.norm(force,axis=-1).mean()), sampled_max_force_N=float(np.linalg.norm(force,axis=-1).max()),
        new_edge_episodes_minimum_persistence=persistence,
        events_count=sum(len(v) for v in events.values()),
        graph_scope='Distance <=2.1R, not fusion; newly formed events shorter than the persistence threshold are right-censored at trajectory end')
    assert legal_pressure
    if r.get('action', r.get('source')) == 'OFF': assert zero_motion==0
    return r, m, (t,x)

def compare(a, b):
    ra, ta, xa = load(a); rb, tb, xb = load(b); p=policy()
    if xa.shape[1:]!=xb.shape[1:] or not np.array_equal(xa[0],xb[0]):
        return dict(admitted=False,reason='INITIAL_STATE_OR_PARTICLE_COUNT_MISMATCH')
    j=np.searchsorted(tb,ta-1e-13)
    complete=len(j)>0 and np.all(j<len(tb)) and np.all(abs(tb[j]-ta)<1e-13)
    if not complete or ra['status']!='PASS' or rb['status']!='PASS' or abs(ta[-1]-tb[-1])>1e-13:
        return dict(admitted=False, reason='INCOMPLETE_MATCHED_WINDOW')
    ma, mb=series(ta,xa),series(ta,xb[j]); dx=xa-xb[j]
    ident=ma['axis_identifiable']&mb['axis_identifiable']
    angle=np.degrees(np.arccos(np.clip(abs(np.einsum('ti,ti->t',ma['principal_axis'],mb['principal_axis'])),0,1)))
    rms=[]
    for q,s in zip(ma['centered_positions_m'], mb['centered_positions_m']):
        cost=cdist(q,s,metric='sqeuclidean'); ii,jj=linear_sum_assignment(cost)
        rms.append(np.sqrt(cost[ii,jj].mean())/R)
    vals=dict(fixed_id_m=float(np.linalg.norm(dx,axis=-1).max()),
        COM_m=float(np.linalg.norm(ma['centroid_m']-mb['centroid_m'],axis=1).max()),
        Rg_m=float(abs(ma['rg_m']-mb['rg_m']).max()),
        Q_relative=float((np.linalg.norm(ma['tensor_m2']-mb['tensor_m2'],axis=(-2,-1))/np.linalg.norm(ma['tensor_m2'],axis=(-2,-1))).max()),
        aspect_relative=float((abs(ma['aspect']-mb['aspect'])/ma['aspect']).max()),
        directional_xy=float(abs(ma['directional_xy']-mb['directional_xy']).max()),
        matched_RMS_over_R=float(max(rms)), axis_angle_deg=float(angle[ident].max()) if ident.any() else None)
    limits=dict(fixed_id_m=p['fixed_id_tolerance_m'], COM_m=p['group_COM_tolerance_m'], Rg_m=p['group_Rg_tolerance_m'],
        Q_relative=p['group_Q_relative'], aspect_relative=p['group_aspect_relative'], directional_xy=p['group_xy_tolerance'],
        matched_RMS_over_R=p['distribution_RMS_over_R'], axis_angle_deg=p['principal_axis_angle_deg'])
    gates={k: vals[k] is None or vals[k]<=v for k,v in limits.items()}
    # Event-time checks use every accepted frame in each trajectory, not only common samples.
    ea=edge_events(ta,ma['edges']); eb=edge_events(tb,series(tb,xb)['edges']); same=set(ea)==set(eb); shifts=[]
    if same:
        for k in ea:
            aa,bb=ea[k],eb[k]
            if len(aa)!=len(bb) or [v[1] for v in aa]!=[v[1] for v in bb]: same=False; break
            shifts.extend(abs(u[0]-v[0]) for u,v in zip(aa,bb))
    event_error=float(max(shifts,default=0)) if same else None
    topology=bool(same and event_error<=p['instant_edge_time_tolerance_s'] and np.array_equal(ma['edges'][-1],mb['edges'][-1]))
    return dict(admitted=all(gates.values()), metric_gates=gates, errors=vals, limits=limits,
        complete_window_s=float(ta[-1]), common_frames=len(ta), event_sequences_match=same,
        max_event_time_difference_s=event_error, event_topology_agrees=topology,
        topology_scope='New dispersed-network events, distinct from survival of a pre-existing connected 32-bubble network')

def main():
    verify(); output=BASE/'analysis'; output.mkdir(exist_ok=True)
    runs={}; paths={}; summaries={}; metrics={}; arrays={}
    for path in sorted(BASE.glob('*/runs/*/RESULT.json')):
        r,m,tx=summarize(path.parent); runs[r['name']]=r; paths[r['name']]=path.parent
        metrics[r['name']]=m; arrays[r['name']]=tx
    checks={}
    for n in [64,128]:
        candidates=[r for r in runs.values() if r['N']==n and r['action']=='PAIR' and r['variant']=='DT1' and r['requested_s']>1e-4]
        if not candidates:
            summaries[str(n)]={'status':'FORMAL_NOT_RUN'}; continue
        r=max(candidates,key=lambda a:a['requested_s']); stem=r['name'][:-3]
        comp=[]
        for suffix in ['DT05','FINE']:
            name=stem+suffix
            comp.append(compare(paths[r['name']],paths[name]) if name in paths else dict(admitted=False,reason='NOT_RUN'))
        checks[str(n)]=comp
        off=next(q for q in runs.values() if q['N']==n and q['action']=='OFF' and q['requested_s']==r['requested_s'])
        summaries[str(n)]=dict(result=r, off=off, numeric_admitted=all(c['admitted'] for c in comp),
            topology_agrees=all(c.get('event_topology_agrees',False) for c in comp),
            Rg_difference_bound_m=max([c['errors']['Rg_m'] for c in comp if 'errors' in c],default=None))
    old=next(HISTORICAL.glob('*/runs/medium_train__PAIR__2000us__DT1/RESULT.json')).parent
    r,m,tx=summarize(old); summaries['32']=dict(result=r, source='Historical F10Z 2ms baseline; current 100us replay independently verified')
    metrics['historical32']=m; arrays['historical32']=tx
    campaign=next(BASE.glob('*/CAMPAIGN.json')); cost=read(campaign)
    write(output/'SUMMARY.json',dict(groups=summaries, runs=runs, comparisons=checks, campaign=cost,
        labels=LABELS, training_ready=False, single_realization=True))
    for name,m in metrics.items():
        np.savez_compressed(output/(name+'.npz'), **{k:v for k,v in m.items() if k not in ['edges','centered_positions_m']})
    from .figures import draw
    draw(summaries,runs,metrics,arrays)
    print({n:{k:v for k,v in s.items() if k not in ['result','off']} for n,s in summaries.items()})

if __name__=='__main__': main()
