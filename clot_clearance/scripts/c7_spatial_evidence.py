"""Physical-space diagnostics and reduced residual Jv gates; no Hessians."""
import csv,time
import numpy as np
from scripts import c7_early_branch as a
from scripts import c7_reduced_vi as reduced


def region(p,ids,amplitude=None):
    ids=np.asarray(ids,int);f=p.family;L=a.old.v2.L
    if not len(ids):return dict(count=0,fraction=0.,centroid_m=None,RMS_radius_m=None,nearest_clamp_percentiles_m=None,equivalent_volume_m3=0.)
    end=f.X[f.pairs[ids]];mid=end.mean(axis=1);ori=f.xi[ids]/f.length[ids,None]
    amp=np.ones(len(ids)) if amplitude is None else np.asarray(amplitude,float)
    w=p.weights[ids]*amp;center=np.average(mid,weights=w,axis=0);radius=float(np.sqrt(np.average(np.sum((mid-center)**2,axis=1),weights=w)))
    clamp=L/2-abs(mid[:,2]);lateral=L/2-np.max(abs(mid[:,:2]),axis=1);angle=np.degrees(np.arccos(np.clip(abs(ori[:,2]),0,1)))
    values=np.zeros(p.ns);values[ids]=amp;aggregate=a.old.v2.aggregate(f,values);V=f.volume
    effective=float((V@aggregate)**2/(V@(aggregate**2)))
    return dict(count=len(ids),fraction=len(ids)/p.ns,centroid_m=center.tolist(),RMS_radius_m=radius,
        midpoint_bbox_min_m=mid.min(axis=0).tolist(),midpoint_bbox_max_m=mid.max(axis=0).tolist(),midpoint_bbox_dimensions_m=np.ptp(mid,axis=0).tolist(),
        endpoint_bbox_dimensions_m=np.ptp(end.reshape(-1,3),axis=0).tolist(),nearest_clamp_percentiles_m=np.percentile(clamp,[0,25,50,75,100]).tolist(),
        top_distance_percentiles_m=np.percentile(L/2-mid[:,2],[0,50,100]).tolist(),bottom_distance_percentiles_m=np.percentile(L/2+mid[:,2],[0,50,100]).tolist(),
        lateral_distance_percentiles_m=np.percentile(lateral,[0,50,100]).tolist(),boundary_fraction=float(np.mean(clamp<=a.cfg()['localization']['boundary_distance_fraction']*L)),
        orientation_angle_to_axis_deg_edges=[0,15,30,45,60,75,90],orientation_counts=np.histogram(angle,[0,15,30,45,60,75,90])[0].tolist(),
        orientation_abs_mean=np.average(abs(ori),weights=w,axis=0).tolist(),equivalent_volume_m3=effective,equivalent_volume_fraction=effective/L**3,
        endpoint_support_volume_m3=float(V[np.unique(f.pairs[ids])].sum()),
        equivalent_volume_definition='(sum V*a)^2 / sum V*a^2 for particle aggregate a; discrete participation volume, not crack width',
        particle_aggregate=aggregate.tolist())


def pattern(r):
    t=a.cfg()['localization'];L=a.old.v2.L
    if r['count']==0:return 'INCONCLUSIVE'
    if r['boundary_fraction']>=t['boundary_fraction']:return 'BOUNDARY_ASSOCIATED'
    if r['fraction']<=t['localized_fraction'] and r['RMS_radius_m']<=t['localized_RMS_fraction']*L:return 'LOCALIZED'
    if r['fraction']<=t['partially_localized_fraction']:return 'PARTIALLY_LOCALIZED'
    return 'DIFFUSE'


def activation_spatial():
    weak=a.read(a.OUT/'activation/WEAK.json');bracket=a.read(a.OUT/'activation/BRACKET.json');state=a.arrays(a.OUT/f'intact/{weak["label"]}/STATE.npz');p=a.problem(7,weak['strain'])
    ids=np.flatnonzero(abs(state['q_scaled'])<=a.cfg()['q_tol']);r=region(p,ids);r['pattern']=pattern(r)
    f=p.family;end=f.X[f.pairs[ids]];mid=end.mean(axis=1);ori=f.xi[ids]/f.length[ids,None];rows=[]
    for k,b in enumerate(ids):
        rows.append(dict(ID=int(b),endpoint_0_m=end[k,0].tolist(),endpoint_1_m=end[k,1].tolist(),midpoint_m=mid[k].tolist(),orientation=ori[k].tolist(),reference_length_m=float(f.length[b]),Psi_plus=float(state['Y'][b]),q=float(state['q'][b]),top_clamp_distance_m=float(a.old.v2.L/2-mid[k,2]),bottom_clamp_distance_m=float(a.old.v2.L/2+mid[k,2]),lateral_distance_m=float(a.old.v2.L/2-max(abs(mid[k,:2]))),cube_center_distance_m=float(np.linalg.norm(mid[k]))))
    first_violated=bracket['upper']['q_IDs']['VIOLATED_LOWER']
    r.update(strain=weak['strain'],IDs=ids.tolist(),first_weak_IDs=ids.tolist(),first_violated_IDs=first_violated,first_violated_strain=bracket['upper']['strain'],bonds=rows,title='C7 AT1 ACTIVATION SET — NOT CRACK PREDICTION')
    a.save(a.OUT/'activation/SPATIAL.json',r)
    np.savez_compressed(a.RESULTS/'ACTIVATION_SET.npz',IDs=ids,endpoints=end,midpoints=mid,orientation=ori,particle_aggregate=r['particle_aggregate'])
    oldset=next(v for v in a.read(a.old.OUT/'activation/SPATIAL.json') if v['role']=='activation');c5=a.problem(5,oldset['strain']);c5ids=[v['ID'] for v in oldset['candidates']];r5=region(c5,c5ids)
    # Historical C5 positions/orientations are checked directly before comparing.
    for row in oldset['candidates']:
        b=row['ID'];assert np.array_equal(c5.family.X[c5.family.pairs[b]].mean(axis=0),row['midpoint_m']);assert np.allclose(c5.family.xi[b]/c5.family.length[b],row['orientation'],rtol=0,atol=1e-15)
    r5.update(strain=oldset['strain'],IDs=c5ids,pattern=pattern(r5),source=str(a.old.OUT/'activation/SPATIAL.json'))
    shift=weak['strain']-r5['strain'];same_axis=np.allclose(r['orientation_abs_mean'],r5['orientation_abs_mean'],atol=.05)
    relation='CONSISTENT' if r['pattern']==r5['pattern']=='BOUNDARY_ASSOCIATED' and same_axis and np.linalg.norm(np.array(r['centroid_m'])-r5['centroid_m'])<.1*a.old.v2.L else 'INCONCLUSIVE'
    result=dict(C5=r5,C7=r,strain_difference=shift,strain_relative_difference=shift/r5['strain'],strain_relation='SHIFTING' if abs(shift)>1e-7 else 'INCONCLUSIVE',region_relation=relation,scope='Two-grid discrete physical-region comparison only; no convergence or C7 stability claim')
    a.save(a.OUT/'activation/C5_C7_COMPARISON.json',result);return result


def jacobian_gate(p,x,s,ids,pool,label):
    start=time.perf_counter();rp=reduced.ReducedProblem(p,ids);u=rp.compress(p.pack(x,s));J,rss=reduced.jacobian(rp,u,pool,a.cfg()['workers']);rng=np.random.default_rng(a.cfg()['jacobian_seed']);r0=reduced.vi.residual(rp,u);checks=[]
    for k in range(a.cfg()['jacobian_directions']):
        v=rng.normal(size=len(u));v[p.nx:]=abs(v[p.nx:]);v/=np.linalg.norm(v);jv=J@v
        for h in a.cfg()['jacobian_steps']:
            # Legal second-order forward probe avoids extrapolating phase below
            # zero. All derivatives remain those of raw, unprojected q/D.
            fd=(-3*r0+4*reduced.vi.residual(rp,u+h*v)-reduced.vi.residual(rp,u+2*h*v))/(2*h)
            checks.append(dict(direction=k,epsilon=h,relative_error=float(np.linalg.norm(fd-jv)/max(np.linalg.norm(jv),1e-30))))
    result=dict(label=label,rows=checks,max_relative_error=max(v['relative_error'] for v in checks),tolerance=a.cfg()['jacobian_tolerance'],scheme='Second-order feasible forward residual difference',pool_size=len(ids),shape=list(J.shape),worker_peak_RSS_bytes=rss,wall_seconds=time.perf_counter()-start)
    result['passed']=result['max_relative_error']<=result['tolerance'];a.save(a.OUT/f'jacobian/{label}.json',result)
    assert result['passed'],'REDUCED_JACOBIAN_GATE_FAILED';return result


def c5_gate(pool):
    r=a.read(a.OUT/'c5_regression/ALL.json');assert r['passed']
    for label in ['first_A1','same_load_as_original']:
        folder=a.OUT/f'c5_regression/{label}';row=a.read(folder/'SUMMARY.json');s=a.arrays(folder/'STATE.npz');p=a.problem(5,row['strain'],s['lower'])
        jacobian_gate(p,s['x'],s['s'],s['pool'],pool,'C5_'+label)


def branch_spatial(row,data):
    p=a.problem(7,row['strain'],data['lower']);ids=np.flatnonzero(data['s']>a.cfg()['phase_detection_tol'])
    r=region(p,ids,data['s'][ids]);a.save(a.OUT/row['label']/'SPATIAL.json',r);return r


if __name__=='__main__':
    import sys
    if sys.argv[1]=='spatial':activation_spatial()
    else:
        pool=a.prior.vi.make_pool(a.cfg()['workers'])
        try:c5_gate(pool)
        finally:pool.shutdown()
