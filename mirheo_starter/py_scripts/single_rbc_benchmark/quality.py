"""Stored-state numerical checks. No solver imports or launches.

Intersection and containment checks apply to saved frames, not continuous time.
They provide bounded geometric evidence, not a strict permeability proof.
"""
from pathlib import Path
import numpy as np


def _segment_triangle(p,q,t,tol=1e-9):
    direction=q-p;e1=t[:,1]-t[:,0];e2=t[:,2]-t[:,0]
    h=np.cross(direction,e2);det=np.einsum('ij,ij->i',e1,h)
    inv=np.divide(1.,det,out=np.zeros_like(det),where=np.abs(det)>tol)
    s=p-t[:,0];u=np.einsum('ij,ij->i',s,h)*inv
    qq=np.cross(s,e1);v=np.einsum('ij,ij->i',direction,qq)*inv
    along=np.einsum('ij,ij->i',e2,qq)*inv
    return (np.abs(det)>tol)&(u>=-tol)&(v>=-tol)&(u+v<=1+tol)&(along>=-tol)&(along<=1+tol)


def mesh_intersections(vertices,faces):
    """AABB broad phase, nonadjacent triangle edges, coplanar separating axes."""
    tri=np.asarray(vertices)[faces];lo=tri.min(axis=1);hi=tri.max(axis=1)
    overlap=np.all(lo[:,None]<=hi[None,:]+1e-9,axis=2)&np.all(hi[:,None]>=lo[None,:]-1e-9,axis=2)
    i,j=np.where(np.triu(overlap,1))
    independent=~np.any(faces[i,:,None]==faces[j,None,:],axis=(1,2));i=i[independent];j=j[independent]
    a,b=tri[i],tri[j];hit=np.zeros(len(i),bool)
    for k in range(3):
        hit|=_segment_triangle(a[:,k],a[:,(k+1)%3],b)
        hit|=_segment_triangle(b[:,k],b[:,(k+1)%3],a)
    if len(i):
        normal=np.cross(a[:,1]-a[:,0],a[:,2]-a[:,0]);norm=np.linalg.norm(normal,axis=1)
        distance=np.max(np.abs(np.einsum('nkj,nj->nk',b-a[:,:1],normal)),axis=1)
        for k in np.where((distance<=1e-9*np.maximum(norm,1e-30))&~hit)[0]:
            axes=[x for x in range(3) if x!=np.argmax(np.abs(normal[k]))]
            aa=a[k][:,axes];bb=b[k][:,axes];separated=False
            for polygon in (aa,bb):
                for edge in np.roll(polygon,-1,axis=0)-polygon:
                    axis=np.array([-edge[1],edge[0]]);pa=aa@axis;pb=bb@axis
                    if pa.max()<pb.min()-1e-9 or pb.max()<pa.min()-1e-9:separated=True
            hit[k]=not separated
    return dict(nonadjacent_intersections=int(hit.sum()),candidate_pairs=len(i),first_pairs=np.column_stack([i[hit],j[hit]])[:12].tolist(),
                scope='saved frames; pairs sharing a vertex excluded; not a continuous collision proof')


def points_inside(points,vertices,faces):
    """Closed-surface solid-angle winding test, chunked to bound memory."""
    tri=np.asarray(vertices)[faces];answer=[]
    for start in range(0,len(points),32):
        q=tri[None]-np.asarray(points)[start:start+32,None,None,:]
        a,b,c=q[:,:,0],q[:,:,1],q[:,:,2];la=np.linalg.norm(a,axis=2);lb=np.linalg.norm(b,axis=2);lc=np.linalg.norm(c,axis=2)
        numerator=np.einsum('nki,nki->nk',a,np.cross(b,c))
        denominator=la*lb*lc+np.einsum('nki,nki->nk',a,b)*lc+np.einsum('nki,nki->nk',b,c)*la+np.einsum('nki,nki->nk',c,a)*lb
        angle=np.sum(2*np.arctan2(numerator,denominator),axis=1)
        answer.extend(np.abs(angle)>2*np.pi)
    return np.asarray(answer,bool)


def leakage_probes(directory,frames,faces,c,completion):
    rows=[];lookup={(f['phase'],f['step']):f for f in frames}
    for path in sorted(Path(directory).glob('fluid_probes_*.npz')):
        _,_,phase,done=path.stem.split('_');step=int(done)+(completion.get('prep_steps',0) if phase=='shear' else 0)
        frame=lookup.get((phase,step))
        if frame is None:continue
        with np.load(path) as raw:probe=raw['probes']
        vertices=np.array(frame['vertices']);points=probe[:,2:5].copy();center=vertices.mean(axis=0);L=c['geometry']['periodic_length']
        points[:,:2]+=L*np.rint((center[:2]-points[:,:2])/L)
        inside=points_inside(points,vertices,faces);expected=probe[:,0]==1;mismatch=inside!=expected
        rows.append(dict(phase=phase,step=step,strain=frame['strain'],tested=len(probe),mismatches=int(mismatch.sum()),
                         mismatch_ids=probe[mismatch,1].astype(int).tolist(),file=path.name))
    return dict(status=('NO_MISMATCH_IN_SAMPLED_PROBES' if all(x['mismatches']==0 for x in rows) else 'SAMPLED_MEMBERSHIP_MISMATCH') if rows else 'NOT_MEASURED',
                tested_point_frames=sum(x['tested'] for x in rows),mismatches=sum(x['mismatches'] for x in rows),frames=rows,
                limitation='At most 96 fixed IDs per side per coordinator at output times. Native reclassification at the stage handoff is disclosed; this is not strict impermeability verification.')


def fluid_summary(directory,c,reader):
    p=Path(directory)/'profiles.csv'
    if not p.exists() or p.stat().st_size<120:return None
    a=reader(p);shear=a[a['phase']=='shear'];result={}
    if not len(shear):return result
    selected=shear[shear['strain']>=2.-1e-8]
    if not len(selected):selected=shear[shear['step']==shear['step'].max()]
    bins=np.unique(selected['bin']);mean={k:np.array([selected[selected['bin']==b][k].mean() for b in bins]) for k in ('z','ux','uy','uz','density')}
    z=mean['z'];expected=c['protocol']['shear_rate']*(z-c['geometry']['gap']/2)
    rate,intercept=np.polyfit(z,mean['ux'],1)
    mass=[]
    for phase,step in dict.fromkeys(zip(a['phase'].tolist(),a['step'].tolist())):
        q=a[(a['phase']==phase)&(a['step']==step)];mass.append(float(np.average(q['density'],weights=q['count'])) if not (Path(directory)/'moments.csv').exists() else float(q['count'].sum()))
    result.update(mean_profile={k:v.tolist() for k,v in mean.items()},window_strain=[float(selected['strain'].min()),float(selected['strain'].max())],
                  profile_relative_l2=float(np.linalg.norm(mean['ux']-expected)/max(np.linalg.norm(expected),1e-30)),
                  effective_shear=float(rate),effective_shear_relative_error=float(abs(rate/c['protocol']['shear_rate']-1)),
                  fitted_wall_velocities=[float(intercept),float(intercept+rate*c['geometry']['gap'])],
                  near_wall_density_relative_error=float(max(abs(mean['density'][[0,-1]]/(c['dpd']['number_density']*c['dpd']['mass'])-1))),
                  mass_relative_drift=float(max(abs(np.array(mass)/mass[0]-1))))
    mp=Path(directory)/'moments.csv'
    if mp.exists():
        m=reader(mp);ms=m[(m['phase']=='shear')&(m['strain']>=selected['strain'].min()-1e-8)]
        result.update(temperature_mean=float(ms['temperature'].mean()),temperature_relative_error=float(abs(ms['temperature'].mean()/c['dpd']['kBT']-1)),
                      temperature_observed_range=[float(ms['temperature'].min()),float(ms['temperature'].max())],
                      particle_count_range=[int(m['N'].min()),int(m['N'].max())],inner_count_range=[int(m['inner_N'].min()),int(m['inner_N'].max())],
                      wall_crossings_max=int(m['wall_crossings'].max()),max_particle_speed=float(m['max_speed'].max()))
    else:result.update(temperature_status='LBM_WITHOUT_THERMAL_FLUCTUATIONS',wall_crossings_status='NO_PARTICLE_TRACKING_IN_LBM')
    return result


def local_flow_summary(directory,reader,window=(2.,4.)):
    p=Path(directory)/'local_flow.csv'
    if not p.exists():return None
    a=reader(p)
    if not len(a):return None
    a=a[(a['phase']=='shear')&(a['strain']>=window[0]-1e-8)&(a['strain']<=window[1]+1e-8)]
    if not len(a):return None
    keys=list(dict.fromkeys(zip(a['ix'].tolist(),a['iy'].tolist(),a['iz'].tolist())))
    rows=[]
    for ix,iy,iz in keys:
        r=a[(a['ix']==ix)&(a['iy']==iy)&(a['iz']==iz)]
        blocks=[]
        for lo,hi in zip(np.linspace(window[0],window[1],5)[:-1],np.linspace(window[0],window[1],5)[1:]):
            q=r[(r['strain']>=lo)&(r['strain']<(hi if hi<window[1] else hi+1e-8))]
            if len(q):blocks.append([float(q[k].mean()) for k in ('ux','uy','uz')])
        rows.append(dict(cell=[int(ix),int(iy),int(iz)],position=[float(r[k][0]) for k in ('x','y','z')],velocity=[float(r[k].mean()) for k in ('ux','uy','uz')],
                         block_std=np.std(blocks,axis=0,ddof=1).tolist() if len(blocks)>1 else None))
    return dict(window=[float(a['strain'].min()),float(a['strain'].max())],requested_window=list(window),cells=rows,definition='Time mean of native local velocities in shared 8^3 bins; HemoCell node average and DPD particle average differ. Block scatter is descriptive, not a precise confidence interval.')


def feedback_comparison(cell,empty,a):
    if not cell or not empty:return dict(status='NOT_MEASURED')
    if not np.allclose(cell['window'],empty['window'],atol=1e-8):return dict(status='INCOMPATIBLE_SAMPLING_WINDOWS',cell_window=cell['window'],empty_window=empty['window'])
    e={tuple(r['cell']):r for r in empty['cells']};delta=[];noise=[];near=[];plot=[]
    for r in cell['cells']:
        ref=e[tuple(r['cell'])];dv=np.array(r['velocity'])-ref['velocity'];delta.append(dv)
        near.append(np.linalg.norm(np.array(r['position'])-12)<=2*a)
        noise.append(np.array(r['block_std'])**2+np.array(ref['block_std'])**2 if r['block_std'] and ref['block_std'] else [np.nan]*3)
        if r['cell'][1]==4:plot.append(dict(position=r['position'],delta_velocity=dv.tolist()))
    delta=np.array(delta);near=np.array(near);noise=np.array(noise)
    scatter=float(np.sqrt(np.mean(np.sum(noise[near],axis=1))))
    return dict(status='ACTUAL_WITH_CELL_MINUS_EMPTY_FLOW',rms_velocity_difference=float(np.sqrt(np.mean(np.sum(delta**2,axis=1)))),
                near_cell_rms_difference=float(np.sqrt(np.mean(np.sum(delta[near]**2,axis=1)))),
                near_cell_combined_block_scatter=scatter if np.isfinite(scatter) else None,slice_y_index=4,slice=plot,
                interpretation='Actual flow changes support the source-audited feedback path; DPD differences also contain thermal noise. This subtraction alone is not a noise-free causal force measurement.')


def dynamics_window(run,window):
    frames=[f for f in run.get('frames',[]) if f['phase']=='shear' and window[0]-1e-8<=f['strain']<=window[1]+1e-8]
    if not frames:return None
    D=np.array([f['metrics']['D'] for f in frames]);angles=[f['metrics']['theta_deg'] for f in frames if f['metrics']['theta_deg'] is not None]
    theta=float(np.degrees(np.angle(np.mean(np.exp(2j*np.radians(angles)))))/2) if angles else None
    return dict(window=list(window),samples=len(frames),D_mean=float(D.mean()),D_std=float(D.std()),D_range=[float(D.min()),float(D.max())],theta_axial_mean_deg=theta,
                invalid_angle_frames=len(frames)-len(angles),center_start=frames[0]['metrics']['center'],center_end=frames[-1]['metrics']['center'])


def strict_comparison(primary,strict,c,window):
    p=dynamics_window(primary,window);s=dynamics_window(strict,window)
    if p is None or s is None:return dict(status='NOT_MEASURED',primary=p,strict=s)
    expected=np.arange(window[0],window[1]+c['protocol']['sample_strain']/2,c['protocol']['sample_strain'])
    def covers(run):
        actual=np.array([f['strain'] for f in run['frames'] if f['phase']=='shear'])
        return all(np.any(np.abs(actual-g)<1e-7) for g in expected)
    covered=covers(primary) and covers(strict)
    difference=abs(s['D_mean']-p['D_mean']);limit=max(c['screen']['own_strict_D_absolute_noise_floor'],c['screen']['own_strict_D_relative_difference']*abs(p['D_mean']))
    theta=abs((s['theta_axial_mean_deg']-p['theta_axial_mean_deg']+90)%180-90) if s['theta_axial_mean_deg'] is not None and p['theta_axial_mean_deg'] is not None else None
    return dict(status=('PASS' if difference<=limit else 'FAILED') if covered else 'INCOMPLETE_WINDOW',window=list(window),primary=p,strict=s,
                D_absolute_difference=difference,D_relative_difference=difference/max(abs(p['D_mean']),1e-30),frozen_absolute_tolerance=limit,
                theta_mean_difference_deg=theta,scope='Mean D over the same strain window, no time warp or rotation. This tests only sampled temporal sensitivity; spatial convergence and full nonlinear model equivalence remain unverified.')
