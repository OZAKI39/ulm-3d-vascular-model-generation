"""Scientific classification from saved full-physics velocities and fixed IDs."""
import time
import numpy as np,h5py
from scipy.spatial.distance import pdist
from .common import *
from .metrics import spatial,rotation,endpoint

def one(r):
    p=BASE/r['trajectory_path'];out=BASE/'series';out.mkdir(exist_ok=True);numeric=r['status']=='PASS'
    with h5py.File(p) as h:
        t=h['time_s'][:];x=h['positions_m'];f=h['force_samples'];n=x.shape[1];ii,jj=np.triu_indices(n,1);om=[];rigid=[];comv=[];torque=[]
        for lo in range(0,len(t)-1,256):
            hi=min(lo+256,len(t)-1);q=x[lo:hi];v=f['velocity_m_s'][lo:hi];r0=q-q.mean(1)[:,None,:];v0=v-v.mean(1)[:,None,:];omega=np.cross(r0,v0)[:,:,2].sum(1)/np.sum(r0*r0,axis=(1,2));omega2=np.cross(r0,v0)[:,:,2].sum(1)/np.sum(r0[:,:,:2]**2,axis=(1,2));w=np.stack([-omega2[:,None]*r0[:,:,1],omega2[:,None]*r0[:,:,0],np.zeros((len(q),n))],2);energy=np.sum(v0*v0,axis=(1,2));frac=np.divide(energy-np.sum((v0-w)**2,axis=(1,2)),energy,out=np.zeros(len(q)),where=energy>1e-40);om.extend(omega);rigid.extend(frac);comv.extend(np.linalg.norm(v.mean(1),axis=1));force=f['total_force_N'][lo:hi];force-=force.mean(1)[:,None,:];torque.extend(np.cross(r0,force)[:,:,2].sum(1))
        om=np.asarray(om);rigid=np.asarray(rigid);cum=np.r_[0,np.cumsum(om*np.diff(t))];sgn=np.sign(cum[-1]);persist=float(np.mean(np.sign(om)==sgn)) if len(om) and sgn else 0.;idx=np.unique(np.rint(np.linspace(0,len(t)-1,min(101,len(t)))).astype(int));xx=x[idx];stats=[]
        for j,k in enumerate(idx):
            m=spatial(xx[j]);stats.append(m)
        a=stats[0];b=stats[-1];gain=b['stripe_order']-a['stripe_order'];null=read(BASE/'METRIC_NULL_CALIBRATION.json');threshold=next(z['admission_peak'] for z in null['rows'] if z['N']==n)
        stripe=bool(numeric and b['stripe_descriptor_pass'] and b['structure_peak_S']>=threshold and gain>=.15 and not a['stripe_descriptor_pass']);lattice=bool(stripe and b['lattice_descriptor_pass'] and b['second_peak_S']>=threshold);rot=bool(numeric and abs(np.degrees(cum[-1]))>=1 and persist>=.8 and np.mean(rigid)>=.8)
        q=xx[-1,:2]-xx[-1,:2].mean(0);rad=np.linalg.norm(q,axis=1);hist,xe,ye=np.histogram2d(xx[-1,:,0],xx[-1,:,1],bins=30,range=[[-.0005,.0005],[-.0004,.0004]]);rh,re=np.histogram(rad,bins=20,range=(0,.00065))
        np.savez_compressed(out/(r['id']+'.npz'),time_s=t[idx],positions_m=xx,Omega_z_s=om[np.minimum(idx,len(om)-1)],cumulative_angle_deg=np.degrees(cum[idx]),rigid_fraction=rigid[np.minimum(idx,len(rigid)-1)],full_rotation_time_s=t[:-1],full_Omega_z_s=om,full_rigid_fraction=rigid,relative_force_torque_z_N_m=torque,COM_velocity_m_s=comv,structure_peak_S=[z['structure_peak_S'] for z in stats],stripe_order=[z['stripe_order'] for z in stats],stripe_spacing_m=[z['stripe_spacing_m'] for z in stats],stripe_normal_deg=[z['stripe_normal_deg'] for z in stats],Rg_m=[z['Rg_3d_m'] for z in stats],aspect_xy=[z['xy_f10y']['aspect'] for z in stats],local_components=[z['local_components_20um'] for z in stats],covariance_3d_m2=[z['covariance_3d_m2'] for z in stats],center_m=[z['center_m'] for z in stats],histogram_2d=hist,histogram_x_edges_m=xe,histogram_y_edges_m=ye,radial_histogram=rh,radial_edges_m=re,initial_spectrum_S=a['spectrum_S'],final_spectrum_S=b['spectrum_S'],spectrum_freq_m1=b['spectrum_freq_m1'],spectrum_theta_rad=b['spectrum_theta_rad'])
    m=endpoint(xx[0],xx[-1]);return dict(id=r['id'],task=r['task'],status=r['status'],N=n,completed_s=r['completed_s'],wall_s=r['task_wall_s'],initial_peak=a['structure_peak_S'],final_peak=b['structure_peak_S'],stripe_order_gain=gain,periods=b['periods'],spacing_m=b['stripe_spacing_m'],stripe_candidate=stripe,lattice_candidate=lattice,rotation_candidate=rot,cumulative_rotation_deg=float(np.degrees(cum[-1])),rotation_sign_persistence=persist,mean_rigid_fraction=float(np.mean(rigid)),min_rigid_fraction=float(rigid.min()),final_COM_displacement_m=m['COM_displacement_m'],max_fixed_ID_displacement_m=m['max_fixed_id_displacement_m'],relative_fixed_ID_rms_m=m['relative_fixed_ID_rms_m'],final_Rg_m=m['Rg_3d_m'],initial_Rg_m=a['Rg_3d_m'],components=m['proximity_components'],largest_component=m['largest_component'],max_pressure_Pa=r['max_metrics'].get('max_pressure_Pa'),series_path=str((out/(r['id']+'.npz')).relative_to(BASE)),series_sha256=sha(out/(r['id']+'.npz')),raw_sha256=r['trajectory_sha256'],scope='CURRENT_MODEL_DYNAMICS; no experimental observability claim')
def main():
    begin=time.perf_counter();verify();runs=read(BASE/'ALL_HF_RUNS.json');rows=[one(r) for r in runs if r.get('trajectory_path') and r.get('accepted_steps',0)>0]
    reverse=[]
    for a in rows:
        if a['task'].get('action')!='ROT_PLUS' or a['status']!='PASS':continue
        for b in rows:
            if b['task'].get('action')=='ROT_MINUS' and b['status']=='PASS' and all(a['task'].get(k)==b['task'].get(k) for k in ['state','duration','frequency','offset_m','variant']):reverse.append(dict(plus=a['id'],minus=b['id'],angle_sign_reversal=a['cumulative_rotation_deg']*b['cumulative_rotation_deg']<0,both_rotation_confirmed=a['rotation_candidate'] and b['rotation_candidate'],plus_deg=a['cumulative_rotation_deg'],minus_deg=b['cumulative_rotation_deg']))
    pairsets=[read(BASE/k) for k in ['PAIRED_VALIDATION.json','EXTENSION_RESULTS.json','GENERALIZATION_RESULTS.json','LONG_TIME_RESULTS.json'] if (BASE/k).exists()];pairs=[];conv=[]
    for d in pairsets:pairs+=d.get('paired',d.get('rows',[]) if 'numeric_pass' in d else []);conv+=d.get('convergence',[])
    good=[r for r in rows if r['status']=='PASS'];fail=[dict(id=r['id'],reason=r.get('reason',r.get('reasons')),completed_s=r.get('completed_s')) for r in runs if r['status']!='PASS']
    screening=[]
    for state in ['N64_wide_r0','N64_wide_v0','N64_wide_v1','N64_wide_v2']:
        from scipy.stats import spearmanr
        aa=[];bb=[];signs=[]
        for r in runs:
            lp=BASE/'lf'/r['id']/'SUMMARY.json'
            if r['task']['state']!=state or r['task']['duration']!=.002 or r['status']!='PASS' or not lp.exists():continue
            l=read(lp)
            if l['status']!='PASS':continue
            a=r['final']['stripe_order']-r['initial']['stripe_order'];b=l['final']['stripe_order']-l['initial']['stripe_order'];aa.append(a);bb.append(b);signs.append(bool(np.sign(a)==np.sign(b)))
        screening.append(dict(state=state,paired_actions=len(aa),stripe_order_gain_rank=float(spearmanr(aa,bb).statistic),sign_agreement=all(signs),positive_pattern_sensitivity='UNIDENTIFIABLE without positive HF patterns'))
    write(BASE/'SCIENTIFIC_RESULTS.json',dict(screening=screening,rows=rows,reverse_phase=reverse,paired=pairs,convergence=conv,HF_attempts=len(runs),HF_complete=len(good),failures=fail,max_HF_time_s=max(r['completed_s'] for r in good),max_fixed_ID_m=max(r['max_fixed_ID_displacement_m'] for r in good),stripe_candidates=[r['id'] for r in good if r['stripe_candidate']],lattice_candidates=[r['id'] for r in good if r['lattice_candidate']],rotation_candidates=[r['id'] for r in good if r['rotation_candidate']],screening_false_negative_sensitivity='INCONCLUSIVE_NO_HF_POSITIVE_PATTERNS' if not any(r['stripe_candidate'] or r['lattice_candidate'] or r['rotation_candidate'] for r in good) else 'REQUIRES_INDEPENDENT_REFINEMENT',wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
