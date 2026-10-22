"""Attach original solver column meanings and all-step pressure diagnostics."""
import h5py,numpy as np
from mevo_mb.fast_sim.f10g_contact import resistance
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import *
def main():
 runs=read(BASE/'analysis/RUN_ANALYSIS.json')
 for name,r in runs.items():
  raw=ROOT/r['path']/'trajectory.h5';path=BASE/'dataset'/(name+'.h5')
  with h5py.File(raw) as source,h5py.File(path,'a') as target:
   for k in ['step_metrics','step_time_duration_segment']:
    for a,v in source[k].attrs.items():target[k].attrs[a]=v
   target['step_time_duration_segment'].attrs['columns']='step_start_s,actual_dt_s,segment_id'
   if 'extra_metrics' not in target:target.create_dataset('extra_metrics',data=source['extra_metrics'][:],compression='gzip')
   for a,v in source['extra_metrics'].attrs.items():target['extra_metrics'].attrs[a]=v
   target.attrs['training_ready']=False;target.attrs['labels_json']=__import__('json').dumps(LABELS)
   target.attrs['raw_path_relative_to_project']=str(raw.relative_to(ROOT))
   fs=target['force_samples']
   physics=F2Physics(radius_m=R,frequency_Hz=2.5e6);zeta=6*np.pi*physics.viscosity_Pa_s*R
   target.attrs['translation_viscosity_Pa_s']=physics.viscosity_Pa_s;target.attrs['radial_effective_viscosity_Pa_s']=physics.radial_effective_viscosity_Pa_s
   times=target['time_s'][:];x=target['positions_m'][:];drag=[];ts=[];balance=[]
   for i,tau in enumerate(fs['time_s'][:]):
    j=int(np.searchsorted(times,tau-1e-13))
    if j>=len(times)-1:continue
    v=(x[j+1]-x[j])/(times[j+1]-times[j]);f=-zeta*(resistance(x[j],R)['A']@v.ravel()).reshape(32,3);drag.append(f);ts.append(tau);balance.append(float(np.linalg.norm(f+fs['total_force_N'][i]+fs['contact_force_N'][i],axis=1).max()))
   for key,values in [('derived_drag_N',drag),('derived_drag_time_s',ts)]:
    if key in fs:del fs[key]
    fs.create_dataset(key,data=values,compression='gzip')
   r['force_balance_sample_max_N']=max(balance,default=None)
   assert r['force_balance_sample_max_N']<=r['max_force_residual_N']+1e-18, (name,r['force_balance_sample_max_N'])
   fs.attrs['contact_force_semantics']='Original solver reaction for the step starting at force_samples/time_s'
   fs.attrs['derived_drag_semantics']='Original -zeta A(x) v; v is the actual accepted interval average; derived_drag_time_s is separate'
  r['dataset_sha256']=sha(path)
 write(BASE/'analysis/RUN_ANALYSIS.json',runs)
 print('Preserved original numerical diagnostic columns and complete pressure-range evidence.')
if __name__=='__main__':main()

