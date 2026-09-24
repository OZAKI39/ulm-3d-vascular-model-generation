#!/usr/bin/env python3
"""Independent complete-ledger conservation/distribution audit; no re-simulation."""
from pathlib import Path
import sys,csv,gzip,json,math
from collections import Counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy.stats import kstest,ks_2samp
from particle_3d.particle7_cases import REPO,SONOVUE
from particle_3d.particle3_cases import write_json,write_rows
from particle_3d.rbc_distribution import sample_rbc_geometries
from particle_3d.sonovue_adapter import read_sonovue


def main():
 folder=REPO/'particle_3d/reports/particle7/data'; summary=json.loads((folder/'11_long_summary.json').read_text()); real=json.loads((folder/'13_real_smoke.json').read_text())
 columns=['D_um','V_fL','r','jeffery_lambda','mb_diameter_um']; values={k:[] for k in columns}; volumes=[]; rbc_times=[]; rbc_exits=[]; mb_times=[]; counts=Counter(); modes=Counter(); last_id=0; last_time=0.; wall_min=math.inf
 with gzip.open(folder/'11_long_lifecycle.csv.gz','rt',newline='') as stream:
  for row in csv.DictReader(stream):
   tag=int(row['particle_id']); t=float(row['scheduled_time_s']); sp=row['species']; status=row['final_status']
   assert tag==last_id+1 and t>=last_time and row['scheduled_time_s']==row['admitted_time_s']
   assert float(row['z_m'])==0 and status in ['ACTIVE','EXITED']
   if status=='EXITED': assert t<float(row['exit_time_s'])<=summary['horizon_s'] and row['outlet_role']=='OUTLET_01'
   else: assert row['exit_time_s']==row['outlet_role']==''
   margin=float(row['wall_path_lower_bound_m']); assert margin>=0
   if sp=='MB': assert margin>=max(2e-9,.001*float(row['mb_diameter_um'])*.5e-6)
   wall_min=min(wall_min,margin); counts[sp]+=1; counts[sp+'_'+status]+=1; modes[row['shape_mode']]+=1
   assert abs(sum(float(row[q])**2 for q in ['qw','qx','qy','qz'])-1)<2e-14
   if sp=='RBC':
    volumes.append(float(row['volume_m3'])); rbc_times.append(t); rbc_exits.append(float(row['exit_time_s']) if row['exit_time_s'] else math.inf)
   else: mb_times.append(t)
   for key in columns:
    if row[key]!='': values[key].append(float(row[key]))
   last_id=tag; last_time=t
 for sp in ['MB','RBC']:
  key=sp.lower(); assert counts[sp]==summary['final']['scheduled_'+key+'_count']==summary['final']['admitted_'+key+'_count']
  for state in ['ACTIVE','EXITED']: assert counts[sp+'_'+state]==summary['final'][state.lower()+'_'+key+'_count']
 total=math.fsum(volumes); assert abs(total-summary['final']['scheduled_rbc_volume_m3'])<1e-25
 prefix=np.cumsum(np.array(volumes,dtype=np.longdouble)); error=float(np.max(np.abs(np.array(rbc_times,dtype=np.longdouble)*np.longdouble(.45*summary['Q_m3_s'])-prefix)))
 assert error<1e-24
 expected=np.arange(1,len(mb_times)+1)/(8.5e12*summary['Q_m3_s']); assert np.array_equal(mb_times,expected)
 _,dist,_=read_sonovue(SONOVUE); mbks=kstest(values['mb_diameter_um'],dist.cdf)
 reference=sample_rbc_geometries(20000,seed=2026092117).samples
 tests={key:ks_2samp(values[key],reference[{'jeffery_lambda':'jeffery_lambda'}.get(key,key)]) for key in columns if key!='mb_diameter_um'}
 assert mbks.pvalue>1e-5 and all(r.pvalue>1e-5 for r in tests.values())
 edges=dict(D_um=np.linspace(4,9.58,45),V_fL=np.linspace(26.6324,69.1676,45),r=np.linspace(0,1,45),jeffery_lambda=np.linspace(-1,0,45),mb_diameter_um=np.linspace(.0,12.,61))
 out=dict(long={},real={},independent_original_sampler_comparison=dict(MB=dict(statistic=mbks.statistic,p_value=mbks.pvalue),
   RBC={k:dict(statistic=v.statistic,p_value=v.pvalue,reference_count=20000) for k,v in tests.items()}),
   scheduled_equals_admitted_exact_in_long=True,real_admission_bias='INSUFFICIENT_ADMITTED_SAMPLE_FOR_POPULATION_INFERENCE; KEEP_SCHEDULED_ADMITTED_PENDING_SEPARATE')
 for key in columns:
  bins=edges[key]; hist=np.histogram(values[key],bins=bins)[0]; assert hist.sum()==len(values[key])
  out['long'][key]=dict(edges=bins,scheduled=hist,admitted=hist.copy(),pending=np.zeros_like(hist),mean=float(np.mean(values[key])))
  groups=dict(scheduled=real['events'],admitted=list(real['births'].values()),pending=real['pending']['RBC']+real['pending']['MB']); group_hist={}
  for name,events in groups.items():
   sample=[]
   for e in events:
    if key=='mb_diameter_um':
     if e['species']=='MB': sample.append(e['diameter_um'])
    elif e['species']=='RBC':
     geom=e['geometry']; ratio=geom['c_m']/geom['a_m']
     sample.append(dict(D_um=e['provenance']['D_um'],V_fL=e['provenance']['V_fL'],r=ratio,jeffery_lambda=(ratio*ratio-1)/(ratio*ratio+1))[key])
   group_hist[name]=np.histogram(sample,bins=bins)[0]; assert sum(group_hist[name])==len(sample)
  out['real'][key]=dict(edges=bins,**group_hist)
 times=np.linspace(0,summary['horizon_s'],1201); births=np.searchsorted(rbc_times,times,side='right'); exits=np.searchsorted(rbc_exits,times,side='right')
 cumulative=np.r_[np.longdouble(0),prefix]; volume_at_time=cumulative[births]-cumulative[exits]; control=summary['width_m']**2*summary['length_m']
 hct=np.asarray(volume_at_time/control,dtype=float)
 write_rows(folder/'14_uniform_time_tube_hct.csv',[dict(time_s=t,active_rbc_count=int(i-j),active_rbc_volume_m3=float(v),control_volume_m3=control,observed_tube_hct=h) for t,i,j,v,h in zip(times,births,exits,volume_at_time,hct)])
 residence=np.minimum(summary['length_m']/summary['velocity_m_s'][2],summary['horizon_s']-np.array(rbc_times))
 time_mean=float(np.sum(np.array(volumes,dtype=np.longdouble)*residence)/(control*summary['horizon_s']))
 write_json(folder/'14_tube_hct_statistics.json',dict(sampling='UNIFORM_PHYSICAL_TIME_NOT_BIRTH_CONDITIONED',uniform_samples=len(times),
  exact_residence_time_weighted_mean=time_mean,uniform_snapshot_mean=float(hct.mean()),uniform_snapshot_min=float(hct.min()),uniform_snapshot_max=float(hct.max()),
  final_snapshot=summary['final']['observed_tube_hct'],real_lumen_snapshot=real['observed_local_tube_hct'],feed_H_D=.45,enforced=False,
  volume_definition='CENTER_ASSIGNED_WHOLE_RBC_VOLUME; NOT_SHAPE_LUMEN_INTERSECTION_VOLUME',interpretation='COMMON_PLUG_VALIDATION_ONLY_NOT_PHYSIOLOGICAL_TUBE_HCT'))
 write_json(folder/'12_distribution_diagnostics.json',out)
 write_json(folder/'11_independent_ledger_audit.json',dict(status='PASS',complete_csv_rows=last_id,counts=dict(counts),shape_modes=dict(modes),
   duplicate_ids=0,lost_particles=0,all_birth_centers_on_inlet=True,wall_path_minimum_lower_bound_m=wall_min,
   rbc_volume_fsum_m3=total,max_rbc_threshold_reconstruction_error_m3=error,mb_times_exact=True,
   volume_balance_error_m3=total-summary['final']['scheduled_rbc_volume_m3']))
 print('Complete ledger audit PASS:',last_id,'rows',counts,flush=True)
if __name__=='__main__':main()
