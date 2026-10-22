"""Read-only raw-trajectory analysis, original independent KKT checks, explicit noise."""
import argparse,time,json,h5py,numpy as np
from scripts.straight_tube_exp1.f10aa.analyze import independent,extra_compare
from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
from .metrics import series,measure
from .common import *

def describe(m):
 # Keep variables with distinct units; principal angle only when identifiable.
 c=m['center_3d_m'][-1]-m['center_3d_m'][0]
 q=m['covariance_3d_m2'];g=m['Rg_3d_m']
 return dict(com_um=c*1e6,rg_change_um=float((g[-1]-g[0])*1e6),rg_fraction=float(g[-1]/g[0]-1),
   aspect_change=float(m['aspect_3d'][-1]-m['aspect_3d'][0]),aspect_xy_change=float(m['aspect_xy'][-1]-m['aspect_xy'][0]),
   aspect_fraction=float(m['aspect_3d'][-1]/m['aspect_3d'][0]-1),
   rigid_rotation_deg=float(np.degrees(np.arccos(np.clip((np.trace(m['rigid_rotation_3d'][-1])-1)/2,-1,1)))),
   similarity_residual_nm=float(m['similarity_rmsd_3d_m'][-1]*1e9),kabsch_rmsd_nm=float(m['kabsch_rmsd_3d_m'][-1]*1e9),
   nonaffine_nm=float(m['nonaffine_residual_3d_m'][-1]*1e9),deviatoric=float(m['deviatoric_stretch_3d'][-1]),
   signed_xy=float(m['directional_xy_stretch_3d'][-1]),shear=float(m['shear_norm_3d'][-1]),
   isotropic_stretch=float(m['isotropic_stretch_3d'][-1]),components=[int(m['proximity_components'][k]) for k in [0,-1]],
   largest=[int(m['largest_component'][k]) for k in [0,-1]],contact_peak=int(m['hard_contact_pairs'].max()),
   affine_valid=bool(m['affine_valid_3d'].all()),orientation_identifiable=bool(m['orientation_valid_3d'][-1,0]))

def analyze_stage(stage):
 start=time.perf_counter();out=BASE/'analysis'/stage;out.mkdir(parents=True,exist_ok=True);summaries=[]
 for record in read(BASE/stage/'INDEX.json'):
  if 'trajectory_path' not in record:continue
  raw=BASE/record['trajectory_path'];folder=raw.parent
  with h5py.File(raw) as h:
   times=h['time_s'][:];x=h['positions_m'][:];m=series(times,x,stride=10)
   dataset=BASE/'dataset'/(record['run_id']+'.h5');dataset.parent.mkdir(exist_ok=True)
   # Derived files identify original arrays through a relative external link.
   import os
   with h5py.File(dataset,'w') as d:
    d['raw']=h5py.ExternalLink(os.path.relpath(raw,dataset.parent),'/')
    for key in ['particle_id']:
     d[key]=h5py.ExternalLink(os.path.relpath(raw,dataset.parent),'/'+key)
    d['raw_positions_m']=h5py.ExternalLink(os.path.relpath(raw,dataset.parent),'/positions_m')
    d.create_dataset('positions_m',data=x[m['sample_index']],compression='gzip')
    d['solver_time_s']=h5py.ExternalLink(os.path.relpath(raw,dataset.parent),'/time_s')
    for key in ['pressure_Pa','primary_force_N','secondary_force_N','total_force_N','contact_force_N']:
     d[key]=h5py.ExternalLink(os.path.relpath(raw,dataset.parent),'/force_samples/'+key)
    d['force_time_s']=h5py.ExternalLink(os.path.relpath(raw,dataset.parent),'/force_samples/time_s')
    gate_rows=read(folder/'SOLVER_GATES.json');names=sorted(gate_rows[0]) if gate_rows else []
    d['solver_gates']=np.asarray([[r[k] for k in names] for r in gate_rows],bool)
    d['solver_gates'].attrs['columns']=json.dumps(names)
    d.attrs['random_seed']=str(read(BASE/'INITIAL_STATE_AUDIT.json')[record['state_id']]['seed'])
    for key,val in m.items():d.create_dataset(key,data=val,compression='gzip' if np.ndim(val)>0 else None)
    for alias,key in {'principal_axes':'principal_axes_3d','rigid_rotation':'rigid_rotation_3d','kabsch_rmsd_m':'kabsch_rmsd_3d_m','affine_deformation':'affine_deformation_3d','nonaffine_residual':'nonaffine_residual_3d_m','pair_distance_statistics':'pair_distance_statistics_m'}.items():d[alias]=h5py.SoftLink('/'+key)
    ii,jj=np.triu_indices(x.shape[1],1);d['proximity_pair_ids']=np.c_[ii,jj]
    sampled=x[m['sample_index']];d.create_dataset('proximity_adjacency',data=np.linalg.norm(sampled[:,ii]-sampled[:,jj],axis=-1)<=2.1*R,compression='gzip')
    d['contact_pair_reaction_N']=h5py.ExternalLink(os.path.relpath(raw,dataset.parent),'/force_samples/pair_contact_reaction_N')
    d['gate_time_s']=h5py.ExternalLink(os.path.relpath(raw,dataset.parent),'/force_samples/time_s')
    d.attrs.update(run_id=record['run_id'],initial_state_id=record['state_id'],particle_count=x.shape[1],particle_radius_m=R,
       model_id=LABELS[1],backend_id='F10AA_Warp_SciPy_CPU1',source_commit=record['task'].get('source_commit',read(folder/'INPUT.json')['identity']['execution_commit']),
       trajectory_sha256=sha(raw),scopes='SIMULATION_3D_TRUTH;SIMULATION_XY_PROJECTION',training_ready=False,
       validity='Completed physical window' if record['status']=='PASS' else 'Accepted prefix only; see raw failure',
       metrics_sampling='Every 10 accepted steps plus endpoint; raw positions and gates at every step',
       units='positions, COM, Rg, RMSD: m; covariance: m^2; stretches, aspect: dimensionless; times: s',
       action_schedule=json.dumps(read(folder/'INPUT.json')['schedule']),failure_reason=json.dumps(record.get('reasons',[])))
  audit_path=out/(record['run_id']+'_AUDIT.json')
  if audit_path.exists() and read(audit_path).get('trajectory_sha256')==sha(raw):audit=read(audit_path)
  elif record['status']=='PASS':audit=independent(folder);write(audit_path,audit)
  else:
   audit=dict(PASS=False,checked_steps=0,trajectory_sha256=sha(raw),reason='Incomplete trajectory; original failure retained. Independent full-window audit not admissible.')
   write(audit_path,audit)
  row={k:record[k] for k in ['run_id','state_id','stage','status','action','variant','wall_s','completed_s','trajectory_path']}
  row.update(metrics=describe(m),independent_audit=audit,trajectory_sha256=sha(raw),dataset=str(dataset.relative_to(BASE)))
  summaries.append(row);write(out/'SUMMARY.json',summaries)
  print('AUDIT',record['run_id'],audit['PASS'],round(row['metrics']['similarity_residual_nm'],5),flush=True)
 write(out/'COMPLETE.json',dict(wall_s=time.perf_counter()-start,identity=identity(),runs=len(summaries),checked_steps=sum(r['independent_audit']['checked_steps'] for r in summaries)))

def numerical(a,b):
 pa=BASE/a['trajectory_path'];pb=BASE/b['trajectory_path'];basic=compare(pa.parent,pb.parent)
 with h5py.File(pa) as ha,h5py.File(pb) as hb:
  ta=ha['time_s'][:];tb=hb['time_s'][:];ia=np.arange(0,len(ta),10);ia=np.unique(np.r_[ia,len(ta)-1]);ib=np.searchsorted(tb,ta[ia]-1e-13)
  assert np.allclose(tb[ib],ta[ia],atol=1e-13,rtol=0)
  ma=series(ta[ia],ha['positions_m'][:][ia]);mb=series(tb[ib],hb['positions_m'][:][ib]);dif={}
  for key in ['similarity_rmsd_3d_m','kabsch_rmsd_3d_m','nonaffine_residual_3d_m','deviatoric_stretch_3d','shear_norm_3d','directional_xy_stretch_3d','isotropic_stretch_3d','aspect_3d','aspect_xy']:
   d=abs(ma[key]-mb[key]);dif[key]=float(np.nanmax(d)) if np.isfinite(d).any() else None
 return dict(reference=a['run_id'],comparison=b['run_id'],original=basic,shape_uncertainty=dif,extra=extra_compare(pa.parent,pb.parent))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage');a=ap.parse_args();analyze_stage(a.stage)
if __name__=='__main__':main()
