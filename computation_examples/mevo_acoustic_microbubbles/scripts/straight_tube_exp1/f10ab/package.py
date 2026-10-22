"""Self-contained sampled HDF5; full every-step originals remain immutable."""
import json,h5py,numpy as np,time,importlib.metadata,platform
from .common import *

def main():
 out=BASE/'portable_dataset';out.mkdir(exist_ok=True);records=[];raws=[]
 for summary in sorted((BASE/'analysis').glob('*/SUMMARY.json')):
  for row in read(summary):
   raw=BASE/row['trajectory_path'];derived=BASE/row['dataset'];dest=out/(row['run_id']+'.h5')
   with h5py.File(raw) as h,h5py.File(derived) as d,h5py.File(dest,'w') as p:
    p.attrs.update(dict(d.attrs));ix=d['sample_index'][:];times=h['time_s'][:];px=h['positions_m'][ix]
    skip={'raw','positions_m','raw_positions_m','solver_time_s','particle_id','pressure_Pa','primary_force_N','secondary_force_N','total_force_N','contact_force_N','force_time_s','contact_pair_reaction_N','gate_time_s'}
    for key in d:
     if key in skip:continue
     if isinstance(d[key],h5py.Dataset):
      obj=p.create_dataset(key,data=d[key][:],compression='gzip');obj.attrs.update(dict(d[key].attrs))
    p.create_dataset('positions_m',data=px,compression='gzip');p['particle_id']=h['particle_id'][:]
    p['solver_time_s']=times
    group=h['force_samples'];ft=group['time_s'][:];fi=ix[ix<len(ft)]
    p['force_time_s']=ft[fi];p['gate_time_s']=ft
    for key in ['pressure_Pa','primary_force_N','secondary_force_N','total_force_N','contact_force_N','velocity_m_s','pair_contact_reaction_N']:
     destkey='contact_pair_reaction_N' if key=='pair_contact_reaction_N' else key
     p.create_dataset(destkey,data=group[key][fi],compression='gzip')
    p['force_sample_index']=fi
    if 'failure_proposal_m' in h:p['failure_proposal_m']=h['failure_proposal_m'][:]
    for alias,key in {'principal_axes':'principal_axes_3d','rigid_rotation':'rigid_rotation_3d','kabsch_rmsd_m':'kabsch_rmsd_3d_m','affine_deformation':'affine_deformation_3d','nonaffine_residual':'nonaffine_residual_3d_m','pair_distance_statistics':'pair_distance_statistics_m'}.items():
     if alias not in p:p[alias]=h5py.SoftLink('/'+key)
    p.attrs.update(portable=True,position_sampling='Every 10 accepted solver steps plus endpoint; no interpolation',
      force_sampling='Actual force evaluations at force_time_s; endpoint has no next-step force unless a rejected proposal was attempted',
      full_raw_relative_path=row['trajectory_path'],full_raw_server_root='/workspace/mevo_acoustic_microbubbles/f10ab_shape_20261018/results/straight_tube_exp1/f10ab',
      f10ab_preregistration_sha256=sha(BASE/'PREREGISTRATION.json'),f10ab_initial_states_sha256=sha(BASE/'INITIAL_STATES.h5'),
      f10ab_stage_plan_sha256=sha(BASE/'plans'/(row['stage']+'.json')),
      parent_input_preregistration_scope='Original inherited F10AA INPUT.json preregistration_sha256 identifies the parent F10AA protocol; F10AB identities are explicit separate attributes here')
    assert len(p['time_s'])==len(p['positions_m']) and np.array_equal(p['time_s'][:],times[ix])
    assert np.array_equal(p['positions_m'][:],px) and p['particle_id'].shape==(px.shape[1],)
    for key in p:assert not isinstance(p.get(key,getlink=True),h5py.ExternalLink)
   rec=dict(run_id=row['run_id'],file=str(dest.relative_to(BASE)),sha256=sha(dest),bytes=dest.stat().st_size,
    raw_path=row['trajectory_path'],raw_sha256=sha(raw),raw_bytes=raw.stat().st_size,status=row['status'],sampled=True)
   records.append(rec);raws.append({k:v for k,v in rec.items() if k.startswith('raw') or k in ['run_id','status']})
 write(out/'INDEX.json',dict(trajectories=records,identity=identity(),training_ready=False,complete_raw_on_server=True))
 write(BASE/'RAW_DATA_MANIFEST.json',dict(records=raws,server_root='/workspace/mevo_acoustic_microbubbles/f10ab_shape_20261018/results/straight_tube_exp1/f10ab',identity=identity()))
 write(BASE/'DATASET_SCHEMA.json',dict(position_time='time_s matches positions_m, derived metrics and proximity adjacency',
  force_time='force_time_s matches pressure and F1/F2/contact forces; evaluation occurs before the corresponding attempted step',
  gate_time='gate_time_s matches full solver_gates, including rejected last proposal if any',
  full_time='solver_time_s contains all accepted position times; full original positions in raw_path',
  aliases='principal_axes, rigid_rotation, kabsch_rmsd_m, affine_deformation, nonaffine_residual refer to 3D definitions; _xy fields retain projection',
  distance_quantiles=[0,.1,.5,.9,1],proximity_rule='center distance <=2.1R; not adhesion or fusion',
  scope=['SIMULATION_3D_TRUTH','SIMULATION_XY_PROJECTION'],training_ready=False))
 packages=['numpy','scipy','h5py','matplotlib','pillow','av','pymupdf']
 write(BASE/'evidence/POSTPROCESS_ENVIRONMENT.json',dict(python=platform.python_version(),packages={k:importlib.metadata.version(k) for k in packages},identity=identity()))
 print('Portable trajectories',len(records),'bytes',sum(r['bytes'] for r in records))
if __name__=='__main__':main()
