"""Compact 50-us review copy; full 10-us export and every-step raw stay unchanged."""
import sys,json,hashlib
from pathlib import Path
import h5py,numpy as np

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def indices(t):
 selected=[0]
 for i in range(1,len(t)):
  if t[i]-t[selected[-1]]>=50e-6-1e-13:selected.append(i)
 return np.unique(selected+[len(t)-1])
def main():
 b=Path(sys.argv[1]);out=b/'review_dataset';out.mkdir(exist_ok=True);meta=json.loads((b/'portable_dataset/INDEX.json').read_text());records=[]
 aliases={'principal_axes':'principal_axes_3d','rigid_rotation':'rigid_rotation_3d','kabsch_rmsd_m':'kabsch_rmsd_3d_m','affine_deformation':'affine_deformation_3d','nonaffine_residual':'nonaffine_residual_3d_m','pair_distance_statistics':'pair_distance_statistics_m'}
 force={'pressure_Pa','primary_force_N','secondary_force_N','total_force_N','contact_force_N','velocity_m_s','contact_pair_reaction_N','force_sample_index','force_time_s'}
 static={'particle_id','proximity_pair_ids','failure_proposal_m','solver_gates','solver_time_s','gate_time_s'}
 for r in meta['trajectories']:
  src=b/r['file'];dest=out/src.name
  with h5py.File(src) as h,h5py.File(dest,'w') as d:
   t=h['time_s'][:];ix=indices(t);ft=h['force_time_s'][:];fi=indices(ft) if len(ft) else np.array([],int)
   for key in h:
    if key in aliases:continue
    data=h[key][:]
    if key in force:data=data[fi]
    elif key not in static and len(data)==len(t):data=data[ix]
    obj=d.create_dataset(key,data=data,compression='gzip');obj.attrs.update(dict(h[key].attrs))
   for alias,target in aliases.items():d[alias]=h5py.SoftLink('/'+target)
   d.attrs.update(dict(h.attrs));d.attrs.update(position_sampling='Actual saved states >=50us apart plus endpoint; no interpolation; full solver gates retained',
    source_portable_sha256=r['sha256'],metrics_sampling='50us selected actual states plus endpoint; do not use this review copy for event-time convergence checks')
   assert len(d['time_s'])==len(d['positions_m']) and np.array_equal(d['positions_m'][:],h['positions_m'][ix])
   assert np.array_equal(d['solver_gates'][:],h['solver_gates'][:])
  records.append(dict(run_id=r['run_id'],file=str(dest.relative_to(b)),sha256=sha(dest),bytes=dest.stat().st_size,source_portable=r['file'],source_portable_sha256=r['sha256'],raw_path=r['raw_path'],raw_sha256=r['raw_sha256'],status=r['status']))
 (out/'INDEX.json').write_text(json.dumps(dict(trajectories=records,sampling_s=50e-6,all_solver_gates_retained=True,training_ready=False,script_sha256=sha(Path(__file__))),indent=2)+'\n')
 print(len(records),'review files',sum(r['bytes'] for r in records),'bytes')
if __name__=='__main__':main()
