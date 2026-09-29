"""Real BraVa P1 interpolation/identity check, including shared mesh boundaries."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import torch

ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);args=ap.parse_args()
sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.gpu_field import GPUField
from particle_3d.population_v1.provenance import write_json
root=Path(__file__).resolve().parents[1]
spec=json.loads((root/'data/environment.json').read_text());started=time.perf_counter()
env=load_environment(spec);field=env.field
print('ENVIRONMENT_LOADED',time.perf_counter()-started,flush=True)
gpu=GPUField(field,cache_path=root/'data/tetra_acceleration.npz')
print('GPU_FIELD_READY',time.perf_counter()-started,flush=True)
rng=np.random.default_rng(2026092905);cells=rng.integers(0,len(field.tetra),20000)
w=rng.dirichlet(np.ones(4),len(cells));xyz=np.einsum('ni,nij->nj',w,field.points_m[field.tetra[cells]])
nodes=rng.integers(0,len(field.points_m),1000)
tri=field.points_m[field.tetra[cells[:1000],:3]].mean(1)
edges=field.points_m[field.tetra[cells[1000:1500],:2]].mean(1)
outside=np.array([field.points_m.min(0)-.01,field.points_m.max(0)+.01])
positions=np.vstack([xyz,field.points_m[nodes],tri,edges,outside]);types=['interior']*len(xyz)+['vertex']*len(nodes)+['face']*len(tri)+['edge']*len(edges)+['outside']*2
cpu_start=time.perf_counter();reference=field.sample_many(positions);cpu_s=time.perf_counter()-cpu_start
x=torch.as_tensor(positions,device='cuda',dtype=torch.float64)
torch.cuda.synchronize();start=time.perf_counter();ids,value,stats=gpu.sample(x);torch.cuda.synchronize();gpu_s=time.perf_counter()-start
data=value.cpu().numpy();ids_host=ids.cpu().numpy();mismatch=np.flatnonzero(ids_host!=reference.tetra_id)
inside=reference.inside_lumen;errors={}
for name,cols,ref in [('velocity_m_s',slice(0,3),reference.velocity_m_s),('pressure_Pa',3,reference.pressure_pa),
    ('gradient_s_inv',slice(4,13),reference.velocity_gradient_s_inv.reshape(-1,9)),
    ('vorticity_s_inv',slice(13,16),reference.vorticity_s_inv),('strain_s_inv',slice(16,25),reference.strain_rate_s_inv.reshape(-1,9))]:
    errors[name]=float(np.max(np.abs(data[inside,cols]-ref[inside]),initial=0.))
# Correct cached IDs and deliberately stale hints exercise walking/recovery.
hint=torch.as_tensor(np.r_[cells,np.full(len(positions)-len(cells),-1)],device='cuda')
walk_ids,walk_value,walk_stats=gpu.sample(x,hint)
wrong_hint=torch.zeros(len(positions),device='cuda',dtype=torch.int64)
recovery_ids,_,_=gpu.sample(x,wrong_hint)
ok=not len(mismatch) and np.array_equal(walk_ids.cpu().numpy(),reference.tetra_id) and np.array_equal(recovery_ids.cpu().numpy(),reference.tetra_id)
ok=ok and errors['velocity_m_s']<1e-12 and errors['pressure_Pa']<1e-9 and errors['gradient_s_inv']==0 and errors['vorticity_s_inv']==0
ok=ok and np.isnan(data[~inside]).all()
record=dict(PASS=bool(ok),phase=5,flow_sha256=env.flow_sha256,point_count=len(positions),
    fixture_types={k:types.count(k) for k in sorted(set(types))},maximum_absolute_errors=errors,
    tetra_id_mismatch_count=len(mismatch),mismatches=[dict(index=int(i),type=types[i],position_m=positions[i].tolist(),
        cpu_id=int(reference.tetra_id[i]),gpu_id=int(ids_host[i])) for i in mismatch[:20]],
    cached_and_stale_hint_parity=bool(np.array_equal(walk_ids.cpu().numpy(),reference.tetra_id) and np.array_equal(recovery_ids.cpu().numpy(),reference.tetra_id)),
    cpu_sample_s=cpu_s,gpu_sample_s=gpu_s,gpu=torch.cuda.get_device_name(),torch_version=torch.__version__,
    maximum_walk_steps=int(walk_stats[:,0].max()),cached_global_recovery_count=int(walk_stats[:,1].sum()),
    geometry='CANONICAL_CONFORMING_NONOVERLAPPING_TETRA_MESH',outside_policy='ID_MINUS_ONE_ALL_PHYSICAL_FIELDS_NAN',
    compiler='FLOAT64_NO_FAST_MATH_FMA_DISABLED',elapsed_s=time.perf_counter()-started)
write_json(root/'validation/phase5_real_fem_gpu.json',record)
print(json.dumps(record,indent=2),flush=True)
if not ok:raise SystemExit(1)
