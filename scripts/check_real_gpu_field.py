"""Small read-only BraVa CPU/GPU field check; original tolerances, no integration."""
from pathlib import Path
import time, json
import numpy as np
import torch
from human_vessel_flow.paths import ROOT, fresh_directory, write_json
from human_vessel_flow.particles import load_environment
from particle_3d.population_v1.gpu_field import GPUField

def main():
    output=fresh_directory(ROOT/'validation/real_gpu_field')
    start=time.perf_counter();env=load_environment();field=env.field
    gpu=GPUField(field,cache_path=output/'tetra_acceleration.npz')
    rng=np.random.default_rng(2026092905);cells=rng.integers(0,len(field.tetra),2048)
    w=rng.dirichlet(np.ones(4),len(cells));xyz=np.einsum('ni,nij->nj',w,field.points_m[field.tetra[cells]])
    nodes=rng.integers(0,len(field.points_m),32)
    tri=field.points_m[field.tetra[cells[:32],:3]].mean(1)
    edges=field.points_m[field.tetra[cells[32:64],:2]].mean(1)
    outside=np.array([field.points_m.min(0)-.01,field.points_m.max(0)+.01])
    positions=np.vstack([xyz,field.points_m[nodes],tri,edges,outside]);reference=field.sample_many(positions)
    x=torch.as_tensor(positions,device='cuda',dtype=torch.float64)
    ids,value,stats=gpu.sample(x);data=value.cpu().numpy();inside=reference.inside_lumen
    errors={}
    for name,cols,ref in [('velocity_m_s',slice(0,3),reference.velocity_m_s),('pressure_Pa',3,reference.pressure_pa),
        ('gradient_s_inv',slice(4,13),reference.velocity_gradient_s_inv.reshape(-1,9)),
        ('vorticity_s_inv',slice(13,16),reference.vorticity_s_inv),('strain_s_inv',slice(16,25),reference.strain_rate_s_inv.reshape(-1,9))]:
        errors[name]=float(np.max(np.abs(data[inside,cols]-ref[inside]),initial=0.))
    hint=torch.as_tensor(np.r_[cells,np.full(len(positions)-len(cells),-1)],device='cuda')
    walk,_,_=gpu.sample(x,hint);stale,_,_=gpu.sample(x,torch.zeros(len(x),device='cuda',dtype=torch.int64))
    ok=all(np.array_equal(a.cpu().numpy(),reference.tetra_id) for a in [ids,walk,stale])
    ok=ok and errors['velocity_m_s']<1e-12 and errors['pressure_Pa']<1e-9 and errors['gradient_s_inv']==0 and errors['vorticity_s_inv']==0
    ok=ok and np.isnan(data[~inside]).all()
    result=dict(status='PASS' if ok else 'FAIL',flow_sha256=env.flow_sha256,queries=len(x),maximum_absolute_errors=errors,
        cached_and_stale_hint_parity=all(np.array_equal(a.cpu().numpy(),reference.tetra_id) for a in [walk,stale]),
        torch_version=torch.__version__,cuda=torch.version.cuda,gpu=torch.cuda.get_device_name(),
        elapsed_s=time.perf_counter()-start,new_trajectory_integrations=0,physical_accuracy_certification=False,
        original_tolerances='u<1e-12 m/s, p<1e-9 Pa, G/curl bit-identical, IDs identical, outside NaN')
    write_json(output/'RESULT.json',result);print(json.dumps(result,indent=2))
    if not ok:raise SystemExit(1)

if __name__=='__main__':main()
