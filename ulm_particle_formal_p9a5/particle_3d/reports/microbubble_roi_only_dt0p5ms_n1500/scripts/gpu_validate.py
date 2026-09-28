"""Real CUDA/float64 computation; never substitutes point tracers for microbubbles.

Prepare: independently calculate all tetra gradients and wall scale on CUDA.
Tracks: compute trajectory displacements/speeds/gaps/time differences in batches.
CPU geometry and finite-size event certificates remain the authoritative checks.
"""
from pathlib import Path
import argparse, json, time
import numpy as np
import torch

HERE=Path(__file__).resolve().parents[1]


def tensor(a):
    return torch.as_tensor(a,device='cuda',dtype=torch.float64)


def main(mode):
    assert torch.cuda.is_available()
    torch.set_num_threads(1)
    start=time.time();torch.cuda.reset_peak_memory_stats()
    event0=torch.cuda.Event(enable_timing=True);event1=torch.cuda.Event(enable_timing=True)
    event0.record()
    info=dict(device=torch.cuda.get_device_name(),torch=torch.__version__,cuda=torch.version.cuda,
              dtype='float64',role='INDEPENDENT_BULK_DIAGNOSTICS_NOT_DYNAMICS_KERNEL')
    if mode=='mesh':
        a=np.load(HERE/'data/gpu_mesh_input.npz')
        nodes=tensor(a['velocity'])[torch.as_tensor(a['tetra'],device='cuda')]
        inverse=tensor(a['inverse'])
        du=nodes[:,1:]-nodes[:,:1]
        grad=torch.einsum('nai,naj->nij',du,inverse)
        cpu=tensor(a['gradients'])
        error=(grad-cpu).abs()
        # Absolute error scaled by each contraction's magnitude; near-zero
        # tensor entries should not be judged with relative error alone.
        bound=64*np.finfo(float).eps*torch.einsum('nai,naj->nij',du.abs(),inverse.abs()).clamp_min(1.)
        passed=bool(torch.all(error<=bound).item())
        strain=.5*(grad+grad.transpose(1,2))
        shear=torch.sqrt(2*(strain*strain).sum(dim=(1,2)))
        curl=torch.stack((grad[:,2,1]-grad[:,1,2],grad[:,0,2]-grad[:,2,0],grad[:,1,0]-grad[:,0,1]),dim=1)
        wall_scale=tensor(a['wall_triangles']).abs().max().item()
        assert wall_scale==float(np.max(np.abs(a['wall_triangles'])))
        np.savez_compressed(HERE/'data/gpu_derived_field.npz',
            shear_rate_s_inv=shear.cpu().numpy(),vorticity_s_inv=curl.cpu().numpy())
        info.update(tetra_count=len(grad),gradient_max_absolute_error_s_inv=error.max().item(),
                    gradient_max_budget_ratio=(error/bound).max().item(),
                    wall_abs_coordinate_max_m=wall_scale,CPU_CUDA_gradient_parity=passed,
                    all_finite=bool(torch.isfinite(grad).all().item()))
        assert passed and info['all_finite']
    else:
        rows=[]
        for folder in sorted((HERE/'tracks').glob('mb_*')):
            if not (folder/'COMPLETE.json').exists():continue
            s=np.load(folder/'trajectory.npz')['samples'];v=tensor(s)
            delta=v[1:,1:4]-v[:-1,1:4];dt=v[1:,0]-v[:-1,0]
            length=torch.linalg.vector_norm(delta,dim=1).sum().item()
            maximum=torch.linalg.vector_norm(v[:,4:7],dim=1).max().item()
            cpu=json.loads((folder/'metrics.json').read_text())
            # Saved velocity is held on its accepted interval, not a separately
            # interpolated flow speed. This check reports data integrity only.
            rows.append(dict(particle_id=cpu['particle_id'],samples=len(s),
                path_length_m=length,maximum_speed_m_s=maximum,
                path_error_m=abs(length-cpu['path_length_m']),
                speed_error_m_s=abs(maximum-cpu['maximum_speed_m_s']),
                dt_min_s=dt.min().item(),dt_max_s=dt.max().item(),
                nonfinite=int((~torch.isfinite(v)).sum().item()),
                minimum_saved_wall_gap_m=v[:,14].min().item()))
        info.update(count=len(rows),sample_count=sum(r['samples'] for r in rows),rows=rows)
        info['PASS']=len(rows)==1500 and all(r['nonfinite']==0 and r['dt_min_s']>0 and r['dt_max_s']<=.0005+1e-12 and r['path_error_m']<1e-15 and r['speed_error_m_s']<1e-15 for r in rows)
    event1.record();torch.cuda.synchronize()
    info.update(cuda_event_elapsed_ms=event0.elapsed_time(event1),
        peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated(),wall_seconds=time.time()-start)
    (HERE/f'data/gpu_{mode}_validation.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps({k:v for k,v in info.items() if k!='rows'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['mesh','tracks']);main(p.parse_args().mode)
