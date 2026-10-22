"""Static canonical P1 field on GPU, adjacency walk plus exact BVH recovery."""
from pathlib import Path
import numpy as np
from .bvh import build_bvh,tetra_adjacency
from .gpu_extension import extension
from .provenance import content_sha


class GPUField:
    def __init__(self,field,*,device='cuda',cache_path=None):
        import torch
        self.device=torch.device(device)
        if self.device.index is None:self.device=torch.device('cuda',torch.cuda.current_device())
        self.flow_field=field
        identity=content_sha(dict(tetra_shape=list(field.tetra.shape),point_shape=list(field.points_m.shape),
            points_sha=__import__('hashlib').sha256(field.points_m.tobytes()).hexdigest(),
            tetra_sha=__import__('hashlib').sha256(field.tetra.tobytes()).hexdigest(),schema=1))
        path=None if cache_path is None else Path(cache_path)
        if path is not None and path.exists():
            with np.load(path) as cache:
                if str(cache['identity'])!=identity:raise ValueError('GPU_FIELD_ACCELERATOR_CACHE_MISMATCH')
                adjacency=cache['adjacency'];bounds=cache['bounds'];links=cache['links'];order=cache['order']
        else:
            adjacency=tetra_adjacency(field.tetra);xyz=field.points_m[field.tetra]
            bounds,links,order=build_bvh(xyz.min(1),xyz.max(1))
            if path is not None:
                path.parent.mkdir(parents=True,exist_ok=True)
                np.savez_compressed(path,identity=identity,adjacency=adjacency,bounds=bounds,links=links,order=order)
        shape_norm=np.linalg.norm(field.geometry.shape_gradients_m_inv,axis=2)
        arrays=[field.tetra,field.geometry.origins,field.geometry.inverse,field.geometry.weight_tolerance,
            adjacency,shape_norm,bounds,links,order,field.velocity_nodes_m_s,field.pressure_nodes_pa,
            field.gradients_s_inv,field.curl_s_inv,field.strain_s_inv]
        self.tensors=[torch.as_tensor(np.array(a,copy=True),device=self.device).contiguous() for a in arrays]
        self.padding=field.geometry.candidate_padding_m
        self.extension=extension()

    def sample(self,positions,tetra_ids=None):
        import torch
        if positions.dtype!=torch.float64 or positions.device!=self.device:
            raise ValueError('GPU FEM requires device-resident FP64 positions')
        if positions.ndim!=2 or positions.shape[1]!=3:raise ValueError('N by 3 positions required')
        if not bool(torch.isfinite(positions).all()):raise ValueError('Finite positions required')
        if tetra_ids is None:tetra_ids=torch.full((len(positions),),-1,dtype=torch.int64,device=self.device)
        return self.extension.fem_sample([positions.contiguous(),tetra_ids.contiguous(),*self.tensors],self.padding)
