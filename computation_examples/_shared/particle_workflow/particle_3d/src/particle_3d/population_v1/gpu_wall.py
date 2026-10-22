"""Real finite-triangle WALL queries and held-path certificates on CUDA."""
import numpy as np
import torch
from .bvh import build_bvh
from .gpu_extension import extension


class GPUWall:
    def __init__(self,wall,device='cuda'):
        self.wall=wall;self.roundoff=wall.roundoff_m;tri=wall.triangles
        boxes,links,order=build_bvh(tri.min(1),tri.max(1))
        rim=np.zeros((len(tri),6),dtype=np.int64);topology=getattr(wall,'open_boundary_topology',None)
        if topology is not None:
            for i,gids in enumerate(wall.global_node_ids):
                for k in range(3):
                    rim[i,k]=int(int(gids[k]) in topology['vertices'])
                    edge=tuple(sorted(int(gids[j]) for j in range(3) if j!=k))
                    rim[i,3+k]=int(edge in topology['edges'])
        self.arrays=[torch.tensor(np.array(a,copy=True),device=device) for a in [tri,wall.normal_in,boxes,links,order,rim]]
        self.ext=extension()

    def query(self,x,radius):
        return self.ext.wall_query([x.contiguous(),radius.contiguous(),*self.arrays],self.roundoff)

    def path(self,x,velocity,radius,dt):
        if not isinstance(dt,torch.Tensor):dt=torch.full((len(x),),float(dt),dtype=x.dtype,device=x.device)
        return self.ext.wall_path([x.contiguous(),velocity.contiguous(),radius.contiguous(),dt.contiguous(),*self.arrays],self.roundoff)
