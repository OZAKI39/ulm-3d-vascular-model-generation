"""GPU hash/cell list with explicit realization isolation and exact spheres."""
import torch
from .gpu_extension import extension


def neighbors(x,radius,realization,*,delta=None):
    """Return conservative canonical pairs, exact start gaps/normals, diagnostics.

    Swept queries enclose complete paths in balls, as CPU candidates(). Only
    scalar allocation metadata crosses to the host; no population state does.
    """
    n=len(x)
    if x.dtype!=torch.float64 or radius.dtype!=torch.float64:raise ValueError('Float64 required')
    if not n:
        return torch.empty((0,2),dtype=torch.int64,device=x.device),torch.empty((0,4),dtype=x.dtype,device=x.device),dict(candidate_count=0)
    if delta is None:delta=torch.zeros_like(x)
    center=x+.5*delta
    pad=512*torch.finfo(torch.float64).eps*torch.maximum(torch.maximum(x.abs().max(),(x+delta).abs().max()),radius.max())
    bound=radius+.5*torch.linalg.vector_norm(delta,dim=1)+torch.maximum(.05*radius,torch.full_like(radius,2e-9))+pad
    cell_size=float(2*bound.max())
    cells=torch.floor(center/cell_size).to(torch.int64)
    keys=((cells[:,0]*73856093)^(cells[:,1]*19349663)^(cells[:,2]*83492791)^(realization*2654435761))&0x7fffffffffffffff
    order=torch.argsort(keys,stable=True);sorted_keys=keys[order].contiguous()
    args=[center.contiguous(),radius.contiguous(),bound.contiguous(),realization.contiguous(),cells.contiguous(),sorted_keys,order.contiguous()]
    ext=extension();counts=ext.neighbor_count(args,cell_size)
    cumulative=counts.cumsum(0);total=int(cumulative[-1]);offsets=(cumulative-counts).contiguous()
    pairs,geometry=ext.neighbor_fill(args,cell_size,offsets,total)
    if total:
        perm=torch.argsort(pairs[:,0]*n+pairs[:,1],stable=True);pairs=pairs[perm].contiguous()
        # Swept balls choose candidates; geometry always describes actual start spheres.
        d=x[pairs[:,0]]-x[pairs[:,1]];length=torch.linalg.vector_norm(d,dim=1)
        normal=d/length[:,None];coincident=length==0
        normal[coincident]=torch.tensor([1.,0.,0.],device=x.device,dtype=x.dtype)
        geometry=torch.cat([(length-radius[pairs[:,0]]-radius[pairs[:,1]])[:,None],normal],dim=1)
    return pairs,geometry,dict(candidate_count=total,cell_size_m=cell_size,all_pairs_used=False,
        realization_isolation=True,swept=bool(torch.any(delta)),host_metadata_only=True)
