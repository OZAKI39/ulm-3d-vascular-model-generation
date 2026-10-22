"""FP64 batched P6.5 normal + P9-A.1 planar resistance and wall contact.

The formulas are literal tensor counterparts of the frozen Python modules.
Contact enumerates all independent active bases of at most three normals for
an isolated sphere. This preserves all inequalities, including redundant or
opposing normals; it does not silently discard linearly dependent constraints.
"""
from itertools import combinations
import math
import torch


def assemble_isolated(radius,mu,wall,gradient,free):
    n=len(radius);a=radius;h=wall[:,0];normal=wall[:,7:10]
    normal=normal/torch.linalg.vector_norm(normal,dim=1)[:,None]
    zt=6*math.pi*mu*a;zr=8*math.pi*mu*a**3
    diag=torch.cat([zt[:,None].expand(-1,3),zr[:,None].expand(-1,3)],dim=1)
    matrix=torch.diag_embed(diag);rhs=diag*free
    chi=h/a;nf=((chi-.01)/.04).clamp(0.,1.)
    weight=(1-nf)**2*(1+2*nf);lower=torch.maximum(torch.full_like(a,2e-9),.001*a)
    coefficient=weight*6*math.pi*mu*a**2/torch.maximum(h,lower)
    matrix[:,:3,:3]+=coefficient[:,None,None]*(normal[:,:,None]*normal[:,None,:])
    xi=(h/a).clamp(.001,.1);log=torch.log(xi)
    ft=(8/15)*log-.9588;tt=.1*log-.1895;fr=(2/15)*log-.2526;tr=.4*log-.3187
    determinant=tt*fr-ft*tr;mt=tr/determinant;mc=.5*(.75*fr/determinant+tt/determinant);mr=.75*ft/determinant
    near=torch.zeros((n,4,4),device=a.device,dtype=a.dtype)
    near[:,0,0]=mt;near[:,1,1]=mt;near[:,2,2]=mr;near[:,3,3]=mr
    near[:,0,3]=-mc;near[:,1,2]=mc;near[:,2,1]=mc;near[:,3,0]=-mc
    s=((h/a-.1)/.9).clamp(0.,1.);w=1.-3*s*s+2*s*s*s
    # Topology-only open-rim exclusion, identical to the original P9-A.1 gate.
    w=torch.where(wall[:,13]!=0,torch.zeros_like(w),w)
    eye4=torch.eye(4,device=a.device,dtype=a.dtype).expand(n,-1,-1)
    fdiag=torch.tensor([1.,1.,.75,.75],device=a.device,dtype=a.dtype)
    free_mobility=torch.diag_embed(fdiag.expand(n,-1))/zt[:,None,None]
    near=near/zt[:,None,None];effective=(1-w[:,None,None])*free_mobility+w[:,None,None]*near
    reff,info=torch.linalg.solve_ex(effective,eye4)
    rfree=torch.diag_embed((zt[:,None]/fdiag[None,:]))
    axis=torch.nn.functional.one_hot(torch.argmin(normal.abs(),dim=1),3).to(a.dtype)
    t1=axis-normal*(axis*normal).sum(1)[:,None];t1/=torch.linalg.vector_norm(t1,dim=1)[:,None]
    t2=torch.linalg.cross(normal,t1);tangent=torch.stack([t1,t2],dim=2)
    lift=torch.zeros((n,4,6),device=a.device,dtype=a.dtype)
    lift[:,:2,:3]=tangent.transpose(1,2);lift[:,2:,3:]=a[:,None,None]*tangent.transpose(1,2)
    gn=(gradient@normal[:,:,None]).squeeze(2)
    shear=gn-normal*(normal*gn).sum(1)[:,None]
    height=a+h.clamp(min=0.)
    force=(1.7005*zt*height)[:,None]*shear
    torque=(.9440*4*math.pi*mu*a**3)[:,None]*torch.linalg.cross(shear,normal)
    local=lambda v:(tangent.transpose(1,2)@v[:,:,None]).squeeze(2)
    load=torch.cat([local(force),local(torque)/a[:,None]],dim=1)
    linear=torch.cat([local(height[:,None]*shear),a[:,None]*local(.5*torch.linalg.cross(normal,shear))],dim=1)
    reference=(1-w[:,None])*linear+w[:,None]*(near@load[:,:,None]).squeeze(2)
    dr4=reff-rfree;db4=(reff@reference[:,:,None]-rfree@linear[:,:,None]).squeeze(2)
    dr4=torch.where((w==0)[:,None,None],torch.zeros_like(dr4),dr4)
    db4=torch.where((w==0)[:,None],torch.zeros_like(db4),db4)
    dr=lift.transpose(1,2)@dr4@lift;db=(lift.transpose(1,2)@db4[:,:,None]).squeeze(2)
    return matrix+dr,rhs+db,diag,dict(affine_solve_info=info,wall_weight=w,normal_coefficient=coefficient)


def solve_isolated(matrix,rhs,diagonal,free,contact_normals,contact_count):
    n=len(rhs);d=1/torch.sqrt(diagonal);a=matrix*d[:,:,None]*d[:,None,:]
    chol,info=torch.linalg.cholesky_ex(a);scaled=torch.cholesky_solve((d*rhs)[:,:,None],chol).squeeze(2)
    hydro=d*scaled;velocity=hydro.clone()
    eigen=torch.linalg.eigvalsh(a);condition=eigen[:,-1]/eigen[:,0]
    eps=torch.finfo(rhs.dtype).eps
    failed=(info!=0)|(~torch.isfinite(hydro).all(1))|(eigen[:,0]<=0)|(eps*6*condition>=1/512)
    multipliers=torch.zeros((n,32),device=rhs.device,dtype=rhs.dtype)
    budgets=torch.zeros(n,device=rhs.device,dtype=rhs.dtype)
    rank=torch.zeros(n,device=rhs.device,dtype=torch.int64)
    # Large or unusually ill-conditioned contact sets use the preserved CPU
    # component authority; this is an explicit subset fallback, not a new model.
    fallback=failed|(contact_count>8)
    for k in range(1,9):
        ix=torch.nonzero((contact_count==k)&(~failed)).flatten()
        if not len(ix):continue
        j=torch.zeros((len(ix),k,6),device=rhs.device,dtype=rhs.dtype);j[:,:,:3]=contact_normals[ix,:k]
        jd=j*d[ix,None,:]
        white=torch.linalg.solve_triangular(chol[ix],jd.transpose(1,2),upper=False).transpose(1,2)
        gram=white@white.transpose(1,2);sv=torch.linalg.svdvals(white)
        significant=sv>max(k,6)*eps*sv[:,:1]
        # The preserved CPU dual resolves W W^T at sqrt(n*eps), not eps.
        # Nearly dependent rows need its deterministic positive-cone reduction;
        # never inflate the admissible normal velocity by their spurious rank.
        unresolved=(significant & (sv <= (max(k,6)*eps)**.5*sv[:,:1])).any(1)
        fallback[ix]|=unresolved
        rank[ix]=significant.sum(1)
        minimum=torch.where(significant,sv,torch.full_like(sv,float('inf'))).min(1).values
        cw=(sv[:,0]/minimum).clamp(min=1.)
        scale=torch.maximum(hydro[ix,:3].abs().max(1).values,free[ix,:3].abs().max(1).values).clamp(min=torch.finfo(rhs.dtype).tiny)
        budget=512*eps*max(k,6)*cw*scale;budgets[ix]=budget
        h=(j@hydro[ix,:,None]).squeeze(2)
        dual_budget=budget/gram.diagonal(dim1=1,dim2=2).max(1).values
        best=torch.where((h>=-budget[:,None]).all(1),torch.zeros_like(budget),torch.full_like(budget,float('inf')))
        chosen=torch.zeros((len(ix),k),device=rhs.device,dtype=rhs.dtype)
        for size in range(1,min(k,3)+1):
            for subset in combinations(range(k),size):
                sub=torch.tensor(subset,device=rhs.device,dtype=torch.int64)
                gs=gram[:,sub[:,None],sub[None,:]];hs=h[:,sub]
                lam,linfo=torch.linalg.solve_ex(gs,-hs[:,:,None]);lam=lam.squeeze(2)
                full=torch.zeros_like(chosen);full[:,sub]=lam
                speed=h+(gram@full[:,:,None]).squeeze(2)
                objective=.5*(lam*(gs@lam[:,:,None]).squeeze(2)).sum(1)
                comp=(full*speed).abs().max(1).values
                comp_budget=budget*torch.maximum(scale,full.abs().max(1).values)
                valid=(linfo==0)&torch.isfinite(lam).all(1)&(lam>=-dual_budget[:,None]).all(1)&(speed>=-budget[:,None]).all(1)&(comp<=comp_budget)
                valid&=objective<best
                best=torch.where(valid,objective,best);chosen=torch.where(valid[:,None],full,chosen)
        chosen=chosen.clamp(min=0.)
        correction=d[ix]*torch.cholesky_solve((d[ix]*(j.transpose(1,2)@chosen[:,:,None]).squeeze(2))[:,:,None],chol[ix]).squeeze(2)
        corrected=hydro[ix]+correction
        speed=(j@corrected[:,:,None]).squeeze(2)
        valid=torch.isfinite(best)&(speed>=-budget[:,None]).all(1)&torch.isfinite(corrected).all(1)
        velocity[ix]=corrected;multipliers[ix,:k]=chosen;fallback[ix]|=~valid
    residual=(a@(velocity/d)[:,:,None]).squeeze(2)-d*rhs
    # Contact stationarity includes J^T lambda, with every original normal.
    jt=(contact_normals*multipliers[:,:,None]).sum(1)
    residual[:,:3]-=d[:,:3]*jt
    backward=torch.linalg.vector_norm(residual,dim=1)/(torch.linalg.matrix_norm(a,ord=2)*torch.linalg.vector_norm(velocity/d,dim=1)+torch.linalg.vector_norm(d*rhs,dim=1)).clamp(min=torch.finfo(rhs.dtype).tiny)
    fallback|=(backward>512*eps*6)|(~torch.isfinite(velocity).all(1))
    return velocity,hydro,dict(fallback=fallback,condition=condition,backward_residual=backward,
        multipliers=multipliers,velocity_budget=budgets,contact_rank=rank,
        solver='CUDA_BATCHED_SPHERE_RESISTANCE_AND_ALL_INEQUALITY_CONTACT_BASES')
