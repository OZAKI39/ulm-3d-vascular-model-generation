import os
from pathlib import Path
import numpy as np
import pytest
torch=pytest.importorskip('torch')
pytestmark=pytest.mark.skipif(not torch.cuda.is_available(),reason='Requires CUDA')
from scipy import sparse
from particle_3d.particle_shapes import Sphere
from particle_3d.hydrodynamic_resistance import PhysicalNearField,sphere_self_diagonal
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block
from particle_3d.resistance_assembly import ResistanceSystem
from particle_3d.resistance_solver import solve_resistance
from particle_3d.kinematic_contact import ContactConstraint
from particle_3d.population_v1.gpu_resistance import assemble_isolated,solve_isolated
from particle_3d.population_v1.provenance import write_json


def test_cuda_affine_resistance_and_all_contact_inequalities():
    rng=np.random.default_rng(2026092971);n=96;mu=.00345312
    radius=rng.uniform(.375e-6,2e-6,n);normal=rng.normal(size=(n,3));normal/=np.linalg.norm(normal,axis=1)[:,None]
    gap=radius*10**rng.uniform(-3,1,n);gap=np.maximum(gap,2e-9)
    gradient=rng.normal(0,100,(n,3,3));free=rng.normal(size=(n,6))*np.array([.02,.02,.02,100,100,100])
    wall=np.zeros((n,14));wall[:,0]=gap;wall[:,7:10]=normal;wall[::9,13]=1
    cn=np.zeros((n,32,3));cc=np.arange(n)%5
    for i,k in enumerate(cc):
        cn[i,:k]=rng.normal(size=(k,3));cn[i,:k]/=np.linalg.norm(cn[i,:k],axis=1)[:,None]
        if k==4:cn[i,3]=cn[i,0] # exact duplicate: keep all inequalities
        if i%7==0 and k>=2:cn[i,1]=-cn[i,0] # opposing normals must not be dropped
    tensor=lambda v:torch.tensor(v,device='cuda')
    mat,rhs,diag,ad=assemble_isolated(tensor(radius),mu,tensor(wall),tensor(gradient),tensor(free))
    velocity,hydro,audit=solve_isolated(mat,rhs,diag,tensor(free),tensor(cn),tensor(cc))
    ma=mat.cpu().numpy();ra=rhs.cpu().numpy();va=velocity.cpu().numpy();fallback=audit['fallback'].cpu().numpy()
    errors=dict(relative_matrix=0.,relative_rhs=0.,velocity_m_s=0.,omega_s_inv=0.);rows=[]
    for i in range(n):
        shape=Sphere([0,0,0],radius[i]);di=sphere_self_diagonal(shape,mu)
        coefficient,nn,record=NearFieldRegularizationV1().evaluate(PhysicalNearField(1,None,gap[i],normal[i],0.),{1:shape},mu)
        j=np.r_[nn,np.zeros(3)];matrix=np.diag(di)+coefficient*np.outer(j,j);b=di*free[i]
        dr,db,_=planar_wall_affine_block(radius[i],mu,gap[i],normal[i],gradient[i],free[i])
        if not wall[i,13]:matrix+=dr;b+=db
        system=ResistanceSystem((1,),sparse.csr_matrix(matrix),b,di,free[i],[],[])
        constraints=[ContactConstraint(('TEST',k),1,None,cn[i,k],np.zeros(3)) for k in range(cc[i])]
        ref=solve_resistance(system,particles={1:shape},constraints=constraints)
        errors['relative_matrix']=max(errors['relative_matrix'],float(np.linalg.norm(ma[i]-matrix)/np.linalg.norm(matrix)))
        errors['relative_rhs']=max(errors['relative_rhs'],float(np.linalg.norm(ra[i]-b)/np.linalg.norm(b)))
        errors['velocity_m_s']=max(errors['velocity_m_s'],float(np.max(abs(va[i,:3]-ref.velocity[:3]))))
        errors['omega_s_inv']=max(errors['omega_s_inv'],float(np.max(abs(va[i,3:]-ref.velocity[3:]))))
        rows.append(dict(index=i,contact_count=int(cc[i]),gpu_fallback=bool(fallback[i]),
            max_velocity_error_m_s=float(np.max(abs(va[i,:3]-ref.velocity[:3])))))
    result=dict(PASS=not np.any(fallback) and errors['velocity_m_s']<1e-10 and errors['omega_s_inv']<1e-4,
        phase=7,scope='ISOLATED_GPU_RESISTANCE_AND_WALL_CONTACT_ALGEBRA',n=n,errors=errors,
        fallback_count=int(fallback.sum()),fixtures=rows,
        includes_exact_redundant_and_opposing_constraints=True)
    write_json(Path(os.environ.get('POPULATION_VALIDATION_DIR','/tmp/population_v1_validation'))/'phase7_gpu_resistance.json',result)
    assert result['PASS'],result
