"""P1 surface-work conjugate dead loads. Never reads an internal force/solution."""
import numpy as np

BARY=np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]])

def traction(normals,p,tau):
    N=np.asarray(normals);ex=np.array([1.,0.,0.])
    return -p*N+tau*(ex-N[:,0,None]*N)

def assemble(X,triangles,areas,normals,wall_faces,p,tau):
    ids=np.flatnonzero(~wall_faces);faces=triangles[ids]
    T=traction(normals[ids],p,tau);weights=np.repeat(areas[ids,None]/3,3,axis=1)
    points=np.einsum('qa,fac->fqc',BARY,X[faces])
    f=np.zeros_like(X)
    for a in range(3):
        np.add.at(f,faces[:,a],np.sum(weights*BARY[None,:,a],axis=1)[:,None]*T)
    return f,dict(face_ids=ids,face_nodes=faces,quadrature_points=points,
                  quadrature_weights=weights,traction_Pa=T,shape_values=BARY,
                  normals=normals[ids])

def external_work(u,data):
    uq=np.einsum('qa,fac->fqc',data['shape_values'],u[data['face_nodes']])
    return float(np.einsum('fq,fc,fqc->',data['quadrature_weights'],data['traction_Pa'],uq))

def ideal_integrals(R,p,tau):
    F=np.array([4*np.pi/3*tau*R**2,0.,-p*np.pi*R**2])
    M=np.array([0.,np.pi*tau*R**3,0.])
    H=-(2*np.pi/3)*p*R**3*np.eye(3)
    H[0,2]+=3*np.pi/4*tau*R**3;H[2,0]-=np.pi/4*tau*R**3
    return F,M,H

def independent_ideal_quadrature(R,p,tau):
    z,w=np.polynomial.legendre.leggauss(32);z=(z+1)/2;w=w/2
    phi=2*np.pi*(np.arange(128)+.5)/128
    zz,pp=np.meshgrid(z,phi,indexing='ij');rr=np.sqrt(1-zz*zz)
    N=np.stack([rr*np.cos(pp),rr*np.sin(pp),zz],axis=-1).reshape(-1,3)
    weights=np.broadcast_to(w[:,None]*R**2*2*np.pi/128,(32,128)).ravel()
    T=traction(N,p,tau);X=R*N
    return np.einsum('n,ni->i',weights,T),np.einsum('n,ni->i',weights,np.cross(X,T)),np.einsum('n,ni,nj->ij',weights,T,X)

def load_audit(X,triangles,areas,normals,wall_faces,fixed,R,p,tau):
    f,data=assemble(X,triangles,areas,normals,wall_faces,p,tau)
    ids=data['face_ids'];centers=X[triangles[ids]].mean(1);T=data['traction_Pa'];a=areas[ids]
    # Independent flat-face integration, no nodal forces in these three quantities.
    F=np.einsum('n,ni->i',a,T);M=np.einsum('n,ni->i',a,np.cross(centers,T));H=np.einsum('n,ni,nj->ij',a,T,centers)
    Fi,Mi,Hi=ideal_integrals(R,p,tau);Fq,Mq,Hq=independent_ideal_quadrature(R,p,tau)
    fs=3*np.pi*R**2;ms=fs*R
    v=np.column_stack([X[:,2],X[:,0]+.2*X[:,2],-.3*X[:,1]])
    u=np.column_stack([.003*X[:,0],.002*X[:,1],.001*X[:,2]])
    eps=.001;fd=(external_work(u+eps*v,data)-external_work(u-eps*v,data))/(2*eps)
    B=data['shape_values'];coords=np.einsum('qa,fac->fqc',B,X[data['face_nodes']])
    tangent=traction(normals[ids],0,1)
    algebra=max(np.linalg.norm(f.sum(0)-F)/fs,np.linalg.norm(np.cross(X,f).sum(0)-M)/ms,
                np.linalg.norm(np.einsum('ni,nj->ij',f,X)-H)/ms,abs(fd-np.sum(f*v))/ms)
    r=dict(p_Pa=p,tau_Pa=tau,mode='reference-flat-face normal; nominal dead traction',
           units=dict(traction='Pa',weights='m^2',nodal_force='N',work='J'),
           ideal_resultant_N=Fi,face_integral_resultant_N=F,nodal_resultant_N=f.sum(0),
           ideal_reference_moment_Nm=Mi,face_integral_reference_moment_Nm=M,nodal_reference_moment_Nm=np.cross(X,f).sum(0),
           ideal_affine_work_tensor_J=Hi,face_affine_work_tensor_J=H,
           resultant_error=np.linalg.norm(F-Fi)/max(np.linalg.norm(Fi),fs),moment_error=np.linalg.norm(M-Mi)/max(np.linalg.norm(Mi),ms),
           affine_work_max_error=np.max(np.abs(H-Hi))/ms,
           independent_ideal_quadrature_error=max(np.linalg.norm(Fq-Fi)/fs,np.linalg.norm(Mq-Mi)/ms,np.linalg.norm(Hq-Hi)/ms),
           discrete_assembly_error=algebra,normal_length_error=np.max(np.abs(np.linalg.norm(normals,axis=1)-1)),
           tangent_normal_dot_error=np.max(np.abs(np.einsum('ij,ij->i',tangent,normals[ids]))),
           shape_constant_error=np.max(np.abs(B.sum(1)-1)),affine_coordinate_error=np.max(np.abs(coords-data['quadrature_points'])),
           known_wall_external_resultant_N=f[fixed].sum(0),fixed_edge_loads_retained=True,
           force_reference_N=fs,moment_work_reference_Nm=ms,
           no_displacement_or_internal_force_input=True)
    r['PASS']=bool(algebra<=1e-10 and max(r['resultant_error'],r['moment_error'],r['affine_work_max_error'])<=.01)
    return r,f,data

def reactions(X,x,internal,external,fixed,R,force_floor):
    residual=internal+external;wall=np.zeros_like(residual);wall[fixed]=-residual[fixed]
    total=external+wall;fs=3*np.pi*R**2
    return wall,dict(free_residual_L2_N=np.linalg.norm(residual[~fixed]),
        free_residual_relative=np.linalg.norm(residual[~fixed])/max(np.linalg.norm(external[~fixed]),force_floor),
        wall_on_clot_resultant_N=wall.sum(0),clot_on_wall_resultant_N=-wall.sum(0),
        wall_known_external_resultant_N=external[fixed].sum(0),
        total_force_imbalance=np.linalg.norm(total.sum(0))/fs,
        total_current_moment_Nm=np.cross(x,total).sum(0),
        total_current_moment_imbalance=np.linalg.norm(np.cross(x,total).sum(0))/(fs*R),
        wall_current_moment_Nm=np.cross(x,wall).sum(0),
        internal_resultant_N=internal.sum(0),internal_current_moment_Nm=np.cross(x,internal).sum(0))
