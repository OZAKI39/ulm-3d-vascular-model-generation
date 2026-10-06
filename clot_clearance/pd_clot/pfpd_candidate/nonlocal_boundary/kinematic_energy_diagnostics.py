"""Read-only quadratic moment/projection diagnostics, never a correction operator."""
import numpy as np
from ..energy import constitutive
from ..kernel import radial
from ..svk_spectral_split import green_lagrange,svk_split_energy


def tensor_norm(a):return np.sqrt(np.sum(np.asarray(a)**2,axis=tuple(range(1,np.ndim(a)))))


def moment_chunk(args):
    """Independent scalar component sums; reverse bonds are odd for M3."""
    n,i,j,xi,w,Vi,Vj,h3=args
    kp=np.zeros((n,3,3));ku=np.zeros_like(kp);mp=np.zeros((n,3,3,3));mu=np.zeros_like(mp)
    def scatter(a,b):return np.bincount(i,weights=a,minlength=n)+np.bincount(j,weights=b,minlength=n)
    for c in range(3):
        for b in range(3):
            v=w*xi[:,c]*xi[:,b]
            kp[:,c,b]=scatter(v*Vj,v*Vi);ku[:,c,b]=scatter(v*h3,v*h3)
            for d in range(3):
                t=v*xi[:,d]
                mp[:,c,d,b]=scatter(t*Vj,-t*Vi);mu[:,c,d,b]=scatter(t*h3,-t*h3)
    return kp,mp,ku,mu


def quadratic_prediction(M3,invK,H):
    return .5*np.einsum('acd,icdb,ibn->ian',H,M3,invK)


def continuous_ball(delta):
    """Independent polar Gauss reference plus radial analytic beta integral."""
    t,wr=np.polynomial.legendre.leggauss(8);r=(t+1)*delta/2;wr=wr*delta/2
    z,wz=np.polynomial.legendre.leggauss(12);phi=np.arange(24)*2*np.pi/24
    ir,iz,ip=np.indices((8,12,24)).reshape(3,-1)
    direction=np.column_stack([np.sqrt(1-z[iz]**2)*np.cos(phi[ip]),np.sqrt(1-z[iz]**2)*np.sin(phi[ip]),z[iz]])
    xi=r[ir,None]*direction;weights=wr[ir]*wz[iz]*(2*np.pi/24)*r[ir]**2*radial(r[ir]/delta)
    Kq=np.einsum('i,ic,ib->cb',weights,xi,xi)
    Mq=np.einsum('i,ic,id,ib->cdb',weights,xi,xi,xi)
    K=4*np.pi*delta**5/315*np.eye(3)
    return dict(K_analytic=K,M3_analytic=np.zeros((3,3,3)),K_numerical=Kq,M3_numerical=Mq,
        K_relative_error=float(np.linalg.norm(Kq-K)/np.linalg.norm(K)),M3_numerical_norm=float(np.linalg.norm(Mq)),
        M3_numerical_scaled=float(np.linalg.norm(Mq)/(np.linalg.norm(K)*delta)),
        omega0_analytic=2*np.pi*delta**3/15,omega0_numerical=float(weights.sum()),
        K_formula='K=(4*pi/3)*integral_0^delta r^4*(1-r/delta)^2 dr*I = 4*pi*delta^5/315*I',
        M3_proof='Radial weight and full ball invariant under xi -> -xi; every rank-3 monomial changes sign, so each integral equals its negative.',
        selected_analytic_components=dict(M3_zzx=0.,M3_zzy=0.,M3_zzz=0.),quadrature='8 radial Gauss x 12 cos(theta) Gauss x 24 azimuth trapezoid points')


def lattice_reconstruction(args):
    """Quadratic full-family displacement reconstruction at actual observation X."""
    X,xi,w,h3,invK,field,row=args
    Y=X[:,None,:]+xi[None,:,:];yi=field.position(X);yj=field.position(Y)
    du=(yj-Y)[:,:,row]-(yi-X)[:,None,row]
    A=np.einsum('ij,j,jb->ib',du,w*h3,xi)
    F=np.broadcast_to(np.eye(3),(len(X),3,3)).copy();F[:,row,:]+=A@invK
    return F-field.response(X[:,2])['F']


def psi(F,material):return svk_split_energy(green_lagrange(F),material.lame,material.mu)[0]


def node_diagnostics(F,exact,prediction,volume,material,gradient_scale,psi_scale):
    actual=F-exact;residual=actual-prediction
    actual_psi=psi(F,material);exact_psi=psi(exact,material);predicted_psi=psi(exact+prediction,material)
    return dict(actual_delta_F=actual,predicted_delta_F=prediction,prediction_residual=residual,
        actual_F_error=tensor_norm(actual),predicted_F_error=tensor_norm(prediction),F_prediction_residual=tensor_norm(residual),
        Psi_actual=actual_psi,Psi_exact=exact_psi,Psi_predicted=predicted_psi,
        Psi_actual_bias=actual_psi-exact_psi,Psi_predicted_bias=predicted_psi-exact_psi,
        Psi_causal_residual=actual_psi-predicted_psi)


def bond_diagnostics(args):
    """Block task: 3x3 tensors only, with full nonlinear spectral evaluations."""
    xi,Ft,exact,Ei,Ej,Pi,Pj,weights,W0,material,ratio_floor=args
    n=xi/np.linalg.norm(xi,axis=1)[:,None]
    projector=np.eye(3)-n[:,:,None]*n[:,None,:]
    D=.5*(Ei+Ej);pred=D@projector;moment_pred=(.5*(Pi+Pj))@projector
    actual=Ft-exact;den=tensor_norm(D);pn=tensor_norm(pred)
    ratio=np.full(len(D),np.nan);np.divide(pn,den,out=ratio,where=den>ratio_floor)
    pa=psi(Ft,material);pe=psi(exact,material);pp=psi(exact+moment_pred,material)
    W1=constitutive(exact,material)[0]
    return dict(actual_F_error=tensor_norm(actual),predicted_F_error=pn,moment_predicted_F_error=tensor_norm(moment_pred),
        projection_residual=tensor_norm(actual-pred),moment_chain_residual=tensor_norm(actual-moment_pred),D_norm=den,
        projection_reduction_ratio=ratio,Psi_actual=pa,Psi_exact=pe,Psi_predicted=pp,
        Psi_actual_bias=pa-pe,Psi_predicted_bias=pp-pe,Psi_causal_residual=pa-pp,
        U0_bond_J=weights*W0,U1_bond_J=weights*W1,Uplus0_bond_J=weights*pa,Uplus1_bond_J=weights*pe)


def weighted_metrics(actual,predicted,residual,weights,scale):
    w=weights/np.sum(weights)
    return dict(actual_max=float(np.max(np.abs(actual))),predicted_max=float(np.max(np.abs(predicted))),
        actual_RMS=float(np.sqrt(w@(actual*actual))),predicted_RMS=float(np.sqrt(w@(predicted*predicted))),
        residual_max=float(np.max(np.abs(residual))),residual_RMS=float(np.sqrt(w@(residual*residual))),
        actual_max_scaled=float(np.max(np.abs(actual))/scale),predicted_max_scaled=float(np.max(np.abs(predicted))/scale),
        residual_max_scaled=float(np.max(np.abs(residual))/scale),residual_RMS_scaled=float(np.sqrt(w@(residual*residual))/scale))


def continuum_inner_energy(field,material,length,depth,order=32):
    """Nested centered cube, same global polynomial. Returned intact and Psi+."""
    width=max(length-2*depth,0.)
    if width==0:return np.zeros(2)
    z,w=np.polynomial.legendre.leggauss(order);z*=width/2;w*=width/2
    r=field.response(z)
    return width**2*np.array([np.dot(w,r['W']),np.dot(w,psi(r['F'],material))])


def energy_summary(bonds,field,material,length):
    u0=float(bonds['U0_bond_J'].sum());u1=float(bonds['U1_bond_J'].sum())
    p0=float(bonds['Uplus0_bond_J'].sum());p1=float(bonds['Uplus1_bond_J'].sum())
    c=continuum_inner_energy(field,material,length,0,32);c64=continuum_inner_energy(field,material,length,0,64)
    def row(a,b,c):return dict(U0_J=a,U1_J=b,U2_J=float(c),kinematic_bias_J=a-b,weight_quadrature_bias_J=b-float(c),
        total_bias_J=a-float(c),kinematic_relative_to_U2=float((a-b)/c),weight_quadrature_relative_to_U2=float((b-c)/c),
        total_relative_to_U2=float((a-c)/c),sum_closure_J=float(((a-b)+(b-c))-(a-c)))
    return dict(intact=row(u0,u1,c[0]),positive=row(p0,p1,c[1]),Gauss32_vs64_relative=(abs(c-c64)/abs(c)).tolist())
