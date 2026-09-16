"""Independent SI reference. See REFERENCE_CONTRACT and THEORY_DERIVATION.

Scalar source: Jeffrey & Onishi (1984) §§3–7; Townsend (2023), (10–12),
(16), (19), (24), (33–34a). Transverse block: Radhakrishnan (2017),
§6 resistance functions, DOI 10.5281/zenodo.1137305. The finite-gap
objective completion is a declared project model, not an exact Stokes pair.
"""
import itertools
import numpy as np
MU=1e-3

def coefficients(a,b,h,regularize=True):
    a,b,h=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float),np.asarray(h,float))
    s=a+b; re=a*b/s; he=np.maximum(h,1e-3*re) if regularize else h
    xi=2*he/s; lam=b/a; t=1+lam; L=np.log(1/xi)
    g1=2*lam**2/t**3
    g2=lam*(1+7*lam+lam**2)/(5*t**3)
    g3=(1+18*lam-29*lam**2+18*lam**3+lam**4)/(42*t**3)
    X=6*np.pi*MU*a*(g1/xi+g2*L+g3*xi*L)
    ya2=4*lam*(2+lam+2*lam**2)/(15*t**3)
    ya3=2*(16-45*lam+58*lam**2-45*lam**3+16*lam**4)/(375*t**3)
    Y0=6*np.pi*MU*a*ya2*L; Y=Y0+6*np.pi*MU*a*ya3*xi*L
    Bi=4*np.pi*MU*a*a*lam*(4+lam)/(5*t*t)*L
    Bj=4*np.pi*MU*b*b*(1/lam)*(4+1/lam)/(5*(1+1/lam)**2)*L
    Cii=8*np.pi*MU*a**3*2*lam/(5*t)*L
    Cij=8*np.pi*MU*a**3*lam**2/(10*t)*L
    # Work-conjugate objective completion: lever sum is actual separation r.
    # Ratio Bi/Y0 is not the particle radius for unequal spheres.
    r=s+h; li=r*Bi/(Bi+Bj); lj=r*Bj/(Bi+Bj)
    pump=(Cii-Bi**2/Y0)*(r/s)**2
    return np.stack([X,Y,pump,np.zeros_like(X),li,lj,Y0,Bi,Bj,Cii,Cij],axis=-1)

def crossmat(n):
    x,y,z=n
    return np.array([[0,-z,y],[z,0,-x],[-y,x,0.]])

def pair_matrix(a,b,xij):
    r=np.linalg.norm(xij);n=xij/r;h=r-a-b;c=coefficients(a,b,h);X,Y,P,_,li,lj=c[:6]
    T=np.eye(3)-np.outer(n,n); A=X*np.outer(n,n)+Y*T; C=crossmat(n)
    # Explicit 3x3 blocks, independent of production's modal B^T D B assembly.
    M=np.zeros((12,12)); vi=slice(0,3);wi=slice(3,6);vj=slice(6,9);wj=slice(9,12)
    M[vi,vi]=A;M[vj,vj]=A;M[vi,vj]=-A;M[vj,vi]=-A
    M[vi,wi]=Y*li*C; M[vi,wj]=Y*lj*C; M[vj,wi]=-Y*li*C;M[vj,wj]=-Y*lj*C
    M[wi,vi]=M[vi,wi].T;M[wj,vi]=M[vi,wj].T;M[wi,vj]=M[vj,wi].T;M[wj,vj]=M[vj,wj].T
    M[wi,wi]=(Y*li**2+P)*T;M[wj,wj]=(Y*lj**2+P)*T
    M[wi,wj]=(Y*li*lj-P)*T;M[wj,wi]=M[wi,wj].T
    return M,c

def assemble(x,a,U,Omega,E,pairs):
    N=len(a); q0=np.column_stack([U,Omega]).reshape(-1)
    drag=np.column_stack([np.repeat((6*np.pi*MU*a)[:,None],3,axis=1),np.repeat((8*np.pi*MU*a**3)[:,None],3,axis=1)]).reshape(-1)
    R=np.diag(drag);rhs=drag*q0;active=[]
    for i,j in pairs:
        rvec=x[i]-x[j];r=np.linalg.norm(rvec);h=r-a[i]-a[j]
        if h>=.2*a[i]*a[j]/(a[i]+a[j]):continue
        M,c=pair_matrix(a[i],a[j],rvec); ids=np.r_[np.arange(6*i,6*i+6),np.arange(6*j,6*j+6)]
        n=rvec/r;A=c[0]*np.outer(n,n)+c[1]*(np.eye(3)-np.outer(n,n));correction=(a[i]+a[j])*.5*(E[i]+E[j])@n
        f=A@correction; torque_i=-c[4]*np.cross(n,f);torque_j=-c[5]*np.cross(n,f)
        strain_rhs=np.r_[f,torque_i,-f,torque_j]
        R[np.ix_(ids,ids)]+=M;rhs[ids]+=M@q0[ids]-strain_rhs
        active.append((i,j))
    return R,rhs,active

def dense_solve(x,a,U,Omega,E,pairs,base=None,dt=0,driving=None):
    R,rhs,active=assemble(x,a,U,Omega,E,pairs)
    if driving is not None:
        drag=np.column_stack([np.repeat((6*np.pi*MU*a)[:,None],3,axis=1),np.repeat((8*np.pi*MU*a**3)[:,None],3,axis=1)]).reshape(-1)
        rhs+=drag*np.asarray(driving).reshape(-1)
    scale=np.column_stack([np.ones((len(a),3)),1/np.repeat(a[:,None],3,axis=1)]).reshape(-1)
    A=R*scale[:,None]*scale[None,:]/1e-8;b=rhs*scale/1e-8
    y=np.linalg.solve(A,b)
    if base is not None and dt>0:
        J=[];bounds=[]
        for i,j in pairs:
            d=base[i]-base[j];r=np.linalg.norm(d);h=r-a[i]-a[j]
            if h>.2*a[i]*a[j]/(a[i]+a[j]):continue
            row=np.zeros(6*len(a));row[6*i:6*i+3]=d/r;row[6*j:6*j+3]=-d/r;J.append(row);bounds.append(-max(h,0)/dt)
        if J:
            J=np.array(J);bounds=np.array(bounds);m=len(J)
            if np.min(J@y-bounds)<-1e-14:
                # Independently enumerate the small validation systems' active faces.
                if m>12:raise RuntimeError('Independent constraint enumerator capped at 12 active candidates')
                found=False
                for count in range(1,m+1):
                    for face in itertools.combinations(range(m),count):
                        B=J[list(face)];K=np.block([[A,-B.T],[B,np.zeros((count,count))]])
                        try:z=np.linalg.solve(K,np.r_[b,bounds[list(face)]])
                        except np.linalg.LinAlgError:continue
                        if np.min(z[len(b):])>=-1e-12 and np.min(J@z[:len(b)]-bounds)>=-1e-12:
                            y=z[:len(b)];found=True;break
                    if found:break
                if not found:raise RuntimeError('No admissible KKT face')
    return (y*scale).reshape(-1,6)

class FlowGradientReference:
    def __init__(self,path):
        import h5py
        with h5py.File(path,'r') as f:
            self.dims=f['dims'][...];self.origin=f['origin_m'][...];self.dx=float(f['dx_m'][()]);self.ids=f['linear_index'][...];self.fluid=f['fluid_linear_index'][...];self.u=f['Velocity_m_s'][...]
    def query(self,x):
        x=np.atleast_2d(x);q=(x-self.origin)/self.dx
        if not np.isfinite(q).all() or (q<0).any() or (q>=self.dims-1).any():raise ValueError('OUTSIDE_OR_NONFINITE')
        base=np.floor(q).astype(int);t=q-base; vel=np.zeros((len(x),3));grad=np.zeros((len(x),3,3))
        corners=list(itertools.product((0,1),repeat=3))
        values=[]
        for corner in corners:
            idx=base+corner;ids=idx[:,0]+self.dims[0]*(idx[:,1]+self.dims[1]*idx[:,2]);slots=np.searchsorted(self.ids,ids)
            if np.any(slots>=len(self.ids)) or np.any(self.ids[np.minimum(slots,len(self.ids)-1)]!=ids):raise ValueError('MISSING_OR_SOLID')
            vals=self.u[slots];values.append(vals);factors=np.where(np.array(corner),t,1-t);vel+=np.prod(factors,axis=1)[:,None]*vals
        # Derivative as interpolated edge differences rather than summation of signed corner weights.
        table={corner:value for corner,value in zip(corners,values)}
        for axis in range(3):
            other=[k for k in range(3) if k!=axis]
            for bits in itertools.product((0,1),repeat=2):
                lo=[0,0,0];hi=[0,0,0];hi[axis]=1;weight=np.ones(len(x))
                for k,bit in zip(other,bits):lo[k]=hi[k]=bit;weight*=t[:,k] if bit else 1-t[:,k]
                grad[:,:,axis]+=weight[:,None]*(table[tuple(hi)]-table[tuple(lo)])/self.dx
        omega=.5*np.column_stack([grad[:,2,1]-grad[:,1,2],grad[:,0,2]-grad[:,2,0],grad[:,1,0]-grad[:,0,1]])
        return vel,omega,.5*(grad+grad.transpose(0,2,1)),grad
