"""Independent common-geometry reference; P2 CFD remains scientifically unqualified."""
from pathlib import Path
import numpy as np

EDGES=np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])
DL=np.vstack([np.eye(3),-np.ones(3)])
RN=np.vstack([np.eye(3),np.zeros((1,3))])
RN10=np.vstack([RN,RN[EDGES].mean(1)])
HP2=np.concatenate([4*np.einsum('ai,aj->aij',DL,DL),
    np.array([4*(np.outer(DL[i],DL[j])+np.outer(DL[j],DL[i])) for i,j in EDGES])])
PACK=[(0,0),(1,1),(2,2),(0,1),(1,2),(0,2)]


def basis(xi,n):
    L=np.r_[xi,1-np.sum(xi)]
    if n==4:return L,DL,np.zeros((4,3,3)),L
    N=np.r_[L*(2*L-1),4*L[EDGES[:,0]]*L[EDGES[:,1]]]
    D=np.vstack([(4*L-1)[:,None]*DL,4*(L[EDGES[:,0],None]*DL[EDGES[:,1]]+L[EDGES[:,1],None]*DL[EDGES[:,0]])])
    return N,D,HP2,L

def geometry(xi,X):
    n=len(X);N,D,H,P=basis(xi,n);Y=X-X[3]
    J=Y.T@D;inv=np.linalg.inv(J);det=np.linalg.det(J)
    Dx=D@inv;Dp=DL@inv
    geomH=np.einsum('ai,ajk->ijk',Y,H)
    Hx=np.array([inv.T@(H[a]-np.einsum('i,ijk->jk',Dx[a],geomH))@inv for a in range(n)])
    return dict(N=N,P=P,D=Dx,Dp=Dp,H=np.array([[h[i,j] for i,j in PACK] for h in Hx]),
        J=J,inv=inv,det=det,metric=inv.T@inv,x=N@Y+X[3])

def unpack(rows,n):
    rows=np.atleast_2d(rows);k=9
    out=dict(opt=rows[:,0].astype(int),g=rows[:,1].astype(int),wref=rows[:,4],det=rows[:,5],xi=rows[:,6:9])
    for name,shape in [('N',(n,)),('P',(4,)),('D',(n,3)),('Dp',(4,3)),('H',(n,6)),('metric',(3,3))]:
        count=int(np.prod(shape));out[name]=rows[:,k:k+count].reshape((len(rows),)+shape);k+=count
    assert k==rows.shape[1]
    return out

def read_RK(prefix,n):
    lines=Path(str(prefix)+'_RK.txt').read_text().splitlines()
    return np.fromstring(lines[0],sep=' ').reshape(4,n),np.fromstring(lines[1],sep=' ').reshape(16,n,n)

def fields(X,label,H):
    n=len(X);d=X-X[3];y=np.zeros((n,4));a=np.zeros((n,4));bf=np.zeros((n,3))
    y[:,3]=.3+d[:,2]/H
    if label=='A':y[:,:3]=d
    if label=='B':y[:,0]=d[:,0];y[:,1]=-d[:,1]
    if label=='C':y[:,2]=d[:,2]**2
    if label=='P':y[:,3]=.3+d@np.array([1.,2.,3.])/H
    if label=='R':
        rng=np.random.default_rng(61006+n);y[:,:3]=rng.normal(size=(n,3))*H
        y[:,3]=rng.normal(size=n);a[:,:3]=rng.normal(size=(n,3))*.2;bf=rng.normal(size=(n,3))*.1
    return y,a,bf
