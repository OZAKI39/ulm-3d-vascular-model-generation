"""Original P1 particle-field diagnostics in SI."""
import numpy as np

def curl(g):
    return np.stack([g[:,2,1]-g[:,1,2],g[:,0,2]-g[:,2,0],g[:,1,0]-g[:,0,1]],axis=1)

def shear(g,n):
    gn=np.einsum('nij,nj->ni',g,n)
    return gn-n*np.sum(n*gn,axis=1)[:,None]

def errors(actual,exact):
    delta=(actual-exact).reshape(len(actual),-1)
    norm=float(np.linalg.norm(delta));den=float(np.linalg.norm(exact))
    return dict(relative_L2=norm/den if den else None,absolute_L2=norm,
                absolute_RMS=norm/np.sqrt(len(actual)),
                max_absolute_error=float(np.linalg.norm(delta,axis=1).max()),
                reference_L2=den)
