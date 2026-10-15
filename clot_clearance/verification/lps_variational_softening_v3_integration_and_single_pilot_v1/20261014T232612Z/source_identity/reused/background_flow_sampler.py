"""Native P1 field interpolation and strict tetra membership; no extrapolation."""
import numpy as np
import vtk
from .y_bifurcation_sources import YSources

class BackgroundFlowSampler:
    def __init__(self,sources=None):
        self.sources=sources or YSources();self.queries=0;self.rejected=0;self.fallbacks=0

    def locate(self,positions):
        s=self.sources;X=np.atleast_2d(np.asarray(positions,float));self.queries+=len(X)
        ids=np.array([s.locator.FindCell(p) for p in X],int)
        def bary(p,ii):
            w=np.einsum('nij,nj->ni',s.inverse[ii],p-s.origins[ii]);return np.c_[1-w.sum(1),w]
        weights=bary(X,np.maximum(ids,0));bad=(ids<0)|(weights.min(1)<-1e-10)|(weights.max(1)>1+1e-10)
        for i in np.flatnonzero(bad):
            found=False;candidate=vtk.vtkIdList();eps=1e-13
            s.locator.FindCellsWithinBounds(np.c_[X[i]-eps,X[i]+eps].ravel(),candidate)
            cells=np.array([candidate.GetId(j) for j in range(candidate.GetNumberOfIds())],int)
            if len(cells):
                w=bary(X[i:i+1],cells);good=(w.min(1)>=-1e-10)&(w.max(1)<=1+1e-10)
                if good.any():j=np.flatnonzero(good)[0];ids[i]=cells[j];weights[i]=w[j];found=True;self.fallbacks+=1
            if not found:ids[i]=-1;weights[i]=np.nan
        self.rejected+=int((ids<0).sum());return ids,weights

    def sample(self,positions):
        ids,w=self.locate(positions);valid=ids>=0;s=self.sources
        u=np.full((len(ids),3),np.nan);p=np.full(len(ids),np.nan)
        u[valid]=np.einsum('ni,nij->nj',w[valid],s.u[s.tet[ids[valid]]]);p[valid]=np.einsum('ni,ni->n',w[valid],s.p[s.tet[ids[valid]]])
        return dict(velocity_m_s=u,pressure_Pa=p,valid=valid,tetra_id=ids,barycentric_weights=w)
