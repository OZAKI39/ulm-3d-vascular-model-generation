"""Conservative cKDTree skin list; exact original sorted narrow-phase decisions."""
from time import perf_counter
import numpy as np
from scipy.spatial import cKDTree


class ContactNeighborCache:
    def __init__(self, skin_m, maximum_gap_m):
        if not np.isfinite(skin_m) or skin_m<=0 or maximum_gap_m<0:
            raise ValueError("Positive finite skin and nonnegative maximum gap required")
        self.skin=float(skin_m)
        self.maximum_gap=float(maximum_gap_m)
        self.anchor=self.radii=None
        self.pairs=np.empty((0,2),int)
        self.calls=self.rebuilds=self.reuses=0
        self.maximum_candidates=0

    def candidates(self, x, radii, gap_m=0., timings=None):
        begin=perf_counter()
        x,radii=np.asarray(x),np.asarray(radii)
        if len(x)!=len(radii) or not np.isfinite(x).all() or not np.isfinite(radii).all():
            raise ValueError("Finite matching neighbor coordinates and radii required")
        rebuild=(self.anchor is None or self.anchor.shape!=x.shape
                 or not np.array_equal(radii,self.radii) or gap_m>self.maximum_gap)
        # Numerical guard enlarges/rebuilds broad phase only; narrow phase unchanged.
        guard=32*np.finfo(float).eps*(float(np.max(np.abs(x)))+float(radii.max())+self.skin) if len(x) else 0.
        if not rebuild and len(x):
            displacement=float(np.linalg.norm(x-self.anchor,axis=1).max())
            rebuild=2*displacement>=self.skin-guard
        if rebuild:
            self.maximum_gap=max(self.maximum_gap,float(gap_m))
            cutoff=2*float(radii.max())+self.maximum_gap+self.skin+guard if len(radii) else 0.
            self.pairs=cKDTree(x).query_pairs(np.nextafter(cutoff,np.inf),output_type="ndarray")
            if len(self.pairs):
                self.pairs=self.pairs[np.lexsort((self.pairs[:,1],self.pairs[:,0]))]
            self.anchor,self.radii=x.copy(),radii.copy()
            self.rebuilds+=1
            self.maximum_candidates=max(self.maximum_candidates,len(self.pairs))
        else:
            self.reuses+=1
        pairs=self.pairs
        if len(pairs):
            i,j=pairs.T
            distances=np.linalg.norm(x[j]-x[i],axis=1)
            pairs=pairs[distances<=radii[i]+radii[j]+gap_m]
        self.calls+=1
        if timings is not None:
            timings["contact_search_s"]+=perf_counter()-begin
            timings["contact_search_calls"]+=1
        return pairs.reshape(-1,2)
