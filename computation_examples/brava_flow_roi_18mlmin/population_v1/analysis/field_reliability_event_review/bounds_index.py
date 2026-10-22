"""Conservative scipy bounding-box candidate index, never a field interpolator.

If two AABBs intersect, their centers are at most the sum of their half-diagonal
lengths apart. A global maximum cell half diagonal therefore makes the tree
ball query conservative; the subsequent exact AABB predicate removes extras.
The existing exact barycentric segment clipping remains authoritative.
"""
import numpy as np
from scipy.spatial import cKDTree

class BoundsIndex:
    def __init__(self,points,tetra):
        xyz=np.asarray(points,float)[tetra]
        self.lower=xyz.min(1);self.upper=xyz.max(1)
        self.centers=(self.lower+self.upper)/2
        self.radius=float(np.linalg.norm((self.upper-self.lower)/2,axis=1).max())
        self.roundoff=128*np.finfo(float).eps*float(np.abs(xyz).max())
        self.tree=cKDTree(self.centers)

    def candidates(self,bounds):
        b=np.asarray(bounds,float).reshape(3,2);lo=b[:,0];hi=b[:,1]
        ids=np.asarray(self.tree.query_ball_point((lo+hi)/2,np.linalg.norm((hi-lo)/2)+self.radius+self.roundoff),int)
        return np.sort(ids[np.all(self.lower[ids]<=hi,axis=1)&np.all(self.upper[ids]>=lo,axis=1)])

    def FindCellsWithinBounds(self,bounds,result):
        result.Reset()
        for cid in self.candidates(bounds):result.InsertNextId(int(cid))
