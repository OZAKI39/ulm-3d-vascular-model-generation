"""Static exact-geometry acceleration data; CPU preprocessing only."""
import numpy as np


def build_bvh(lower,upper,leaf_size=16):
    lower=np.asarray(lower,dtype=float);upper=np.asarray(upper,dtype=float)
    centers=.5*(lower+upper);order=np.arange(len(lower),dtype=np.int64);bounds=[];links=[]
    def build(start,end):
        node=len(links);links.append([-1,-1,start,end-start])
        ids=order[start:end];lo=lower[ids].min(0);hi=upper[ids].max(0);bounds.append(np.r_[lo,hi])
        if end-start>leaf_size:
            axis=int(np.argmax(hi-lo));mid=(start+end)//2
            permutation=np.argpartition(centers[ids,axis],mid-start)
            order[start:end]=ids[permutation]
            left=build(start,mid);right=build(mid,end);links[node]=[left,right,0,0]
        return node
    build(0,len(order))
    return np.asarray(bounds),np.asarray(links,dtype=np.int64),order


def tetra_adjacency(tetra):
    tetra=np.asarray(tetra,dtype=np.int64);n=len(tetra)
    faces=np.stack([np.sort(tetra[:,[j for j in range(4) if j!=k]],axis=1) for k in range(4)],axis=1).reshape(-1,3)
    order=np.lexsort((faces[:,2],faces[:,1],faces[:,0]));sorted_faces=faces[order]
    same=np.all(sorted_faces[1:]==sorted_faces[:-1],axis=1);where=np.flatnonzero(same)
    if np.any(np.diff(where)==1):raise ValueError('NONMANIFOLD_TETRA_FACE')
    a,b=order[where],order[where+1]
    result=np.full(n*4,-1,dtype=np.int64);result[a]=b//4;result[b]=a//4
    return result.reshape(n,4)
