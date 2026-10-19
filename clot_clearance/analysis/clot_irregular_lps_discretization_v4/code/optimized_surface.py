"""Exact caches for fixed reference queries and current-state containment data."""
import numpy as np
class ReferenceQueryCache:
    def __init__(self,tree):self.tree=tree;self.cache={};self.hits=0;self.misses=0
    def query_ball_point(self,points,radii):
        answers=[]
        for point,radius in zip(points,np.broadcast_to(radii,(len(points),))):
            key=(*point,float(radius))
            if key not in self.cache:
                self.cache[key]=self.tree.query_ball_point(point,radius);self.misses+=1
            else:self.hits+=1
            answers.append(self.cache[key])
        return answers
    def __getattr__(self,name):return getattr(self.tree,name)
def install(module):
    original_init=module.CurrentSurface.__init__
    def init(self,*a,**k):
        original_init(self,*a,**k);self.reference_tree=ReferenceQueryCache(self.reference_tree)
    module.CurrentSurface.__init__=init
    original_inside=module.OldSurface.occlusion_checker
    def inside(self,x,solid):
        previous=getattr(self,'_inside_exact_cache',None)
        if previous is not None and np.array_equal(previous[0],x) and np.array_equal(previous[1],solid):return previous[2]
        result=original_inside(self,x,solid);self._inside_exact_cache=(x.copy(),solid.copy(),result);return result
    module.OldSurface.occlusion_checker=inside
def optimize_lineage_source(code):
    old="""   for old,oldset in previous_free.items():
    children=[cid for _,cid,ids in current if oldset.intersection((ids+1).tolist())]"""
    new="""   # A partition can only split after irreversible edge removal. Map
   # old members directly to current labels instead of all component pairs.
   current_cid=np.zeros(len(sizes),np.int64)
   for label,cid,ids in current:current_cid[label]=cid
   for old,oldset in previous_free.items():
    child_labels=np.unique(lab[np.fromiter(oldset,dtype=np.int64)-1])
    children=current_cid[child_labels].tolist()"""
    assert code.count(old)==1
    return code.replace(old,new)
