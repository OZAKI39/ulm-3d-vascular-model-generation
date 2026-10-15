from common import *
import importlib.util
spec=importlib.util.spec_from_file_location('frozen_softening_reference',R/'analysis/softening_feasibility_v1/model.py')
reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
class TypedLPS(reference.SofteningLPS):
 def __init__(self,ref):
  super().__init__(ref,K=K,G=G);self.on=S0[ref['bond_types']];self.complete=SC[ref['bond_types']]
 def update(self,x,kappa):
  e,_=self.geometry(x);new=np.maximum(kappa,np.maximum(e/self.r,0.))
  d=np.zeros(self.nb);part=(new>self.on)&(new<self.complete);full=new>=self.complete
  d[part]=self.complete[part]*(new[part]-self.on[part])/(new[part]*(self.complete[part]-self.on[part]));d[full]=1
  return new,d
