"""Change only the direction operand in the frozen historical load formula."""
from common import *
_historical=module_from_file('direction_historical_external',HIST_A/'external.py')
HistoricalExternal=_historical.External
FIELDS=_historical.FIELDS
CONTACT=_historical.CONTACT
class External(HistoricalExternal):
 def __init__(self,ns,restart_state=None,query=None):
  super().__init__(ns,restart_state,query)
  normal=load_npz(O/'REFERENCE_NORMAL_FIELD.npz');ids=normal['atom_ids']-1
  assert np.array_equal(ids,np.flatnonzero(ns.ref['areas']>0))
  direction=ns.ref['tangent'].copy()
  direction[ids]=normal['normal_ref']*np.linalg.norm(ns.ref['tangent'][ids],axis=1)[:,None]
  self.ref=dict(ns.ref,tangent=direction)
  # amp(), compute(), callback(), accept(), drag and wall contact are inherited
  # unchanged. ns.ref and its original tangent/areas/geometry remain immutable.
