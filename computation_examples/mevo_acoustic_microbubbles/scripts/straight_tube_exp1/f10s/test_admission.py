"""Meaningful regression: legacy 2ms identity and false partial 3ms rejection."""
import copy
import unittest
from .common import R, read, policy
from .admission import layered
from scripts.straight_tube_exp1.f10p.analyze import layered as old
class Coverage(unittest.TestCase):
 def test_legacy_and_partial(self):
  p=policy(); record=read(R/'navigation_v1/runs/X_PAIR_DT1/RESULT.json');record['source']='X_PAIR'
  comps=read(R/'analysis/COMPARISONS.json')
  # Select actual saved X pair numerical comparisons, retain original metric values.
  c=[v for k,v in comps.items() if k in ['X_PAIR_DT1__DT05','X_PAIR_DT1__FINE']]
  self.assertEqual(len(c),2)
  t=[{'initial_network_persists_entire_2ms':True}]*3
  legacy=old(record,*c,t,p)
  cc=[dict(v,complete_window=v['complete_2ms']) for v in c]
  tt=[{'initial_network_persists_entire_window':True}]*3
  generalized=layered(record,*cc,tt,p,.002)
  legacy.pop('topology_scope');generalized.pop('topology_scope')
  self.assertEqual(legacy,generalized)
  partial=layered(record,*cc,tt,p,.003)
  self.assertEqual(partial['BASIC_NUMERICS_PASS'],'FAIL')
  rr=dict(record,completed_s=.003)
  for v in cc:v['complete_window']=False
  incomplete=layered(rr,*cc,tt,p,.003)
  self.assertEqual(incomplete['GROUP_POSITION_NUMERICALLY_ADMITTED'],'FAIL')
  self.assertEqual(incomplete['PARTICLE_ID_NUMERICALLY_ADMITTED'],'FAIL')
if __name__=='__main__':unittest.main()
