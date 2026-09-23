from stage016_helpers import *
import pytest
from fem3d.planar_port import volume_gates

def test_malicious_good_quality_tetra_over_200000_rejected():
 v,q,p=good_volume();p['N_tetra']=200001;r=volume_gates(v,q,p,POLICY)
 assert r['status']=='FAIL' and not r['checks']['tetra_budget'] and r['checks']['P1']

def test_malicious_excellent_quality_P2_over_1_35_rejected():
 v,q,p=good_volume();p['N_P2_velocity_proxy']=int(POLICY['cost']['baseline']['N_P2_velocity_proxy']*1.35)+1;r=volume_gates(v,q,p,POLICY)
 assert r['status']=='FAIL' and not r['checks']['P2_budget'] and r['checks']['median']

def test_good_quality_and_budget_allowed():
 v,q,p=good_volume();assert volume_gates(v,q,p,POLICY)['status']=='PASS'
