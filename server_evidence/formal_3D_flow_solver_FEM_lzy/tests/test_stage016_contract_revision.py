from stage016_helpers import *
import pytest
from fem3d.audit import sha256
from fem3d.planar_port import verify_plane

def test_derived_contract_references_unchanged_stage0():
 assert CONTRACT['version']==2 and not CONTRACT['legacy_scalar_area_is_gate']
 assert sha256(ROOT/'reports/stage00/source_contract.json')==CONTRACT['source_contract_sha256']
 old=read(ROOT/'reports/stage00/source_contract.json')
 for p in old['inlets']+old['outlets']:
  n=CONTRACT['ports'][p['name']]
  assert n['plane_origin_m']==p['plane_origin_m'] and n['outward_normal']==p['outward_normal']
 assert CONTRACT['facet_names']=={'1':'WALL','2':'OUTLET_03','3':'OUTLET_01','4':'INLET','5':'OUTLET_02'}

@pytest.mark.parametrize('field',['origin','normal','basis'])
def test_malicious_plane_identity_change_rejected(field):
 p=CONTRACT['ports']['inlet']; o=np.array(p['plane_origin_m']);n=np.array(p['outward_normal']);b=np.array(p['basis'])
 if field=='origin':o[0]+=1e-15
 if field=='normal':n[0]+=1e-6
 if field=='basis':b[0,0]+=1e-6
 with pytest.raises(ValueError,match='Frozen plane'):verify_plane(p,o,n,b)
