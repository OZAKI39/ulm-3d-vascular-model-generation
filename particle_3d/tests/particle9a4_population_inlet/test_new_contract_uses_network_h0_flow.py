from .conftest import ROOT
from particle_3d.population_inlet_p9a4 import make_contract,make_source
from particle_3d.continuous_infusion import NEW_FLOW_SHA256,sha256

def test_new_contract_and_concentration(real_env):
 c=make_contract(ROOT,real_env)
 assert c['flow_SHA']==NEW_FLOW_SHA256==sha256(real_env.new_flow_path)
 assert c['Q_in_m3_s']==1.5513591604402322e-14
 assert c['C_modeled_D_le_4um_m3']==c['C_total_m3']*c['F_original_4um']
 assert c['C_modeled_D_le_4um_m3']<c['C_total_m3']
 assert c['source_distribution']=='SONOVUE_D_LE_4UM_CONDITIONAL'
 assert c['concentration_provenance']['measured'] is False
