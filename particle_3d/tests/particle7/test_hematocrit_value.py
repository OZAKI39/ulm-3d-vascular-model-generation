from particle_3d.injection_population import H_D
def test_discharge(contract):
 assert H_D==contract['RBC_DISCHARGE_HEMATOCRIT']==.45
 assert contract['hematocrit_role']=='FEED_DISCHARGE_VOLUME_FLUX'
