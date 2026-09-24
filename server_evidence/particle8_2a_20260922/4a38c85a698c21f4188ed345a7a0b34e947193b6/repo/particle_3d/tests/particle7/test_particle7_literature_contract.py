def test_sources(contract):
 a,b,c=contract['sources']
 assert (a['male_hct_percent'],a['female_hct_percent'])==(45.7,44.8)
 assert b['iv_bubble_dose']==20000000 and b['body_mass_g']==32.6
 assert c['mouse_blood_volume_ml_kg']==72
 assert abs(b['iv_bubble_dose']/(b['body_mass_g']/1000*c['mouse_blood_volume_ml_kg'])/8.5e6-1)<.003
