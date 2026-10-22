import json

def test_saved_real_radius_exact(real_saved):
 a=real_saved['initialization']['original_particle']['radius_m']
 assert all(p['radius_m']==a for p in real_saved['positions'])
 assert real_saved['positions'][0]['center_m']==real_saved['initialization']['original_particle']['center_m']
