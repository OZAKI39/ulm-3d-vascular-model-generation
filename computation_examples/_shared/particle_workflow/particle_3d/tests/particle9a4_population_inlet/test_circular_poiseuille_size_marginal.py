import numpy as np

def test_flux_weighted_size_probability(poiseuille):
 d=np.array([r['diameter_m'] for r in poiseuille.rows if r['particle_id'] is not None]);x=(2e-6-np.array([.6e-6,2.4e-6])/2-2e-9)/2e-6
 weights=np.array([.4,.6])*(2*x*x-x**4);expected=weights[0]/weights.sum()
 assert abs(np.mean(d==.6e-6)-expected)<.018
