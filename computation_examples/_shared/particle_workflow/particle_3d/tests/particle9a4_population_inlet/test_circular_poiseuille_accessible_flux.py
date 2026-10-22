import numpy as np

def test_analytic_flux_fraction_by_size(poiseuille):
 for d in [.6e-6,2.4e-6]:
  rows=[r for r in poiseuille.rows if r['diameter_m']==d];x=(2e-6-d/2-2e-9)/2e-6;p=2*x*x-x**4
  assert abs(np.mean([r['particle_id'] is not None for r in rows])-p)<.018
