import numpy as np

def test_finite_size_changes_marginal(poiseuille):
 d=np.array([r['diameter_m'] for r in poiseuille.rows]);keep=np.array([r['particle_id'] is not None for r in poiseuille.rows])
 assert d[keep].mean()<d.mean()-.1e-6
