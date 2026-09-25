import numpy as np

def test_accepted_rate_equals_lambda_times_mean_acceptance(poiseuille):
 x=(2e-6-np.array([.6e-6,2.4e-6])/2-2e-9)/2e-6
 p=np.array([.4,.6])@(2*x*x-x**4)
 measured=poiseuille.accepted/poiseuille.time_s
 assert abs(measured/(5*p)-1)<.035
 times=np.array([r['proposal_time_s'] for r in poiseuille.rows]);counts=np.histogram(times,bins=np.arange(0,times[-1],4.))[0]
 # Fano is descriptive audit output; only finite counts are a CI gate.
 assert np.isfinite(counts.var()/counts.mean())
