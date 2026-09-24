import numpy as np
from scipy.stats import chi2
def test_actual(real,report):
 s=real[5]['INLET']; x=np.loadtxt(report/'data/02_actual_inlet_samples.csv',delimiter=',',skiprows=1)
 assert len(x)>=100000
 counts=np.bincount(x[:,3].astype(int),minlength=len(s.weights)); expected=s.weights/s.weights.sum()*len(x)
 select=expected>=5; statistic=np.sum((counts[select]-expected[select])**2/expected[select])
 assert chi2.sf(statistic,select.sum()-1)>1e-5
 assert np.linalg.norm(x[:,:3].mean(0)-s.expectation())<2e-8
