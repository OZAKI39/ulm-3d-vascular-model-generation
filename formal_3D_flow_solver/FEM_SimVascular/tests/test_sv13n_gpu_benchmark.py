from sv13n_support import *
def test_gpu_median_recomputed():
 d=accepted('gpu_benchmark');assert timing_gate([read(n+'_acceptance') for n in d['names']])==d['median_s']
