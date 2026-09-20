from sv13n_support import *
def test_fallback_justified():
 d=accepted('cuda_fallback_policy');assert not d['main_used']
 if d['fallback_used']:assert d['G1_failure_evidence'] and d['release_fix_audit']
