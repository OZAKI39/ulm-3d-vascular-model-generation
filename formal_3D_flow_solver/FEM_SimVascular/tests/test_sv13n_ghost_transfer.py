from sv13n_support import *
def test_transfers_are_measured_or_explicitly_unknown():
 d=accepted('ghost_transfer');assert d['evidence']
 for k in ('H2D','D2H','ghost_update','per_KSP_iteration'):assert k in d['measurements']
