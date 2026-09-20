from sv13m_support import *
from sv_validation.sv13m import *
def test_transfer_attribution_explicit():
 d=accepted('ghost_transfer');assert all(k in d for k in ('ghost','MatMult','PC','growth_relation','limitations'))
