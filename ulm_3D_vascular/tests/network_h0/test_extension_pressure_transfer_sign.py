import numpy as np
import pytest
from network_1d0d.boundary_conditions import cap_pressures

@pytest.mark.parametrize('q',[3e-15,-3e-15,0.])
def test_transfer_respects_outward_signed_flow(q):
    # Positive flow loses pressure toward the distal cap; reversal increases it.
    p=float(cap_pressures(50.,q,2e16))
    assert p==50.-2e16*q
    if q>0: assert p<50.
    if q<0: assert p>50.

def test_real_o3_terminal_is_not_an_invented_downstream_resistance(bc):
    o3=bc['ports'][2]
    assert o3['pressure_realcut_Pa']==0
    assert o3['pressure_cap_raw_Pa']==-o3['extension_deltaP_Pa']
