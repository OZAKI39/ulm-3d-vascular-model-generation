from sv12_support import accepted
from sv_validation.sv11 import measured_mass_gate

def test_accepted_final_mass():
    final=accepted();assert measured_mass_gate(final)
    assert final['native_flux_difference_over_Q']<=1e-6
