from sv_validation.validation import mass_balance

def test_actual_field_mass_conservation(measured):
    q=measured[3]
    mass_balance(q["Q_target_m3_s"],-q["Q_in_m3_s"],q["outlet_flows_m3_s"])
