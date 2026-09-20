from sv13m_support import *
from sv_validation.sv13m import *
def test_observed_repeatability():
 d=accepted('benchmark')
 for runs in d['runs'].values():timing_gate(runs)
 assert d['speedup']==min(d['medians']['CPU1'],d['medians']['CPU4'])/d['medians']['GPU1']
def test_single_sample_no_conclusion():
 with pytest.raises(GateError):timing_gate([dict(wall_time_s=1)])
def test_profile_timing_rejected():
 row=dict(accepted=True,steps_completed=20,initial_state='t=0',wall_time_s=1,profiling=True)
 with pytest.raises(GateError):timing_gate([row,row])
def test_variability_requires_third_and_median_uses_all():
 row=dict(accepted=True,steps_completed=20,initial_state='t=0',profiling=False)
 with pytest.raises(GateError):timing_gate([dict(row,wall_time_s=10),dict(row,wall_time_s=12)])
 assert timing_gate([dict(row,wall_time_s=10),dict(row,wall_time_s=12),dict(row,wall_time_s=11)])==11
