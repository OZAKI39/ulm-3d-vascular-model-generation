from sv13m_support import *
from sv_validation.sv13m import *
def test_reverse_including_nonzero_remote_contributions():
 d=ghost();ghost_run_gate(d)
 assert any(ph=='reverse' and float(a)>100 for ph,rank,i,a,e in d['values'])
def test_forward_pass_reverse_fail_rejected():
 d=ghost();d['checks'][0][2]='1'
 with pytest.raises(GateError):ghost_run_gate(d)
def test_repair02_error_suppression_not_accepted():
 d=next(r for r in read('ghost_probe_repair_02')['runs'] if r['variant']=='CUDA_MPI')
 with pytest.raises(GateError):ghost_run_gate(d)
