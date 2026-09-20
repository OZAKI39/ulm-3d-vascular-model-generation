from sv13m_support import *
from sv_validation.sv13m import *
def test_three_actual_runs():
 rows=[r for r in read('ghost_probe_after')['runs'] if r['variant']=='CPU_seq']
 assert len(rows)==3
 for r in rows:ghost_run_gate(r)
def test_cpu_regression_must_be_rejected():
 d=ghost('CPU_MPI');d['values'][0][-2]='-100'
 with pytest.raises(GateError):ghost_run_gate(d)
