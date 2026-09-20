from sv13m_support import *
from sv_validation.sv13m import *
def test_once_actual_cuda_and_ghost_error():
 d=accepted('baseline_ghost_failure');assert d['run_count']==1 and d['reproduced']
 assert 'seqcuda' in d['vector_types'] and 'seqaijcusparse' in d['matrix_types']
 assert any('VecGhostUpdateBegin' in l for l in d['stack'])
def test_unpatched_cpu_pass_cuda_fail():
 d=read('ghost_probe_before')
 for r in d['runs']:
  if r['variant'].startswith('CPU'):ghost_run_gate(r)
  else:
   with pytest.raises(GateError):ghost_run_gate(r)
