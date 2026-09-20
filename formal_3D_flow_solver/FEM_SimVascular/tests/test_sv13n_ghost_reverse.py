from sv13n_support import *
def test_actual_reverse_values():
 for v in ('CPU_seq','CUDA_seq','CPU_MPI'):ghost_run_gate(ghost(v))
 r=ghost('CUDA_MPI')
 assert any(float(a)!=float(e) for phase,rank,index,a,e in r['values'] if phase=='reverse')
 assert not r['accepted']  # Optional unsupported MPI CUDA failure cannot become PASS.
def test_wrong_reverse_rejected():
 r=synthetic_ghost();row=next(v for v in r['values'] if v[0]=='reverse');row[3]=row[4]=str(float(row[3])+1)
 with pytest.raises(GateError):ghost_run_gate(r)
