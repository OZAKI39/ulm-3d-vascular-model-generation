from sv13n_support import *
@pytest.mark.parametrize('variant',['CPU_seq','CUDA_seq'])
def test_sequential_ghost(variant):ghost_run_gate(ghost(variant))
