from sv13n_support import *
def test_requested_one_rank_scope():
 d=accepted('gpu_ghost_target');assert d['ranks']==1 and not d['MPI_CUDA_supported']
 for v in d['target_variants']:ghost_run_gate(ghost(v))
