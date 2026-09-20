from sv13n_support import *
def test_residency_is_evidenced():
 d=accepted('gpu_residency');assert d['evidence']
 for k in ('matrix','vectors','ASM','ILU_factor','triangular_solve','ghost_storage'):assert k in d['residency']
