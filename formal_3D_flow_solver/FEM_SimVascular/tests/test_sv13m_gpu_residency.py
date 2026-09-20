from sv13m_support import *
from sv_validation.sv13m import *
def test_actual_profile_separate_from_benchmark():
 d=accepted('gpu_residency');assert d['profiling_case'].startswith('PROFILE_') and d['evidence']
 assert d['peak_VRAM_MiB']>0
