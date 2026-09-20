from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_official_fluid_gpu_smoke():
    d=actual('svmp_gpu_smoke');assert d['exit_code']==0 and d['linear_failures']==0
    assert d['velocity_finite'] and d['pressure_finite'] and d['VTU_sha256']
    cuda_types_gate(d['mat_type'],d['vec_type'])
