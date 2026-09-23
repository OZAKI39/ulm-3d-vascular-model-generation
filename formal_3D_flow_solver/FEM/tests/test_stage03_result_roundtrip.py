from stage03_helpers import *
from fem3d.audit import sha256

def test_new_process_restores_solution_and_reintegrates_without_second_solve():
    require_solution()
    r=read('outputs/stage03/reference/qc/roundtrip.json')
    assert r['status']=='PASS' and r['new_process'] and r['pid']!=r['solve_pid']
    assert not r['pde_solve_called'] and r['exact_owned_coefficients_restored']
    assert r['checkpoint_sha256']==sha256(BASE/'checkpoints/primary.npz')
    assert {'inlet','outlet_01','outlet_02','outlet_03','lambda','velocity_L2_norm_m_pow_2p5_s','pressure_integral_pa_m3'}<=r['comparisons'].keys()
    assert all(v['pass'] for v in r['comparisons'].values())
    assert len([k for k in r['comparisons'] if k.startswith('residual_')])>=18
