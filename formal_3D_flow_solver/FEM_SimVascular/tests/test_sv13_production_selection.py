from sv_validation.sv13 import gpu_gate,select_production
from test_sv13_gpu_transfer_audit import good
from sv13_support import artifact,ROOT
import yaml
def cpu():
    return dict(science_equivalent=True,mass_pass=True,steady_pass=True,checkpoint_complete=True,linear_failures=0,nonlinear_failures=0)
def test_faster_gpu_with_bad_mass_rejected():
    d=good();d.update(mass_pass=False,name='GPU_A',full_run_pass=True,full_run_speedup=3)
    assert gpu_gate(d)=='GPU_REJECTED_NUMERICAL_DIFFERENCE'
    assert select_production(cpu(),d)=='CPU_EARLY_STOP_PRODUCTION'
def test_gpu_without_full_validation_not_selected():
    d=good();d.update(name='GPU_A',full_run_speedup=3)
    assert select_production(cpu(),d)=='CPU_EARLY_STOP_PRODUCTION'
def test_actual_production_config():
    d=artifact('production_selection');c=yaml.safe_load((ROOT/'configs/sv1_3/production_performance.yaml').read_text())
    assert d['status']=='PASS' and d['selected'] in ('CPU_EARLY_STOP_PRODUCTION','GPU_A_PRODUCTION','GPU_B_PRODUCTION')
    assert c['designation']=='FAST_PRODUCTION_CONFIGURATION' and c['validation_reference']['stage']=='SV1.2'
    assert c['expected_termination_mechanism']=='native STOP_SIM at a saved timestep'
