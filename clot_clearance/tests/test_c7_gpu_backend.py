from scripts import c7_early_branch as a


def test_gpu_gate_is_explicit_and_cpu_remains_canonical():
    r=a.read(a.OUT/'GPU_BACKEND.json');assert r['CPU_scientific_canonical']
    assert len(r['CPU_full_bond_batch_seconds'])>=10 and r['CPU_batch_median_seconds']>0
    if r['used']:
        assert r['energy_relative']<=1e-11 and r['q_scaled_error']<=1e-10
        assert r['classification_IDs_equal'] and r['speedup']>=1.3
    else:
        assert r['status'] in ['RUNTIME_UNAVAILABLE','AVAILABLE_NOT_BENEFICIAL']
        assert r['reason'] and r['GPU_peak_VRAM_bytes'] is None
