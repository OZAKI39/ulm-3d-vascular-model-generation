from sv13_support import artifact
def test_live_environment_evidence_exists():
    e=artifact('gpu_environment')
    assert e['timestamp'] and e['probes']
    assert e['status'] in ('AVAILABLE','BLOCKED')
    if e['status']=='AVAILABLE':assert '4090' in e['gpu_model'] and e['gpu_memory_bytes']>0
    else:assert e['reason']
